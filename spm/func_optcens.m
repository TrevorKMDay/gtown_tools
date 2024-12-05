function [out_cens, out_int, out_both] = ...
    func_optcens(ospms, approach, docalc, o_minrem, o_maxrem, zeropad, adapt)

% Inputs:
%
%          ospms        a char array with the name(s) of one or several
%                       estimated spm12 SPM.mat files;
%                       (can also be one directory which will be searched,
%                       including all (!) subdirectories); alternatively, a
%                       cell array of images pathnames can be passed in
%                       column 1, with a matching mask in the second column
%                       (optional, if not provided, one will be generated);
%                       this will only use delta-D-var and STS (a.k.a.
%                       approach 2); required, will prompt if empty
%
%          approach     approch to take, different for task-based (1,
%                       default) or for resting state images (2);
%                       effectively defines thresholds for delta-D-var and
%                       STS and whether to calculate R2 and AIC (yes for 1,
%                       no for 2); required, but may not prompt if empty as
%                       the value can also be set in the defaults, below;
%                       also see the settings section, below, for the
%                       specific threshold settings
%
%          docalc       specify how you want the identified outlier volumes
%                       to be handled: docalc = 1 indicates
%                       that outlying datapoints are censored in the design
%                       matrix, while docalc = 2 will do a linear
%                       interpolation of outlying datapints in the raw data
%                       instead; docalc = 3 [default] will do both, alone
%                       and in combination; additionally, this parameter
%                       allows to set the interpolation method. If the
%                       inpaints approach is requested instead of the
%                       default (art_repair), use -2 and -3, respectively;
%                       required, but may not prompt if empty as the value
%                       can also be set in the defaults, below
%
%          o_minrem     minimum number of images to remove (in percent); if
%                       specified, this percentage of
%                       images will always be removed, even if not formally
%                       designated as outliers; required, but may not
%                       prompt if empty as the value can also be set in the
%                       defaults, below
%
%          o_maxrem     maximum number of images to remove (in percent);
%                       even if designated as outliers,
%                       no more images will be removed; required, but may
%                       not prompt if empty as the value can also be set in
%                       the defaults, below
%
%          zeropad      do (1) or do not (0) make sure that, if several
%                       files are supplied, the same number
%                       of datapoints is censored in each; this is achieved
%                       by adding dummy regressors to those analyses where
%                       fewer outliers are present so that over all
%                       analyses, identical degrees of freedom are ensured;
%                       integer; required, but may not prompt if empty as
%                       the value can also be set in the defaults, below
%
%          adapt        do (1) or do not (0) zeropad in an adaptive manner:
%                       if several files are supplied (cf. zeropad),
%                       this will only zeropad to the minimum number
%                       required across the studies supplied; required, but
%                       may not prompt if empty as the value can also be
%                       set in the defaults, below
%
% Outputs:8
%
%          out_cens     a cell array with the name(s) of the optimally
%                       censored, modified spm12 SPM.mat files;
%                       will be empty if docalc is set to 2
%
%          out_int      a cell array with the name(s) of the optimally
%                       interpolated, modified spm12 SPM.mat files;
%                       will be empty if docalc is set to 1
%
%          out_both     a cell array with the name(s) of the optimally
%                       interpolated & censored, modified spm12 SPM.mat 
%                       files; will be empty if docalc is set to 1 or 2
%
% Idea and implementation by Marko Wilke, see below for version
% information.
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details. You should have received a copy of the
% GNU General Public License along with this program. If not, see
% <http://www.gnu.org/licenses/>.
%


%% =========================================================================================================
%                                          Preludes: settings, inputs, etc.
% ==========================================================================================================


% version information v1.0   for initial public distribution; 2019-03-08
% v1.1   replacing calls to the statistics and the curve fitting toolbox;
% 2019-05-10 v1.11  minor bugfixes, mainly for assessing moved datasets;
% 2019-05-17 v1.12  included half IQR approach, updated documentation;
% 2019-05-20 v1.13  fixed bugs with older Matlab versions and when running
% approach 2; 
% 2019-05-23 v1.14  allow for processing of SPM.mat files with
% no conditions in approach 2; 
% 2019-05-24 v1.2   allow for processing of
% image sessions without associated SPM.mat files in approach 2; 
% 2019-07-22 v1.21  code prettifications and better cooperation with the toolbox
% config file; 
% 2020-01-10 v1.22  removed resetting aic_min in case of no
% improvement due to downstream error; 
% 2021-11-09 v1.23  corrected syntax
% errors in part for parametric modulation; some code prettifications;
% 2022-08-29 v1.24  added an option to select rp-files based on name
% instead of order; 2022-09-03


% default settings
approach_d  = 1;                   % default value for approach; if set here, will not ask below (else set to [])
o_maxrem_d  = 50;                  % same for o_maxrem
o_minrem_d  = 0;                   % same for o_minrem
zeropad_d   = 0;                   % same for zeropad
adapt_d     = 0;                   % same for adapt
docalc_d    = 3;                   % same for docalc (0 will skip, 1 will censor in design matrix, 2 will interpolate timeseries, and 3 will do both)


% settings dependent on approach
thr_s_1     = [0.3 1.5];           % defaults for scan-to-scan displacement in approach 1 (task based fMRI)
thr_s_2     = [0   0.3];           % defaults for scan-to-scan displacement in approach 2 (resting state fMRI)
% NB: first value is the lower threshold below which no image will be
% censored; second value is the upper threshold above which all images will
% be censored;
thr_d_1     = [0.05 50 15];        % delta%D-var defaults for approach 1 (task based fMRI)
thr_d_2     = [0.05 50  5];        % delta%D-var defaults for approach 2 (resting state fMRI)
% NB: first value is statistical threshold, second value is percentage of
% slices required above this level third value is the absolute percent
% deviation (fallback) factor


% debugging, display and processing settings
skip        = 0;                   % if set, skip analyses if previous results detected
rp_1        = 1;                   % if set (1), use first of > 1 rp_files, otherwise (0), select based on name
davg        = 65;                  % average cortical distance (in mm); cf. Wilke, P1, 2014
f_tukey     = 1.5;                 % factor for initial Tukey outlier criterion (1.5 = outlier; 3 = far outlier)
f_tukey2    = [25 1.5];            % factors for second Tukey outlier criterion, to be applied if more than f_tukey2(1) outliers are removed
dmask       = 1;                   % delta-D-var: mask timeseries by existing mask?
dd          = 1/3;                 % delta%D-var: power of transformation
f_aic       = 2;                   % factor for AIC (allow a max of "this much worse" AIC value)
perc_warn   = 40;                  % warn user if percentage of removed datapoints exceeds this value; to disable, set to 100
pval        = [NaN NaN NaN];       % p-values for quality control overlap maps for FWE, FDR, or none (set to NaN to skip)


% set half IQR approach and check for statistics toolbox (if still
% necessary)
hIQR = 1;
if exist('bootstrp.m', 'file') ~= 0 && hIQR == 0,  doboot = 1;  else,  doboot = 0;  end
% doboot = 0;  % enable to disable bootstrapping manually


% warnings, variables...
warning off MATLAB:MKDIR:DirectoryExists
warning off MATLAB:Figure:FigureSavedToMATFile
warning off MATLAB:singularMatrix


% say hello
ori= pwd;
clc;
dirname = mfilename;
try

    disp(' ');
    disp(['... Hello ' spm('GetUser') ', welcome to ' dirname '...']);
    disp(' ');
    dirname = 'optcens';
    % dirname = dirname(4:end)

catch

    disp('   ... there was an error calling SPM, please make sure it is installed; aborting...');
    return;

end
myout = [dirname '_qc.txt'];


% check that spm is running and that we have the correct version
if isempty(spm('Show')) && isempty(getCurrentTask)

    disp('   ... SPM needs to be up and running while executing this function; aborting...');
    return;

end
ver = spm('ver');
if str2double(ver(4:5)) < 12

    disp(['   ... you need to use at least SPM12 to run this script, not ' ver '; aborting...']);
    return;

end


% inputs: SPM mat file OR directory OR images only
if nargin == 0 || isempty(ospms)

    ospms = spm_select(Inf, 'mat',  'Select estimated SPM.mat file(s) to assess (or done to select a directory)', {}, pwd, 'SPM.mat');
    if isempty(ospms)

        ospms = spm_select([0 1], 'dir',  'Select directory to search for estimated SPM.mat files', {}, pwd);

    end
    if isempty(ospms)

        temp1 = spm_select(Inf, 'image',  'Select images belonging to one session', {}, pwd);
        if ~isempty(temp1)

            temp2 = spm_select([0 1], 'image',  'Optional: select mask image belonging to this session', {}, pwd);
            ospms{1,1} = temp1;
            if ~isempty(temp2),  ospms{1,2} = temp2;  end

        end

    end
    if isempty(ospms),  disp('... nothing selected, aborting...');  return;  end

end


% check inpput (may only be images) or find SPM.mat file; if directory,
% search all (!) subdirectories
if iscell(ospms)

    % adapt settings
    nmats       = size(ospms,1);
    approach_d  = 2;
    approach    = 2;
    docalc_d    = 2;
    docalc      = 2;
    zeropad_d   = 0;
    zeropad     = 0;
    adapt       = 0;
    pval(1:end) = NaN;


    % looks like images only were passed; verify
    try

        for j = 1:nmats

            spm_check_orientations(spm_vol(char(ospms{j,1})));

        end

    catch

        disp(['      ... sorry, some of the selected images (likely in session ' num2str(j) ') were not in identical orientation, aborting...']);

    end

else

    if isfolder(ospms(1,:))  && size(ospms,1) > 1,  error('Input error: you may pass more than one file, but only one directory,');  end
    if isfolder(ospms)

        % search subdirs
        disp(['   ... searching ' ospms ', please wait...']);
        ospms = spm_select('FPListRec', ospms, 'SPM.mat');
        disp(['   ... found ' num2str(size(ospms,1)) ' SPM.mat files...']);

    end
    nmats = size(ospms,1);

end


% inputs: approach?
if nargin < 2 || isempty(approach)

    approach = approach_d;
    if isempty(approach)

        approach = spm_input('Use approach...','+1','m','  ... for task-based fMRI|  ... for resting state fMRI',[1 2],1);

    end

end
if approach == 1

    thr_s = thr_s_1;
    thr_d = thr_d_1;
    disp('   ... using task-based fMRI approach ...');

elseif approach == 2

    thr_s = thr_s_2;
    thr_d = thr_d_2;
    disp('   ... using resting-state fMRI approach ...');

else

    error(['Unknown data processing approach approach (' num2str(approach) ')!']);

end
dmax        = thr_d(3)*2;          % set delta-%D-var maximum value for color scaling


% inputs: do which calculations?
if nargin < 3 || isempty(docalc)

    docalc = docalc_d;
    if isempty(docalc)

        docalc = spm_input('Censor, interpolate, or both?','+1','b', 'C | I | B', [1 2 3], 0);

    end

end
if docalc < 0,  docalc = abs(docalc);  useart = 0;
else,           useart = 1;            end


% inputs: fixed minimum percentage to remove?
if nargin < 4 || isempty(o_minrem)

    o_minrem = o_minrem_d;
    if isempty(o_minrem)

        o_minrem = spm_input('Min percentage of datapoints to remove?','+1','n', o_minrem);

    end

end


% inputs: maximum percentage of datapoints to remove
if nargin < 5 || isempty(o_maxrem)

    o_maxrem = o_maxrem_d;
    if isempty(o_maxrem)

        o_maxrem = spm_input('Max percentage of datapoints to remove?','+1','n', o_maxrem);

    end

end
disp(['   ... removing ' num2str(o_minrem) '-' num2str(o_maxrem) '% (MIN-MAX) of images...']);


% inputs: zeropad regressors?
if nargin < 6 || isempty(zeropad)

    zeropad = zeropad_d;
    if isempty(zeropad)

        zeropad = spm_input(['Zeropad regressors (up to ' num2str(o_maxrem) ')?'],'+1','b', 'Yes | No', [1 0], 0);

    end

end


% if several files: adaptive zeropadding (to joint number of covariates)?
if nargin < 7 || isempty(adapt)

    adapt = adapt_d;
    if isempty(adapt) && nmats > 1

        adapt = spm_input('Apply adaptive zeropadding across supplied sessions?','+1','b', 'Yes | No', [1 0], 0);

    end

end
if nmats > 1 && adapt == 1,  zeropad = -1;  end
if zeropad == 0

    disp('   ... no zeropadding of regressors will be applied...');

elseif zeropad > 0

    disp(['   ... will apply zeropadding up to ' num2str(zeropad) ' regressors...']);

elseif zeropad < 0

    disp('   ... applying adaptive zeropadding...');

end


% inform user
disp('   ... thank you, now gathering information from input files...');


% temporarily change some spm defaults to avoid annoying feedback
ocmd = spm_get_defaults('cmdline');
spm_get_defaults('cmdline', 1);


%% =========================================================================================================
%                                          Loop over input data
% ==========================================================================================================



% create storage, get going
bestcoms = zeros(1, nmats);
goon = zeros(1, nmats);
tic;
for j = 1:nmats


    % find current file(s)
    if ~iscell(ospms)

        ospm = deblank(ospms(j,:));
        [ospm_p, ~, ~, ~] = spm_fileparts(ospm);
        if isempty(ospm_p),  ospm_p = pwd;  end
        disp(['      ... analyzing model ' num2str(j) '/' num2str(nmats) ' from ' ospm_p '...']);


        % double check SPM.mat: still available ? estimated? One of ours?
        % Correct version?
        if exist(ospm, 'file') ~= 2

            disp('         ... ERROR: file not found anymore; skipping...');
            continue;

        end
        try

            ospm = load(ospm);
            ospm.SPM.xY.VY;
            ospm.SPM.swd;

        catch

            disp('         ... ERROR: this model has not been estimated yet, or is not readable; skipping...');
            continue;

        end
        if ~isempty(strfind(ospm.SPM.SPMid, 'mw_'))

            disp('         ... ERROR: this model already is an optimally censored model; skipping...');
            continue;

        end
        if isempty(strfind(ospm.SPM.SPMid, 'SPM12'))

            disp('         ... ERROR: this model was not specified using SPM12; skipping...');
            continue;

        end


        % check if a condition/contrast is specified (necessary only for
        % task-based approach)
        try

            if size(ospm.SPM.xCon,2) == 0 && approach == 1

                disp('         ... WARNING: this model has no contrast/condition set, trying resting state approach...');
                approach = 2;
                thr_s = thr_s_2;
                thr_d = thr_d_2;

            end

        catch

            if approach == 1

                disp('         ... WARNING: this model has no contrast/condition set, trying resting state approach...');
                approach = 2;
                thr_s = thr_s_2;
                thr_d = thr_d_2;

            end

        end


        % check if first level model to begin with
        if isfield(ospm.SPM, 'Sess') == 0

            disp('         ... ERROR: this model does not seem to be a first level model; skipping...');
            continue;

        end


        % check if this is a single session model
        if numel(ospm.SPM.nscan) > 1

            disp('         ... ERROR: this model is for more than one session, which is not supported; skipping...');
            continue;

        else

            nscans = ospm.SPM.nscan;

        end


        % check for previous analyses
        if exist([ospm_p filesep dirname filesep dirname '.mat'], 'file') == 2

            if skip == 1

                disp(['         ... previous results detected in ' ospm_p filesep dirname ' ; skipping calculations...']);
                goon(1,j) = 1;
                continue;

            else

                movefile([ospm_p filesep dirname filesep dirname '.mat'], [ospm_p filesep dirname filesep dirname '_old.mat']);

            end

        end


        % convert percentage to absolute number of scans
        maxrem = round(nscans * o_maxrem/100);
        minrem = round(nscans * o_minrem/100);


        % start harvesting information from original model: general ...
        clear matlabbatch matlabbatch_o;
        matlabbatch{1}.spm.stats.fmri_spec.timing.units   = ospm.SPM.xBF.UNITS;
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = ospm.SPM.xX.K.RT;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = ospm.SPM.xBF.T;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = ospm.SPM.xBF.T0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf       = ospm.SPM.xX.K.HParam;
        matlabbatch{1}.spm.stats.fmri_spec.volt           = ospm.SPM.xBF.Volterra;
        matlabbatch{1}.spm.stats.fmri_spec.global         = ospm.SPM.xGX.iGXcalc;


        % ... scans (make sure they are accessible) ...
        try

            temp = spm_vol(ospm.SPM.xY.P);
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(ospm.SPM.xY.P);

        catch

            disp('      ... sorry, the scans from this model were not found at the original location...');
            disp(['          (' spm_str_manip(ospm.SPM.xY.P(1,:),'h') ')']);
            disp('      ... you now have the option to respecify, but BE SURE TO SELECT THE CORRECT IMAGES!');
            temp = spm_select(nscans,'image','Images not found in original location, select now or close to skip...',[],pwd,'.*');
            if isempty(temp)

                disp('      ... okay, skipping this session...');
                continue;

            end
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(temp);
            try

                spm_check_orientations(spm_vol(temp));

            catch

                disp('      ... sorry, selected images were not in identical orientation, skipping this session...');
                continue;
            end

        end


        % ... conditions ...
        for i = 1:size(ospm.SPM.Sess.U,2)

            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).name     = char(ospm.SPM.Sess.U(i).name);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).onset    = ospm.SPM.Sess.U(i).ons;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).duration = ospm.SPM.Sess.U(i).dur;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).tmod     = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).orth     = ospm.SPM.Sess.U(i).orth;

            if strcmp(ospm.SPM.Sess.U(i).P.name, 'time')

                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).tmod = ospm.SPM.Sess.U(i).P(1).h;

            elseif strcmp(ospm.SPM.Sess.U(i).P.name, 'none')

                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).pmod = struct('name', {}, 'param', {}, 'poly', {});

            else

                for ii = 1:numel(ospm.SPM.Sess.U(i).P)

                    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).pmod(ii).name  = ospm.SPM.Sess.U(i).P(ii).name;
                    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).pmod(ii).param = ospm.SPM.Sess.U(i).P(ii).P;
                    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).pmod(ii).poly  = ospm.SPM.Sess.U(i).P(ii).h;

                end

            end

        end


        % ... regressors ...
        if isempty(ospm.SPM.Sess.C.name)

            matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});

        else

            for ii = 1:numel(ospm.SPM.Sess.C.name)

                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(ii).name = ospm.SPM.Sess.C.name{ii};
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(ii).val  = ospm.SPM.Sess.C.C(:,ii);

            end

        end


        % ... response function ...
        if ~isempty(strfind(ospm.SPM.xBF.name, 'hrf'))

            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];

        end
        if ~isempty(strfind(ospm.SPM.xBF.name, 'time'))

            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs(1) = 1;

        end
        if ~isempty(strfind(ospm.SPM.xBF.name, 'dispersion'))

            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs(2) = 1;

        end
        if ~isempty(strfind(ospm.SPM.xBF.name, 'Fourier')) && isempty(strfind(ospm.SPM.xBF.name, 'Hanning'))

            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier.length = ospm.SPM.xBF.length;
            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier.order  = ospm.SPM.xBF.order;

        end
        if ~isempty(strfind(ospm.SPM.xBF.name, 'Hanning'))

            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier_han.length = ospm.SPM.xBF.length;
            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier_han.order  = ospm.SPM.xBF.order;

        end
        if ~isempty(strfind(ospm.SPM.xBF.name, 'Gamma'))

            matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.length = ospm.SPM.xBF.length;
            matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.order  = ospm.SPM.xBF.order;

        end
        if ~isempty(strfind(ospm.SPM.xBF.name, 'Finite'))

            matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length = ospm.SPM.xBF.length;
            matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order  = ospm.SPM.xBF.order;

        end


        % ... masking ...
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = ospm.SPM.xM.gMT;
        if isempty(ospm.SPM.xM.VM)

            matlabbatch{1}.spm.stats.fmri_spec.mask = {''};

        else

            matlabbatch{1}.spm.stats.fmri_spec.mask = cellstr([ospm_p filesep ospm.SPM.VM.fname]);

        end


        % ... AR options (somewhat convoluted due to an inconsistency in
        % the naming convention, see Mail from Torben, July 30, 2016)
        if strcmp(ospm.SPM.xVi.form, 'AR(0.2)')

            matlabbatch{1}.spm.stats.fmri_spec.cvi  = 'AR(1)';

        elseif strcmp(ospm.SPM.xVi.form,'FAST')

            matlabbatch{1}.spm.stats.fmri_spec.cvi  = 'FAST';

        else

            matlabbatch{1}.spm.stats.fmri_spec.cvi  = 'none';

        end


        % ... and finally, directory
        matlabbatch{1}.spm.stats.fmri_spec.dir = {[ospm_p filesep dirname]};
        try

            if exist([ospm_p filesep dirname],'dir') == 7 && skip == 0,  rmdir([ospm_p filesep dirname], 's');  end

        catch

            error(['Sorry, old results were found in ' [ospm_p filesep dirname] ' and could not be removed; please repair and try again']);

        end
        try

            mkdir([ospm_p filesep dirname]);

        catch

            error(['Sorry, the directory ' [ospm_p filesep dirname] ' could not be created; please repair and try again']);
        end
        cd([ospm_p filesep dirname]);


        % for completeness' sake
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi     = {''};
        matlabbatch{1}.spm.stats.fmri_spec.fact           = struct('name', {}, 'levels', {});


        % in case only images were passed
    else


        % harvesting information from passed images (need to be in the same
        % directory for this to work)
        ospm = deblank(ospms{j,1}(1,:));
        [ospm_p, ~, ~, ~] = spm_fileparts(ospm);
        if isempty(ospm_p),  ospm_p = pwd;  end
        disp(['      ... analyzing dataset ' num2str(j) '/' num2str(nmats) ' from ' ospm_p '...']);


        % directory
        try

            if exist([ospm_p filesep dirname],'dir') == 7 && skip == 0,  rmdir([ospm_p filesep dirname], 's');  end

        catch

            error(['Sorry, old results were found in ' [ospm_p filesep dirname] ' and could not be removed; please repair and try again']);

        end
        try

            mkdir([ospm_p filesep dirname]);

        catch

            error(['Sorry, the directory ' [ospm_p filesep dirname] ' could not be created; please repair and try again']);
        end
        cd([ospm_p filesep dirname]);


        % if not passed, generate mask
        try

            mask = deblank(ospms{j,2});
            mask = spm_read_vols(spm_vol(mask));

        catch

            mask = mw_brainmask2(deblank(ospms{j,1}), 1);
            mask = spm_read_vols(spm_vol(mask));

        end


        % gather, prepare required information
        nscans       = size(deblank(ospms{j,1}),1);
        matlabbatch  = [];
        ospm         = [];


        % convert percentage to absolute number of scans
        maxrem = round(nscans * o_maxrem/100);
        minrem = round(nscans * o_minrem/100);

    end


    %% =========================================================================================================
    %                                          Calculations
    % ==========================================================================================================


    % storage (depending on approach)
    if approach == 1

        store_t       = zeros([ospm.SPM.VM.dim maxrem size(ospm.SPM.xCon,2)]);
        store_aic     = zeros([ospm.SPM.VM.dim maxrem]);

    else

        store_t       = [];
        store_aic     = [];

    end
    thrs          = zeros(maxrem, 3);
    matlabbatch_o = matlabbatch;


    % load data, using existing mask (if necessary)
    if ~exist('mask', 'var')

        mask = spm_read_vols(spm_vol([ospm_p filesep ospm.SPM.VM.fname]));

    end
    curr = find(mask>0);
    Y = zeros(nscans, numel(curr));
    for i = 1:nscans

        % Y(i,:) = spm_data_read(spm_vol(ospm.SPM.xY.P(i,:)),curr);
        try

            Y(i,:) = spm_data_read(spm_vol(char(matlabbatch{1}.spm.stats.fmri_spec.sess.scans(i,:))),curr);

        catch

            Y(i,:) = spm_data_read(spm_vol(char(ospms{j,1}(i,:))),curr);

        end

    end
    Y_o = Y;


    %% ====== calculate % change in D-var according to Afyouni & Nichols ===============
    disp(['         ... assessing data quality (' char(916) '%D-var), please wait...']);


    % mask timeseries by existing mask? If not, need to read in again
    if dmask == 0

        curr = 1:numel(mask);
        Y = zeros(nscans, numel(curr));
        for i = 1:nscans

            try

                Y(i,:) = spm_data_read(spm_vol(char(matlabbatch{1}.spm.stats.fmri_spec.sess.scans(i,:))),curr);

            catch

                Y(i,:) = spm_data_read(spm_vol(char(ospms{j,1}(i,:))),curr);

            end

        end

    end


    % flip, remove empty and NaN timeseries
    Y = Y';
    nan_idx    = find(isnan(sum(Y,2)));
    zeros_idx  = find(sum(Y,2)==0);
    Y([nan_idx;zeros_idx],:) = [];
    curr([nan_idx;zeros_idx]) = [];


    % de-mean, then difference in the time domain
    Y  = Y-repmat(mean(Y,2),[1,size(Y,2)]);
    DY = diff(Y,1,2);


    % slicewise calculation: average timecourses, prepare for calculation
    Y_MS   = mean(Y.^2, 2);
    DVARS2 = DY.^2;
    [~,~,slices] = ind2sub(size(mask),curr);
    temp   = zeros(max(slices), size(DVARS2,2));
    temp2  = zeros(max(slices), 1);
    div    = zeros(max(slices),1);
    for  i = 1:numel(slices)

        temp(slices(i),:)  = temp(slices(i),:)  + DVARS2(i,:);
        temp2(slices(i),:) = temp2(slices(i),:) + Y_MS(i,:);
        div(slices(i),1)   = div(slices(i),1) + 1;

    end
    div(div == 0) = 1;
    DVARS2        = temp ./ repmat(div, [1 size(temp,2)]);
    Y_MS          = temp2 ./ div;


    % calculate, clean Delta-D-var DeltapDvar =
    % (DVARS2-median(DVARS2,2))./(4*Y_MS)*100; % may lead to errors in old
    % Matlab releases
    DeltapDvar = (DVARS2-repmat(median(DVARS2, 2), 1, size(DVARS2, 2)))./(4*repmat(Y_MS, 1, size(DVARS2, 2)))*100;
    DeltapDvar(isnan(DeltapDvar)) = 0;  % NaN may occurr in case of missing data


    % now on to statistical inference
    Mn  = median(DVARS2,2);
    Z   = DVARS2.^dd;
    M_Z = median(Z,2);


    % use publically available replacement to avoid dependency on the
    % statistics toolbox... replace NaN with Inf so that results correspond
    % to Matlab's quantile function
    Z(isnan(Z)) = Inf;
    H_IQRsd = (ah_quantile(Z',0.5)-ah_quantile(Z',0.25))/1.349*2;
    dvars_p = NaN(max(slices), size(DVARS2,2));
    for i = 1:max(slices)

        Va        = sqrt((1/dd*M_Z(i)^(1/dd-1)*H_IQRsd(i))^2);
        if Va == 0,  continue;  end  % may happen for empty slices
        arg1      = 2*Mn(i) / Va^2 * DVARS2(i,:);
        arg2      = 2*Mn(i)^2 / Va^2;

        % another new file to avoid chi2cdf, kindly provided by Guillaume
        dvars_p(i,:) = gf_spm_Gcdf(arg1,arg2/2,1/2,'upper');

    end


    % actually find significant delta%D-var outliers (in at least thr_d(2)%
    % of slices!)
    corrfact = round(max(slices) .* thr_d(2) / 100);
    ind_d = sum((dvars_p * max(slices)) < thr_d(1))';  % Bonferroni-correct for # slices
    ind_d = [0 ind_d'];


    % fallback position: use "delta%D-var fallback value" if this is
    % stricter
    temp = [0 sum(DeltapDvar > thr_d(3))];
    if sum(temp >= corrfact) > sum(ind_d >= corrfact)

        ind_d = temp;
        dfb   = 1;

    else

        dfb   = 0;

    end
    outl_d = sum(ind_d >= corrfact);
    [~, rem_d] = sort(ind_d(:),'descend');


    %% ====== now STS: scan-to-scan displacement at davg ========================
    try

        [odat, oname, ~, ~] = spm_fileparts(char(matlabbatch{1}.spm.stats.fmri_spec.sess.scans(1,:)));

    catch

        [odat, oname, ~, ~] = spm_fileparts(char(ospms{j,1}(1,:)));

    end
    ind_s = zeros(1, nscans);

    try

        rps = spm_select('FPlist', odat, '^rp_.*.txt');
        if ~isempty(rps)

            disp('         ... assessing data quality (STS), please wait...');

        else

            %error(['Realignment parameter file (^rp_.*.txt) not found in '
            %odat '...']); AG 2023: Instead of throwing an error, if not
            %found in odat, search one directory up
            oneup = fileparts(odat);
            rps = spm_select('FPlist', oneup, '^rp_.*.txt');

        end

        % this may be an issue so include a workaround
        if size(rps, 1) > 1

            revert = 0;
            if rp_1 == 0

                temp = [odat filesep 'rp_' oname '.txt'];
                if exist(temp, 'file') == 2

                    disp(['             ... WARNING: more than one rp_.*.txt file found, proceeding with closest name match (rp_' oname '.txt)...']);
                    rps = temp;

                else

                    disp(['             ... WARNING: rp_file selection based on name did not work, reverting to selecting first file...']);
                    revert = 1;

                end

            end
            if rp_1 == 1 || revert == 1

                [~, temp, ~,~] = spm_fileparts(rps(1,:));
                disp(['             ... WARNING: more than one rp_.*.txt file found, proceeding with first file (' temp '.txt)...']);
                rps = rps(1,:);

            end

        end

        pr = load(rps);
        for i = 2:size(pr,1)

            dx = (pr(i,1) - pr(i-1,1)) + (pr(i,4) .* davg - pr(i-1,4) .* davg);
            dy = (pr(i,2) - pr(i-1,2)) + (pr(i,5) .* davg - pr(i-1,5) .* davg);
            dz = (pr(i,3) - pr(i-1,3)) + (pr(i,6) .* davg - pr(i-1,6) .* davg);
            ind_s(1,i) = sqrt(dx.^2 + dy.^2 + dz.^2);

        end


        % find STS outliers  according to Tukey (bootstrapping if possible)
        if doboot == 1

            q3 = mean(bootstrp(10000, @quantile, ind_s, 0.75));
            ubound_s = q3 + f_tukey * mean(bootstrp(10000, @iqr, ind_s));

        else

            % potentially branch further to use half IQR
            if hIQR == 1

                myIQR = (ah_quantile(ind_s, 0.50)-ah_quantile(ind_s, 0.25)) .* 2;

            else

                myIQR = ah_quantile(ind_s, 0.75)-ah_quantile(ind_s, 0.25);

            end

            q3 = ah_quantile(ind_s, 0.75);
            ubound_s = q3 + f_tukey * myIQR;

        end


        % check if relative or absolute bounds
        if any(thr_s > 9)


            vs = mean(sqrt(sum(ospm.SPM.VM.mat(1:3,1:3).^2)));
            if thr_s(1) > 9,  thr_s(1) = thr_s(1) * vs / 100;  end
            if thr_s(2) > 9,  thr_s(2) = thr_s(2) * vs / 100;  end

        end


        % respect pre-specified lower, enforce pre-specified upper
        % threshold
        ubound_s = max([ubound_s thr_s(1)]);
        ubound_s = min([ubound_s thr_s(2)]);


        % apply
        outl_s = sum(ind_s > ubound_s);
        [~, rem_s] = sort(ind_s(:),'descend');


        % potentially do second round if many timepoints were removed
        if ceil(outl_s/nscans*100) > f_tukey2(1)

            disp(['            ... removing > ' num2str(round(f_tukey2(1))) ' datapoints, adapting Tukey threshold for STS, please wait...']);
            temp = ind_s;
            temp(temp > ubound_s) = [];
            if doboot == 1

                q3 = mean(bootstrp(10000, @quantile, temp, 0.75));
                ubound_s = q3 + f_tukey2(2) * mean(bootstrp(10000, @iqr, temp));

            else

                if hIQR == 1

                    myIQR = (ah_quantile(temp, 0.50)-ah_quantile(temp, 0.25)) .* 2;

                else

                    myIQR = ah_quantile(temp, 0.75)-ah_quantile(temp, 0.25);

                end

                q3 = ah_quantile(temp, 0.75);
                ubound_s = q3 + f_tukey * myIQR;

            end
            outl_s = sum(ind_s > ubound_s);

        end

    catch

        disp(['            ... realignment parameters not found in ' odat ', proceeding without STS...']);
        outl_s = 0;
        rem_s = ones(1, nscans)';
        ind_s = zeros(1, nscans);
        ubound_s = 'N/A';

    end


    %% ====== calculate R2 to assess each scan's individual contribution ==========
    % this may not be necessary
    if approach == 1


        % allow to find out the hard way
        try


            % run original model for reference
            disp('         ... assessing explained variance (R2), please wait...');
            X = ospm.SPM.xX.X;
            Y = Y_o;


            % estimate model and get residuals (code from
            % www.sbirc.ed.ac.uk/cyril/glm/in_depth/cp_glm1.html)
            B = inv(X'*X)*X'*Y;
            Yhat = X*B;
            Res  = Y - Yhat;
            try

                Y_dm = Y-mean(Y);

            catch

                Y_dm = Y-repmat(mean(Y), size(Y, 1), 1);  % another workaround for older Matlab versions

            end


            % calculate original R2
            SStotal  = (sum(abs(Y_dm).^2).^(1/2)) .^2;  % get around norm limitation on matrices
            SSerror  = (sum(abs(Res).^2).^(1/2))  .^2;  % dito
            ori_r2    = 1 - SSerror./SStotal;


            % double check
            if sum(isnan(ori_r2)) == numel(ori_r2), error('error in calculating original r2...');  end


            % calculate original corrected AIC
            try

                vol = spm_read_vols(spm_vol([ospm_p filesep 'ResMS.img']));

            catch

                vol = spm_read_vols(spm_vol([ospm_p filesep 'ResMS.nii']));

            end
            k   = repmat(size(ospm.SPM.xX.X,2)-1, ospm.SPM.VM.dim);
            n   = repmat(nscans, ospm.SPM.VM.dim);
            df  = repmat(ospm.SPM.xX.erdf, ospm.SPM.VM.dim);
            ori_aic = 2*k + n .* (log((vol .* df) ./ n)) + (2.*k.*(k+1) ./ (n-k-1));
            ori_aic(isnan(ori_aic)) = [];
            % ori_aic = sum(sum(sum(ori_aic))) / nnz(ori_aic);


            % now estimate each scan's contribution
            spm_progress_bar('Init',nscans, 'Explained variance', ['Iterations 1:' num2str(nscans)])
            inds = eye(nscans);
            ind_r2 = zeros(1,nscans);
            for i = 1:nscans

                % inform user
                spm_progress_bar('set', i);


                % get design, include regressor to current image
                X = ospm.SPM.xX.X;  X(:,end) = inds(:,i);  X = [X ones(nscans,1)];


                % estimate model and get residuals
                B = inv(X'*X)*X'*Y;
                Yhat = X*B;
                Res  = Y - Yhat;
                try

                    Y_dm = Y-mean(Y);

                catch

                    Y_dm = Y-repmat(mean(Y), size(Y, 1), 1);

                end


                % calculate new R2 from current model
                SStotal  = (sum(abs(Y_dm).^2).^(1/2)) .^2;
                SSerror  = (sum(abs(Res).^2).^(1/2))  .^2;
                ind_r2(1,i) = mean(1 - SSerror./SStotal);

            end
            ind_r2 = ind_r2 ./ mean(ori_r2);
            spm_progress_bar('Clear');


            % find R2 outliers according to Tukey (using a botstrap if
            % possible)
            if doboot == 1

                q3 = mean(bootstrp(10000, @quantile, ind_r2, 0.75));
                ubound_r2 = q3 + f_tukey * mean(bootstrp(10000, @iqr, ind_r2));

            else

                if hIQR == 1

                    myIQR = (ah_quantile(ind_r2, 0.50)-ah_quantile(ind_r2, 0.25)) .* 2;

                else

                    myIQR = ah_quantile(ind_r2, 0.75)-ah_quantile(ind_r2, 0.25);

                end

                q3 = ah_quantile(ind_r2, 0.75);
                ubound_r2 = q3 + f_tukey * myIQR;

            end
            outl_r2 = sum(ind_r2 > ubound_r2);
            [~, rem_r2] = sort(ind_r2(:),'descend');


            % potentially do second round if many timepoints were removed
            if ceil(outl_r2/nscans*100) > f_tukey2(1)

                disp(['            ... removing > ' num2str(round(f_tukey2(1))) ' datapoints, adapting Tukey threshold for R2, please wait...']);
                temp = ind_r2;
                temp(temp > ubound_r2) = [];
                if doboot == 1

                    q3 = mean(bootstrp(10000, @quantile, temp, 0.75));
                    ubound_r2 = q3 + f_tukey2(2) * mean(bootstrp(10000, @iqr, temp));

                else

                    if hIQR == 1

                        myIQR = (ah_quantile(temp, 0.50)-ah_quantile(temp, 0.25)) .* 2;

                    else

                        myIQR = ah_quantile(temp, 0.75)-ah_quantile(temp, 0.25);

                    end

                    q3 = ah_quantile(temp, 0.75);
                    ubound_r2 = q3 + f_tukey2(2) * myIQR;

                end
                outl_r2 = sum(ind_r2 > ubound_r2);

            end

        catch


            % looks like we failed
            spm_progress_bar('Clear');
            disp('             ... assessing explained variance FAILED, reverting to resting state approach, please wait...');
            approach  = 2;
            ori_r2    = NaN;
            outl_r2   = NaN;
            ind_r2    = NaN;
            rem_r2    = NaN;
            ubound_r2 = NaN;

        end

    else

        ori_r2    = NaN;
        outl_r2   = NaN;
        ind_r2    = NaN;
        rem_r2    = NaN;
        ubound_r2 = NaN;

    end

    %% ====== now combine outlier information from all outlier parameters
    % depending on approach, consider R2 or not
    if approach == 1


        % make sure all "hard removals" are actually removed in order of
        % respective percentile
        perc_r2 = zeros(1, nscans);
        perc_d = zeros(1, nscans);
        perc_s = zeros(1, nscans);
        for i = 1:nscans

            perc_r2(1,i)  = (sum(ind_r2 < ind_r2(rem_r2(i))) + 0.5) / numel(ind_r2) * 100;
            perc_d(1,i)  = (sum(ind_d < ind_d(rem_d(i))) + 0.5) / numel(ind_d) * 100;
            perc_s(1,i)  = (sum(ind_s < ind_s(rem_s(i))) + 0.5) / numel(ind_s) * 100;

        end
        sortout = [[perc_r2(1:outl_r2); rem_r2(1:outl_r2)']' ; [perc_d(1:outl_d); rem_d(1:outl_d)']' ; [perc_s(1:outl_s); rem_s(1:outl_s)']'];
        sortout = sortrows(sortout, -1);
        sortout = unique(sortout(:,2), 'stable')';
        outl_a = numel(sortout);


        % add other "deviators" successively, avoiding repeats
        temp = [[perc_r2(outl_r2+1:end); rem_r2(outl_r2+1:end)']' ; [perc_d(outl_d+1:end); rem_d(outl_d+1:end)']' ; [perc_s(outl_s+1:end); rem_s(outl_s+1:end)']'];
        temp = sortrows(temp, -1);
        temp = unique(temp(:,2), 'stable')';
        temp(ismember(temp, sortout)) = [];
        sortout = [sortout temp];
        sortout(sortout==0) = [];

    else


        % make sure all "hard removals" from D-var and STS are actually
        % removed in order of respective percentile
        perc_d = zeros(1, nscans);
        perc_s = zeros(1, nscans);
        for i = 1:nscans

            perc_d(1,i)  = (sum(ind_d < ind_d(rem_d(i))) + 0.5) / numel(ind_d) * 100;
            perc_s(1,i)  = (sum(ind_s < ind_s(rem_s(i))) + 0.5) / numel(ind_s) * 100;

        end
        sortout = [[perc_d(1:outl_d); rem_d(1:outl_d)']' ; [perc_s(1:outl_s); rem_s(1:outl_s)']'];
        sortout = sortrows(sortout, -1);
        sortout = unique(sortout(:,2), 'stable')';
        outl_a = numel(sortout);


        % add other "deviators" successively, avoiding repeats
        temp = [[perc_d(outl_d+1:end); rem_d(outl_d+1:end)']' ; [perc_s(outl_s+1:end); rem_s(outl_s+1:end)']'];
        temp = sortrows(temp, -1);
        temp = unique(temp(:,2), 'stable')';
        temp(ismember(temp, sortout)) = [];
        sortout = [sortout temp];
        sortout(sortout==0) = [];

    end


    %% ====== generate AIC; depending on approach, this is possible - or not
    if approach == 1


        spm_progress_bar('Init',maxrem, 'Model complexity', ['Iterations 1:' num2str(maxrem)])
        disp('         ... assessing model complexity (AIC), please wait...');
        for i = 1:maxrem


            % clean slate
            matlabbatch = matlabbatch_o;


            % inform user
            spm_progress_bar('set', i);


            % index to current image
            for ii = 1:i

                ind = zeros(nscans,1);
                ind(sortout(ii),1) = 1;
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1, size(matlabbatch{1}.spm.stats.fmri_spec.sess.regress,2)+1).name = ['RN' num2str(ii)];
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1, size(matlabbatch{1}.spm.stats.fmri_spec.sess.regress,2)).val = ind;

            end


            % set up model, reload
            evalc('mw_spm_run_fmri_spec(matlabbatch{1}.spm.stats.fmri_spec);');
            load([pwd filesep 'SPM.mat']);

            % designate nuisance variables as such (see mail from Torben
            % from 2016-12-20)
            temp = 1:(SPM.xX.iB-1);
            temp(ismember(temp, ospm.SPM.xX.iC)) = [];
            SPM.xX.iC = ospm.SPM.xX.iC;    % original columns of interest
            SPM.xX.iG = temp;              % our regressors = columns of no interest


            % and estimate
            if i == 1,  Y = [];  cmask = [];  [~, SPM, Y, cmask] = evalc('mw_spm_spm(SPM, Y, cmask);');
            else,                             [~, SPM, ~, ~]     = evalc('mw_spm_spm(SPM, Y, cmask);');  end

            % disp("estimate succeeeded")

            % calculate corrected AIC
            try

                vol = spm_read_vols(spm_vol([pwd filesep 'ResMS.img']));

            catch

                vol = spm_read_vols(spm_vol([pwd filesep 'ResMS.nii']));

            end

            k   = repmat(size(SPM.xX.X,2)-1, ospm.SPM.VM.dim);
            n   = repmat(SPM.nscan, ospm.SPM.VM.dim);
            df  = repmat(SPM.xX.erdf, ospm.SPM.VM.dim);
            store_aic(:,:,:,i) = 2*k + n .* (log((vol .* df) ./ n)) + (2.*k.*(k+1) ./ (n-k-1));


            % loop over contrasts
            for ii = 1:size(ospm.SPM.xCon,2)


                % get original contrast, add zeros for our regressor (or
                % skip entirely)
                if any(~isnan(pval)) == 0,  continue;  end
                contrast = ii;
                c = [ospm.SPM.xCon(contrast).c; zeros(size(SPM.xX.xKXs.X,2)-numel(ospm.SPM.xCon(contrast).c),1)];
                cname = ospm.SPM.xCon(contrast).name;


                % add and evaluate original contrast
                if ii == 1

                    evalc('SPM.xCon = spm_FcUtil(''Set'',cname,ospm.SPM.xCon(contrast).STAT,''c'',c,SPM.xX.xKXs);');

                else

                    evalc(['SPM.xCon(' num2str(ii) ') = spm_FcUtil(''Set'',cname,ospm.SPM.xCon(contrast).STAT,''c'',c,SPM.xX.xKXs);']);

                end
                evalc('spm_contrasts(SPM);');


                % read in results: original contrast
                try

                    temp = spm_read_vols(spm_vol([pwd filesep 'spm' ospm.SPM.xCon(contrast).STAT '_' sprintf('%04.f', contrast) '.img']));

                catch

                    temp = spm_read_vols(spm_vol([pwd filesep 'spm' ospm.SPM.xCon(contrast).STAT '_' sprintf('%04.f', contrast) '.nii']));

                end
                store_t(:,:,:,i,ii) = temp;


                % get T-values corresponding to requested threshold
                % (original and new)
                if ~isnan(pval(1)),  thrs(i,1) = spm_uc(pval(1),[1 SPM.xX.erdf],SPM.xCon(contrast).STAT,SPM.xVol.R,1,SPM.xVol.S);  end
                if ~isnan(pval(2))

                    try

                        thrs(i,2) = spm_uc_FDR(pval(2),[1 SPM.xX.erdf],SPM.xCon(contrast).STAT,1,spm_vol([SPM.swd filesep 'spmT_' sprintf('%04.f', contrast) '.img']),0);

                    catch

                        thrs(i,2) = spm_uc_FDR(pval(2),[1 SPM.xX.erdf],SPM.xCon(contrast).STAT,1,spm_vol([SPM.swd filesep 'spmT_' sprintf('%04.f', contrast) '.nii']),0);

                    end

                end
                if ~isnan(pval(3)),  thrs(i,3) = spm_u(pval(3),[1 SPM.xX.erdf],SPM.xCon(contrast).STAT);  end

            end

        end
        spm_progress_bar('Clear');


        % find aic
        aic = zeros(1,maxrem);
        for i = 1:maxrem

            temp = store_aic(:,:,:,i);
            temp(isnan(temp)) = [];
            aic(1,i) = sum(sum(sum(temp))) / nnz(temp);

        end


        % relate to original, apply moving average
        aic = aic ./ (sum(sum(sum(ori_aic))) / nnz(ori_aic));
        tempsmo = to_fastsmooth(aic, ceil(maxrem/10),1,1);


        % fit, find global minimum (beware curve fitting requirements...)
        try

            [cc, ~] = fit((1:maxrem)', tempsmo', 'linearinterp');
            [temp, aic_min] = min(cc(1:maxrem));

        catch

            cc = tempsmo;
            [temp, aic_min] = min(cc);

        end


        % potentially enforce aic? Only if necessary...
        f_aic_curr = NaN;
        aic_ind = NaN;
        if outl_a > aic_min

            try

                % relate currently suggested to actually selected aic ...
                f_aic_curr = cc(outl_a)/cc(aic_min);

            catch

                % ... which sometimes fails if outl_a is too far out
                f_aic_curr = cc(maxrem)/cc(aic_min);

            end
            if f_aic_curr > f_aic

                % deviation too strong, find last acceptable value
                aic_ind = find(cc(aic_min:maxrem) <= f_aic*cc(aic_min),1 ,'last');
                aic_ind = aic_min + aic_ind -1;
                disp(['            ... removing ' num2str(outl_a) ' outliers worsens AIC more than allowed (' sprintf('%02.1f%', f_aic_curr) ' > ' sprintf('%02.1f%', f_aic) '); correcting to ' num2str(aic_ind) ', please wait...']);
                outl_a = aic_ind;
                %else disp(['            ... accepting AIC deviation factor
                %(' sprintf('%02.1f%', f_aic_curr) ' < ' sprintf('%02.1f%',
                %f_aic) ')...']);

            end

        end


        % other way round...
        if outl_a < aic_min

            % also don't allow aic to dominate exclusion
            disp(['            ... removing ' num2str(aic_min) ' outliers as suggested by AIC exceeds combined number (' num2str(outl_a) '), correcting, please wait...']);
            aic_min = outl_a;

        end


        % now compare and get the maximum, not forgetting minrem (and
        % maxrem)
        bestcom = max([aic_min, outl_a, minrem]);
        if bestcom > maxrem,  bestcom = maxrem;  end


    elseif approach == 2


        % get number of outliers without AIC, not forgetting minrem (and
        % maxrem)
        bestcom = max([outl_a, minrem]);
        if bestcom > maxrem,  bestcom = maxrem;  end
        aic_min    = NaN;
        aic        = NaN;
        aic_ind    = NaN;
        cc         = NaN;
        f_aic_curr = NaN;
        ori_aic    = NaN;

    end


    %% potentially warn user
    if (bestcom/nscans*100) > perc_warn

        disp(['            ... warning: suggested number of datapoints (' num2str(bestcom) ') exceeds suggested maximum (' num2str(perc_warn) '%), please check data quality carefully!']);

    end


    % clean main directory
    files = {'^mask\..{3}$','^ResMS\..{3}$','^RPV\..{3}$','^beta_.{4}\..{3}$','^con_.{4}\..{3}$','^ResI_.{4}\..{3}$','^ess_.{4}\..{3}$', '^spm\w{1}_.{4}\..{3}$', '^SPM.mat'};
    for i=1:length(files)

        junk = spm_select('List',pwd,files{i});
        for k=1:size(junk,1),  spm_unlink(deblank(junk(k,:)));  end

    end


    % save everything of interest
    if approach == 1

        disp(['      ... done, censoring ' num2str(bestcom) ' datapoints, based on AIC (' num2str(aic_min) '), R2 (' num2str(outl_r2) '), DVARS (' num2str(outl_d) ') and STS (' num2str(outl_s) '), namely']);

    else

        disp(['      ... done, censoring ' num2str(bestcom) ' datapoints, based on DVARS (' num2str(outl_d) ') and STS (' num2str(outl_s) '), namely']);

    end
    disp(['          #' regexprep(num2str(sort(sortout(1:bestcom))), '\s*', ', ') '...']);
    goon(1,j) = 1;
    bestcoms(1,j) = bestcom;
    if verLessThan('matlab','7.3')

        save([dirname '.mat'], 'aic', 'aic_ind', 'aic_min', 'approach', 'bestcom', 'cc', 'corrfact', 'dvars_p', 'dfb', 'DeltapDvar', 'f_aic', 'f_aic_curr', 'f_tukey', 'hIQR', 'ind_d', 'ind_r2', 'ind_s', 'matlabbatch_o', 'nscans', 'outl_a', 'outl_d', 'outl_r2', 'outl_s', 'ospm', 'ori_r2', 'ori_aic', 'ospm_p', 'o_maxrem', 'o_minrem', 'minrem', 'maxrem', 'rem_r2', 'rem_d', 'rem_s', 'thrs', 'thr_s', 'sortout', 'store_t', 'useart', 'ubound_r2', 'ubound_s');

    else

        save([dirname '.mat'], 'aic', 'aic_ind', 'aic_min', 'approach', 'bestcom', 'cc', 'corrfact', 'dvars_p', 'dfb', 'DeltapDvar', 'f_aic', 'f_aic_curr', 'f_tukey', 'hIQR', 'ind_d', 'ind_r2', 'ind_s', 'matlabbatch_o', 'nscans', 'outl_a', 'outl_d', 'outl_r2', 'outl_s', 'ospm', 'ori_r2', 'ori_aic', 'ospm_p', 'o_maxrem', 'o_minrem', 'minrem', 'maxrem', 'rem_r2', 'rem_d', 'rem_s', 'thrs', 'thr_s', 'sortout', 'store_t', 'useart', 'ubound_r2', 'ubound_s', '-v7.3');

    end

end


% reconsider zero-padding?
disp(['   ... done with first stage analyses of ' num2str(nmats) ' models, continuing ...']);
if zeropad < 0

    zeropad = max(bestcoms);
    disp(['      ... setting adaptive zeropadding to ' num2str(zeropad) ' regressors ...']);

end


% inform about individual contributions for possible quality control
if nmats > 1 && sum(goon(:)) > 0


    % additionally write to text file, but watch out for permissions
    try

        % inform on screen, sort by number of outliers
        disp(['      ... quality control report (see also ' myout ' in ' ori '):']);

        new = exist([ori filesep myout],'file');
        fid = fopen([ori filesep myout],'At+');
        if new == 0,  fprintf(fid, ['\t' 'Outlying datapoints' '\t' 'Input file' '\n']);  end


        for i = 1:nmats


            % maybe not necessary
            if goon(1,i) == 0,  continue;  end


            % inform on screen and write to file
            try

                fprintf(fid, ['\t' num2str(bestcoms(i)) '\t' strrep(deblank(ospms(i,:)), '\', '/') '\n']);
                disp(['          ... identified ' num2str(bestcoms(i)) ' outlying datapoints in ' deblank(ospms(i,:)) '...']);

            catch

                [temp, ~,~,~] = spm_fileparts(ospms{i,1}(1,:));
                fprintf(fid, ['\t' num2str(bestcoms(i)) '\t' strrep(deblank(temp), '\', '/') '\n']);
                disp(['          ... identified ' num2str(bestcoms(i)) ' outlying datapoints in the session from ' odat '...']);

            end

        end
        fclose(fid);

    catch

        disp(['      ... sorry, writing quality control report  to ' myout ' in ' ori ' did not work, proceeding...']);

    end

end


%% =========================================================================================================
%                                          Begin stage 2...
% ==========================================================================================================


% storage
out_cens = cell(nmats,1);
out_int  = cell(nmats,1);
out_both = cell(nmats,1);


% get going (again...)
for j = 1:nmats


    % do we need to do this at all?
    if docalc == 0,  disp('      ... skipping final analysis');  continue;  end


    % find current file and previous results
    if goon(1,j) == 0,  continue;  end
    if ~iscell(ospms)

        ospm = deblank(ospms(j,:));
        [ospm_p, ~, ~, ~] = spm_fileparts(ospm);

    else

        [ospm_p, ~, ~, ~] = spm_fileparts(ospms{j,1}(1,:));

    end
    if isempty(ospm_p),  ospm_p = pwd;  end
    cd([ospm_p filesep dirname]);
    old = pwd;
    load([dirname '.mat']);
    addto = [pwd filesep dirname '.mat'];


    % create storage and evaluate effects on DeltapDVAR and STS
    pchange      = zeros(3,4);
    curr         = [zeros(size(DeltapDvar,1),1) DeltapDvar];
    temp         = sum(sum(abs(curr(:)))) - sum(sum(abs(curr(:,sortout(1:bestcom)))));
    pchange(:,1) = temp ./ sum(abs(curr(:))) .* 100;
    temp         = sum(ind_s) - sum(ind_s(sortout(1:bestcom)));
    pchange(:,2) = temp ./ sum(ind_s) .* 100;


    % prettify so that timepoints are removed in temporal order
    sortout(1:bestcom) = sort(sortout(1:bestcom));


    % final analyses 1 (censoring)
    if docalc == 1 || docalc == 3


        % inform about approach
        disp(['      ... running final analysis on model ' num2str(j) '/' num2str(nmats) ' (censoring ' num2str(bestcom) ' images)...']);
        matlabbatch = matlabbatch_o;


        % change directory
        mkdir([pwd filesep '01_censored'])
        matlabbatch{1}.spm.stats.fmri_spec.dir = {[pwd filesep '01_censored']}
        cd([pwd filesep '01_censored']);
        out_cens{j,1} = [pwd filesep 'SPM.mat'];


        % index to best combination (and potentially zerofill)
        ind = zeros(nscans,maxrem);
        for ii = 1:maxrem

            if ii <= bestcom

                ind(sortout(ii),ii) = 1;
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1, size(matlabbatch{1}.spm.stats.fmri_spec.sess.regress,2)+1).name = ['RI' num2str(ii)];
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1, size(matlabbatch{1}.spm.stats.fmri_spec.sess.regress,2)).val = ind(:,ii);

            else

                if ii <= zeropad

                    ind(end,ii) = 1;
                    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1, size(matlabbatch{1}.spm.stats.fmri_spec.sess.regress,2)+1).name = ['ZP' num2str(ii-bestcom)];
                    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1, size(matlabbatch{1}.spm.stats.fmri_spec.sess.regress,2)).val = ind(:,ii);

                end

            end

        end


        % set up model, reload
        myback_cens = matlabbatch;
        evalc('mw_spm_run_fmri_spec(matlabbatch{1}.spm.stats.fmri_spec);');
        load([pwd filesep 'SPM.mat']);


        % designate nuisance variables as such and estimate
        temp = 1:(SPM.xX.iB-1);
        temp(ismember(temp, ospm.SPM.xX.iC)) = [];
        SPM.xX.iC = ospm.SPM.xX.iC;    % original columns of interest
        SPM.xX.iG = temp;              % our regressors = columns of no interest
        SPM = mw_spm_spm(SPM, [], []);


        % assess original contrast(s)
        for ii = 1:size(ospm.SPM.xCon,2)


            % get original contrast, add zeros for our regressor
            contrast = ii;
            c = [ospm.SPM.xCon(contrast).c; zeros(size(SPM.xX.xKXs.X,2)-numel(ospm.SPM.xCon(contrast).c),1)];
            if zeropad == 0

                cname = [ospm.SPM.xCon(contrast).name ' (censoring ' num2str(bestcom) ', no zeropadding)'];

            else

                cname = [ospm.SPM.xCon(contrast).name ' (censoring ' num2str(bestcom) ', zeropadded to ' num2str(zeropad) ')'];

            end


            % ...and evaluate
            if ii == 1

                evalc('SPM.xCon = spm_FcUtil(''Set'',cname,ospm.SPM.xCon(contrast).STAT,''c'',c,SPM.xX.xKXs);');

            else

                evalc(['SPM.xCon(' num2str(ii) ') = spm_FcUtil(''Set'',cname,ospm.SPM.xCon(contrast).STAT,''c'',c,SPM.xX.xKXs);']);

            end

        end
        evalc('spm_contrasts(SPM);');


        % re-run f-statistics on our optimized model (re-using old mask to
        % ensure comparability)
        mask = spm_read_vols(spm_vol([ospm_p filesep ospm.SPM.VM.fname]));
        curr = find(mask>0);
        Y = zeros(SPM.nscan, numel(curr));
        for i = 1:SPM.nscan

            Y(i,:) = spm_data_read(spm_vol(SPM.xY.P(i,:)),curr);

        end


        % get data from new model and remove dummy regressors to enable
        % comparability
        X = SPM.xX.X;
        if zeropad > bestcom

            temp = size(X,2) - (zeropad-bestcom);
            X(:,temp:end-1) = [];

        end


        % estimate model and get residuals
        B = inv(X'*X)*X'*Y;
        Yhat = X*B;
        Res  = Y - Yhat;
        try

            Y_dm = Y-mean(Y);

        catch

            Y_dm = Y-repmat(mean(Y), size(Y, 1), 1);

        end


        % calculate new R2
        SSerror  = (sum(abs(Res).^2).^(1/2))  .^2;
        SStotal  = (sum(abs(Y_dm).^2).^(1/2)) .^2;
        new_r2   = 1 - SSerror./SStotal;


        % calculate percent change F-values due to new model
        myend = min(numel(new_r2), numel(ori_r2));
        temp  = (new_r2(1:myend) / ori_r2(1:myend) .* 100);
        pchange(1,3) = mean(temp(~isnan(temp)));


        % calculate new corrected AIC
        try

            vol = spm_read_vols(spm_vol([pwd filesep 'ResMS.img']));

        catch

            vol = spm_read_vols(spm_vol([pwd filesep 'ResMS.nii']));

        end
        k    = repmat(size(SPM.xX.X,2)-1, ospm.SPM.VM.dim);
        n    = repmat(SPM.nscan, ospm.SPM.VM.dim);
        df   = repmat(SPM.xX.erdf, ospm.SPM.VM.dim);
        temp = 2*k + n .* (log((vol .* df) ./ n)) + (2.*k.*(k+1) ./ (n-k-1));
        temp(isnan(temp)) = [];
        myend = min(numel(temp), numel(ori_aic));
        temp = (temp(1:myend) / ori_aic(1:myend) .* 100);
        pchange(1,4) = mean(temp(~isnan(temp)));


    end
    cd(old);


    % final analyses 2 (interpolation)
    if docalc == 2 || docalc == 3


        % inform about approach
        if useart == 1

            if ~iscell(ospms)

                disp(['      ... running final analysis on model ' num2str(j) '/' num2str(nmats) ' (interpolating (a) ' num2str(bestcom) ' images)...']);

            else

                disp(['      ... running final analysis on dataset ' num2str(j) '/' num2str(nmats) ' (interpolating (a) ' num2str(bestcom) ' images)...']);

            end

        else

            if ~iscell(ospms)

                disp(['      ... running final analysis on model ' num2str(j) '/' num2str(nmats) ' (interpolating (i) ' num2str(bestcom) ' images)...']);

            else

                disp(['      ... running final analysis on dataset ' num2str(j) '/' num2str(nmats) ' (interpolating (i) ' num2str(bestcom) ' images)...']);

            end

        end
        matlabbatch = matlabbatch_o;


        % change directory
        mkdir([pwd filesep '02_interp']);
        matlabbatch{1}.spm.stats.fmri_spec.dir = {[pwd filesep '02_interp']};
        cd([pwd filesep '02_interp']);
        out_int{j,1} = [pwd filesep 'SPM.mat'];
        mkdir([pwd filesep 'data']);


        % index to images, remove selected frames
        remove = sortout(1:bestcom);
        if useart == 0


            % use inpaints_nan for interpolation (use all voxels!)
            try

                curr = spm_read_vols(spm_vol(char(matlabbatch{1}.spm.stats.fmri_spec.sess.scans)));

            catch

                curr = spm_read_vols(spm_vol(char(ospms{j,1})));

            end
            sc   = size(curr);
            curr = reshape(curr, [sc(1)*sc(2)*sc(3) sc(4)]);
            ind  = sum(curr,2);
            curr(:,remove) = NaN;
            new = zeros(size(curr));
            try

                parfor i = 1:size(curr,1)

                    if (ind(i) == 0 || ~any(diff(curr(i,:)))),  continue;  end
                    new(i,:) = inpaint_nans(curr(i,:), 3);

                end

            catch

                for i = 1:size(curr,1)

                    if (ind(i) == 0 || ~any(diff(curr(i,:)))),  continue;  end
                    new(i,:) = inpaint_nans(curr(i,:), 3);

                end

            end
            new = reshape(new, sc);


            % write out and adapt design
            try

                vols = char(matlabbatch{1}.spm.stats.fmri_spec.sess.scans);

            catch

                vols = char(ospms{j,1});

            end
            newvols = '';
            for i = 1:size(vols,1)

                V = spm_vol(vols(i,:));
                [~, nm, e, v] = spm_fileparts(V.fname);
                if ismember(i, remove)

                    V.fname = [pwd filesep 'data' filesep nm '_interp' e v];

                else

                    V.fname = [pwd filesep 'data' filesep nm e v];

                end
                spm_write_vol(V, new(:,:,:,i));
                newvols = char(newvols, V.fname);

            end
            newvols = newvols(2:end,:);

        else

            % use our modified art_repairvol for interpolation
            try

                newvols = mw_art_repairvol(char(matlabbatch{1}.spm.stats.fmri_spec.sess.scans), remove, [pwd filesep 'data']);

            catch

                newvols = mw_art_repairvol(char(ospms{j,1}), remove, [pwd filesep 'data']);

            end

        end


        % set up model, reload and estimate (only if model was passed)
        if ~iscell(ospms)

            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(newvols);
            myback_int = matlabbatch;
            evalc('mw_spm_run_fmri_spec(matlabbatch{1}.spm.stats.fmri_spec);');
            load([pwd filesep 'SPM.mat']);
            SPM = mw_spm_spm(SPM, [], []);


            % assess original contrast(s)
            for ii = 1:size(ospm.SPM.xCon,2)

                % get original contrast, add zeros for our regressor
                contrast = ii;
                c = [ospm.SPM.xCon(contrast).c; zeros(size(SPM.xX.xKXs.X,2)-numel(ospm.SPM.xCon(contrast).c),1)];
                cname = [ospm.SPM.xCon(contrast).name ' (interpolating ' num2str(bestcom) ')'];


                % ...and evaluate
                if ii == 1

                    evalc('SPM.xCon = spm_FcUtil(''Set'',cname,ospm.SPM.xCon(contrast).STAT,''c'',c,SPM.xX.xKXs);');

                else

                    evalc(['SPM.xCon(' num2str(ii) ') = spm_FcUtil(''Set'',cname,ospm.SPM.xCon(contrast).STAT,''c'',c,SPM.xX.xKXs);']);

                end

            end
            evalc('spm_contrasts(SPM);');


            % re-run f-statistics on our optimized model (re-using old mask
            % to ensure comparability)
            mask = spm_read_vols(spm_vol([ospm_p filesep ospm.SPM.VM.fname]));
            curr = find(mask>0);
            Y = zeros(SPM.nscan, numel(curr));
            for i = 1:SPM.nscan

                Y(i,:) = spm_data_read(spm_vol(SPM.xY.P(i,:)),curr);

            end


            % specify, estimate model and get residuals
            X = SPM.xX.X;
            B = inv(X'*X)*X'*Y;
            Yhat = X*B;
            Res  = Y - Yhat;
            try

                Y_dm = Y-mean(Y);

            catch

                Y_dm = Y-repmat(mean(Y), size(Y, 1), 1);

            end


            % calculate new R2
            SSerror  = (sum(abs(Res).^2).^(1/2))  .^2;
            SStotal  = (sum(abs(Y_dm).^2).^(1/2)) .^2;
            new_r2   = 1 - SSerror./SStotal;


            % calculate percent change R2-values due to new model
            myend = min(numel(new_r2), numel(ori_r2));
            temp = new_r2(1:myend) / ori_r2(1:myend) .* 100;
            pchange(2,3) = mean(temp(~isnan(temp)));


            % calculate new corrected AIC
            try

                vol = spm_read_vols(spm_vol([pwd filesep 'ResMS.img']));

            catch

                vol = spm_read_vols(spm_vol([pwd filesep 'ResMS.nii']));

            end
            k   = repmat(size(SPM.xX.X,2)-1, ospm.SPM.VM.dim);
            n   = repmat(SPM.nscan, ospm.SPM.VM.dim);
            df  = repmat(SPM.xX.erdf, ospm.SPM.VM.dim);
            temp = 2*k + n .* (log((vol .* df) ./ n)) + (2.*k.*(k+1) ./ (n-k-1));
            temp(isnan(temp)) = [];
            myend = min(numel(temp), numel(ori_aic));
            temp = temp(1:myend) / ori_aic(1:myend) .* 100;
            pchange(2,4) = mean(temp(~isnan(temp)));

        else

            % special case if only images were passed
            if iscell(ospms)

                out_int{j,1} = spm_select('FPList', [pwd filesep 'data'], '.*(img|nii)');

            end

        end
        cd(old);

    end


    % final analyses 3 (interpolation & censoring)
    if docalc == 3


        % inform about approach
        if useart == 1

            disp(['      ... running final analysis on model ' num2str(j) '/' num2str(nmats) ' (interpolating (a) & censoring ' num2str(bestcom) ' images)...']);

        else

            disp(['      ... running final analysis on model ' num2str(j) '/' num2str(nmats) ' (interpolating (i) & censoring ' num2str(bestcom) ' images)...']);

        end


        % change directory
        mkdir([pwd filesep '03_both']);
        cd([pwd filesep '03_both']);
        out_both{j,1} = [pwd filesep 'SPM.mat'];


        % reuse aspects from models above
        matlabbatch = myback_cens;
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(newvols);
        matlabbatch{1}.spm.stats.fmri_spec.dir = {pwd};


        % set up model, reload
        evalc('mw_spm_run_fmri_spec(matlabbatch{1}.spm.stats.fmri_spec);');
        load([pwd filesep 'SPM.mat']);


        % designate nuisance variables as such and estimate
        temp = 1:(SPM.xX.iB-1);
        temp(ismember(temp, ospm.SPM.xX.iC)) = [];
        SPM.xX.iC = ospm.SPM.xX.iC;    % original columns of interest
        SPM.xX.iG = temp;              % our regressors = columns of no interest
        SPM = mw_spm_spm(SPM, [], []);


        % assess original contrast(s)
        for ii = 1:size(ospm.SPM.xCon,2)


            % get original contrast, add zeros for our regressor
            contrast = ii;
            c = [ospm.SPM.xCon(contrast).c; zeros(size(SPM.xX.xKXs.X,2)-numel(ospm.SPM.xCon(contrast).c),1)];
            if zeropad == 0

                cname = [ospm.SPM.xCon(contrast).name ' (interpolating & censoring ' num2str(bestcom) ', no zeropadding)'];

            else

                cname = [ospm.SPM.xCon(contrast).name ' (interpolating & censoring ' num2str(bestcom) ', zeropadded to ' num2str(zeropad) ')'];

            end


            % ...and evaluate
            if ii == 1

                evalc('SPM.xCon = spm_FcUtil(''Set'',cname,ospm.SPM.xCon(contrast).STAT,''c'',c,SPM.xX.xKXs);');

            else

                evalc(['SPM.xCon(' num2str(ii) ') = spm_FcUtil(''Set'',cname,ospm.SPM.xCon(contrast).STAT,''c'',c,SPM.xX.xKXs);']);

            end

        end
        evalc('spm_contrasts(SPM);');


        % re-run f-statistics on our optimized model (re-using old mask to
        % ensure comparability)
        mask = spm_read_vols(spm_vol([ospm_p filesep ospm.SPM.VM.fname]));
        curr = find(mask>0);
        Y = zeros(SPM.nscan, numel(curr));
        for i = 1:SPM.nscan

            Y(i,:) = spm_data_read(spm_vol(SPM.xY.P(i,:)),curr);

        end


        % specify, estimate model and get residuals
        X = SPM.xX.X;
        B = inv(X'*X)*X'*Y;
        Yhat = X*B;
        Res  = Y - Yhat;
        try

            Y_dm = Y-mean(Y);

        catch

            Y_dm = Y-repmat(mean(Y), size(Y, 1), 1);

        end


        % calculate new R2
        SSerror  = (sum(abs(Res).^2).^(1/2))  .^2;
        SStotal  = (sum(abs(Y_dm).^2).^(1/2)) .^2;
        new_r2    = 1 - SSerror./SStotal;


        % calculate percent change R2-values due to new model
        myend = min(numel(new_r2), numel(ori_r2));
        temp = new_r2(1:myend) / ori_r2(1:myend) .* 100;
        pchange(3,3) = mean(temp(~isnan(temp)));


        % calculate new corrected AIC
        try

            vol = spm_read_vols(spm_vol([pwd filesep 'ResMS.img']));

        catch

            vol = spm_read_vols(spm_vol([pwd filesep 'ResMS.nii']));

        end
        k   = repmat(size(SPM.xX.X,2)-1, ospm.SPM.VM.dim);
        n   = repmat(SPM.nscan, ospm.SPM.VM.dim);
        df  = repmat(SPM.xX.erdf, ospm.SPM.VM.dim);
        temp = 2*k + n .* (log((vol .* df) ./ n)) + (2.*k.*(k+1) ./ (n-k-1));
        temp(isnan(temp)) = [];
        myend = min(numel(temp), numel(ori_aic));
        temp = temp(1:myend) / ori_aic(1:myend) .* 100;
        pchange(3,4) = mean(temp(~isnan(temp)));


    else

        % define variable needed below
        new_r2 = NaN;

    end
    cd(old);


    % create a figure (for some strange reason, this sometimes fails,
    % so...)
    try

        pos = spm('WinSize','Graphics');
        if isempty(pos) || any(pos(3:4)<200)

            pos = [500 20 550 800];

        end
        if approach == 1

            figname = [ dirname ' results plot (== task-based == fMRI approach)'];

        else

            figname = [ dirname ' results plot (== resting-state == fMRI approach)'];

        end
        fig = figure('Name', figname, 'position', pos, 'Color', 'w', 'Toolbar', 'none', 'visible', 'on');

        figure(fig); subplot(4,2,[1 2], 'replace');
        dvarsd = [DeltapDvar'; ones(1, size(DeltapDvar,1)) .* dmax];
        dvarsd(dvarsd > dmax) = dmax;
        surf(dvarsd, 'linestyle', 'none');
        set(gca, 'view', [90 90], 'XDir', 'reverse');
        ylim([1 size(DeltapDvar,2)]);
        xlim([1 size(DeltapDvar,1)]);
        title('\Delta%D-var (over slices and time)');
        xlabel('Slices', 'FontAngle', 'italic');
        ylabel('Timepoints', 'FontAngle', 'italic');

        figure(fig); subplot(4,2,3, 'replace');
        if dfb == 0,  dvarsd = [0 sum((dvars_p*corrfact)< thr_d(1))];
        else,         dvarsd = [0 sum(DeltapDvar > thr_d(3))];  end
        bar(dvarsd, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k', 'LineWidth', 0.25);
        xlim([1 numel(ind_d)]);
        ylim([0 size(DeltapDvar,1)]);
        if dfb == 0,  title('\Delta%D-var: significance');  else,  title('Excessive \Delta%D-var');  end
        xlabel('Timepoints', 'FontAngle', 'italic');
        ylabel('# slices', 'FontAngle', 'italic');
        plotme = NaN(1, numel(ind_d));
        plotme(rem_d(1:outl_d)) = dvarsd(rem_d(1:outl_d));
        hold on; plot(plotme, 'sr')
        line([0 numel(ind_d)], [corrfact corrfact], 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
        if dfb == 0,  msg = char(['p < ' sprintf('%0.2f', thr_d(1)) ' in \geq ' num2str(corrfact) ' slices (n = ' num2str(outl_d) ')']);
        else,         msg = char(['> ' sprintf('%1.0f', thr_d(3)) '% in \geq ' num2str(corrfact) ' slices (n = ' num2str(outl_d) ')']);  end
        drawnow;
        text(0.1*max(get(gca, 'xlim')), 0.9*diff(get(gca, 'YLim'))+min(get(gca, 'ylim')), msg, 'FontSize',8);

        figure(fig); subplot(4,2,4, 'replace');
        bar(ind_s, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k', 'LineWidth', 0.25);
        xlim([0 numel(ind_s)]);
        title(['STS @ ' num2str(davg) ' mm'], 'Interpreter', 'none');
        xlabel('Timepoints', 'FontAngle', 'italic');
        ylabel('mm', 'FontAngle', 'italic');
        hold on;
        plotme = NaN(1, numel(ind_s));
        plotme(rem_s(1:outl_s)) = ind_s(rem_s(1:outl_s));
        plot(plotme, 'sr')
        if sum(ind_s(:)) == 0

            msg = char('STS not available');

        else

            line([0 numel(ind_s)], [ubound_s ubound_s], 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
            msg = char(['STS > ' sprintf('%0.1f', ubound_s) ' mm (n = ' num2str(outl_s) ')']);

        end
        drawnow;
        text(0.1*max(get(gca, 'xlim')), 0.9*diff(get(gca, 'YLim'))+min(get(gca, 'ylim')), msg, 'FontSize',8);

        figure(fig); subplot(4,2,5, 'replace');
        bar(ind_r2, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k', 'LineWidth', 0.25);
        xlim([0 numel(ind_r2)]);
        title('Explained variance (R2)', 'Interpreter', 'none');
        ylabel('Ratio', 'FontAngle', 'italic');
        hold on;
        if approach == 1

            ylim([min(ind_r2)*0.9 max(ind_r2)*1.1]);
            xlabel('Timepoints', 'FontAngle', 'italic');
            plotme = NaN(1, numel(ind_r2));
            plotme(rem_r2(1:outl_r2)) = ind_r2(rem_r2(1:outl_r2));
            plot(plotme, 'sr')
            line([0 numel(ind_r2)], [ubound_r2 ubound_r2], 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
            msg = char(['Tukey_R_2 (n = ' num2str(outl_r2) ')' ]);

        else

            msg = char('R2 not available');

        end
        drawnow;
        text(0.1*max(get(gca, 'xlim')), 0.9*diff(get(gca, 'YLim'))+min(get(gca, 'ylim')), msg, 'FontSize',8);

        figure(fig); subplot(4,2,6, 'replace');
        plot(aic, 'Color', [0.5 0.5 0.5]);
        xlim([0 maxrem]);
        ylim([min(get(gca, 'ylim')) max(get(gca, 'ylim'))]);
        xlabel('Removed images', 'FontAngle', 'italic');
        ylabel('Ratio', 'FontAngle', 'italic');
        title(['AIC_C (-1:' num2str(maxrem) ' images)']);
        if approach == 1

            hold on; plot(cc(1:maxrem), 'r.');
            if isnan(aic_ind)

                msg = char(['AIC_C (suggested = ' num2str(aic_min) ')' ]);

            else

                msg = char(['AIC_C (suggested = ' num2str(aic_min) ', acceptable < ' num2str(aic_ind) ' )' ]);
                hold on; line([0 max(get(gca, 'xlim'))], [cc(aic_min)*f_aic cc(aic_min)*f_aic], 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);

            end
            hold on; line([aic_min aic_min], [min(get(gca, 'ylim')) max(get(gca, 'ylim'))], 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
            drawnow;

        else

            msg = char('AIC not available');

        end
        text(0.1*max(get(gca, 'xlim')), 0.9*diff(get(gca, 'YLim'))+min(get(gca, 'ylim')), msg, 'FontSize',8);

        figure(fig); subplot(4,2,[7 8], 'replace');
        plotme = NaN(1, nscans);
        plotme(sortout(1:bestcom)) = 1;
        bar(plotme, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k', 'LineWidth', 0.1);
        ylim([0 1]);

        try, ax = gca; for ii = 1:numel(ax.YTickLabel), ax.YTickLabel{ii} = ['\color{white}' ax.YTickLabel{ii}]; end; end
        xlim([0 nscans]);
        xlabel('Timepoints', 'FontAngle', 'italic');
        ylabel('Removed images', 'FontAngle', 'italic');
        title('Optimized censoring');
        if (aic_min < minrem) && (outl_a < minrem)

            msg = char(['Final: remove ' num2str(bestcom) ' images (owing to MinRem) = ' sprintf('%0.1f', bestcom/nscans*100) '% of datapoints']);

        elseif max([aic_min, outl_a, minrem]) > bestcom

            msg = char(['Final: remove ' num2str(bestcom) ' images (owing to MaxRem) = ' sprintf('%0.1f', bestcom/nscans*100) '% of datapoints']);

        else

            msg = char(['Final: remove ' num2str(bestcom) ' images = ' sprintf('%0.1f', bestcom/nscans*100) '% of datapoints']);

        end
        drawnow;
        text(0.5*nscans, 0.5, msg, 'FontSize',8, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', [0.5 0.5 0.5], 'HorizontalAlignment', 'center');

        axes('Position',[0.01,0.01,0.1,0.1],'Visible','off');
        if length(ospm_p) > 60,  ospm_p_d = ['...' ospm_p(end-55:end)];  else,  ospm_p_d = ospm_p;  end
        if zeropad == 0

            text(1,0.2,{['Results from ' ospm_p_d ', no zeropadding, removing'] ['#' regexprep(num2str(sort(sortout(1:bestcom))), '\s*', ', ') '.']},'FontSize',7,'interpreter','none');

        else

            text(1,0.2,{['Results from ' ospm_p_d ', zeropadded to ' num2str(zeropad) ' regressors, removing'] ['#' regexprep(num2str(sort(sortout(1:bestcom))), '\s*', ', ') '.']},'FontSize',7,'interpreter','none');

        end
        drawnow;
        set(fig,'PaperPositionMode','auto');
        try

            saveas(fig, [pwd filesep dirname '.png'], 'png');

        catch

            print(fig, '-dpng', '-noui', '-r300', [pwd filesep dirname '.png']);

        end

    catch

        fig = figure('Name', ['SOMETHING WENT WRONG WITH ' dirname ' FIGURE CREATION!'], 'position', [500 20 550 800], 'Toolbar', 'none', 'visible', 'on');

    end


    % add new results (including figure) to results file
    disp(['   ... saving results in ' ospm_p filesep dirname '...']);
    set(fig, 'visible', 'on');


    % potentially consider Octave users
    try

        myfigure = handle2struct(fig);

    catch

        myfigure = hdl2struct(fig);

    end
    eval(['save(''' addto ''', ''myfigure'', ''zeropad'', ''new_r2'', ''pchange'', ''-append'');']);
    if nmats > 1 || nargin > 1 || nargout > 0,  close(fig);  end


    %% =========================================================================================================
    %                                          Further outputs
    % ==========================================================================================================


    % loop over t-map store and generate overlap maps at the threshold
    % specified above
    if any(~isnan(pval)) == 1

        for ii = 1:size(ospm.SPM.xCon,2);

            disp(['      ... done, preparing further outputs (contrast ' num2str(ii) '/' num2str(size(ospm.SPM.xCon,2)) ')...']);
            curr_c = squeeze(store_t(:,:,:,:,ii));
            for iii = 1:numel(pval)

                if ~isnan(pval(iii))

                    % create storage
                    overlap = zeros(ospm.SPM.VM.dim);
                    for iiii = 1:maxrem

                        overlap = overlap + (curr_c(:,:,:,iiii) > thrs(iiii,iii));

                    end
                    overlap = overlap .* 100 ./ maxrem;


                    % get a handle
                    try

                        V = spm_vol([ospm_p filesep 'spm' ospm.SPM.xCon(ii).STAT '_' sprintf('%04.f', ii) '.img']);

                    catch

                        V = spm_vol([ospm_p filesep 'spm' ospm.SPM.xCon(ii).STAT '_' sprintf('%04.f', ii) '.nii']);

                    end


                    % store results as percent overlap map
                    [~, ~, e, ~] = spm_fileparts(V.fname);
                    V.pinfo = [1 0 0]';
                    if iii == 1,  mc_d = 'FWE';  elseif iii == 2,  mc_d = 'FDR';  else,  mc_d = 'NONE';  end
                    V.fname   = [pwd filesep dirname '_pom_c' num2str(ii) '_' mc_d e];
                    V.descrip = ['POM, n=' num2str(maxrem) ', p=' sprintf('%0.3f', pval) ' (' mc_d ')'];
                    spm_write_vol(V, overlap);

                end

            end

        end

    end

end


%% =========================================================================================================
%                                          Housekeeping
% ==========================================================================================================


% undo changes to spm defaults
spm_get_defaults('cmdline', ocmd);


% inform user
ende = toc;
tt(1) = floor(ende/3600);
tt(2) = floor((ende - tt(1)*3600) / 60);
tt(3) = ceil(ende - (tt(1)*3600) - (tt(2)*60));
disp(' ');
if ende >= 3600

    disp(['... Done - this run took ' num2str(tt(1)) ' hours, ' num2str(tt(2)) ' minutes and ' num2str(tt(3)) ' seconds;']);

elseif ende >= 60

    disp(['... Done - this run took ' num2str(tt(2)) ' minutes and ' num2str(tt(3)) ' seconds;']);

else

    disp(['... Done - this run took ' num2str(tt(3)) ' seconds;']);

end
disp('... Thank you for using this script and have a nice day :)');
disp(' ');
cd(ori);


% outputs?
if nargout == 0

    clear out_cens out_int out_both;

else

    if ~(docalc == 1 || docalc == 3),  out_cens = [];  end
    if ~(docalc == 2 || docalc == 3),  out_int  = [];  end
    if docalc ~= 3,                    out_both = [];  end

end

return;


%% =========================================================================================================
%                                        Embedded functions
% ==========================================================================================================
function mask = mw_brainmask2(P, silent);

% little function to generate a brainmask from a series of images, usually
% an fMRI timeseries; it uses the correlation approach described by Gerard
% Ridgeway in NeuroImage 2009, see
%  http://dx.doi.org/10.1016/j.neuroimage.2008.08.045
% respective bits of code are used here with kind permission!
%
% Otherwise, implementation is by Marko Wilke, USE AT YOUR OWN RISK!
%

% ======================================================================================
%                      Preludes
% ======================================================================================


% settings
fillin = 100;    % fill holes in the mask up to this size


% check inputs
if nargin == 0

    P = spm_select(Inf,'image','Select timeseries images to analyze',[],pwd,'.*');

elseif ~ischar(P)

    error('... sorry, the first input array has to be a char array!');

end


% silent mode?
if nargin < 2 || isempty(silent)

    silent = 0;

end


% ======================================================================================
%                      computations
% ======================================================================================


% check for 4D, read images
if size(P,1) == 1,  [p, nm, e, ~] = spm_fileparts(P);  P = [p filesep nm e];  end
V = spm_vol(P);


% only proceed if enough images are found
if size(V,1) < 10

    if silent == 0

        error('   Cannot work on only 10 images or less, please reconsider!')

    else

        return;

    end

else

    if silent == 0,  disp('   ... generating brain mask, please wait...');  end

end


% hopefully this works...
try

    evalc('vols = spm_read_vols(V);');

catch

    if silent == 0,  disp('   ... reading images in one go did not work, trying successive approach...');  end
    vols = zeros([V(1).dim size(V,1)]);
    for i = 1:size(V,1)

        try

            vols(:,:,:,i) = spm_read_vols(V(i));

        catch

            if silent == 0,  error('   ... sorry, successive reading of images also did not work, my give up!');  end

        end

    end

end


% clean data, get optimum threshold, apply
vols(~isfinite(vols)) = 0;
vols(isnan(vols)) = 0;
thr = opt_thr_corr(vols);
vols = mean(vols,4) > thr;


% only use largest cluster
[temp2, num] = spm_bwlabel(double(vols > 0.5),18);
store = zeros(num,1);
for i = 1:num,  store(i) = sum(sum(sum(temp2==i)));  end
targ = find(store == max(store));
vols = double(temp2 == targ);


% fill holes by inverting and filling small clusters
[temp2,num] = spm_bwlabel(double(vols < 1),18);
targ = zeros(size(temp2));
for i = 1:num

    temp = sum(sum(sum(temp2==i)));
    if temp < fillin

        targ(temp2 == i) = 1;

    end

end
vols = double((vols + targ) > 0).*255;
vols = vols > 0;


% write out
V = V(1);
[p, nm, e, ~] = spm_fileparts(V.fname);
V.fname = [p filesep nm '_mask' e];
V.dt = [spm_type('uint8') spm_platform('bigend')];
V.private = [];
V.pinfo = [1 0 0]';
V.descrip = 'Brainmask, created with mw_brainmask2';
spm_write_vol(V, vols);


% that's it folks
if silent == 0,  disp(['   ... done, mask was saved to ' p '.']);  end
mask = V.fname;
return


% ======================================================================================
%                       Embedded functions in mw_brainmask2
% ======================================================================================

% #1: opt_thr_corr, used with kind permission from GR
function thr = opt_thr_corr(vols)
costfunc = @(thr) -correlation(vols, vols > thr);
[thr, ncc] = fminbnd(costfunc, min(vols(:)), max(vols(:)));
return;


% #2: correlation, used with kind permission from GR
function c = correlation(x, y)
cs = corrcoef(x, double(y));
c = cs(1, 2);
return;


%% =========================================================================================================

function  q = ah_quantile(x,p,method)
%QUANTILE Empirical quantile (percentile).
%
%         q = ah_quantile(x,p)
%
%	  The empirical quantile is computed by one of three ways determined by
%	  a third input argument (with default 1).
%
%	  1. Interpolation so that F(X_(k)) == (k-0.5)/n. 2. Interpolation so
%	  that F(X_(k)) == k/(n+1). 3. Based on the empirical distribution.
%
%	  If input  x  is a matrix then the quantile is computed for every
%	  column. Input  p  may be vector also. It even accepts  x  being a
%	  vector and  p  a matrix!

%     Copyright (c) Anders Holtsberg, 1995, 1998 File from
%     http://www.maths.lth.se/matstat/stixbox/Contents.html

if nargin<3, method=1; end
if min(size(x)) == 1
    x = x(:);
    q = zeros(size(p));
else
    if min(size(p)) > 1
        error('Not both matrix x and matrix p input')
    end
    q = zeros(length(p),size(x,2));
end
if any(any((p>1|p<0)))
    error('Input p is not probability')
end

x = sort(x);
p = p(:);
n = size(x,1);
if method == 3
    qq1 = x(ceil(max(1,p*n)),:);
    qq2 = x(floor(min(p*n+1,n)),:);
    qq = (qq1+qq2)/2;
else
    x = [x(1,:); x; x(n,:)];
    if method == 2
        % This method is from Hjort's "Computer intensive statistical
        % methods" page 102
        i = p*(n+1)+1;
    else % Metod 1
        i = p*n+1.5;
    end
    iu = ceil(i);
    il = floor(i);
    d = (i-il)*ones(1,size(x,2));
    qq = x(il,:).*(1-d)+x(iu,:).*d;
end

q(:) = qq;
return;

%% ==========================================================================================================

function F = gf_spm_Gcdf(x,h,l,tail)
% Cumulative Distribution Function (CDF) of Gamma distribution FORMAT F =
% spm_Gcdf(x,h,l,tail)
%
% x    - Gamma-variate   (Gamma has range [0,Inf) ) h    - Shape parameter
% (h>0) l    - Scale parameter (l>0) tail - if 'upper', return the upper
% tail probability of the Gamma
%        distribution [Default: 'lower']
% F    - CDF of Gamma-distribution with shape & scale parameters h & l
%__________________________________________________________________________
%
% spm_Gcdf implements the Cumulative Distribution of the
% Gamma-distribution.
%
% Definition:
%--------------------------------------------------------------------------
% The CDF F(x) of the Gamma distribution with shape parameter h and scale l
% is the probability that a realisation of a Gamma random variable X has
% value less than x F(x)=Pr{X<x} for X~G(h,l). The Gamma distribution is
% defined for h>0 & l>0 and for x in [0,Inf) (See Evans et al., Ch18, but
% note that this reference uses the alternative parameterisation of the
% Gamma with scale parameter c=1/l)
%
% Variate relationships: (Evans et al., Ch18 & Ch8)
%--------------------------------------------------------------------------
% For natural (strictly +ve integer) shape h this is an Erlang
% distribution.
%
% The Standard Gamma distribution has a single parameter, the shape h. The
% scale taken as l=1.
%
% The Chi-squared distribution with v degrees of freedom is equivalent to
% the Gamma distribution with scale parameter 1/2 and shape parameter v/2.
%
% Algorithm:
%--------------------------------------------------------------------------
% The CDF of the Gamma distribution with scale parameter l and shape h is
% related to the incomplete Gamma function by
%
%       F(x) = gammainc(l*x,h)
%
% See Abramowitz & Stegun, 6.5.1; Press et al., Sec6.2 for definitions of
% the incomplete Gamma function. The relationship is easily verified by
% substituting for t/c in the integral, where c=1/l.
%
% MATLAB's implementation of the incomplete gamma function is used.
%
% References:
%--------------------------------------------------------------------------
% Evans M, Hastings N, Peacock B (1993)
%       "Statistical Distributions"
%        2nd Ed. Wiley, New York
%
% Abramowitz M, Stegun IA, (1964)
%       "Handbook of Mathematical Functions"
%        US Government Printing Office
%
% Press WH, Teukolsky SA, Vetterling AT, Flannery BP (1992)
%       "Numerical Recipes in C"
%        Cambridge
%__________________________________________________________________________
% Copyright (C) 1992-2019 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes $Id: spm_Gcdf.m 7582 2019-05-01 15:37:16Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<3, error('Insufficient arguments.'), end

ad = [ndims(x);ndims(h);ndims(l)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
    [size(h),ones(1,rd-ad(2))];...
    [size(l),ones(1,rd-ad(3))]];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 && any(any(diff(as(xa,:)),1))
    error('Non-scalar arguments must match in size.');
end

if nargin<4, tail='lower'; end
if ~strcmpi(tail,'lower') && ~strcmpi(tail,'upper')
    error('Tail must be ''lower'' or ''upper''.');
end

%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
F = zeros(rs);

%-Only defined for strictly positive h & l. Return NaN if undefined.
md = ( ones(size(x))  &  h>0  &  l>0 );
if any(~md(:))
    F(~md) = NaN;
    % warning('Returning NaN for out of range arguments.'); disabled by
    % Marko, may happen for empty timeseries
end

%-Non-zero where defined and x>0
Q  = find( md  &  x>0 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qh=Q; else Qh=1; end
if xa(3), Ql=Q; else Ql=1; end

%-Compute
F(Q) = gammainc(l(Ql).*x(Qx),h(Qh),tail);
return;

%% ==========================================================================================================

function SmoothY=to_fastsmooth(Y,w,type,ends)
% Version 3.0, October 2016, taken from
% https://de.mathworks.com/matlabcentral/fileexchange/19998-fast-smoothing-function
%
% Copyright (c) 2012, Thomas C. O'Haver
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.

if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
switch type
    case 1
        SmoothY=sa(Y,w,ends);
    case 2
        SmoothY=sa(sa(Y,w,ends),w,ends);
    case 3
        SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends);
    case 4
        SmoothY=sa(sa(sa(sa(Y,w,ends),w,ends),w,ends),w,ends);
    case 5
        SmoothY=sa(sa(sa(sa(Y,round(1.6*w),ends),round(1.4*w),ends),round(1.2*w),ends),w,ends);
end
function SmoothY=sa(Y,smoothwidth,ends)
w=round(smoothwidth);
SumPoints=sum(Y(1:w));
s=zeros(size(Y));
halfw=round(w/2);
L=length(Y);
for k=1:L-w,
    s(k+halfw-1)=SumPoints;
    SumPoints=SumPoints-Y(k);
    SumPoints=SumPoints+Y(k+w);
end
s(k+halfw)=sum(Y(L-w+1:L));
SmoothY=s./w;
% Taper the ends of the signal if ends=1.
if ends==1,
    startpoint=(smoothwidth + 1)/2;
    SmoothY(1)=(Y(1)+Y(2))./2;
    for k=2:startpoint,
        SmoothY(k)=mean(Y(1:(2*k-1)));
        SmoothY(L-k+1)=mean(Y(L-2*k+2:L));
    end
    SmoothY(L)=(Y(L)+Y(L-1))./2;
end
return

%% ==========================================================================================================

function [SPM, Y, cmask] = mw_spm_spm(SPM, Y, cmask);
%  this is a very stripped down version of Karl & Guillaume's spm_spm.m
%  (6678), removing everything unnecessary and implementing some ideas from
%  spm_spm_WB from the PPI toolbox

try
    cd(SPM.swd);
catch
    SPM.swd = pwd;
end

VY = SPM.xY.VY;
for i = 1:numel(VY)
    if ~spm_existfile(VY(i).fname)
        error('File not found: %s',VY(i).fname);
    end
    if ~spm_mesh_detect(VY)
        VY(i).private.dat.scl_slope = VY(i).pinfo(1);
        VY(i).private.dat.scl_inter = VY(i).pinfo(2);
    end
end

M       = VY(1).mat;
DIM     = VY(1).dim;
YNaNrep = spm_type(VY(1).dt(1),'nanrep');
if spm_mesh_detect(VY)
    file_ext = '.gii';
    g        = VY(1).private;
    metadata = {g.private.metadata(1).name, g.private.metadata(1).value};
else
    file_ext = spm_file_ext;
    metadata = {};
end

xX             = SPM.xX;
[nScan, nBeta] = size(xX.X);

if isfield(SPM,'xM')
    xM         = SPM.xM;
else
    xM         = -Inf(nScan,1);
end
if ~isstruct(xM)
    xM         = struct(...
        'T',  [],...
        'TH', xM,...
        'I',  0,...
        'VM', {[]},...
        'xs', struct('Masking','analysis threshold'));
end

mask           = true(DIM);

if ~isfield(xX,'K')
    xX.K       = 1;
end

if isfield(SPM,'xVi')
    xVi        = SPM.xVi;
else
    xVi        = struct('form', 'i.i.d.',...
        'V',    speye(nScan,nScan));
end

if ~isfield(xVi,'V')
    SPM.xY.VY  = VY;
    SPM.xM     = xM;
    SPM.xX.K   = xX.K;
    [xVi, am]  = mw_spm_est_non_sphericity(SPM, Y, cmask);
    % evalc('[xVi, am]  = spm_est_non_sphericity(SPM);');
    mask       = mask & am;
end

if isfield(xX,'W')
    W          = xX.W;
else
    W          = spm_sqrtm(spm_inv(xVi.V));
    W          = W.*(abs(W) > 1e-6);
    xX.W       = sparse(W);
end

xX.xKXs        = spm_sp('Set',spm_filter(xX.K,W*xX.X));    % KWX
xX.xKXs.X      = full(xX.xKXs.X);
xX.pKX         = spm_sp('x-',xX.xKXs);                     % Projector
erdf           = spm_SpUtil('trRV',xX.xKXs);               % error df

xX.V           = spm_filter(xX.K,spm_filter(xX.K,W*xVi.V*W')'); % KWVW'K'
[trRV, trRVRV] = spm_SpUtil('trRV',xX.xKXs,xX.V);          % trRV (for X)
xX.trRV        = trRV;                                     % <R'*y'*y*R>
xX.trRVRV      = trRVRV;                                   %-Satterthwaite
xX.erdf        = trRV^2/trRVRV;                            % approximation
xX.Bcov        = xX.pKX*xX.V*xX.pKX';                      % Cov(beta)

VM = struct(...
    'fname',   ['mask' file_ext],...
    'dim',     DIM,...
    'dt',      [spm_type('uint8') spm_platform('bigend')],...
    'mat',     M,...
    'pinfo',   [1 0 0]',...
    'descrip', 'spm_spm:resultant analysis mask',...
    metadata{:});
VM = spm_data_hdr_write(VM);

Vbeta(1:nBeta) = deal(struct(...
    'fname',   [],...
    'dim',     DIM,...
    'dt',      [spm_type('float32') spm_platform('bigend')],...
    'mat',     M,...
    'pinfo',   [1 0 0]',...
    'descrip', 'spm_spm:beta',...
    metadata{:}));

for i = 1:nBeta
    Vbeta(i).fname   = [sprintf('beta_%04d',i) file_ext];
    Vbeta(i).descrip = sprintf('spm_spm:beta (%04d) - %s',i,xX.name{i});
end
Vbeta = spm_data_hdr_write(Vbeta);

VResMS = struct(...
    'fname',   ['ResMS' file_ext],...
    'dim',     DIM,...
    'dt',      [spm_type('float64') spm_platform('bigend')],...
    'mat',     M,...
    'pinfo',   [1 0 0]',...
    'descrip', 'spm_spm:Residual sum-of-squares',...
    metadata{:});
VResMS = spm_data_hdr_write(VResMS);

nSres = min(nScan, spm_get_defaults('stats.maxres'));
resInMem = spm_get_defaults('stats.resmem');
VResI(1:nSres) = deal(struct(...
    'fname',   [],...
    'dim',     DIM,...
    'dt',      [spm_type('float64') spm_platform('bigend')],...
    'mat',     M,...
    'pinfo',   [1 0 0]',...
    'descrip', 'spm_spm:StandardisedResiduals',...
    metadata{:}));
if resInMem, for i=1:nSres, VResI(i).dat = zeros(VResI(i).dim); end; end

for i = 1:nSres
    VResI(i).fname   = [sprintf('ResI_%04d', i) file_ext];
    VResI(i).descrip = sprintf('spm_spm:ResI (%04d)', i);
end
VResI = spm_data_hdr_write(VResI);

iRes = round(linspace(1,nScan,nSres))';

for i = 1:numel(xM.VM)
    C = spm_bsplinc(xM.VM(i), [0 0 0 0 0 0]');
    v = true(DIM);
    [x1,x2] = ndgrid(1:DIM(1),1:DIM(2));
    for x3 = 1:DIM(3)
        M2  = inv(M\xM.VM(i).mat);
        y1 = M2(1,1)*x1+M2(1,2)*x2+(M2(1,3)*x3+M2(1,4));
        y2 = M2(2,1)*x1+M2(2,2)*x2+(M2(2,3)*x3+M2(2,4));
        y3 = M2(3,1)*x1+M2(3,2)*x2+(M2(3,3)*x3+M2(3,4));
        v(:,:,x3) = spm_bsplins(C, y1,y2,y3, [0 0 0 0 0 0]') > 0;
    end
    mask = mask & v;
    clear C v x1 x2 x3 M2 y1 y2 y3
end

%-Do not split data into chunks (edit by Marko)
chunk = [1:prod(DIM)];
%-Get data & construct analysis mask (or use data from memory, by Marko)
if nargin == 1 || (isempty(Y) && isempty(cmask))

    cmask = mask(chunk);
    Y     = zeros(nScan,numel(chunk));
    for j=1:nScan
        if ~any(cmask), break, end

        Y(j,cmask) = spm_data_read(VY(j),chunk(cmask));

        cmask(cmask) = Y(j,cmask) > xM.TH(j);
        if xM.I && ~YNaNrep && xM.TH(j) < 0
            cmask(cmask) = abs(Y(j,cmask)) > eps;
        end
    end
    cmask(cmask) = any(diff(Y(:,cmask),1));
    Y            = Y(:,cmask);
end

KWY          = spm_filter(xX.K,W*Y);

beta         = xX.pKX*KWY;
if any(cmask)
    res      = spm_sp('r',xX.xKXs,KWY);
else
    res      = zeros(nScan,0);
end
ResSS        = sum(res.^2);
res          = res(iRes,:);

c            = NaN(numel(chunk),1);
mask(chunk)  = cmask;
VM           = spm_data_write(VM, cmask', chunk);
for j=1:nBeta
    c(cmask) = beta(j,:);
    Vbeta(j) = spm_data_write(Vbeta(j), c, chunk);
end
c(cmask)     = ResSS / xX.trRV;
VResMS       = spm_data_write(VResMS, c, chunk);
for j=1:nSres
    c(cmask) = res(j,:)./sqrt(ResSS/erdf); % or xX.erdf
    VResI(j) = spm_data_write(VResI(j), c, chunk);
end

try
    if ~strcmpi(spm_get_defaults('modality'),'fmri')
        ResMS  = spm_data_read(VResMS);
        ResMS  = ResMS + 1e-3 * max(ResMS(isfinite(ResMS)));
        VResMS = spm_data_write(VResMS, ResMS);
        clear ResMS
    end
end

if ~spm_mesh_detect(VY)
    [FWHM,VRpv,R] = mw_spm_est_smoothness(VResI,VM,[nScan erdf]);
else
    VRpv = struct(...
        'fname',   ['RPV' file_ext],...
        'dim',     DIM,...
        'dt',      [spm_type('float64') spm_platform('bigend')],...
        'mat',     M,...
        'pinfo',   [1 0 0]',...
        'descrip', 'spm_spm: resels per voxel',...
        metadata{:});
    VRpv = spm_data_hdr_write(VRpv);
    ResI = zeros(prod(DIM),numel(VResI));
    for i=1:numel(VResI)
        ResI(:,i) = spm_data_read(VResI(i));
    end
    g = metadata{2};
    if isempty(spm_file(g,'path'))
        g = fullfile(spm_file(VY(1).fname,'path'),g);
    end
    [R, RPV] = spm_mesh_resels(gifti(g),mask,ResI,[nScan erdf]);
    RPV(~mask) = NaN;
    VRpv = spm_data_write(VRpv,RPV);
    FWHM = [1 1 1] * (1/mean(RPV(mask))).^(1/3);
end

fres = cellstr(spm_select('FPList',SPM.swd,'^ResI_.{4}\..{3}$'));
for i=1:numel(fres)
    spm_unlink(fres{i});
end

xX.nKX         = spm_DesMtx('sca',xX.xKXs.X,xX.name);
[x,y,z]        = ind2sub(DIM,find(mask));
XYZ            = [x y z]';

SPM.xVol.XYZ   = XYZ;
SPM.xVol.M     = M;
SPM.xVol.iM    = inv(M);
SPM.xVol.DIM   = DIM';
SPM.xVol.FWHM  = FWHM;
SPM.xVol.R     = R;
SPM.xVol.S     = nnz(mask);
SPM.xVol.VRpv  = VRpv;
if spm_mesh_detect(VY)
    SPM.xVol.G = g;
end
SPM.Vbeta      = Vbeta;
SPM.VResMS     = VResMS;
SPM.VM         = VM;
SPM.xVi        = xVi;
SPM.xX         = xX;
SPM.xM         = xM;
SPM.xCon       = struct([]);
SPM.SPMid      = 'SPM12: mw_spm_spm (v6678) from mw_optcens';
SPM.swd        = pwd;
save('SPM.mat','SPM', spm_get_defaults('mat.format'));
return;

%% ==========================================================================================================

function [xVi, mask] = mw_spm_est_non_sphericity(SPM, Y, cmask);
% Marko's stripped down version of spm_est_non_sphericity (6015), including
% an option to read stuff from memory

VY           = SPM.xY.VY;
DIM          = VY(1).dim;
YNaNrep      = spm_type(VY(1).dt(1),'nanrep');
xX           = SPM.xX;
nScan        = size(xX.X,1);
if ~isfield(xX,'W')
    xX.W     = speye(nScan,nScan);
end
xM           = SPM.xM;
xVi          = SPM.xVi;

xX.xKXs      = spm_sp('Set',spm_filter(xX.K,xX.W*xX.X));
xX.xKXs.X    = full(xX.xKXs.X);
xX.pKX       = spm_sp('x-',xX.xKXs);
if isfield(xVi,'Fcontrast')
    Fcname   = 'User-specified contrast';
    xCon     = spm_FcUtil('Set',Fcname,'F','c',xVi.Fcontrast,xX.xKXs);
else
    Fcname   = 'effects of interest';
    iX0      = [xX.iB xX.iG];
    xCon     = spm_FcUtil('Set',Fcname,'F','iX0',iX0,xX.xKXs);
end
if ~isempty(xCon(1).c)
    X1o      = spm_FcUtil('X1o', xCon(1),xX.xKXs);
    Hsqr     = spm_FcUtil('Hsqr',xCon(1),xX.xKXs);
    trMV     = spm_SpUtil('trMV',X1o);
else
    trMV     = 1;
    Hsqr     = Inf;
end
trRV         = spm_SpUtil('trRV',xX.xKXs);

%try
%    modality = lower(spm_get_defaults('modality')); UFp      =
%    spm_get_defaults(['stats.' modality '.ufp']);
%catch
%    UFp      = 0.001;
%end
UFp          = 0.5;  % adapted by Marko
xVi.UFp      = UFp;
UF           = spm_invFcdf(1 - UFp,[trMV,trRV]);

mask = true(DIM);
for i = 1:numel(xM.VM)
    C = spm_bsplinc(xM.VM(i), [0 0 0 0 0 0]');
    v = true(DIM);
    [x1,x2] = ndgrid(1:DIM(1),1:DIM(2));
    for x3 = 1:DIM(3)
        M2 = inv(VY(1).mat\xM.VM(i).mat);
        y1 = M2(1,1)*x1+M2(1,2)*x2+(M2(1,3)*x3+M2(1,4));
        y2 = M2(2,1)*x1+M2(2,2)*x2+(M2(2,3)*x3+M2(2,4));
        y3 = M2(3,1)*x1+M2(3,2)*x2+(M2(3,3)*x3+M2(3,4));
        v(:,:,x3) = spm_bsplins(C, y1,y2,y3, [0 0 0 0 0 0]') > 0;
    end
    mask = mask & v;
    clear C v x1 x2 x3 M2 y1 y2 y3
end

Cy        = 0;
Cm        = mask;

chunk = 1:prod(DIM);

if nargin == 1 || (isempty(Y) && isempty(cmask))

    Y       = zeros(nScan,numel(chunk));
    cmask   = mask(chunk);
    for j=1:nScan
        if ~any(cmask), break, end

        Y(j,cmask) = spm_data_read(VY(j),chunk(cmask));

        cmask(cmask) = Y(j,cmask) > xM.TH(j);
        if xM.I && ~YNaNrep && xM.TH(j) < 0
            cmask(cmask) = abs(Y(j,cmask)) > eps;
        end
    end
    cmask(cmask) = any(diff(Y(:,cmask),1));
    mask(chunk)  = cmask;
    Cm(chunk)    = cmask;
    if ~any(cmask), return, end
    Y       = Y(:,cmask);

end

KWY      = spm_filter(xX.K,xX.W*Y);

beta    = xX.pKX*KWY;
if any(cmask)
    res = spm_sp('r',xX.xKXs,KWY);
else
    res = zeros(nScan,0);
end
ResSS   = sum(res.^2);
clear res

j       = sum((Hsqr*beta).^2,1)/trMV > UF*ResSS/trRV;
Cm(chunk(cmask)) = j;
q       = nnz(j);
if q
    q   = spdiags(sqrt(trRV./ResSS(j)'),0,q,q);
    Y   = Y(:,j)*q;
    Cy  = Cy + Y*Y';
end

s = nnz(Cm);
nsv = 0;
if ~s
    % error('Please check your data: There are no significant voxels.');
    nsv = 1;  % potentially use later, but don't break things here
    s = 1;
end
Cy = Cy / s;

try

    if isstruct(xX.K)
        m     = length(xVi.Vi);
        h     = zeros(m,1);
        V     = sparse(nScan,nScan);
        for i = 1:length(xX.K)

            q     = xX.K(i).row;
            p     = [];
            Qp    = {};
            for j = 1:m
                if nnz(xVi.Vi{j}(q,q))
                    Qp{end + 1} = xVi.Vi{j}(q,q);
                    p           = [p j];
                end
            end

            Xp     = xX.X(q,:);
            try
                Xp = [Xp xX.K(i).X0];
            end

            evalc('[Vp,hp] = spm_reml(Cy(q,q),Xp,Qp);');
            V(q,q)  = V(q,q) + Vp;
            h(p)    = hp;
        end

    else
        evalc('[V,h] = spm_reml(Cy,xX.X,xVi.Vi);');
    end

catch
    evalc('[V,h] = spm_reml(Cy,xX.X,xVi.Vi);');

end

V      = V*nScan/trace(V);
xVi.h  = h;
xVi.V  = V;
xVi.Cy = Cy;
return;

%% ==========================================================================================================

function [FWHM,VRpv,R] = mw_spm_est_smoothness(V,VM,ndf)
% Marko's stripped version of spm_est_smoothness (6157)

if nargin < 3, ndf = [NaN NaN]; end

if ~isstruct(V)
    V  = spm_vol(V);
end
if ~isstruct(VM)
    VM = spm_vol(VM);
end
if any(isnan(ndf))
    ndf     = [length(V) length(V)];
end
n_full = ndf(1);
edf    = ndf(2);

VRpv  = struct( 'fname',['RPV' spm_file_ext],...
    'dim',      VM.dim(1:3),...
    'dt',       [spm_type('float64') spm_platform('bigend')],...
    'mat',      VM.mat,...
    'pinfo',    [1 0 0]',...
    'descrip',  'spm_spm: resels per voxel');
VRpv   = spm_create_vol(VRpv);

D        = 3 - sum(VM.dim(1:3) == 1);
if D == 0
    FWHM = [Inf Inf Inf];
    R    = [0 0 0];
    return;
end

d          = spm_read_vols(VM);
[Ix,Iy,Iz] = ndgrid(1:VM.dim(1),1:VM.dim(2),1:VM.dim(3));
Ix = Ix(d~=0); Iy = Iy(d~=0); Iz = Iz(d~=0);

L     = zeros(size(Ix,1),D,D);
ssq   = zeros(size(Ix,1),1);
for i = 1:length(V)

    [d,dx,dy,dz] = spm_sample_vol(V(i),Ix,Iy,Iz,1);

    ssq  = ssq + d.^2;

    if D >= 1
        L(:,1,1) = L(:,1,1) + dx.*dx;
    end
    if D >= 2
        L(:,1,2) = L(:,1,2) + dx.*dy;
        L(:,2,2) = L(:,2,2) + dy.*dy;
    end
    if D >= 3
        L(:,1,3) = L(:,1,3) + dx.*dz;
        L(:,2,3) = L(:,2,3) + dy.*dz;
        L(:,3,3) = L(:,3,3) + dz.*dz;
    end

end

L  = L/length(V);
L  = L*(n_full/edf);

if D == 1
    resel_xyz = L;
    resel_img = L;
end
if D == 2
    resel_xyz = [L(:,1,1) L(:,2,2)];
    resel_img = L(:,1,1).*L(:,2,2) - ...
        L(:,1,2).*L(:,1,2);
end
if D == 3
    resel_xyz = [L(:,1,1) L(:,2,2)  L(:,3,3)];
    resel_img = L(:,1,1).*L(:,2,2).*L(:,3,3) + ...
        L(:,1,2).*L(:,2,3).*L(:,1,3)*2 - ...
        L(:,1,1).*L(:,2,3).*L(:,2,3) - ...
        L(:,1,2).*L(:,1,2).*L(:,3,3) - ...
        L(:,1,3).*L(:,2,2).*L(:,1,3);
end
resel_img(resel_img<0) = 0;
resel_img = sqrt(resel_img/(4*log(2))^D);
resel_xyz = sqrt(resel_xyz/(4*log(2)));

if spm_get_defaults('stats.rft.nonstat')
    fwhm_vox = 3;
else
    fwhm_vox = 0;
end
if any(fwhm_vox)
    if length(fwhm_vox) == 1, fwhm_vox = fwhm_vox([1 1 1]);  end

    mask = spm_read_vols(VM) > 0;
    RPV = zeros(size(mask));
    RPV(mask) = resel_img;
    smask = double(mask & isfinite(RPV));
    spm_smooth(RPV, RPV, fwhm_vox);
    spm_smooth(smask, smask, fwhm_vox);
    infer = smask > 1e-3;
    RPV( infer) = RPV(infer) ./ smask(infer);
    RPV(~infer) = NaN;
    resel_img = RPV(mask);
end

for i = 1:VM.dim(3)
    d = NaN(VM.dim(1:2));
    I = find(Iz == i);
    if ~isempty(I)
        d(sub2ind(VM.dim(1:2), Ix(I), Iy(I))) = resel_img(I);
    end
    VRpv = spm_write_plane(VRpv, d, i);
end

i     = isnan(ssq) | ssq < sqrt(eps);
resel_img = mean(resel_img(~i,:));
resel_xyz = mean(resel_xyz(~i,:));

RESEL = resel_img^(1/D)*(resel_xyz/prod(resel_xyz)^(1/D));
FWHM  = full(sparse(1,1:D,1./RESEL,1,3));
FWHM(isnan(FWHM)) = 0;
FWHM(~FWHM) = 1;

R      = spm_resels_vol(VM,FWHM)';
return;

%% ==========================================================================================================

function out = mw_spm_run_fmri_spec(job)
% Marko's stripped version of spm_run_fmri_spec.m (6562)


files = {'^mask\..{3}$','^ResMS\..{3}$','^RPV\..{3}$',...
    '^beta_.{4}\..{3}$','^con_.{4}\..{3}$','^ResI_.{4}\..{3}$',...
    '^ess_.{4}\..{3}$', '^spm\w{1}_.{4}\..{3}$'};

for i=1:length(files)
    j = spm_select('List',pwd,files{i});
    for k=1:size(j,1)
        spm_unlink(deblank(j(k,:)));
    end
end

SPM.xY.RT = job.timing.RT;
SPM.xBF.UNITS = job.timing.units;
SPM.xBF.T     = job.timing.fmri_t;
SPM.xBF.T0    = job.timing.fmri_t0;

bf = char(fieldnames(job.bases));
if strcmp(bf,'none')
    SPM.xBF.name = 'NONE';
elseif strcmp(bf,'hrf')
    if all(job.bases.hrf.derivs == [0 0])
        SPM.xBF.name = 'hrf';
    elseif all(job.bases.hrf.derivs == [1 0])
        SPM.xBF.name = 'hrf (with time derivative)';
    elseif all(job.bases.hrf.derivs == [1 1])
        SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
    else
        error('Unknown HRF derivative choices.');
    end
else
    switch bf
        case 'fourier'
            SPM.xBF.name = 'Fourier set';
        case 'fourier_han'
            SPM.xBF.name = 'Fourier set (Hanning)';
        case 'gamma'
            SPM.xBF.name = 'Gamma functions';
        case 'fir'
            SPM.xBF.name = 'Finite Impulse Response';
        otherwise
            error('Unknown basis functions.');
    end
    SPM.xBF.length = job.bases.(bf).length;
    SPM.xBF.order  = job.bases.(bf).order;
end

SPM.xBF.Volterra = job.volt;

design_only = ~isfield(job,'mask');

if ~design_only
    SPM.xY.P = [];
end

for i = 1:numel(job.sess)

    sess = job.sess(i);

    if design_only
        SPM.nscan(i) = sess.nscan;
    else
        SPM.nscan(i) = numel(sess.scans);
        if SPM.nscan(i) == 1
            sess.scans   = cellstr(spm_select('Expand',sess.scans{1}));
            SPM.nscan(i) = numel(sess.scans);
        end
        SPM.xY.P     = strvcat(SPM.xY.P,sess.scans{:});
    end
    if SPM.nscan(i) == 1
        error('Not enough scans in session %d.',i);
    end

    if ~isempty(sess.multi{1})

        try
            multicond = load(sess.multi{1});
        catch
            error('Cannot load %s',sess.multi{1});
        end

        if ~all(isfield(multicond, {'names','onsets','durations'})) || ...
                ~iscell(multicond.names) || ...
                ~iscell(multicond.onsets) || ...
                ~iscell(multicond.durations) || ...
                ~isequal(numel(multicond.names), numel(multicond.onsets), ...
                numel(multicond.durations))
            error(['Multiple conditions MAT-file ''%s'' is invalid:\n',...
                'File must contain names, onsets, and durations '...
                'cell arrays of equal length.\n'],sess.multi{1});
        end

        for j=1:numel(multicond.onsets)

            cond.name     = multicond.names{j};
            if isempty(cond.name)
                error('MultiCond file: sess %d cond %d has no name.',i,j);
            end
            cond.onset    = multicond.onsets{j};
            if isempty(cond.onset)
                error('MultiCond file: sess %d cond %d has no onset.',i,j);
            end
            cond.duration = multicond.durations{j};
            if isempty(cond.onset)
                error('MultiCond file: sess %d cond %d has no duration.',i,j);
            end

            if ~isfield(multicond,'tmod');
                cond.tmod = 0;
            else
                try
                    cond.tmod = multicond.tmod{j};
                catch
                    error('Error specifying time modulation.');
                end
            end

            cond.pmod = [];
            if isfield(multicond,'pmod')
                try

                    if (j <= numel(multicond.pmod)) && ...
                            ~isempty(multicond.pmod(j).name)
                        for ii = 1:numel(multicond.pmod(j).name)
                            cond.pmod(ii).name  = multicond.pmod(j).name{ii};
                            cond.pmod(ii).param = multicond.pmod(j).param{ii};
                            cond.pmod(ii).poly  = multicond.pmod(j).poly{ii};
                        end
                    end
                catch
                    warning('Error specifying parametric modulation.');
                    rethrow(lasterror);
                end
            end

            if isfield(multicond,'orth') && (j <= numel(multicond.orth))
                cond.orth    = multicond.orth{j};
            else
                cond.orth    = true;
            end

            sess.cond(end+1) = cond;
        end
    end

    U = [];

    if~isfield(sess, 'cond'),  sess.cond = [];  end
    for j = 1:numel(sess.cond)

        cond      = sess.cond(j);
        U(j).name = {cond.name};
        U(j).ons  = cond.onset(:);
        U(j).dur  = cond.duration(:);
        U(j).orth = cond.orth;
        if isempty(U(j).orth), U(j).orth = true; end
        if length(U(j).dur) == 1
            U(j).dur = repmat(U(j).dur,size(U(j).ons));
        elseif numel(U(j).dur) ~= numel(U(j).ons)
            error('Mismatch between number of onset and number of durations.');
        end

        P  = [];
        q1 = 0;

        switch job.timing.units
            case 'secs'
                sf    = 1 / 60;
            case 'scans'
                sf    = job.timing.RT / 60;
            otherwise
                error('Unknown unit "%s".',job.timing.units);
        end
        if cond.tmod > 0
            P(1).name = 'time';
            P(1).P    = U(j).ons * sf;
            P(1).h    = cond.tmod;
            q1        = 1;
        end
        if ~isempty(cond.pmod)
            for q = 1:numel(cond.pmod)
                q1 = q1 + 1;
                P(q1).name = cond.pmod(q).name;
                P(q1).P    = cond.pmod(q).param(:);
                P(q1).h    = cond.pmod(q).poly;
            end
        end
        if isempty(P)
            P.name = 'none';
            P.h    = 0;
        end

        U(j).P = P;

    end

    SPM.Sess(i).U = U;

    C     = [];
    Cname = cell(1,numel(sess.regress));
    for q = 1:numel(sess.regress)
        Cname{q} = sess.regress(q).name;
        if numel(sess.regress(q).val(:)) ~= SPM.nscan(i)
            error('Length of regressor is not commensurate with data points.');
        end
        C        = [C, sess.regress(q).val(:)];
    end

    if ~isempty(sess.multi_reg{1})
        for q=1:numel(sess.multi_reg)
            tmp   = load(sess.multi_reg{q});
            names = {};
            if isstruct(tmp) % .mat
                if isfield(tmp,'R')
                    R = tmp.R;
                    if isfield(tmp,'names')
                        names = tmp.names;
                    end
                elseif isfield(tmp,'xY');
                    R = tmp.xY.u;
                    names = {tmp.xY.name};
                elseif isfield(tmp,'PPI')
                    R    = [tmp.PPI.ppi tmp.PPI.Y tmp.PPI.P];
                    names = {...
                        ['PPI Interaction: ' tmp.PPI.name],...
                        ['Main Effect: ' tmp.PPI.xY.name ' BOLD'],...
                        'Main Effect: "psych" condition'};
                else
                    error(['Variable ''R'' not found in multiple ' ...
                        'regressor file ''%s''.'], sess.multi_reg{q});
                end
            elseif isnumeric(tmp) % .txt
                R     = tmp;
            end

            if size(R,1) ~= SPM.nscan(i)
                error('Length of regressor is not commensurate with data points.');
            end
            if ~isempty(names) && numel(names) ~= size(R,2)
                warning('Mismatch between number of regressors and their names.');
                names = {};
            end
            C  = [C, R];
            if numel(sess.multi_reg) == 1
                nb_mult = '';
            else
                nb_mult = sprintf('_%d',q);
            end
            for j=1:size(R,2)
                if isempty(names)
                    Cname{end+1} = sprintf('R%d%s',j,nb_mult);
                else
                    Cname{end+1} = sprintf('%s%s',names{j},nb_mult);
                end
            end
        end
    end

    SPM.Sess(i).C.C    = C;
    SPM.Sess(i).C.name = Cname;

end

if ~isempty(job.fact)
    for i=1:numel(job.fact)
        SPM.factor(i).name   = job.fact(i).name;
        SPM.factor(i).levels = job.fact(i).levels;
    end
    if prod([SPM.factor.levels]) ~= numel(SPM.Sess(1).U)
        error('Factors do not match conditions');
    end
else
    SPM.factor = [];
end

SPM.xGX.iGXcalc = job.global;

SPM.xM.gMT = job.mthresh;

for i = 1:numel(job.sess),
    SPM.xX.K(i).HParam = job.sess(i).hpf;
end

SPM.xVi.form = job.cvi;


SPM = mw_spm_fmri_spm_ui(SPM);

if ~design_only
    if ~isempty(job.mask{1})
        SPM.xM.VM         = spm_data_hdr_read(job.mask{1});
        SPM.xM.xs.Masking = [SPM.xM.xs.Masking, '+explicit mask'];
    end
end

save('SPM.mat','SPM', spm_get_defaults('mat.format'));
out.spmmat{1} = fullfile(pwd, 'SPM.mat');

return;

%% ==========================================================================================================

function [SPM] = mw_spm_fmri_spm_ui(SPM)
% Marko's stripped version of spm_fmri_spm_ui.m (6088)

SPM   = spm_fMRI_design(SPM);
nscan = SPM.nscan;
nsess = length(nscan);


try
    HParam     = [SPM.xX.K(:).HParam];
    if length(HParam) == 1
        HParam = repmat(HParam,1,nsess);
    end
catch
    error('High-pass filter not specified.');
end

for  i = 1:nsess
    K(i) = struct('HParam', HParam(i),...
        'row',    SPM.Sess(i).row,...
        'RT',     SPM.xY.RT);
end
SPM.xX.K = spm_filter(K);


try
    cVi  = SPM.xVi.form;
catch
    error('Serial correlations not specified.');
end


if ~ischar(cVi)
    SPM.xVi.Vi = spm_Ce(nscan,cVi(1));
    cVi        = ['AR( ' sprintf('%0.1f ',cVi) ')'];

else
    switch lower(cVi)

        case {'i.i.d', 'none'}
            SPM.xVi.V  = speye(sum(nscan));
            cVi        = 'i.i.d';

        case 'fast'
            dt = SPM.xY.RT;
            Q  = {};
            l  = sum(nscan);
            k  = 0;
            for m=1:length(nscan)
                T     = (0:(nscan(m) - 1))*dt;
                d     = 2.^(floor(log2(dt/4)):log2(64));
                for i = 1:length(d)
                    for j = 0:2
                        QQ = toeplitz((T.^j).*exp(-T/d(i)));
                        [x,y,q] = find(QQ);
                        x = x + k;
                        y = y + k;
                        Q{end + 1} = sparse(x,y,q,l,l);
                    end
                end
                k = k + nscan(m);
            end
            SPM.xVi.Vi = Q;
            cVi        = upper(cVi);

        otherwise
            SPM.xVi.Vi = spm_Ce(nscan,0.2);
            cVi        = 'AR(0.2)';
    end
end

SPM.xVi.form = cVi;


for i     = 1:nsess, ntr(i) = length(SPM.Sess(i).U); end
Fstr      = sprintf('[min] Cutoff: %d {s}',min([SPM.xX.K(:).HParam]));
SPM.xsDes = struct(...
    'Basis_functions',      SPM.xBF.name,...
    'Number_of_sessions',   sprintf('%d',nsess),...
    'Trials_per_session',   sprintf('%-3d',ntr),...
    'Interscan_interval',   sprintf('%0.2f {s}',SPM.xY.RT),...
    'High_pass_Filter',     Fstr);


try, SPM.xY.P; catch, return; end

VY    = spm_data_hdr_read(SPM.xY.P);
SPM.xY.VY = VY;


GM    = 100;
q     = length(VY);
g     = zeros(q,1);
if spm_mesh_detect(VY)
    for i = 1:q
        dat = spm_data_read(VY(i));
        g(i) = mean(dat(~isnan(dat)));
    end
else
    for i = 1:q
        g(i) = spm_global(VY(i));
    end
end

gSF   = GM./g;
if strcmpi(SPM.xGX.iGXcalc,'none')
    for i = 1:nsess
        gSF(SPM.Sess(i).row) = GM./mean(g(SPM.Sess(i).row));
    end
end

for i = 1:q
    SPM.xY.VY(i).pinfo(1:2,:) = SPM.xY.VY(i).pinfo(1:2,:) * gSF(i);
    if spm_mesh_detect(VY)
        SPM.xY.VY(i).private.private.data{1}.data.scl_slope = ...
            SPM.xY.VY(i).private.private.data{1}.data.scl_slope * gSF(i);
        SPM.xY.VY(i).private.private.data{1}.data.scl_inter = ...
            SPM.xY.VY(i).private.private.data{1}.data.scl_inter * gSF(i);
    else
        SPM.xY.VY(i).private.dat.scl_slope = ...
            SPM.xY.VY(i).private.dat.scl_slope * gSF(i);
        SPM.xY.VY(i).private.dat.scl_inter = ...
            SPM.xY.VY(i).private.dat.scl_inter * gSF(i);
    end
end

SPM.xGX.sGXcalc = 'mean voxel value';
SPM.xGX.sGMsca  = 'session specific';
SPM.xGX.rg      = g;
SPM.xGX.GM      = GM;
SPM.xGX.gSF     = gSF;

try
    gMT = SPM.xM.gMT;
catch
    gMT = spm_get_defaults('mask.thresh');
end
TH = g.*gSF*gMT;

SPM.xM = struct(...
    'T',   ones(q,1),...
    'TH',  TH,...
    'gMT', gMT,...
    'I',   0,...
    'VM',  {[]},...
    'xs',  struct('Masking','analysis threshold'));


xs = struct(...
    'Global_calculation',   SPM.xGX.sGXcalc,...
    'Grand_mean_scaling',   SPM.xGX.sGMsca,...
    'Global_normalisation', SPM.xGX.iGXcalc);
for fn=(fieldnames(xs))', SPM.xsDes.(fn{1}) = xs.(fn{1}); end

save('SPM.mat', 'SPM', spm_get_defaults('mat.format'));

return;

%% ==========================================================================================================

function B=inpaint_nans(A,method)
% INPAINT_NANS: in-paints over nans in an array usage: B=INPAINT_NANS(A)
% % default method usage: B=INPAINT_NANS(A,method)   % specify method used
%
% Solves approximation to one of several pdes to interpolate and
% extrapolate holes in an array
%
% arguments (input):
%   A - nxm array with some NaNs to be filled in
%
%   method - (OPTIONAL) scalar numeric flag - specifies
%       which approach (or physical metaphor to use for the interpolation.)
%       All methods are capable of extrapolation, some are better than
%       others. There are also speed differences, as well as accuracy
%       differences for smooth surfaces.
%
%       methods {0,1,2} use a simple plate metaphor. method  3 uses a
%       better plate equation,
%                 but may be much slower and uses more memory.
%       method  4 uses a spring metaphor. method  5 is an 8 neighbor
%       average, with no
%                 rationale behind it compared to the other methods. I do
%                 not recommend its use.
%
%       method == 0 --> (DEFAULT) see method 1, but
%         this method does not build as large of a linear system in the
%         case of only a few NaNs in a large array. Extrapolation behavior
%         is linear.
%
%       method == 1 --> simple approach, applies del^2
%         over the entire array, then drops those parts of the array which
%         do not have any contact with NaNs. Uses a least squares approach,
%         but it does not modify known values. In the case of small arrays,
%         this method is quite fast as it does very little extra work.
%         Extrapolation behavior is linear.
%
%       method == 2 --> uses del^2, but solving a direct
%         linear system of equations for nan elements. This method will be
%         the fastest possible for large systems since it uses the sparsest
%         possible system of equations. Not a least squares approach, so it
%         may be least robust to noise on the boundaries of any holes. This
%         method will also be least able to interpolate accurately for
%         smooth surfaces. Extrapolation behavior is linear.
%
%         Note: method 2 has problems in 1-d, so this method is disabled
%         for vector inputs.
%
%       method == 3 --+ See method 0, but uses del^4 for
%         the interpolating operator. This may result in more accurate
%         interpolations, at some cost in speed.
%
%       method == 4 --+ Uses a spring metaphor. Assumes
%         springs (with a nominal length of zero) connect each node with
%         every neighbor (horizontally, vertically and diagonally) Since
%         each node tries to be like its neighbors, extrapolation is as a
%         constant function where this is consistent with the neighboring
%         nodes.
%
%       method == 5 --+ See method 2, but use an average
%         of the 8 nearest neighbors to any element. This method is NOT
%         recommended for use.
%
%
% arguments (output):
%   B - nxm array with NaNs replaced
%
%
% Example:
%  [x,y] = meshgrid(0:.01:1); z0 = exp(x+y); znan = z0; znan(20:50,40:70) =
%  NaN; znan(30:90,5:10) = NaN; znan(70:75,40:90) = NaN;
%
%  z = inpaint_nans(znan);
%
%
% See also: griddata, interp1
%
% Author: John D'Errico e-mail address: woodchips@rochester.rr.com Release:
% 2 Release date: 4/15/06


% I always need to know which elements are NaN, and what size the array is
% for any method
[n,m]=size(A);
A=A(:);
nm=n*m;
k=isnan(A(:));

% list the nodes which are known, and which will be interpolated
nan_list=find(k);
known_list=find(~k);

% how many nans overall
nan_count=length(nan_list);

% convert NaN indices to (r,c) form nan_list==find(k) are the unrolled
% (linear) indices (row,column) form
[nr,nc]=ind2sub([n,m],nan_list);

% both forms of index in one array: column 1 == unrolled index column 2 ==
% row index column 3 == column index
nan_list=[nan_list,nr,nc];

% supply default method
if (nargin<2) || isempty(method)
    method = 0;
elseif ~ismember(method,0:5)
    error 'If supplied, method must be one of: {0,1,2,3,4,5}.'
end

% for different methods
switch method
    case 0
        % The same as method == 1, except only work on those elements which
        % are NaN, or at least touch a NaN.

        % is it 1-d or 2-d?
        if (m == 1) || (n == 1)
            % really a 1-d case
            work_list = nan_list(:,1);
            work_list = unique([work_list;work_list - 1;work_list + 1]);
            work_list(work_list <= 1) = [];
            work_list(work_list >= nm) = [];
            nw = numel(work_list);

            u = (1:nw)';
            fda = sparse(repmat(u,1,3),bsxfun(@plus,work_list,-1:1), ...
                repmat([1 -2 1],nw,1),nw,nm);
        else
            % a 2-d case

            % horizontal and vertical neighbors only
            talks_to = [-1 0;0 -1;1 0;0 1];
            neighbors_list=identify_neighbors(n,m,nan_list,talks_to);

            % list of all nodes we have identified
            all_list=[nan_list;neighbors_list];

            % generate sparse array with second partials on row variable
            % for each element in either list, but only for those nodes
            % which have a row index > 1 or < n
            L = find((all_list(:,2) > 1) & (all_list(:,2) < n));
            nl=length(L);
            if nl>0
                fda=sparse(repmat(all_list(L,1),1,3), ...
                    repmat(all_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
                    repmat([1 -2 1],nl,1),nm,nm);
            else
                fda=spalloc(n*m,n*m,size(all_list,1)*5);
            end

            % 2nd partials on column index
            L = find((all_list(:,3) > 1) & (all_list(:,3) < m));
            nl=length(L);
            if nl>0
                fda=fda+sparse(repmat(all_list(L,1),1,3), ...
                    repmat(all_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
                    repmat([1 -2 1],nl,1),nm,nm);
            end
        end

        % eliminate knowns
        rhs=-fda(:,known_list)*A(known_list);
        k=find(any(fda(:,nan_list(:,1)),2));

        % and solve...
        B=A;
        B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);

    case 1
        % least squares approach with del^2. Build system for every array
        % element as an unknown, and then eliminate those which are knowns.

        % Build sparse matrix approximating del^2 for every element in A.

        % is it 1-d or 2-d?
        if (m == 1) || (n == 1)
            % a 1-d case
            u = (1:(nm-2))';
            fda = sparse(repmat(u,1,3),bsxfun(@plus,u,0:2), ...
                repmat([1 -2 1],nm-2,1),nm-2,nm);
        else
            % a 2-d case

            % Compute finite difference for second partials on row variable
            % first
            [i,j]=ndgrid(2:(n-1),1:m);
            ind=i(:)+(j(:)-1)*n;
            np=(n-2)*m;
            fda=sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
                repmat([1 -2 1],np,1),n*m,n*m);

            % now second partials on column variable
            [i,j]=ndgrid(1:n,2:(m-1));
            ind=i(:)+(j(:)-1)*n;
            np=n*(m-2);
            fda=fda+sparse(repmat(ind,1,3),[ind-n,ind,ind+n], ...
                repmat([1 -2 1],np,1),nm,nm);
        end

        % eliminate knowns
        rhs=-fda(:,known_list)*A(known_list);
        k=find(any(fda(:,nan_list),2));

        % and solve...
        B=A;
        B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);

    case 2
        % Direct solve for del^2 BVP across holes generate sparse array
        % with second partials on row variable for each nan element, only
        % for those nodes which have a row index > 1 or < n

        % is it 1-d or 2-d?
        if (m == 1) || (n == 1)
            % really just a 1-d case
            error('Method 2 has problems for vector input. Please use another method.')

        else
            % a 2-d case
            L = find((nan_list(:,2) > 1) & (nan_list(:,2) < n));
            nl=length(L);
            if nl>0
                fda=sparse(repmat(nan_list(L,1),1,3), ...
                    repmat(nan_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
                    repmat([1 -2 1],nl,1),n*m,n*m);
            else
                fda=spalloc(n*m,n*m,size(nan_list,1)*5);
            end

            % 2nd partials on column index
            L = find((nan_list(:,3) > 1) & (nan_list(:,3) < m));
            nl=length(L);
            if nl>0
                fda=fda+sparse(repmat(nan_list(L,1),1,3), ...
                    repmat(nan_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
                    repmat([1 -2 1],nl,1),n*m,n*m);
            end

            % fix boundary conditions at extreme corners of the array in
            % case there were nans there
            if ismember(1,nan_list(:,1))
                fda(1,[1 2 n+1])=[-2 1 1];
            end
            if ismember(n,nan_list(:,1))
                fda(n,[n, n-1,n+n])=[-2 1 1];
            end
            if ismember(nm-n+1,nan_list(:,1))
                fda(nm-n+1,[nm-n+1,nm-n+2,nm-n])=[-2 1 1];
            end
            if ismember(nm,nan_list(:,1))
                fda(nm,[nm,nm-1,nm-n])=[-2 1 1];
            end

            % eliminate knowns
            rhs=-fda(:,known_list)*A(known_list);

            % and solve...
            B=A;
            k=nan_list(:,1);
            B(k)=fda(k,k)\rhs(k);

        end

    case 3
        % The same as method == 0, except uses del^4 as the interpolating
        % operator.

        % del^4 template of neighbors
        talks_to = [-2 0;-1 -1;-1 0;-1 1;0 -2;0 -1; ...
            0 1;0 2;1 -1;1 0;1 1;2 0];
        neighbors_list=identify_neighbors(n,m,nan_list,talks_to);

        % list of all nodes we have identified
        all_list=[nan_list;neighbors_list];

        % generate sparse array with del^4, but only for those nodes which
        % have a row & column index >= 3 or <= n-2
        L = find( (all_list(:,2) >= 3) & ...
            (all_list(:,2) <= (n-2)) & ...
            (all_list(:,3) >= 3) & ...
            (all_list(:,3) <= (m-2)));
        nl=length(L);
        if nl>0
            % do the entire template at once
            fda=sparse(repmat(all_list(L,1),1,13), ...
                repmat(all_list(L,1),1,13) + ...
                repmat([-2*n,-n-1,-n,-n+1,-2,-1,0,1,2,n-1,n,n+1,2*n],nl,1), ...
                repmat([1 2 -8 2 1 -8 20 -8 1 2 -8 2 1],nl,1),nm,nm);
        else
            fda=spalloc(n*m,n*m,size(all_list,1)*5);
        end

        % on the boundaries, reduce the order around the edges
        L = find((((all_list(:,2) == 2) | ...
            (all_list(:,2) == (n-1))) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1))) | ...
            (((all_list(:,3) == 2) | ...
            (all_list(:,3) == (m-1))) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1))));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(all_list(L,1),1,5), ...
                repmat(all_list(L,1),1,5) + ...
                repmat([-n,-1,0,+1,n],nl,1), ...
                repmat([1 1 -4 1 1],nl,1),nm,nm);
        end

        L = find( ((all_list(:,2) == 1) | ...
            (all_list(:,2) == n)) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1)));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(all_list(L,1),1,3), ...
                repmat(all_list(L,1),1,3) + ...
                repmat([-n,0,n],nl,1), ...
                repmat([1 -2 1],nl,1),nm,nm);
        end

        L = find( ((all_list(:,3) == 1) | ...
            (all_list(:,3) == m)) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1)));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(all_list(L,1),1,3), ...
                repmat(all_list(L,1),1,3) + ...
                repmat([-1,0,1],nl,1), ...
                repmat([1 -2 1],nl,1),nm,nm);
        end

        % eliminate knowns
        rhs=-fda(:,known_list)*A(known_list);
        k=find(any(fda(:,nan_list(:,1)),2));

        % and solve...
        B=A;
        B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);

    case 4
        % Spring analogy interpolating operator.

        % list of all springs between a node and a horizontal or vertical
        % neighbor
        hv_list=[-1 -1 0;1 1 0;-n 0 -1;n 0 1];
        hv_springs=[];
        for i=1:4
            hvs=nan_list+repmat(hv_list(i,:),nan_count,1);
            k=(hvs(:,2)>=1) & (hvs(:,2)<=n) & (hvs(:,3)>=1) & (hvs(:,3)<=m);
            hv_springs=[hv_springs;[nan_list(k,1),hvs(k,1)]];
        end

        % delete replicate springs
        hv_springs=unique(sort(hv_springs,2),'rows');

        % build sparse matrix of connections, springs connecting diagonal
        % neighbors are weaker than the horizontal and vertical springs
        nhv=size(hv_springs,1);
        springs=sparse(repmat((1:nhv)',1,2),hv_springs, ...
            repmat([1 -1],nhv,1),nhv,nm);

        % eliminate knowns
        rhs=-springs(:,known_list)*A(known_list);

        % and solve...
        B=A;
        B(nan_list(:,1))=springs(:,nan_list(:,1))\rhs;

    case 5
        % Average of 8 nearest neighbors

        % generate sparse array to average 8 nearest neighbors for each nan
        % element, be careful around edges
        fda=spalloc(n*m,n*m,size(nan_list,1)*9);

        % -1,-1
        L = find((nan_list(:,2) > 1) & (nan_list(:,3) > 1));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([-n-1, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end

        % 0,-1
        L = find(nan_list(:,3) > 1);
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([-n, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end

        % +1,-1
        L = find((nan_list(:,2) < n) & (nan_list(:,3) > 1));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([-n+1, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end

        % -1,0
        L = find(nan_list(:,2) > 1);
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([-1, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end

        % +1,0
        L = find(nan_list(:,2) < n);
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([1, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end

        % -1,+1
        L = find((nan_list(:,2) > 1) & (nan_list(:,3) < m));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([n-1, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end

        % 0,+1
        L = find(nan_list(:,3) < m);
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([n, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end

        % +1,+1
        L = find((nan_list(:,2) < n) & (nan_list(:,3) < m));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([n+1, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end

        % eliminate knowns
        rhs=-fda(:,known_list)*A(known_list);

        % and solve...
        B=A;
        k=nan_list(:,1);
        B(k)=fda(k,k)\rhs(k);

end

% all done, make sure that B is the same shape as A was when we came in.
B=reshape(B,n,m);


% ====================================================
%      end of main function (inpaint_nans)
% ====================================================
% ====================================================
%      begin subfunctions for inpaint_nans
% ====================================================
function neighbors_list=identify_neighbors(n,m,nan_list,talks_to)
% identify_neighbors: identifies all the neighbors of
%   those nodes in nan_list, not including the nans themselves
%
% arguments (input):
%  n,m - scalar - [n,m]=size(A), where A is the
%      array to be interpolated
%  nan_list - array - list of every nan element in A
%      nan_list(i,1) == linear index of i'th nan element nan_list(i,2) ==
%      row index of i'th nan element nan_list(i,3) == column index of i'th
%      nan element
%  talks_to - px2 array - defines which nodes communicate
%      with each other, i.e., which nodes are neighbors.
%
%      talks_to(i,1) - defines the offset in the row
%                      dimension of a neighbor
%      talks_to(i,2) - defines the offset in the column
%                      dimension of a neighbor
%
%      For example, talks_to = [-1 0;0 -1;1 0;0 1] means that each node
%      talks only to its immediate neighbors horizontally and vertically.
%
% arguments(output):
%  neighbors_list - array - list of all neighbors of
%      all the nodes in nan_list

if ~isempty(nan_list)
    % use the definition of a neighbor in talks_to
    nan_count=size(nan_list,1);
    talk_count=size(talks_to,1);

    nn=zeros(nan_count*talk_count,2);
    j=[1,nan_count];
    for i=1:talk_count
        nn(j(1):j(2),:)=nan_list(:,2:3) + ...
            repmat(talks_to(i,:),nan_count,1);
        j=j+nan_count;
    end

    % drop those nodes which fall outside the bounds of the original array
    L = (nn(:,1)<1)|(nn(:,1)>n)|(nn(:,2)<1)|(nn(:,2)>m);
    nn(L,:)=[];

    % form the same format 3 column array as nan_list
    neighbors_list=[sub2ind([n,m],nn(:,1),nn(:,2)),nn];

    % delete replicates in the neighbors list
    neighbors_list=unique(neighbors_list,'rows');

    % and delete those which are also in the list of NaNs.
    neighbors_list=setdiff(neighbors_list,nan_list,'rows');

else
    neighbors_list=[];
end
%% ==========================================================================================================

function newvols = mw_art_repairvol(P,out_idx,outpth);
% Function mw_art_repairvol(varargin) from art_repair v5.2b this version
% modified for in- and outputs and stripped off some interactions and
% edited to suit Marko's needs used here with kind permission from Paul
% Mazaika

if nargin==2, outpth = [];  end

inputdir = cd;
cd(fileparts(P(1,:)));
inter_method = 3; % interpolation across scans, default

% Find the correct scans
allscan = [ 1: size(P,1) ];
for k = 1:length(out_idx)
    allscan(out_idx(k)) = 0;
end
in_idx = find(allscan>0);

% First copy over the scans that require no repair.
for j = 1:length(in_idx)
    V = spm_vol(P(in_idx(j),:));
    v = V;
    Y = spm_read_vols(V);
    [currpath, currname, currext] = fileparts(V.fname);
    if ~isempty(outpth),  currpath = outpth;  end

    copyname = ['v' currname currext];
    copy2 = fullfile(currpath, copyname);
    v.fname = copy2;   %  same name with a 'v' added
    v.private = [];    %  makes output files read-write.
    spm_write_vol(v,Y);
end

% Repair the outlier scans, and tag them as outliers.
for i = 1:length(out_idx)
    V = spm_vol(P(out_idx(i),:));
    v = V;
    Y = spm_read_vols(V);
    [currpath, currname, currext] = fileparts(V.fname);
    if ~isempty(outpth),  currpath = outpth;  end

    copyname = ['outlier_' currname currext];
    copyname_hdr = ['v' currname '.hdr'];
    copyname_img = ['v' currname currext];
    copy2 = fullfile(currpath, copyname);
    copy_hdr = fullfile(currpath, copyname_hdr);
    copy_img = fullfile(currpath, copyname_img);
    v.fname = copy2;
    v.private = [];  %  makes output files read-write.
    im_map = spm_vol(P(out_idx(i),:));

    if inter_method == 1;  % Replace outlier with mean image.
        v.fname = copy_img;  % was copy_hdr;
        spm_write_vol(v,Ym);
    end

    if inter_method == 2;  % Replace outlier with interpolated image.
        if out_idx(i) == 1 % Extrapolate for first scan
            im_in = spm_vol([P(2,:);P(3,:)]);
        elseif out_idx(i) == size(P,1) % Extrapolate for last scan
            im_in = spm_vol([P(out_idx(i)-2,:);P(out_idx(i)-1,:)]);
        else  %  Interpolate for most scans
            im_in = spm_vol([P(out_idx(i)-1,:);P(out_idx(i)+1,:)]);
        end
        Y1 = spm_read_vols(im_in(1));
        Y2 = spm_read_vols(im_in(2));
        Ym = (Y1 + Y2 )/2;
        v.fname = copy_img;  % was copy_hdr;
        spm_write_vol(v,Ym);
    end
    % New method. Interpolate between nearest non-outlier scans. Provides
    % linear interpolation over extended outliers.
    if inter_method == 3;  % Replace outlier with interpolated image.

        % Find nearest non-outlier scan in each direction
        yyy = find(allscan>out_idx(i));
        if (length(yyy)>0), highside = min(yyy); else, highside = 0; end
        yyy = find(allscan<out_idx(i) & allscan>0);
        if (length(yyy)>0), lowside = max(yyy); else, lowside = 0; end

        if lowside == 0 % Extrapolate from first good scan
            im_in = spm_vol(P(highside,:));
            Y1 = spm_read_vols(im_in(1));
            Ym = Y1;
        elseif highside == 0 % Extrapolate from last good scan
            im_in = spm_vol(P(lowside,:));
            Y1 = spm_read_vols(im_in(1));
            Ym = Y1;
        else  %  Interpolate for most scans
            im_in = spm_vol([P(lowside,:);P(highside,:)]);
            lenint = highside - lowside;
            hiwt = (out_idx(i)-lowside)/lenint;
            lowwt = (highside - out_idx(i))/lenint;
            Y1 = spm_read_vols(im_in(1));
            Y2 = spm_read_vols(im_in(2));
            Ym = Y1*lowwt + Y2*hiwt;
        end
        v.fname = copy_img;  % was copy_hdr;
        spm_write_vol(v,Ym);
    end
end

newvols = spm_select('FPList', currpath, ['^v.*' currext]);

cd(inputdir);
return;

%% ==========================================================================================================

