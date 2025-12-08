
function func_warpsmooth(home, sub, func, warptype, smoothing)

% maxNumCompThreads(1) ;
addpath('~/Documents/MATLAB/spm12') ;

disp(strcat("Starting func_warpsmooth for ", sub, " in ", home, ...
    " using smoothing=", smoothing)) ;

% SPMwarp, SPMCFMwarp, ANTSwarp (for warps whose folders are named that
% way), 'useThis' or some other label marking the folders with the best
% warps if not named consistently

% Default values (no way to set) 
% warptype = 'SPMwarp_lessRegularized';
% smoothat = 6;

smoothat = str2double(smoothing);

% processed = strcat(home, "/sub-", sub);
sub_anat =  strcat(home, "/sub-", sub, "/anat/") ;

% First, find the warp field
% warpfield = strcat(processed, "/anat/y_sub-", sub, "_T1w.nii") ;
warpfield_info = dir(strcat(sub_anat, "y_sub-*_T1w.nii")) ;
warpfield = strcat(warpfield_info.folder, "/", warpfield_info.name) 

if ~exist(warpfield, 'file')
    disp('Cannot find the requested warpfield. Aborting.')
    exit
else
    disp(strcat("Warpfield: ", warpfield))
end

% Identify the files
processed = strcat(home, "/sub-", sub, "/") ;
func_dir = strcat(processed, "/func/", func);
nii_names = struct2table(dir(strcat(func_dir, ...
    "/03_rmoco/rmoco_sub-*.nii"))).name;
nii_files = cellstr(strcat(func_dir, "/03_rmoco/", nii_names)) ;

clear matlabbatch

% Add SPM to path if necessary...
%addSPMifNecessary;

%% Prepare to specify batch from script
spm('defaults','fmri');
spm_jobman('initcfg');

% M-code from SPM (replace file/folder references with variables as
% necessary)

%% Warp co-registered functional images to standard space
% based on deformation field determined for the anatomical
matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(warpfield);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = nii_files;
% default
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = ...
    [-78 -112 -70; 78 76 85];
% Default is [2 2 2], and that's what WFU pickatlas ROIs fit. Original
% resolution was [3 3 3], but after coreg to the MPRAGE it's [1 1 1].
% Higher resolution = bigger files. PROBLEM: Repeatedly resampling data
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
% default
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
% default
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = ...
    '../04_warped/w';

%% Apply spatial smoothing
matlabbatch{2}.spm.spatial.smooth.data = ...
    cfg_dep('Normalise: Write: Normalised Images (Subj 1)', ...
    substruct('.','val', '{}', {1}, ...
    '.','val', '{}', {1}, ...
    '.','val', '{}', {1}, ...
    '.','val', '{}', {1}), ...
    substruct('()',{1}, '.','files'));
% default is [8 8 8]; set prefix below accordingly!
matlabbatch{2}.spm.spatial.smooth.fwhm = ...
    [smoothat smoothat smoothat];
% default
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
% default
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = ...
    strcat('../05_warpsmooth/s', int2str(smoothat));

%% Run the batch...
spm_jobman('run', matlabbatch);

end