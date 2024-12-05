% This function registers the functionals to native space

function func_1stlevel = func_1stlevel(spm_path, home, sub, rundir, ...
                                        expected_vols, TR, contrasts, ...
                                        events_tsv, moco_file)

% Add SPM library to path
addpath(spm_path) ;

% Setup structure
proc_dir = strcat(home, "/sub-", sub, "/func/", rundir) ;
disp(strcat("Looking in ", proc_dir));
output_dir = strcat(proc_dir, "/results") ;

% Check for input files

%   Removed prefix check - too hard to generalize
search_str = strcat(proc_dir, "/", "*.nii") ;
nii_files = dir(search_str) ;

% Let 0 mean skip verification
if expected_vols > 0
    if length(nii_files) ~= str2double(expected_vols)
        disp(strcat("Expected ", expected_vols, " found ", ...
            num2str(length(nii_files))))
        disp(nii_files)
        error("Length of files mismatch")
    end
elseif isempty(nii_files)
    % Finding 0 files, even if expected # is 0 is still wrong.
    disp("Didn't verify # of files, but found 0. Exiting!")
    exit()
else 
    disp("Not verifying input # of files, found: ")
    disp(length(nii_files))
end

% High pass filter cutoff in seconds. 
% If you set it to run_duration (in seconds), it will include a half 
% cosine and a full cosine over the course of the run, 
% removing linear and quadratic trends in the data. If you set it 
% higher, only linear trends will be removed (unless you set it higher 
% than 2*run_duration, in which case there will be no filtering at 
% all). If you set it to a smaller number, cosines with faster 
% frequency will be included. 
% Make sure not to get into the frequency range of your experimental 
% manipulation! 
run_duration = length(nii_files) * str2double(TR) ; 
hpf_cutoff = run_duration ;
disp(strcat("Using HPF cutoff: ", num2str(hpf_cutoff)))

% Calculate onsets

% events_tsv = strcat(home, "/sub-", sub, "/func/", ...
%                     strrep(rundir, "_bold", "_events.tsv")) 

if ~isfile(events_tsv)
    error(strcat("No events file: ", events_tsv))
else
    % read in TSV and turn it into a MATLAB struct
    events_table = readtable(events_tsv, "FileType", "text", ...
        'Delimiter', '\t');
end

trial_types = unique(events_table.trial_type) ;

if ismember("rest", trial_types)
    disp("Removing 'rest' from trial types.")
    trial_types = trial_types(trial_types~="rest") ;
end

%% Create struct to store onsets from files

% reformat data from input
TR = str2double(TR) ;

% Get trial types by calculating from contrasts
% [trial_types1, ~] = strsplit(contrasts, {' ', '>'}, ...
%     "CollapseDelimiters", true) ;

% trial_types = unique(trial_types1) 
contrasts = strsplit(contrasts, ' ')

disp(events_table) 
disp(trial_types)

% create new struct: name/onsets/duration (1 dur value)
conditions.names = trial_types ;

for i = 1:length(trial_types)

    % events_table.trial_type
    % trial_types{i}

    included_trials = strcmp(events_table.trial_type, trial_types{i}) ;
    new_table = events_table(included_trials, :) ;

    % Convert the onsets and durations to volumes/frames/TRs/scans
    onsets = floor(new_table.onset / TR) ;
    durations = ceil(unique(new_table.duration) / TR) ;

    conditions.onsets{i} = onsets ;
    conditions.durations{i} = durations ;

    disp(conditions.names{i})
    disp(conditions.onsets{i})
    disp(conditions.durations{i})

end

%% Setup matlabbatch

% output_dir = strcat(proc_dir, "/results") ;
mkdir(proc_dir, "results");

clear matlabbatch ;

% Prepare to specify batch from script

disp("Configuring SPM ...")
spm('defaults', 'fmri');
spm_jobman('initcfg');

% M-code from SPM (replace file/folder references with variables 
% as necessary)

% Model specification
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(output_dir);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

% Input files - reformat and include full path
nii_names = struct2table(nii_files).name;
nii_files = cellstr(strcat(proc_dir, "/", nii_names)) ;

matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = nii_files;

for i = 1:length(conditions.names) 

    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(i).name = ...
        conditions.names{i};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(i).onset = ...
        conditions.onsets{i};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(i).duration = ...
        conditions.durations{i};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(i).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(i).pmod = ...
        struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(i).orth = 1;

end

% Continue setting up matlabbatch
% IDK if this has to come after the constrasts, but it did originally
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = cellstr('');
matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = ... 
    struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = ...
    {moco_file};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = hpf_cutoff;

matlabbatch{1}.spm.stats.fmri_spec.fact = ...
    struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

% Model estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('fMRI model specification: SPM.mat File', ...
    substruct('.', 'val', '{}', {1}, ...
                '.', 'val', '{}', {1}, ...
                '.', 'val', '{}', {1}), ...
                substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

%% Contrasts
matlabbatch{3}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
            substruct('.','val', '{}', {2}, ...
                        '.','val', '{}', {1}, ...
                        '.','val', '{}',{1}), ...
                        substruct('.', 'spmmat'));

% Each trial>rest, each contrast, all>rest
names = strings(1, length(trial_types) + length(contrasts) + 1) ;

contrast_index = 1 ;
for i = 1:length(trial_types)

    % Create the design matrix so that all conditions are set to 0
    % except the one of interest 
    design = zeros(1, length(trial_types)) ;
    design(1, i) = 1 ;
    name = strcat(trial_types{i}, '>Rest') ;

    names(i) = name ;

    disp(name) 
    disp(design)

    matlabbatch{3}.spm.stats.con.consess{contrast_index}.tcon.name = ...
        name;
    matlabbatch{3}.spm.stats.con.consess{contrast_index}.tcon.weights = ...
        design;
    matlabbatch{3}.spm.stats.con.consess{contrast_index}.tcon.sessrep = ...
        'both';

    contrast_index = contrast_index + 1;

end

for i = 1:length(contrasts) 

    design = zeros(1, length(trial_types)) ;

    % extract the contrasts, lhs>rhs
    contrast = strsplit(contrasts{i}, '>') ;
    lhs = contrast(1) ;
    rhs = contrast(2) ;

    lhs_where = find(ismember(trial_types, lhs)) ;
    rhs_where = find(ismember(trial_types, rhs)) ; 

    design(lhs_where) = 1 ;
    design(rhs_where) = -1 ;

    disp(contrasts{i})
    disp(design)

    name_index = length(trial_types) + i ;
    names(name_index) = contrasts{i} ;

    matlabbatch{3}.spm.stats.con.consess{contrast_index}.tcon.name = ...
        contrasts{i};
    matlabbatch{3}.spm.stats.con.consess{contrast_index}.tcon.weights = ...
        design ;
    matlabbatch{3}.spm.stats.con.consess{contrast_index}.tcon.sessrep = ...
        'both';

    contrast_index = contrast_index + 1;

end

all_index = contrast_index ;
all_design = ones(1, length(trial_types)) ;
matlabbatch{3}.spm.stats.con.consess{all_index}.tcon.name = 'All>Rest';
matlabbatch{3}.spm.stats.con.consess{all_index}.tcon.weights = ...
    all_design ;
matlabbatch{3}.spm.stats.con.consess{all_index}.tcon.sessrep = 'both';

names(all_index) = "All>Rest" ;

disp("All>Rest")
disp(all_design)

% What are we deleting?
matlabbatch{3}.spm.stats.con.delete = 0;

spm_jobman('run', matlabbatch);

end