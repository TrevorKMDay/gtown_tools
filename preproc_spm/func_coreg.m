% This function registers the functionals to native space

% EDIT ME TO TAKE PATH TO ONE FUNC DIR

% Parameter func: name of directory containing 01_orig/, 02_moco/, etc., e.g.:
%   sub-01_task-line_run-1/. Do not include full path or trailing slash

function func_coreg = func_coreg(home, sub, func)

% maxNumCompThreads(1) ;
addpath('~/Documents/MATLAB/spm12') ;

ref_T1w = strcat(home, "/sub-", sub, "/anat/sub-", sub, "_T1w.nii");

if ~exist(ref_T1w, 'file')
    disp('T1 not found. Aborting.')
else
    disp(strcat("T1:", ref_T1w))
end

% Identify the files in the given <sub>/func/<func>/01_orig directory:
processed = strcat(home, "/sub-", sub, "/");
func_dir = strcat(processed, "/func/", func);
nii_names = struct2table(dir(strcat(func_dir, "/01_orig/sub-*.nii"))).name;
nii_files = cellstr(strcat(func_dir, "/01_orig/", nii_names));

disp(nii_files) ;

if height(nii_files) == 0
    disp("No .nii files found!")
end

clear matlabbatch

% Prepare to specify batch from script
spm('defaults', 'fmri');
spm_jobman('initcfg');

% M-code from SPM (replace file/folder references with variables
% as necessary)

% Realign functionals (aka motion correction)
matlabbatch{1}.spm.spatial.realign.estwrite.data = {nii_files} ;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {''};
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = ...
    '../02_moco/moco_';

% Coregister with native-space MPRAGE
matlabbatch{2}.spm.spatial.coreg.estwrite.ref = cellstr(ref_T1w);
matlabbatch{2}.spm.spatial.coreg.estwrite.source(1) = ...
    cfg_dep(cellstr('Realign: Estimate & Reslice: Mean Image'), ...
    substruct('.', 'val', '{}', {1}, ...
    '.', 'val', '{}', {1}, ...
    '.', 'val', '{}', {1}, ...
    '.', 'val', '{}', {1}), ...
    substruct('.', 'rmean'));

matlabbatch{2}.spm.spatial.coreg.estwrite.other(1) = ...
    cfg_dep(cellstr('Realign: Estimate & Reslice: Resliced Images (Sess 1)'), ...
    substruct('.', 'val', '{}', {1}, ...
    '.', 'val', '{}', {1}, ...
    '.', 'val', '{}', {1}, ...
    '.', 'val', '{}', {1}), ...
    substruct('.', 'sess', '()', {1}, ...
    '.', 'rfiles'));
% defaults
matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.tol = ...
    [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];

% not defaults?
matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.prefix = ...
    '../03_rmoco/r';

%% Run the batch...
% disp(matlabbatch)
spm_jobman('run', matlabbatch);

end