
function func_group = func_group(outfolder, contrast_files)

% % User input - change as needed.
% parent_datadir = '/Volumes/POLISBACKUP/KatrinaLineReanalysis2024/MRI_SPM/Subs/'; 
% subIDs = {'TDCh_121','TDCh_123','TDCh_125','TDCh_133','TDCh_134','TDCh_136',...
%     'TDCh_139','TDCh_140','TDCh_142','TDCh_143','TDCh_144','TDCh_145',...
%     'TDCh_150','TDCh_159','TDCh_163','TDCh_165','TDCh_166','TDCh_167',...
%     'TDCh_169','TDCh_170','TDCh_171','TDCh_172','TDCh_173','TDCh_174','TDCh_176','TDCh_177',...
%     'TDCh_178','TDCh_181','TDCh_183','TDCh_184','TDCh_186','TDCh_187','TDCh_188','TDCh_189',...
%     'TDCh_190','TDCh_191'}; % Only the subject you want to include
% resultsfoldername = '/processed/Results/SPMwarp_lessRegularized_s6_acrossRuns_OCTnoMoco/'; % plain, or with _OCT or _OCTnoMoco
% contrastname = 'con_0002.nii';
% outfolder = '/Volumes/POLISBACKUP/KatrinaLineReanalysis2024/MRI_SPM/GroupMaps/KatrinaSubs_ColorVsRest_OCTnoMoco'; % Description needs to match contrast number above
% % corresponding to two lines above, end plain, or with _OCT or _OCTnoMoco
% 
% confiles = cell(length(subIDs),1);
% for subnr = 1:length(subIDs)
%     confiles{subnr,1} = [parent_datadir, subIDs{subnr}, resultsfoldername, contrastname];
% end

confiles = split(contrast_files, " " );

%% Set up SPM

addpath('~/Documents/MATLAB/spm12') ;

clear matlabbatch

spm('defaults', 'FMRI');
spm_jobman('initcfg');

% M-code from SPM (replace file/folder references with variables as necessary)

matlabbatch{1}.spm.stats.factorial_design.dir = {outfolder};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = confiles;
matlabbatch{1}.spm.stats.factorial_design.cov = ...
    struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = ...
    struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1; 

% Model estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, ...
                '.','val', '{}',{1}), ...
                substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Contrast definition
matlabbatch{3}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.', 'val', '{}', {2}, '.',' val', '{}', {1}, ...
                '.','val', '{}', {1}), ...
                substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'TheContrast';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

spm_jobman('run', matlabbatch);

end