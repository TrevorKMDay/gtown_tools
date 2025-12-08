function func_1stlevel_postOCT = func_1stlevel_postOCT(dirs, TR, ...
    regs, output_dir)

% Trevor refactor:
addpath('~/Documents/MATLAB/spm12') ;

% dirs='first_level_fwhm6/sub-kLat003/func/sub-kLat003_task-ADDT_run-1_space-MNI_desc-smoothedmasked first_level_fwhm6/sub-kLat003/func/sub-kLat003_task-ADDT_run-2_space-MNI_desc-smoothedmasked'

disp("dirs:")
dirs_to_combine = strsplit(dirs, " ")

% High pass filter cutoff in seconds.
hpf_cutoff = 400;

% Takes a space-separated list of directories to combine, and pulls out
% the input data (05_warpsmooth) and OCT info.

TR = str2double(TR) ;

data = struct;
data.dir = dirs_to_combine;
n_runs = length(dirs_to_combine);

for i = 1:n_runs

    fdir = dirs_to_combine{i} ;
    disp(strcat("Working in ", fdir));

    data.nii_files{i} = dir(strcat(fdir, "/*.nii"));
    data.spm_mat{i} = load(strcat(fdir, ...
        "/results/optcens/01_censored/SPM.mat"));

    % disp(data.spm_mat{i}.SPM.xsDes)

    % Extract the design matrix
    data.des_mat{i} = data.spm_mat{i}.SPM.xX.X ;

    head(data.des_mat{i})

    %% Extract the nuisance regressors

    % Get number of conditions
    n_conditions = size(data.spm_mat{i}.SPM.Sess.U, 2) ;

    % Extract the motion columns
    % Columns are:
    %   regressors, 6 RP, censoring, constant
    %   In all cases, ignore the constant

    % Get the dimensions
    des_mat_size = size(data.des_mat{i}) ;
    dm_ncols = des_mat_size(2) ;
    dm_last_col = dm_ncols - 1 ;

    disp(strcat("Starting with size: ", num2str(des_mat_size(1)), ", ", ...
        num2str(des_mat_size(2))))

    if regs == "rp"
        % Only use the 6 motion parameters
        cols_to_keep = (n_conditions + 1):(n_conditions + 6);
    elseif regs == "rp+cens"
        cols_to_keep = (n_conditions + 1):dm_last_col ;
    elseif regs == "cens"
        % Skip over the 6 rp columns
        cols_to_keep = (n_conditions + 1 + 6):dm_last_col ;
    else
        disp(strcat("Illegal parameter: ", regs))
    end

    nuisance_regs = data.des_mat{i}(:,cols_to_keep);
    nregs_mat_size = size(nuisance_regs) ;

    disp(strcat("Nuisance reg size: ", num2str(nregs_mat_size(1)), ...
        ", ", num2str(nregs_mat_size(2))))

    data.nuisance_file{i} = strcat(fdir, ...
        "/results/optcens/01_censored/regsfile.txt");

    % Now create output text files
    save(data.nuisance_file{i}, 'nuisance_regs', '-ascii', '-double', ...
        '-tabs');

end

% return

% Create new output folder
mkdir(output_dir);

%% Set up run

% Prepare to specify batch from script
spm('defaults','fmri');
spm_jobman('initcfg');

% Model specification
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(output_dir);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

% We are going to assume the number of conditions is the same for both
% runs. If it's not, something's gone very wrong

for i = 1:n_runs

    scans = cellstr(strcat(dirs_to_combine{i}, "/", ...
        struct2table(data.nii_files{i}).name));

    matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans = scans;

    for j = 1:n_conditions

        details = data.spm_mat{i}.SPM.Sess.U(j);

        % The U data has the information from the previous run
        name = details.name{1} ;
        onsets = details.ons ;
        duration = unique(details.dur) ;

        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(j).name = name;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(j).onset = onsets;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(j).duration = ...
            duration ;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(j).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(j).pmod = ...
            struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(j).orth = 1;

    end

    % Get the nuisance file from the master data
    nuisance_file = data.nuisance_file{i};

    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress = ...
        struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = ...
        cellstr(nuisance_file);
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = hpf_cutoff;

end

% matlabbatch{1}.spm.stats.fmri_spec.timing

matlabbatch{1}.spm.stats.fmri_spec.fact = ...
    struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('fMRI model specification: SPM.mat File', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, ...
    '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

%% Contrasts
matlabbatch{3}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{2}, '.','val', '{}', {1}, ...
    '.','val', '{}',{1}), substruct('.','spmmat'));

% Contrasts info should be identical between the two runs, so use the
% first one
contrasts_info = data.spm_mat{1}.SPM.xCon;
n_contrasts = size(contrasts_info, 2);

for i = 1:n_contrasts

    % Remove junk from name
    name_raw = contrasts_info(i).name;
    name = regexprep(name_raw, " - .*", "");

    % Extract weights from struct based on the known number of
    % conditions
    w1 = contrasts_info(i).c(1:n_conditions);
    weights = transpose(w1) ;

    matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = name;
    matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = weights;
    matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'both';

end

matlabbatch{3}.spm.stats.con.delete = 0;

% return

% Run the batch...
spm_jobman('run', matlabbatch);

end