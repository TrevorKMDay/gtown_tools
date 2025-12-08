% Applies standard SPM segmentation to an t1, creates a smoothed 
% brain mask from the c1 and c2 tissue, applies it to the t1, and also 
% warps the result into standard space. 
% This assumes that data have been read into .nii format already. Note that
% pre-20221027 structural files are read in as files starting with s, and
% post-20221027 files as files starting with MF.

% SETUP

function anat_segnorm = anat_segnorm(src_t1) 

    maxNumCompThreads(1) ;
    
    addpath('~/Documents/MATLAB/spm12') ;
    
    % % home = '/Volumes/TrevorDay/Projects/LBS/ferrara_lbs';
    % % bids = strcat(home, '/bids/') ;
    % 
    % %% prep the output folder
    % outputdir = strcat(output, "/sub-", sub, "/anat/") ;
    % mkdir(outputdir);
    % 
    % % base T1 for processing
    % src_t1 = strcat(outputdir, 'sub-', sub, '_T1w.nii') ;
    % 
    % if ~isfile(src_t1)
    % 
    %     % Read in the t1
    %     t1_file = strcat(bids, '/sub-', sub, '/anat/', 'sub-', sub, ...
    %         '_T1w.nii.gz') ;
    % 
    %     rawt1 = niftiread(t1_file);
    %     rawt1_hdr = niftiinfo(t1_file);
    % 
    % 
    %     % write a copy of this t1 to the new folder
    %     % QUESTION: Why not compress?
    %     niftiwrite(rawt1, src_t1, rawt1_hdr, 'Compressed', false);
    % 
    % else
    % 
    % 
    %     disp("Skip T1 copy")
    % 
    % end

    % disp(src_t1)
    [outputdir, basename, ext] = fileparts(src_t1) ;
    sub = erase(regexp(basename, "sub-[^_.]*", 'match'), 'sub-')

    % Check that input is uncompressed nifti
    if ext ~= ".nii"
        disp('Input T1 must have .nii extension')
        return
    end

    % Change tehse so they can be red properly
    outputdir = strcat(outputdir, "/") ;
    basename = strcat(basename, ext) ;

    %% START PROCESSING

    % If tissue seg outputs don't exist, create them 
    if ~isfile(strcat(outputdir, "c6", basename))
    
        disp('Starting the seg/norm estimation in SPM12...')
        clear matlabbatch;
        spm('defaults', 'fmri');
        spm_jobman('initcfg');
        
        %% M-code from SPM 
        % (replace file/folder references with variables as necessary)
        
        % TPM = tissue probability map, built into SPM, based on MNI152 
        % template
        tpmfile = fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii'); 
        
        %% Segmentation of cost-function masked t1
        matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(src_t1); 
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
        
        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = ...
            {[tpmfile ',1']};
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
        % writes the native space (c*) and the DARTEL-importable (rc*) 
        % version
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1]; 
        % writes modulated (mwc*) and unmodulated (wc*) warped to standard 
        % space versions
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 1]; 
        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = ... 
            {[tpmfile ',2']};
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 1];
        
        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = ...
            {[tpmfile ',3']};
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 1];
        
        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = ...
            {[tpmfile ',4']};
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [1 1];
        
        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = ...
            {[tpmfile ',5']};
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [1 1];
        
        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = ...
            {[tpmfile ',6']};
        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [1 1];
        
        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        % Default is [0 0.001 0.5 0.05 0.2], but that often "cuts" 
        % ventricles in stroke brains. The *0.1 recommendation comes from 
        % John Ashburner, who developed the algorithm, so I'm not sure why 
        % it isn't the default anyway. See here:
        % https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;c64dcab1.1802
        matlabbatch{1}.spm.spatial.preproc.warp.reg = ... 
            [0 0.0001 0.05 0.005 0.02];  
        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
        
        %% Run the segmentation batch...
        spm_jobman('run', matlabbatch);

    else

        disp("Skipping seg/norm")

    end

    %% Produce warped outputs for 
    stripped_t1_file = fullfile(outputdir, ...
        strcat('sub-', sub, '_stripped_t1.nii')) ;
    smoothed_brainmask_file = fullfile(outputdir, ...
        strcat('sub-', sub, '_brain_mask_blur2.nii')) ;

    c1t1_file = fullfile(outputdir, strcat('c1', basename));
    c2t1_file = fullfile(outputdir, strcat('c2', basename));
    c3t1_file = fullfile(outputdir, strcat('c3', basename));
    c4t1_file = fullfile(outputdir, strcat('c4', basename));
    c5t1_file = fullfile(outputdir, strcat('c5', basename));
    c6t1_file = fullfile(outputdir, strcat('c6', basename));

    if ~isfile(stripped_t1_file) || ~isfile(smoothed_brainmask_file)
    
        %% OK, now create a nice brain mask from WM and GM in native space
        % combines c1 and c2, then blurs it slightly, and outputs a 
        % brain mask

        disp('Computing a brain mask.')
        
        %[hdr,img]=read_nifti(t1_file); % read in anatomical
        hdr = niftiinfo(src_t1);
        img = niftiread(src_t1);

        c1img=niftiread(c1t1_file);
        c2img=niftiread(c2t1_file);
        
        in_brain_mask = any(cat(4, c1img, c2img), 4);
        
        % now blur a bit
        sigma = 2;
        in_brain_mask = imgaussfilt3(100*in_brain_mask, sigma);
        
        % binarize
        % mask at 20... so .2 the spread...
        in_brain_mask = in_brain_mask > 40; 
        % zero out outside of mask
        img(~in_brain_mask) = 0; 
       
        % AG: Note that this does NOT include the ventricles but DOES include 
        % the lesion (as it is replaced with the standard tpm in the c1 and c2 
        % images)
        niftiwrite(feval(hdr.Datatype, in_brain_mask), ...
            smoothed_brainmask_file, hdr) 
        
        niftiwrite(feval(hdr.Datatype, img), stripped_t1_file, hdr) 
    
    else

        disp("Brain mask already ready")

    end
    
    % using this deformation field
    deffield_file = fullfile(outputdir, strcat('y_', basename)); 
    
    % Warp c1 (gm) and c2 (wm)
    for i = [src_t1 stripped_t1_file smoothed_brainmask_file ...
                c1t1_file c2t1_file] 

        % I had to move this into a loop; not sure why, I couldn't figure
        % out how to put multiple images into subj.resample

        [i_dir, i_fn, i_sfx] = fileparts(i) ;
        warped_i = strcat(i_dir, "/w", i_fn, i_sfx);
        % disp(warped_i)

        if ~isfile(warped_i)

            clear matlabbatch;
            spm('defaults', 'fmri');
            spm_jobman('initcfg');
    
            matlabbatch{1}.spm.spatial.normalise.write.subj.def = ...
                cellstr(deffield_file);
            
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample = ...
                cellstr(i);
        
                % { [src_t1, ',1']; 
                %   [stripped_t1_file, ',1'] ; 
                %   [smoothed_brainmask_file, ',1'] }
            
            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = ...
                [-78 -112 -70; 78 76 85];
            
            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = ...
                [1 1 1];
            matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = ...
                4;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = ...
                'w';
            
            spm_jobman('run', matlabbatch);

        else 

            disp(strcat("File ", i, " already exists"))

        end

    end

end