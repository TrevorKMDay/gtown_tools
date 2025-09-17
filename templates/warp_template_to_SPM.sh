#!/bin/bash

# Set up SPM references ============

spm_home=~/Documents/MATLAB/spm12/
spm_template=${spm_home}/toolbox/wfu_pickatlas/MNI_atlas_templates/MNI_T1.nii

spm_ref=SPM_MNI_T1w.nii.gz

# flirt works better with extracted brains, and fnirt with whole brains
# See: https://fsl.fmrib.ox.ac.uk/fsl/docs/#/registration/fnirt/user_guide


# Copy SPM template locally
if [ ! -e ${spm_ref} ] ; then

    echo "Copying SPM files ..."
    rsync -u ${spm_template} SPM_MNI_T1w.nii
    gzip SPM_MNI_T1w.nii

fi

# Create SPM brain mask for registration
if [ ! -e SPM_MNI_desc-brain_mask.nii.gz ] ; then

    echo "Extracting SPM brain mask"
    bet2 ${spm_ref} SPM_MNI_desc-brain --mask --nooutput

fi

if [ ! -e SPM_MNI_desc-brain_T1w.nii.gz ] ; then

    echo "Masking SPM brin"
    fslmaths \
        SPM_MNI_T1w.nii.gz \
        -mul SPM_MNI_desc-brain_mask.nii.gz \
        SPM_MNI_desc-brain_T1w.nii.gz

fi

# Space ======

if [ ${#} -eq 1 ] ; then
    space_name=${1}
else
    echo "Supply space label, starting with 'tpl-'"
    exit 1
fi

# Get relevant template from TemplateFlow
template_name=$(echo "${space_name}" | grep -o 'tpl-[^_]*')

# cd templateflow/ || exit
# datalad get -r "${template_name}"

space_dir="${space_name}/"
mkdir -p "${space_dir}"

rsync -Lu \
    templateflow/"${template_name}"/"${space_name}"_{T1w,desc-brain_mask}.nii.gz \
    "${space_dir}"

if [ ! -e "${space_dir}/${space_name}_desc-brain_T1w.nii.gz" ] ; then

    echo "Applying brain mask to ${space_name} T1w"
    fslmaths \
        "${space_dir}/${space_name}_T1w.nii.gz"                     \
        -mul "${space_dir}/${space_name}_desc-brain_mask.nii.gz"    \
        "${space_dir}/${space_name}_desc-brain_T1w.nii.gz"

fi

# Register the MNI template to the SPM template
to_spm_linear="${space_dir}/${space_name}_desc-toSPMlin_T1w.mat"
if [ ! -e "${to_spm_linear}" ] ; then

    echo "Running flirt ..."
    flirt \
        -in     "${space_dir}/${space_name}_desc-brain_T1w.nii.gz"      \
        -ref    SPM_MNI_desc-brain_T1w.nii.gz                           \
        -out    "${space_dir}/${space_name}_desc-toSPMlin_T1w.nii.gz"   \
        -omat   "${to_spm_linear}"

fi

# Apply the xform to the brian mask
brainmask="${space_dir}/${space_name}_desc-braintoSPMlin_mask.nii.gz"
if [ ! -e "${brainmask}" ] ; then

    # Moving mask to SPM space
    echo "Applying warp to brain mask ..."
    flirt \
        -in     "${space_dir}/${space_name}_desc-brain_mask.nii.gz"         \
        -ref    ${spm_ref}                                                  \
        -applyxfm -init "${to_spm_linear}"                                  \
        -interp nearestneighbour                                            \
        -out    "${brainmask}"

fi

# Do the nonlinear warp using brainmasks
echo "Starting fnirt ..."
fnirt \
    --ref="${spm_ref}"                                                      \
    --refmask="SPM_MNI_desc-brain_mask.nii.gz"                              \
    --in="${space_dir}/${space_name}_T1w.nii.gz"                            \
    --inmask="${space_dir}/${space_name}_desc-brain_mask.nii.gz"            \
    --aff="${to_spm_linear}"                                                \
    --cout="${space_dir}/${space_name}_desc-toSPMnlin_coefs"                \
    --iout="${space_dir}/${space_name}_desc-toSPMnlin_T1w"                  \
    --logout="${space_dir}/${space_name}_toSPMnlin.log"

    # --fout="${space_dir}/${space_name}_desc-field"