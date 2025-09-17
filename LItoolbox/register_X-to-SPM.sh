#!/bin/bash

if [ ${$} -ge 2 ] ; then

    src=${1}
    shift
    files=${*}

else

    echo "Usage: ${0} template moving [moving, ..]"
    exit 1

fi

li_home=~/code/LIToolbox/

# SPM destination space
spm_home=~/Documents/MATLAB/spm12/
spm_template=${spm_home}/toolbox/wfu_pickatlas/MNI_atlas_templates/MNI_T1.nii

# Transforms
xforms=~/code/templates/

# fslroi "${src}" src_refvol.nii.gz 0 1

srcfile=${xforms}/${src}/${src}_desc-toSPMnlin_coefs.nii.gz

if [ ! -e "${srcfile}" ] ; then
    echo "ERROR: ${srcfile} doesn't exist!"
    exit 1
fi

# Nonlinear warp

for f in ${files} ; do

    newname=${f//_contrast-/_space-SPM_contrast-}

    if [ "${newname}" -ot "${f}" ] ; then

        applywarp \
            -i "${f}"                \
            -o "${newname}"          \
            -r ${spm_template}       \
            -w "${srcfile}"

    else

        echo "Not re-warping ${f}, new space file already exists"

    fi

done

# xfm_name=${src}_to-SPM

# if [ ! -e "${xfm_name}.mat" ] ; then

#     echo "Starting non-linear transform ..."

#     fnirt   \
#         --in="${srcfile}"          \
#         --ref="${spm_template}"

# fi