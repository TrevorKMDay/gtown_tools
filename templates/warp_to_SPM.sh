#!/bin/bash

src_space=${1}
file=${2}
oname=${3}

home=~/code/templates/

if [ ! -d "${home}/${src_space}" ] ; then
    echo "No transforms directory for ${src_space}!"
    exit 1
fi

ref=${home}/${src_space}/${src_space}_desc-toSPMnlin_T1w.nii.gz
warp=${ref//T1w/coefs}


# https://fsl.fmrib.ox.ac.uk/fsl/docs/#/registration/fnirt/user_guide?id=applywarp

applywarp \
    --ref="${ref}"      \
    --in="${file}"      \
    --out="${oname}"    \
    --warp="${warp}"