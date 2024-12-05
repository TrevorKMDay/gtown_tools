#!/bin/bash


# Instructions:
#   (1) Cluster correction using SPM output
#   https://www.andysbrainblog.com/andysbrainblog/2017/6/21/using-afnis-3dfwhmx-with-spms-residual-data
#   (2) Actually running cluster correction
#   https://andysbrainbook.readthedocs.io/en/stable/fMRI_Short_Course/fMRI_Appendices/Appendix_A_ClusterCorrection.html


spm_dir=${1}
NN=${2}
side=${3}
thresh=${4}

mask=${spm_dir}/mask.nii

# Square the residual image
resms_sq=${spm_dir}/ResMSsqrt.nii
if [ ! -e "${resms_sq}" ] ; then 
    fslmaths "${spm_dir}/ResMS.nii" -nan -sqrt "${resms_sq}"
fi

# Cluster estimation

fwhmx="${spm_dir}/3dFWHMx.txt"
if [ ! -e "${fwhmx}" ] ; then 

    3dFWHMx -acf -mask "${mask}" "${resms_sq}" > "${fwhmx}"
    # mv 3dFWHMx.1D{,.png} "${spm_dir}"

fi

if [ ! -e "${spm_dir}/cluster.NN1_bisided.1D" ] ; then 

    # remove last value (see ref 1)
    values=$(head -n2 "${fwhmx}" | tail -n1 | sed 's/ [0-9]*[.][0-9]*$//')
    echo "Using values: ${values}"

    # shellcheck disable=SC2086
    3dClustSim \
        -nodec                      \
        -mask   "${mask}"           \
        -acf    ${values}           \
        -athr   0.05                \
        -pthr   0.001               \
        -prefix ${spm_dir}/cluster

fi

nvox=$(tail -n1 "${spm_dir}/cluster.NN${NN}_${side}sided.1D" | \
        tr -s ' ' ' ' | \
        awk '{print $2}')

# Make sure not to get any of the other files
spmT_srcfiles=$(find "${spm_dir}/" -regex ".*/spmT_[0-9]*[.]nii")

for spmT in ${spmT_srcfiles}; do 

    # shellcheck disable=SC2001
    corrected=$(echo "${spmT}" | sed "s/.nii/_cc${thresh}.nii/")

    if [ ! -e "${corrected}" ] ; then 

        3dClusterize                                \
            -inset          "${spmT}"               \
            -ithr           0                       \
            -idat           0                       \
            -mask           "${mask}"               \
            -NN             "${NN}"                 \
            -"${side}sided" RIGHT "${thresh}"       \
            -clust_nvox     "${nvox}"               \
            -pref_dat       "${corrected}"

    fi

done