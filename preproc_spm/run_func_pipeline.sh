#!/bin/bash

bids=${1}
dest=${2}
sub=${3}

# This path must be set to be where the files named as exactly 
#   <function_name>.m are saved
export MATLABPATH=~/code/preproc_spm/
export OMP_NUM_THREADS=1

# First step is to unzip funcs for SPM

# For some reason MacOS cp is missing `-u`
# This should really copy to $dest/sub-$sub/func/ but oops 
rsync --update \
    "${bids}/sub-${sub}/func"/*_bold.nii.gz \
    "${dest}/sub-${sub}/" \
    "${dest}/sub-${sub}/func/"

# Copy inherited events files
bolds=$(find "${bids}/sub-${sub}/func/" -name "*_bold.nii.gz" \
                -exec basename {} \; | \
            sed 's/_bold.nii.gz//' | \
            tr '\n' ' ')

for bold in ${bolds} ; do 

    f=${bold}_events.tsv
    cp "${bids}/sub-${sub}/func/${f}" "${dest}/sub-${sub}/func/sub-${sub}_${f}"

done

for f in "${dest}/sub-${sub}/func"/*_bold.nii.gz ; do 

    func_dir=${dest}/sub-${sub}/func/$(basename "${f}" .nii.gz)
    mkdir -p "${func_dir}/01_orig/"

    prefix=${func_dir}/01_orig/$(basename "${f}" _bold.nii.gz)_
    # echo "    prefix: ${prefix}"

    if [ ! -e "${prefix}0000.nii" ] ; then 

        # If the first file doesn't exist, but the NIFTI-2 version does, just
        # gunzip it; else fslsplit from scratch
        if [ -e "${prefix}0000.nii.gz" ] ; then 
            gunzip "${prefix}"*.gz
            mv "${prefix}"*.nii "${func_dir}/01_orig/"
        else 
            echo "Separating $(basename "${f}")"
            fslsplit "${f}" "${prefix}"
            gunzip "${func_dir}"/01_orig/*.gz
            # mv ${prefix}*.nii.gz ${func_dir}/01_orig/
            # gunzip ${func_dir}/01_orig/*.gz
        fi
    fi 

done

echo "Done setting up for SPM, invoking MATLAB ..."
echo

# Loop over each directory m
for func_dir in "${dest}/sub-${sub}/func"/*_bold ; do 

    echo -e "\n${func_dir}"
    fdir=$(basename "${func_dir}")

    mkdir -p "${func_dir}"/{02_moco,03_rmoco,04_warped,05_warpsmooth}/

    # 01: orig; 02: motion corrected; 03: registered and motion-corrected
    n_sourcefiles=$(find "${func_dir}/01_orig" -name "sub-*.nii" | wc -l)
    n_rmocofiles=$(find "${func_dir}/03_rmoco" -name "rmoco_sub-*.nii" | wc -l)

    if [ "${n_sourcefiles}" -eq 0 ] ; then 

        echo "No sourcefiles found in ${func_dir}/01_orig/, skipping this one."
        continue 

    fi

    echo "${n_sourcefiles} ${n_rmocofiles}"

    if [[ "${n_sourcefiles}" != "${n_rmocofiles}" ]] ; then 

        matlab \
            -batch "func_coreg('${dest}', '${sub}', '${fdir}')"

    else

        echo "Motion correction already done"

    fi

    # 'n' for smoothing value
    n_snwrmocofiles=$(find "${func_dir}/05_warpsmooth" \
                        -name "s[0-9]wrmoco_sub-*.nii" | wc -l)

    if [[ "${n_sourcefiles}" != "${n_snwrmocofiles}" ]] ; then 

        matlab \
            -batch "func_warpsmooth('${dest}', '${sub}', '${fdir}', \
                                    'SPMwarp_lessRegularized', '6')"

    else

        echo "Warp and smoothing already done"

    fi

    # break

done