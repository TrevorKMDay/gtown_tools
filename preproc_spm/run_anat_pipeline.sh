#!/bin/bash

if [ ${#} -eq 3 ] ; then

    bids=${1}
    dest=${2}
    sub=${3}

fi

T1s=$(find ${bids}/sub-${sub}/anat -name "*_T1w.nii.gz")
n_T1s=$(echo ${T1s} | wc -w)

echo "Found ${n_T1s} T1s!"


if [ ${n_T1s} -eq 0 ] ; then

    echo "No T1s to work with ... dying ..."
    exit

elif [ ${n_T1s} -eq 1 ] ; then

    mkdir -p ${dest}/sub-${sub}/anat/

    echo ${T1s//.gz/}
    if [ ! -e ${T1s//.gz/} ] ; then
        rsync -u ${T1s} ${dest}/sub-${sub}/anat/
        find ${dest}/sub-${sub}/anat/ -name "*.gz" -exec gunzip -f {} \;
    fi

    # Set this to exactlly this name
    T1="${dest}/sub-${sub}/anat/sub-${sub}_run-1_T1w.nii"

else

    echo "I can't handle multiple T1s right now"
    exit 1

fi

# This path must be set to be where the files named as exactly
#   <function_name>.m are saved
export MATLABPATH=~/code/preproc_spm/
export OMP_NUM_THREADS=1

matlab -batch "anat_segnorm('${T1}')"
