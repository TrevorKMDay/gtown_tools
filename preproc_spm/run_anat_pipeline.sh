#!/bin/bash

if [ ${#} -eq 3 ] ; then

    bids=${1}
    dest=${2}
    sub=${3}

fi

T1s=$(find "${bids}/sub-${sub}/anat" -name "*_T1w.nii.gz")
n_T1s=$(echo "${T1s}" | wc -w)

echo "Participant: ${sub}"
echo "==================="

echo "    Found ${n_T1s} T1s!"

if [ "${n_T1s}" -eq 0 ] ; then

    echo "    INFO: No T1s to work with ... dying ..."
    exit

elif [ "${n_T1s}" -eq 1 ] ; then

    mkdir -p "${dest}/sub-${sub}/anat/"
    src_T1=$(find "${bids}/sub-${sub}/anat/" -name "*_T1w.nii.gz")

else

    t1_to_use="$(cat "${bids}/sub-${sub}/anat/USE")"
    src_T1=${bids}/sub-${sub}/anat/${t1_to_use}

    if [[ "${t1_to_use}" == 'NONE' ]] ; then
        echo "    INFO: Multiple T1s, none good enough ... dying ..."
        exit
    else
        echo "    INFO: Multiple T1s, using ${src_T1}"
        src_T1=${bids}/sub-${sub}/anat/${t1_to_use}
    fi

    mkdir -p "${dest}/sub-${sub}/anat/"

fi

dest_T1=${dest}/sub-${sub}/anat/$(basename "${src_T1}" .gz)
if [ ! -e "${dest_T1}" ] ; then
    rsync -u "${src_T1}" "${dest}/sub-${sub}/anat/"
    find "${dest}/sub-${sub}/anat/" -name "*.gz" -exec gunzip -f {} \;
fi

# Set this to exactly this name
dest_T1=$(find "${dest}/sub-${sub}/anat/" -name "*_T1w.nii")

# exit

# This path must be set to be where the files named as exactly
#   <function_name>.m are saved
export MATLABPATH=~/code/preproc_spm/
export OMP_NUM_THREADS=1

matlab -batch "anat_segnorm('${dest_T1}')"
