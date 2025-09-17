#!/bin/bash

input=${1}
task=${2}

# Setup output directory
i=$(basename "${input}")
for_li=${i}_forLI/
mkdir -p "${for_li}"/

# Set up directory

echo "Setting up files ..."

files=$(find "${input}" \
            -name "sub-*_task-${task}_contrast-*_stat-t_statmap.nii.gz" \
            -exec readlink -f {} \;)
n_files=$(echo "${files}" | wc -w)

echo "Found ${n_files} files"

for f in ${files} ; do

    # Don't update links each time we run so timestamps are tracked

    target="${for_li}/$(basename "${f}")"

    if [ "${target}" -ot "${f}" ] ; then
        ln -s "${f}" "${target}"
    fi

done

# Exclude files already in SPM space
input_files=$(find "${for_li}/" \
                    -name "sub-*_task-${task}_contrast-*.nii.gz" | \
                sort)

#shellcheck disable=SC2086
~/code/LItoolbox/register_X-to-SPM.sh   \
    "tpl-MNI152NLin2009cAsym_res-02"    \
    ${input_files}

spm_files=$(find "${for_li}/" \
                    -name "sub-*_task-${task}_space-SPM_contrast-*.nii.gz" | \
                sort)

for f in ${spm_files} ; do

    ungzed=${f//.gz/}

    if [ ! -e "${ungzed}" ] || [ "${ungzed}" -ot "${f}" ] ; then
        gunzip -k "${f}"
    fi

done

echo "Done creating unGZed NIFTIs in SPM space"