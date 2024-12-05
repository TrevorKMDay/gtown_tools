#!/bin/bash

spm=~/code/spm/
export MATLABPATH=${spm}

if [ ${#} -ge 4 ] ; then 

    result_dir=${1}
    task=${2}
    con=${3}
    shift 2 
    files=${*}

else 

    echo "Usage: ${0} result_dir task-<task> con-<con> files... "
    exit 1

fi 

echo "${task} ${con}"
echo "======================"
echo "    dir:        ${result_dir}"

# Remove extra whitespace from sed output
n_files=$(echo "${files}" | wc -w | sed 's/[^0-9]//g')
echo "    file count: ${n_files}"

if [ "${n_files}" -eq  0 ] ; then 

    echo "Given 0 files! Quitting."
    exit 1

fi

con_dir="${result_dir}/${task}/${con}"
mkdir -p "${con_dir}"

# Run the contrast
matlab -batch "func_group('${con_dir}', '${files}')"

# Rename files
pfx=${result_dir}/${task}/${task}_${con}

# I was symlinking, but then fsleyes only loads the original name, and not
#   the helpful new labeldu 

cp "$(readlink -f "${con_dir}/spmT_0001.nii")" "${pfx}_stat-t_statmap.nii"

cp "$(readlink -f "${con_dir}/beta_0001.nii")" \
    "${pfx}_stat-beta_statmap.nii"

cp "$(readlink -f "${con_dir}/con_0001.nii")" \
    "${pfx}_stat-effect_statmap.nii"

# Gzip final files
gzip -f "${pfx}"_*.nii

# Fix NaNs
for i in "${pfx}"_*.nii.gz ; do 

    fslmaths "${i}" -nan "${i}"

done