#!/bin/bash

# t threshold
t=${1} ; shift
out=${1} ; shift

files=${*}

n_files=$(echo "${files}" | wc -w | tr -d '[:space:]')

working_dir=$(mktemp -d /tmp/penetrance.XXXX)

echo "    ${working_dir}"
echo "    ${n_files} files"

for f in ${files} ; do 

    bn=$(basename "${f}")

    fslmaths                    \
        "${f}"                  \
        -thr    "${t}"          \
        -bin                    \
        "${working_dir}/${bn}"

done

fslmerge -t "${working_dir}/merged.nii.gz" "${working_dir}"/*.nii.gz

fslmaths \
    "${working_dir}/merged.nii.gz"                  \
    -Tmean                                          \
    "${working_dir}/penetrance_proportion.nii.gz"


fslmaths \
    "${working_dir}/penetrance_proportion.nii.gz"    \
    -mul "${n_files}"                                \
    "${working_dir}/penetrance_count.nii.gz"



cp "${working_dir}/penetrance_count.nii.gz"         "${out}_stat-n.nii.gz"
cp "${working_dir}/penetrance_proportion.nii.gz"    "${out}_stat-pct.nii.gz"

echo "    Done!"
echo 