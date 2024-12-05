#!/bin/bash

# Add MATLAB functions to path
# home=~/Projects/LBS/ferrara_lbs
home=~/code/LItoolbox/
export MATLABPATH=${home}

# Mask
# roi=${home}/data/BA-073940_dil-3.nii

if [ ${#} -ge 3 ] ; then

    roi=${1}
    out=${2} 
    thr=${3}
    shift 3
    # skip="0"

elif [ ${#} -eq 2 ] ; then 

    echo "Skipping to result collocation"
    temp_dir=${1}
    out=${2}
    # skip="1"

else 

    echo "Usage: ${0} <ROI> <out.tsv> <thr> t.nii[, t.nii ..]"
    exit 

fi

# Check arguments

if [[ ${out} =~ .*.nii ]] ; then 
    echo "Error: Output file '${out}' contains '.nii'; check and re-run"
    exit 1
fi

if [[ ${thr} =~ .*.nii ]] ; then 
    echo "Error: Threshold '${thr}' contains '.nii'; check and re-run"
    exit 1
fi

# Function to call LI on a single file
calculate_LI_on_file () {

    input_Tstat=${1}
    roi=${2}
    thr=${3}
    output_dir=${4}

    # Create a file with the given name in the temp directory 
    # shellcheck disable=SC2001
    fname=$(echo "${input_Tstat}" | sed 's|/|_|g').tsv
    output_tsv=${output_dir}/${fname}

    echo "${output_tsv}"

    echo "${input_Tstat}"

    matlab -batch "call_LI('${input_Tstat}', '${roi}', '${thr}', \
                           '${output_tsv}')"

}

export -f calculate_LI_on_file 

# shellcheck source=/opt/homebrew/bin/env_parallel.bash
. "$(which env_parallel.bash)"

# =============================================================================

# This creates a temporary directory to save the output so you don't have to 
# recalculate it later, but can be deleted to tidy up.

temp_dir=temp_${out//.tsv/}
mkdir -p "${temp_dir}/"
# all_LI_files_dir=$(mktemp -d)
n_files=$(echo "${@}" | wc -w)

echo "Working in ${temp_dir}"
echo "I was given ${n_files} stat files to operate on"

# Run the calculation in parallel; this will put n LI files into this temp dir
# for concatenation
files=${*} # $(echo ${*} | tr ' ' '\n')
# shellcheck disable=SC2086
rsync --update ${files} ${temp_dir}/

echo 

gz_files=$(find "${temp_dir}/" -name "*.nii.gz" | wc -l)

if [ "${gz_files}" -gt 0 ] ; then

    for file in "${temp_dir}"/*.nii.gz ; do 

        # echo $file
        fslmaths "${file}" -mul 1 "${file}" -odt float 
        yes n | gunzip "${file}"

    done

fi

echo "Starting parallel run ..."
new_files=$(find "${temp_dir}"/ -name "*.nii")
env_parallel -j 10 --env calculate_LI_on_file \
    calculate_LI_on_file {} "${roi}" "${thr}" "${temp_dir}" ::: "${new_files}"

result_files=($(find "${temp_dir}" -name "*.tsv"))
n_result_files=${#result_files[@]}
# echo "${result_files[@]}"

# Get the header from the first file
echo 
echo "Putting results in ${out} (${n_result_files})"
cat "${result_files[0]}" > "${out}"

max=$(( n_result_files - 1 ))

for i in $(seq 1 ${max} ) ; do 

    # echo "${result_files[${i}]}"
    tail -n +2 "${result_files[${i}]}" >> "${out}"

done 

# Delete generated files (I only care about the LI values rn)
rm LI_{boot,masking}.ps LI_r_*.nii