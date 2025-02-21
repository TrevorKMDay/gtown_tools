#!/bin/bash

home=~/code/LItoolbox/
export MATLABPATH=${home}

roi=${1} ; shift
thr1=${1} ; shift
output_tsv=${1} ; shift
files=${*}

n_files=$(echo "${files}" | wc -w | tr -d '[:space:]')

echo "I was given ${n_files} files" 

if [ "${n_files}" -eq 0 ] ; then 
    echo "Usage: ${0} roi thr1 output.tsv f1[, f2 ...]"
    exit 1
fi

if [ -e "${output_tsv}" ] ; then 

    echo "Output file ${output_tsv} already exists, not overwriting!"
    exit 1

fi

# shellcheck disable=SC2086
matlab -batch "call_LI('${files}', '${roi}', '${thr1}', '${output_tsv}')"
