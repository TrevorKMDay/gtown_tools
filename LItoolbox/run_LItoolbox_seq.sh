#!/bin/bash

home=~/code/LItoolbox/
export MATLABPATH=${home}

force="false"
while getopts "f:" flag; do
 case $flag in
   f) force="true" ; shift ;;
   \?)
   # Handle invalid options
   ;;
 esac
done

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

echo "Starting value check ..."

if [[ ${force} == "false" ]] ; then

    for f in ${files} ; do

        for i in {,pix}dim{1,3} ; do

            srcval=$(fslval "${roi}" "${i}")
            fval=$(fslval "${f}" "${i}")

            if [[ "${srcval}" != "${fval}" ]] ; then
                echo "ERROR: ${i} values not equal:"
                echo "    ${roi}:"
                echo "        ${srcval}"
                echo "    ${f}:"
                echo "        ${fval}"
                exit 1
            fi

        done

        break

    done

fi

echo "Finished"

# shellcheck disable=SC2086
matlab -batch "call_LI('${files}', '${roi}', '${thr1}', '${output_tsv}')"
