#!/bin/bash

home=~/code/LItoolbox/
export MATLABPATH=${home}


force="false"
while getopts "h:f:s:" flag; do
    case ${flag} in
        h) home=${OPTARG} ; shift ;;
        s) spm=${OPTARG}  ; shift ;;
        f) force="true"   ; shift ;;
        \?)
        # Handle invalid options
        ;;
    esac
done

roi=${1} ; shift
thr1=${1} ; shift
output_tsv=${1} ; shift
files=${*}

# spm=~/Documents/MATLAB/spm12/

echo "${spm}"

# Count the number of files
n_files=$(echo "${files}" | wc -w | tr -d '[:space:]')
echo "I was given ${n_files} files"

# Can't do anything with no files
if [ "${n_files}" -eq 0 ] ; then
    echo "Usage: ${0} roi thr1 output.tsv f1[, f2 ...]"
    exit 1
fi

# Don't run LI if the output already exists
if [ -e "${output_tsv}" ] ; then

    echo "Output file ${output_tsv} already exists, not overwriting!"
    exit 1

fi

if [ ! -e "${roi}" ] ; then
    echo "Can't find ROI file: ${roi}"
    exit 1
fi

if [[ ${force} == "false" ]] ; then

    # Check that the ROI is the same dimensions as the input data
    echo "Starting value check ..."

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

    echo "Finished with value check"

else

    echo "Skipping value check"

fi

# Actually invoke MATLAB

# shellcheck disable=SC2086
matlab -batch "call_LI('${files}', '${roi}', '${thr1}', '${output_tsv}')"
