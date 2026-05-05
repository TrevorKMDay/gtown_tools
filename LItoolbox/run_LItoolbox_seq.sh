#!/bin/bash

home=~/code/LItoolbox/
export MATLABPATH=${home}

#defaults
force="false"
thresh="unset"

while getopts "h:f:s:t:" flag; do
    case ${flag} in
        # Make sure to shift 2 when OPTARG is assigned
        h) home=${OPTARG} ; shift 2 ;;
        s) spm=${OPTARG}  ; shift 2 ;;
        f) force="true"   ; shift ;;
        # This is the same as 'thr3'
        # This could be changed to allow multiple, but not right now
        t) thresh=${OPTARG} ; shift 2 ;;
        \?)
        # Handle invalid options
        ;;
    esac
done

roi=${1} ; shift
thr1=${1} ; shift
output_tsv=${1} ; shift
files=${*}

if [ ! -d "$(dirname "${output_tsv}")" ] ; then

    echo -n "Requested output driectory $(dirname "${output_tsv}") does not"
    echo    "exist!"
    exit 1

fi

# spm=~/Documents/MATLAB/spm12/

# First, do some checks

if [ "${thr1}" == 1 ] ; then
    echo "Using fixed thresholding."

    if [ "${thresh}" == "unset" ] ; then
        echo "... but the threshold (-t) wasn't set, dying"
        exit 1
    else
        echo "... using threshold: ${thresh}"
        thr1="1:${thresh}"
    fi

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

echo "${spm}"

# Count the number of files only when we know we can start

# Count the number of files
n_files=$(echo "${files}" | wc -w | tr -d '[:space:]')
echo "I was given ${n_files} files"

# Can't do anything with no files
if [ "${n_files}" -eq 0 ] ; then
    echo "Usage: ${0} roi thr1 output.tsv f1[, f2 ...]"
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
