#!/bin/bash

# When given an fMRIPREP derivatives directory and input BIDS, this script
# sets it up for SPM task analysis

# Trevor Day, Oct. 3, 2024

home=$(dirname "${0}")

smooth_sigma=''
smooth_fwhm=''
while getopts s:f: flag ; do
    case "${flag}" in
        # Accept either smoothing sigma or FWHM and convert FWHM to sigma below
        s) smooth_sigma=${OPTARG} ; shift 2 ;;
        f) smooth_fwhm=${OPTARG} ; shift 2 ;;
        *) echo "Bad flag ${flag}!" ; exit 1 ;;
    esac
done

# Get smoothing parameters
if [ "${smooth_sigma}" != "" ] && [ "${smooth_fwhm}" != "" ] ; then
    echo "Both smoothing sigma and FWHM can't be set; pick one."
    exit 1
elif [ "${smooth_sigma}" == "" ] && [ "${smooth_fwhm}" != "" ] ; then 
    # Convert FWHM (mm) to sigma; round to two decimal places
    smooth_sigma=$(printf %.2f "$(echo "${smooth_fwhm} / 2.354 "| bc -l)")
    SMOOTHING="${smooth_sigma}"
elif [ "${smooth_sigma}" == "" ] && [ "${smooth_fwhm}" == "" ] ; then 
    SMOOTHING="false"
else 
    SMOOTHING="${smooth_sigma}"
fi

if [ ${#} -eq 5 ] ; then 

    src=${1}
    bids=${2}
    # Include "_cohort-?" if necessary
    space=${3}
    space_new=${4}
    dest=${5}

else 

    echo "Usage: ${0} <fmriprep> <bids> <space> <new_space_label> <dest>"
    exit 1

fi

echo "Derivatives: ${src}"
echo "BIDS:        ${bids}"
echo "Space:       ${space} ==> ${space_new}"
echo "Destination: ${dest}"
echo "Smoothing:   sigma=${SMOOTHING}"

# Get data ====

echo 

# echo -n "Identifying subs from derivs ..."
# subs=$(find "${src}" -maxdepth 1 -name "sub-*" -type d -exec basename {} \; | \
#         sed 's/sub-//' | \
#         sort)

echo -n "Identifying subs from existing confounds_timeseries.tsv files ..."

subs=$(find "${src}" -name "*_desc-confounds_timeseries.tsv" \
            -exec basename {} \; | \
        grep -o "sub-[^_]*" | \
        sed 's/sub-//' | \
        sort -u)

n_subs=$(echo "${subs}" | wc -w)
echo " found ${n_subs}."
echo

echo "Starting copy (by default, no overwriting)"

for sub in ${subs} ; do 

    echo "    Copying ${sub} ..."

    sub_dest=${dest}/sub-${sub}/
    mkdir -p "${sub_dest}"/{anat,func}

    rsync -u \
        "${src}"/sub-"${sub}"/anat/*_space-"${space}"_desc-preproc_T1w.nii.gz \
        "${sub_dest}/anat/" 

    rsync -u \
        "${src}"/sub-"${sub}"/anat/*_space-"${space}"_desc-brain_mask.nii.gz \
        "${sub_dest}/anat/"


    # Copy confounds file
    rsync -u \
        "${src}/sub-${sub}"/func/*_desc-confounds_timeseries.* \
        "${sub_dest}/func/"

    # Copy native-space 
    #   task-*_run-* : gets the two tasks an all runs
    #   desc-*       : gets preproc_bold and brain_mask
    source_fmri=$(find "${src}/sub-${sub}/func/" \
                    -name "*_task-*_run-*_space-${space}_desc-*.nii.gz")

    for i in ${source_fmri} ; do

        bn=$(basename "${i}")
        # shellcheck disable=SC2001 
        new_name=$(echo "${bn}" | sed "s/_space-${space}/_space-${space_new}/")
        rsync -u "${i}" "${sub_dest}/func/${new_name}"

    done

done

echo "    Copying events information from BIDS directory ..."

for sub in ${subs} ; do 

    # Copy events files from bids
    sub_dest=${dest}/sub-${sub}/
    rsync -u "${bids}/sub-${sub}"/func/*_events.tsv "${sub_dest}/func/"

done

du -hsc "${dest}"

# Mask files 

echo 
echo "Smoothing and masking func files"
echo "    Smoothing at ${smooth_sigma}"
echo 

for i in "${dest}"/sub-*/func/*_desc-preproc_bold.nii.gz ; do 

    smfile=${i//preproc_bold/smoothedmasked}

    if [ ! -e "${smfile}" ] ; then 
 
        fslmaths "${i}" \
            -s "${smooth_sigma}" \
            -mas "${i//preproc_bold/brain_mask}" \
            "${smfile}"

    fi

    echo "    ${smfile}"

done

du -hsc "${dest}"

# Separating for SPM 

echo
echo "Separating for SPM"

for i in "${dest}"/sub-*/func/*_desc-smoothedmasked.nii.gz ; do 

    bn=$(basename "${i}" .nii.gz)
    dir=$(dirname "${i}")/${bn}

    echo "    ${i}"

    if [ ! -e "${dir}/${bn}_0000.nii" ] ; then

        mkdir -p "${dir}"
        fslsplit "${i}" "${dir}/${bn}_" -t
        gunzip "${dir}"/*.nii.gz

    fi

done 

du -hsc "${dest}"

# Create motion file

echo 
echo "Creating motion files"

for i in "${dest}"/sub-*/func/*_desc-smoothedmasked.nii.gz ; do 

    # Remove space label from basename to identify confounds, confounds are
    #   space-indep.
    bn=$(basename "${i}" _desc-smoothedmasked.nii.gz | sed 's/_space-[^_]*//')
    dir=$(dirname "${i}")
    motion_file=${i//.nii.gz/_rp.txt}

    if [ ! -e "${motion_file}" ] ; then

        # This script extracts the base 6 motion parameters from the confounds
        #   file and writes it to an SPM-friendly .txt file. 
        Rscript "${home}"/convert_confounds_to_rp.R \
            "${dir}/${bn}_desc-confounds_timeseries.tsv" \
            "${motion_file}"

    fi

    echo "    ${motion_file}"

done 
