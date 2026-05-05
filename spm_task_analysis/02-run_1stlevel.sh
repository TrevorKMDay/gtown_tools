#!/bin/bash

# Need to pass full path to directory, or else SPM.mat saves only the relative
#   path ...
analysis_dir=$(readlink -f "${1}")
task=${2}
spm_path=${3}
spm_scripts_path=${4}

# Contrasts should only include non >Rest contrasts
contrasts=${5}

# Set mthresh to 0 to do less aggressive masking for problem runs
mthresh=${6}

shift 6
sub_input=${*}
n_sub_input=$(echo "${sub_input}" | wc -w)

if [ "${n_sub_input}" -eq 0 ] ; then

    echo
    echo -n "Identifying subs from ${analysis_dir} ..."
    subs=$(find "${analysis_dir}" -maxdepth 1 -name "sub-*" -type d \
                -exec basename {} \; | \
            sed 's/sub-//' | \
            sort)
    n_subs=$(echo "${subs}" | wc -w)
    echo " found ${n_subs}."
    echo

else

    echo "Using ${n_sub_input} supplied subs"
    echo "${sub_input}"
    subs=${sub_input}

fi

# Invoking MATLAB

echo "Starting 1st level"

# Add scripts to path
export MATLABPATH=${spm_scripts_path}

for sub in ${subs} ; do

    echo "${sub}"

    # Check that there's an anatomical files to use, if not, skip this
    #   participant
    anat=$(find "${analysis_dir}/sub-${sub}/anat/" -name "wsub-*_T1w.nii")
    if [ "$(echo "${anat}" | wc -w)" -eq 0 ] ; then
        echo "    No anatomical file found, skipping"
        continue
    else
        echo "    Found anat file: ${anat}"
    fi

    # Check that there's at least one func directory, if not, skip this
    #   participant
    dirs=$(find "${analysis_dir}/sub-${sub}/func/" -maxdepth 1 -type d \
            -name "sub-${sub}_${task}_*")

    if [ "$(echo "${dirs}" | wc -w)" -eq 0 ] ; then
        echo "    No ${task} dirs found, skipping"
        continue
    else
        echo "    Found ${task} dirs:"
        for d in ${dirs} ; do
            echo "        ${d}"
        done
    fi

    for rdir in ${dirs} ; do

        # Try resetting the MATLAB cache each time to avoid the seg fault
        MCR_CACHE_ROOT=$(mktemp -d /tmp/mcr.XXXXXX)
        export MCR_CACHE_ROOT

        echo "    ${rdir}"
        bn_rdir="$(basename "${rdir}")"

        # Find the
        motion_file=$(find "${rdir}/01_orig/" -name "rp_*.txt")

        # Replace suffix and remove space label; it's not included in the
        #   events.tsv file name
        events_file=$(echo "${rdir}" | \
                        sed -e 's/_bold//'          \
                            -e 's/_space-[^_]*//'   \
                            -e 's/_desc-[^_]*//'    )

        # Add suffix if not present
        [[ ${events_file} == *+_events.tsv ]] || events_file+=_events.tsv

        if [ ! -f "${events_file}" ] ; then
            echo "Can't find missing events file: ${events_file}"
            echo "Fix and proceed."
            continue 2
        fi

        firstlevel_mat=${rdir}/05_warpsmooth/results/SPM.mat
        if [ ! -e "${firstlevel_mat}" ] ; then

            # Get TR from the first file in the input; if there's a problem,
            # skip this subject
            ws00=$(find "${rdir}"/05_warpsmooth/ -maxdepth 1 \
                    -name "*_0000.nii")

            if [[ "${ws00}" == "" ]] ; then
                echo "        No warpsmooth files to work with, skipping"
                break
            fi

            tr=$(fslval "${ws00}" pixdim4)

            matlab \
                -batch "func_1stlevel('${analysis_dir}', \
                            '${sub}', '${bn_rdir}/05_warpsmooth', \
                            0, '${tr}', '${contrasts}', \
                            '${events_file}', '${motion_file}', '${mthresh}')"

        else

            echo "    Skipping 1st level, already done"

        fi

        # echo "PWD: $(pwd)"

        # rp file must match 'rp_.*.txt'
        rsync -u "${motion_file}" "${rdir}"/rp__.txt

        optcens_mat=${rdir}/05_warpsmooth/results/optcens/optcens.mat

        if [ ! -e "${optcens_mat}" ] ; then

            matlab \
                -nodisplay -nosplash \
                -batch "addpath('${spm_path}') ; \
                            spm_get_defaults('cmdline',true); spm ; \
                            echo on ; \
                            call_optcens('${firstlevel_mat}')"

        else

            echo "    OCT already done"

        fi

        # break 2

    done

    # break

done