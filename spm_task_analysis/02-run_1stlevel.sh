#!/bin/bash

# Need to pass full path to directory, or else SPM.mat saves only the relative
#   path ...
analysis_dir=$(readlink -f "${1}")
task=${2}
spm_path=${3}
spm_scripts_path=${4}
contrasts=${5}

echo 
echo -n "Identifying subs from ${analysis_dir} ..."
subs=$(find "${analysis_dir}" -maxdepth 1 -name "sub-*" -type d \
            -exec basename {} \; | \
        sed 's/sub-//' | \
        sort)
n_subs=$(echo "${subs}" | wc -w)
echo " found ${n_subs}."
echo

# Invoking MATLAB

echo "Starting 1st level"

# Add scripts to path
export MATLABPATH=${spm_scripts_path}

for sub in ${subs} ; do 

    echo "${sub}"

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

        motion_file=${rdir}_rp.txt

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
            break 2
        fi

        tr=$(fslval "${rdir}"/*_0000.nii pixdim4)
        echo "        Found TR: ${tr}"

        if [ ! -e "${rdir}/results/SPM.mat" ] ; then 

            matlab \
                -batch "func_1stlevel('${spm_path}', '${analysis_dir}', \
                            '${sub}', '${bn_rdir}', \
                            0, '${tr}', '${contrasts}', \
                            '${events_file}', '${motion_file}')"

        else 

            echo "    Skipping 1st level, already done"

        fi

        # echo "PWD: $(pwd)"

        # rp file must match 'rp_.*.txt'
        cp "${motion_file}" "${rdir}"/rp__.txt

        firstlevel_mat=${rdir}/results/SPM.mat
        optcens_mat=${rdir}/results/optcens/optcens.mat

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