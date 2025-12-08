#!/bin/bash

if [ ${#} -eq 2 ] ; then

    home=${1}
    task=${2}

else

    echo "Usage: ${0} dir task"
    exit 1

fi

# Allows use of func_1stlevel_postOCT()
export MATLABPATH=~/code/preproc_spm/

# Find all the directories with the task in the name, then get the subject
#   labels and get the unique subs

subs=$(find ${home}/sub-*/func/ -maxdepth 1 -type d \
            -name "sub-[^_]*_task-${task}_*" \
            -exec basename {} \; | \
        grep -o 'sub-[^_]*' | \
        sort -u)

echo ${subs}

for s in ${subs} ; do

    dirs="$(find ${home}/${s}/func -maxdepth 1 -type d \
                -name "${s}_task-${task}_*_desc-smoothedmasked" | \
            sort | \
            tr '\n' ' ' | \
            sed 's/ $//')"

    echo $dirs

    odir=${home}/${s}/func/task-${task}_postOCT/

    echo
    echo "sub: ${s}"
    echo "================"

    # Try resetting the MATLAB cache each time to avoid the seg fault
    MCR_CACHE_ROOT=$(mktemp -d /tmp/mcr.XXXXXX)
    export MCR_CACHE_ROOT

    # if [ ! -e ${odir} ] ; then

        matlab -batch "func_1stlevel_postOCT('${dirs}', '3', 'cens', '${odir}')"

    # fi

    # break

done