#!/bin/bash

home=~/Projects/LBS/ferrara_lbs/code/analysis_spm/fmriprep_input
subs=$(find ${home} -maxdepth 1 -type d -name "sub-TDCh*" \
        -exec basename {} \; | sort)


export MATLABPATH=/Users/tkmd/Projects/LBS/ferrara_lbs/code/preproc_spm/matlab

for s in ${subs} ; do 

    dirs="$(find ${home}/${s}/func -maxdepth 1 -type d \
            -name "*_desc-smoothedmasked" | tr '\n' ' ' | sed 's/ $//')"
    odir=${home}/${s}/func/combined_postOCT/

    echo 
    echo "Sub: ${s}" 
    echo "================"


    # Try resetting the MATLAB cache each time to avoid the seg fault 
    MCR_CACHE_ROOT=$(mktemp -d /tmp/mcr.XXXXXX)
    export MCR_CACHE_ROOT

    if [ ! -e ${odir} ] ; then
    
        matlab -batch "func_1stlevel_postOCT('${dirs}', 3, '${odir}')"

    fi

    # break

done