#!/bin/bash

bids=~/Projects/LBS/ferrara_lbs/bids
subs=$(find ${bids} -maxdepth 1 -type d -name "sub-*" -exec basename {} \; | \
        sort)

ptsv=${bids}/participants.tsv

for s in ${subs} ; do 

    sub=${s//sub-/}
    age=$(grep ${sub} ${ptsv} | awk '{print $2}')


    # Technically, there are participants beyond 11 y with the correct
    #   ages (cohort-4: 7.5-13.5 y). Should probably run fMRIPREP with
    #   cohort-1 (4.5 - 18.5 y)

    # The (( )) wrapper flips bc's 0/1 output to T/F for `if`
    if (( $(echo "${age} < 8" | bc -l) )) ; then
        # 4.5 - 8.5 y
        cohort=2
    else
        # 7.0 - 11.0 y
        cohort=3
    fi

    ./01-1stlevel_fP_in_SPM.sh ${sub} ${cohort}

done