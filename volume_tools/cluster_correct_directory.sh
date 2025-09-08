#!/bin/bash

# This is designed to take an output directory from nilearn, which takes the
#   layout dir/{group,sub-1,sub-2,...}

# nilearn calculates one errorts file per run, but the modeling eventually
#   combines the run. We ignore this

dir=${1}

# Settings for cluster; NN faces, corners edges (default: 1, 2 3)
# side: 1/2 sided test (do 1)
# thresh: t thresh (2-3 or so)
NN=${2}
side=${3}
thresh=${4}

# Maybe make these settable later
# Don't go less stringent than p<.001, based on Andy's Brain Book
athr=0.05
pthr=0.001

cat > "${dir}/clustering_settings.json" <<EOF
{
    "NN": ${NN},
    "side": ${side},
    "thresh": ${thresh},
    "datetime": "$(date)"
}
EOF

# 'sub-' prefix skips group/ directory
errorts_files=$(find "${dir}" -name "sub-*_stat-errorts_statmap.nii.gz" | sort)

# shellcheck disable=SC2086
n_errorts_files=$(echo ${errorts_files} |wc -w)
echo "Found ${n_errorts_files} errorts files!"

for errorts in ${errorts_files} ; do

    # shellcheck disable=SC2001
    pfx=$(echo "${errorts}" | sed 's/_stat-errorts_statmap.nii.gz//')

    echo
    echo "Working on cluster correcting: ${pfx}"

    fwhmx="${pfx}_3dFWHMx.txt"

    echo "  Starting 3dFWHM calculation"

    if [ ! -e "${fwhmx}" ] ; then

        3dFWHMx -acf -automask -input "${errorts}" > "${fwhmx}"

    else

        echo "    ${fwhmx} already exists!"

    fi

    echo "  Starting cluster simulation"

    if [ ! -e "${pfx}_cluster.NN1_bisided.1D" ] ; then

        # remove last value (see ref 1)
        values=$(head -n2 "${fwhmx}" | tail -n1 | sed 's/ [0-9]*[.][0-9]*$//')
        echo "    Using values: ${values}"

        # shellcheck disable=SC2086
        3dClustSim \
            -nodec                      \
            -acf    ${values}           \
            -athr   ${athr}             \
            -pthr   ${pthr}             \
            -prefix ${pfx}_cluster

    else

        echo "    ${pfx}_cluster.NN1_bisided.1D already exists!"

    fi

    echo "  Starting cluster correction"

    nvox=$(tail -n1 "${pfx}_cluster.NN${NN}_${side}sided.1D" | \
        tr -s ' ' ' ' | \
        awk '{print $2}')

    # Make sure not to get any of the other files

    # Find by sub name, because run- details get lost
    sub=$(basename "${pfx}" | grep -o "sub-[^_]*")
    con_files=$(find "${dir}/${sub}" \
                    -name "${sub}_*_contrast-*_stat-t_statmap.nii.gz")

    echo "   Found contrast files: ${con_files}"

    for f in ${con_files}; do

        echo "    ${f}"

        # shellcheck disable=SC2001
        corrected=$(echo "${f}" | sed "s/.nii.gz/_cc-${thresh}.nii.gz/")

        if [ ! -e "${corrected}" ] ; then

            3dClusterize                                \
                -inset          "${f}"                  \
                -ithr           0                       \
                -idat           0                       \
                -NN             "${NN}"                 \
                -"${side}sided" RIGHT "${thresh}"       \
                -clust_nvox     "${nvox}"               \
                -pref_dat       "${corrected}"

        fi

    done

    # break

done

# Clean up directory
# find "${dir}" \( -name "*.1D" -not -name "*NN${NN}" \) -delete