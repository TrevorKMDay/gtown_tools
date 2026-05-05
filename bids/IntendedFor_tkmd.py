#!/usr/bin/env python3
import json
import os
from glob import glob
import argparse
# import numpy
# import subprocess
import pprint

# Script originally by Tim Hendrickson (2019)
#   https://github.com/tjhendrickson/BIDS_scripts
# Modified by Trevor Day (add interface, 2023)
#   Upload to Gist 2025


pp = pprint.PrettyPrinter(indent=4)

parser = argparse.ArgumentParser()

parser.add_argument("SES",
                    help="Path(s) to session directories",
                    nargs="+")

parser.add_argument("-v", "--show-shim-values",
                    help="If set, echoes shim values for fmaps/funcs",
                    action="store_true")

args = parser.parse_args()

sessions = args.SES
verbose = args.show_shim_values

# Function definitions


def setup(subject_path, verbose=False):

    print()
    print(f"Working in {subject_path}")

    path_contents = os.listdir(subject_path)

    if len(path_contents) == 0:
        print(" ++ Directory is empty, skipping.")
        # print(subject_path)
    else:

        # if "ses" in os.listdir(subject_path)[0]:
        #     for item in (os.listdir(subject_path)):
        #         session_path = subject_path + '/' + item
        #         IntendedFor(session_path)

        IntendedFor(subject_path, verbose=verbose)

def IntendedFor(data_path, verbose=False):

    # Identify JSONs in directory
    fmaps = sorted(glob(data_path + '/fmap/*.json'))
    funcs = sorted(glob(data_path + '/func/*bold.json'))

    for fmap in fmaps:

        # change fmap file permissions to allow write access
        # os.system("chmod 660 " + fmap)

        with open(fmap, 'r') as f:
            fmap_json = json.load(f)
            f.close()

        if "IntendedFor" in fmap_json:
            del fmap_json["IntendedFor"]

        shim_fmap = fmap_json["ShimSetting"]

        if verbose:
            print()
            print(f"+ INFO: Shim setting for {os.path.basename(fmap)} "
                  f"is: {shim_fmap}")

        # patient_pos_fmap = fmap_json["ImageOrientationPatientDICOM"]
        # acquisition_time_fmap = fmap_json["AcquisitionTime"]
        func_list = []

        for func in funcs:

            with open(func, 'r') as g:
                func_json = json.load(g)
                shim_func = func_json["ShimSetting"]

                if verbose:
                    print(f"++ INFO: Shim setting for "
                          f"{os.path.basename(func)} is: {shim_func}")

                # patient_pos_func = func_json["ImageOrientationPatientDICOM"]
                # acquisition_time_func = func_json["AcquisitionTime"]
                g.close()

            if shim_fmap == shim_func:

                # Find nifti files that correspond to JSON we are looking at
                func_glob = func.replace("json", "nii*")
                # print(f"glob: {func_glob}")

                func_nii_glob = glob(func_glob)

                if len(func_nii_glob) > 0:

                    func_nii = func_nii_glob[0]

                    if "ses" in data_path:
                        func_nii = "/".join(func_nii.split("/")[-3:])
                        func_list.append(func_nii)
                    else:
                        func_nii = "/".join(func_nii.split("/")[-2:])
                        func_list.append(func_nii)

                else:
                    print(" ++ WARNING: glob did not match any results:")
                    print("  ++ " + func_glob)

        if len(func_list) == 0:
            print(f"++ WARNING: Found no compatible bolds for {fmap}")

        entry = {"IntendedFor": func_list}
        fmap_json.update(entry)

        with open(fmap, 'w') as f:
            json.dump(fmap_json, f, indent=4)
            f.close()

        # print(fmap)
        # pp.pprint(fmap_json)
        # print()

        # change file permissions to read only
        # os.system("chmod 444 " + fmap)


# Run main loop

# Check sessions exist
sessions_abs = [os.path.abspath(x) for x in sessions]
sessions_exist = [x for x in sessions_abs if os.path.exists(x)]
sessions_bad = [x for x in sessions_abs if not os.path.exists(x)]

if len(sessions_bad) > 0:
    print("WARNING: Dropped these bad sessions:")
    pp.pprint(sessions_bad)
    print()

# This code checked that each path has a 'ses-' label, which isn't required
#   (can be: sub-X/func/), but this wasn't flexible enough for that

sessions_ses = sessions_exist

# Check sessions contain ses- label
# sessions_ses = [x for x in sessions_exist if "ses-" in x]
# sessions_noses = [x for x in sessions_exist if "ses-" not in x]

# if len(sessions_noses) > 0:
#     print("WARNING: Dropped these sessions that don't include a 'ses-' label:")
#     pp.pprint(sessions_noses)
#     print()

for subject_path in sessions_ses:
    setup(subject_path, verbose=verbose)
