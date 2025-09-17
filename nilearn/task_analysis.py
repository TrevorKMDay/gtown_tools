
import argparse as ap
import json
import pprint as pp
import numpy as np
import glob
import re
import tempfile as tf

import os
from pathlib import Path
from shutil import copy

import nilearn as nil
from nilearn.glm.first_level import first_level_from_bids
from nilearn.interfaces.bids import save_glm_to_bids

print(f"Using nilearn version: {nil.__version__}")

def clean_contrast_name(contrast_name):

    # Stolen from nilearn dev - Sept. 18, 2024

    new_name = contrast_name[:]

    # Some characters translate to words
    new_name = new_name.replace("-", " Minus ")
    new_name = new_name.replace("+", " Plus ")
    new_name = new_name.replace(">", " Gt ")
    new_name = new_name.replace("<", " Lt ")

    # Others translate to spaces
    new_name = new_name.replace("_", " ")

    # Convert to camelCase
    new_name = new_name.split(" ")
    new_name[0] = new_name[0].lower()
    new_name[1:] = [c.title() for c in new_name[1:]]
    new_name = " ".join(new_name)

    # Remove non-alphanumeric characters
    new_name = "".join(ch for ch in new_name if ch.isalnum())

    return new_name


# Get arguments ====

parser = ap.ArgumentParser()

parser.add_argument("bids_dir", type=str,
                    help="BIDS source directory.")
parser.add_argument("derivatives_dir", type=str,
                    help="fMRIPREP outputs directory.")
parser.add_argument("task",
                    help="Name of task to process, without 'task-'.")
parser.add_argument("contrast_file",
                    help="JSON file with contrasts listed.")
parser.add_argument("out_dir", type=str,
                    help="Directory to save the individual maps to.")

# Input options

input = parser.add_mutually_exclusive_group()

input.add_argument("--subs", "-s", nargs="+",
                   help="Subs to analyze. Without 'sub-'.")

input.add_argument("--fmriprep_only", "-F", action="store_true",
                   help="Use only subs with fMRIPREP output for selected "
                        "task (and space).")

# Other options

parser.add_argument("--skip", nargs="+",
                    help="Skip these participants.")

parser.add_argument("--space", default="MNI152NLin2009cAsym",
                    help="Name of output space to work in.")

parser.add_argument("--overwrite", "-o", action="store_true",
                    help="Overwrite results if they exist?")

parser.add_argument("--test", "-1", action="store_true",
                    help="Only run one participant, as a test.")

parser.add_argument("--filter_file", "-f", metavar="FILE",
                    help="JSON file with which runs per to filter. "
                         "If a participant is missing from this file, all "
                         "runs are used.")

parser.add_argument("--stop-on-failure", action="store_true")

model_settings = parser.add_argument_group("Model Settings")

model_settings.add_argument("--hrf", default="spm",
                            choices=["spm", "glover"],
                            help="HRF model to use.")

model_settings.add_argument("--hrf_deriv", action="store_true",
                            help="Add derivative parameter to HRF model.")

model_settings.add_argument("--hrf_disp", action="store_true",
                            help="Add dispersion parameter to HRF model.")

model_settings.add_argument("-smooth_fwhm",
                            help="Smoothing kernel (mm)")

model_settings.add_argument("--strategy", nargs='*',
                            default=["motion", "high_pass", "wm_csf"],
                            choices=["motion", "wm_csf", "global_signal",
                                     "compcor", "ica_aroma", "scrub",
                                     "high_pass"],
                            help="The type of noise components to include.")

model_settings.add_argument("--motion", default="full",
                            choices=["basic", "power2", "derivatives", "full"],
                            help="Type of confounds extracted from head motion"
                                  " estimates.")

model_settings.add_argument("--wm_csf", default="full",
                            choices=["basic", "power2", "derivatives", "full"],
                            help="Type of confounds extracted from masks of "
                                 "white matter and cerebrospinal fluids.")

model_settings.add_argument("--global_signal", default="full",
                            choices=["basic", "power2", "derivatives", "full"],
                            help="Type of confounds extracted from the global "
                                 "signal.")

model_settings.add_argument("--scrub", default=5, metavar="N",
                            help="After accounting for time frames with "
                                 "excessive motion, further remove segments "
                                 "shorter than the given number.")

model_settings.add_argument("--fd", default=0.5, metavar="mm", type=float,
                            help="Framewise displacement threshold for scrub "
                                 "in mm. Default: 0.5")


args = parser.parse_args()

# pp.pprint(args)

# Main loop =====

# Give CLI args better names
bids_dataset = args.bids_dir
derivatives_folder = args.derivatives_dir
task_label = args.task
space_label = args.space
fail_stop = args.stop_on_failure

if args.filter_file is not None:
    filter_file = args.filter_file
    print(f"Using filter file {filter_file}")
else:
    filter_file = None

contrasts_file = args.contrast_file

out_dir = args.out_dir

# HRF settings
hrf_deriv = args.hrf_deriv
hrf_disp = args.hrf_disp

if filter_file is not None:

    # Get list of subj
    with open(filter_file, "r") as f:
        include1 = json.load(f)

    # These are the subs that have a filter rule for the current task

    subs_to_filter = [x["sub"].replace("sub-", "") for x in include1
                      if x["task"][0]["task"] == f"task-{task_label}"]

    print(f"Found {len(subs_to_filter)} subjects to filter runs for task "
          f"{task_label}.")

if args.fmriprep_only:

    # This checks in the supplied derivatives directory to only use those that
    # are processed.

    the_glob = os.path.join(derivatives_folder, "sub-*", "func",
                            f"*_task-{task_label}*{space_label}*"
                             "_desc-preproc_bold.nii.gz")

    # Get matching files; extract subs with matching files
    matching_files = glob.glob(the_glob)
    task_files = [os.path.basename(x) for x in matching_files]
    subs_with_task = set([re.search("sub-[^_]*", x).group(0)
                          for x in task_files])

    sub_labels = [x.replace("sub-", "") for x in subs_with_task]

elif args.subs is not None:

    sub_labels = args.subs

elif args.subs is None:

    dirs = glob.glob(os.path.join(bids_dataset, "sub-*"))
    sub_labels = [re.search("sub-[^_]*", x).group(0).replace("sub-", "")
                  for x in dirs]

# quit()

if len(sub_labels) == 0:

    raise ValueError("No subjects matching criteria found; please check")

sub_labels.sort()

if args.skip is not None:
    sub_labels = [x for x in sub_labels if x not in args.skip]

# Check arguments

if hrf_disp and not hrf_deriv:
    print("HRF dispersion requires derivative, adding to model.")
    hrf_deriv = True

args.strategy.sort()

print("\nStrategy:")
pp.pprint(args.strategy)

## Get contrasts from file

with open(contrasts_file, 'r') as f:
    contrasts = json.load(f)

# Get contrast labels

contrast_names = list(contrasts.keys())
contrast_names2 = [x.replace(" - ", "Minus") for x in contrast_names]

# Calculate the first-level model

if args.test:
    sub_labels = [sub_labels[0]]

# Display what we have so far
print()
print("Contrasts:")
pp.pprint(contrasts)
print()
print("Participants:")
pp.pprint(sub_labels, compact=True)
print()

# Settings

settings = {}

settings["space"] = space_label

settings["hrf"] = "spm" if args.hrf == "spm" else "glover"

if hrf_deriv:
    settings["hrf"] += " + derivative"

if hrf_disp:
    settings["hrf"] += " + dispersion"

strategy = args.strategy
strategy.sort()

settings["confounds_strategy"] = strategy
settings["confounds_motion"] = args.motion
settings["confounds_fd_threshold"] = args.fd

settings_file = f"{out_dir}/settings.json"

print("Checking settings ...")

os.makedirs(out_dir, exist_ok=True)

# This checks that the current output directory was created with the same
#   settings. If the settings don't match, it dies. If it does match, it
#   continues, assuming that you are rerunning because there's new participants

if not os.path.exists(settings_file):

    with open(settings_file, "w") as f:
        json.dump(settings, f, indent=4)

else:

    with open(settings_file, "r") as f:
        settings_already = json.load(f)

    if settings != settings_already:
        print("ERROR: Settings in existing settings file don't match, dying.")
        print("\nRequested settings:")
        pp.pprint(settings)
        print("Settings in directory:")
        pp.pprint(settings_already)
        exit(1)
    else:
        print("Settings OK!\n")

pp.pprint(settings)
print()

def check_dir(dir, contrasts, img=None):

    ok = True

    print()
    print(f"Checking output dir {dir} ...")

    if img is not None:
        for i in img:
            design_svg = i.replace("_space-*.nii.gz", "_design.tsv")
            if os.path.exists(design_svg):
                print(f"    {design_svg} exists")
            else:
                ok = False

    for contrast in contrasts:

        con = clean_contrast_name(contrast)
        matching_files = glob.glob(f"{dir}/"
                                   f"*task-{task_label}_"
                                   f"contrast-{con}_*.nii.gz")

        n_files = len(matching_files)
        print(f"    Found {n_files} for task-{task_label} {con} (need 5)")

        if n_files != 5:
            ok = False

    check="OK" if ok else "NOT OK"
    print(f"  Result of checking {dir} is: {check}")

    return(ok)


def copy_sub(sub, src, temp, filter):

    # src:  tuple of bids/derivs to copy from
    # temp: tuple of temp bids/derivs to copy to

    # Setup directories
    bids_dest = Path(temp[0]) / f"sub-{sub}" / "func"
    derivs_dest = Path(temp[1]) / f"sub-{sub}" / "func"
    dest = (bids_dest, derivs_dest)

    [os.makedirs(x) for x in dest]

    if filter is not None:

        # Only copy files based on list if they're in the subs_to_filter
        #   list

        for run in sub_task_runnames:

            run_regex = re.search("run-[0-9]*", run)
            # pp.pprint(run_regex.group(0))

            if run_regex:
                run_id = run_regex.group(0)
            else:
                # If this fails, something really weird is going on
                print(f"ERROR: No 'run-#' sequence found for {sub}"
                        f"{task_label}!")
                continue

            # Only need event files from BIDS
            bfiles = glob.glob(f"{src[0]}/sub-{sub}/func/"
                               f"sub-{sub}_task-{task_label}_{run_id}"
                               "*_events.tsv")

            # Copy relevant files from derivative - this could probably
            #   be pared down to save space/copy time, but it's pretty fast
            #   as-is.
            dbolds = f"{src[1]}sub-{sub}/func/" + \
                     f"sub-{sub}_task-{task_label}_{run_id}_{space_label}_" + \
                     "*_bold.*"

            dconfounds =f"{src[1]}sub-{sub}/func/" + \
                        f"sub-{sub}_task-{task_label}_{run_id}_" + \
                        "desc-confounds_timeseries.tsv"

            dfiles = glob.glob(dbolds) + glob.glob(dconfounds)

            # pp.pprint(dfiles)

            [copy(x, dest[0]) for x in bfiles]
            [copy(x, dest[1]) for x in dfiles]

    else:

        # Copy everything

        print(f"    {sub} {task_label} all runs")

        # Only need event files from BIDS
        bfiles = glob.glob(f"{src[0]}/sub-{sub}/func/"
                           f"sub-{sub}_task-{task_label}_*_events.tsv")

        dbolds = glob.glob(f"{src[1]}/sub-{sub}/func/"
                           f"sub-{sub}_task-{task_label}*{space_label}*"
                           "desc-preproc_bold.*")

        dconfounds = glob.glob(f"{src[1]}/sub-{sub}/func/*"
                           "_desc-confounds_timeseries.*")

        dfiles = dbolds + dconfounds

        [copy(x, dest[0]) for x in bfiles]
        [copy(x, dest[1]) for x in dfiles]

# Set up copies ====

sub_labels = [s for s in sub_labels
              if not check_dir(f"{out_dir}/sub-{s}/", contrast_names)]

print("\nReduced subjects to the following subjects needing processing:")
pp.pprint(sub_labels, compact=True)

if filter_file is not None:

    print("\nCreating temporary directories to exclude runs ...")

    bids_dir = tf.TemporaryDirectory(delete = False, suffix="_bids")
    derivs_dir = tf.TemporaryDirectory(delete = False, suffix="_derivs")

    temp_dirs = (bids_dir.name, derivs_dir.name)

    print(f"    BIDS: {temp_dirs[0]}")
    print(f"    Derivatives: {temp_dirs[1]}")
    print()

    for sub in sub_labels:

        if sub in subs_to_filter:

            print(f"Setting up {sub} (filtering) ...")

            # The the sub was in the filter list, then do the subset copy

            # Get the tasks for this subject
            sub_tasks = [x["task"] for x in include1
                        if x["sub"] == f"sub-{sub}"]

            # Filter down to the current task
            sub_task_runs = [x[0]["runs"] for x in sub_tasks
                            if x[0]["task"] == f"task-{task_label}"]

            # And extract the actual run name from the object (also stores
            #   total_pct_outlier)
            sub_task_runnames = [x["run_name"] for x in sub_task_runs[0]]
            sub_task_runnames.sort()

            print(f"    {sub} {task_label} {sub_task_runnames}")

            copy_sub(sub=sub, src=(bids_dataset, derivatives_folder),
                     temp=temp_dirs, filter=sub_task_runnames)

        else:

            print(f"Setting up {sub} (copy all) ...")

            # Otherwise, we still have to copy everything to one place, but
            # we can just copy everything.

            copy_sub(sub=sub, src=(bids_dataset, derivatives_folder),
                     temp=temp_dirs, filter=None)

    # Reset the values for the BIDS/derivatives to the new copied location
    bids_dataset, derivatives_folder = temp_dirs

# CREATE MODELS ====

space_search = re.search("space-[^_]*", space_label)
res_search = re.search("res-[0-9]*", space_label)
cohort_search = re.search("cohort-[0-9]*", space_label)

if space_search:
    space_only = space_search.group(0).replace("space-", "")

img_filter_list = []

if res_search:
    res_label = res_search.group(0)
    img_filter_list += [("res", res_label.replace("res-", ""))]

if cohort_search:
    cohort_label = cohort_search.group(0)
    img_filter_list += [("cohort", cohort_label.replace("cohort-", ""))]

if len(img_filter_list) == 0:
    img_filter_list = None
else:
    print("\nUsing the following filters:")
    pp.pprint(img_filter_list)

(

    models,
    models_run_imgs,
    models_events,
    models_confounds,

) = first_level_from_bids(

    bids_dataset,
    task_label=task_label,
    space_label=space_only,
    sub_labels=sub_labels,
    derivatives_folder=derivatives_folder,
    img_filters=img_filter_list,

    smoothing_fwhm=5.0,
    drift_model=None,

    # Add more things to regression strategy
    # Set motion threshold to 0.5 mm for task
    confounds_strategy=settings["confounds_strategy"],
    confounds_motion=settings["confounds_motion"],
    confounds_fd_threshold=settings["confounds_fd_threshold"],

    # Running settings
    verbose=3,
    minimize_memory=False,

    # Set these values to None to infer
    t_r=None,
    slice_time_ref=None,
    hrf_model = settings["hrf"],

    # job control
    # n_jobs=-2

)

# pp.pprint(models_confounds[0])

# The models are not being generated in alphabetical order, fix that here

models_orig = [m.subject_label for m in models]
new_order = np.argsort(models_orig).tolist()

models_ord = [models[i] for i in new_order]
models_run_imgs_ord = [models_run_imgs[i] for i in new_order]
models_events_ord = [models_events[i] for i in new_order]
models_confounds_ord = [models_confounds[i] for i in new_order]

# Run the models

for subject, model, imgs, event, confound in zip(
    sub_labels, models_ord, models_run_imgs_ord, models_events_ord,
    models_confounds_ord
):

    out_sub = f"{out_dir}/sub-{subject}/"
    os.makedirs(out_sub, exist_ok=True)

    # If there is only one run, then we can't loop over a single pandas DF;
    #   so put this in a list to loop over it in the confound-checking stage.
    if type(confound) is not list:
        confound = [confound]

    for i in confound:

        # Check that there's not more predictors than observations (i.e.
        #   frames). If there is, matrix error will occur.

        n_frames, n_confounds = i.shape
        print(f"Design matrix is {n_frames}x{n_confounds}")
        print(f"    {n_confounds} should be < {n_frames}")

    try:

        m = model.fit(imgs, events=event, confounds=confound)

    except ValueError:

        total = len(event)
        for i, x in enumerate(event):

            print(f"Event table {i + 1}/{total}")
            print(event[i])

            print(f"Confounds table {i + 1}/{total}")
            print(confound[i])

        if not fail_stop:
            print(f"Continuing past ValueError for {subject} ... ")
            continue

    # print(out_sub)

    save_glm_to_bids(
        m,
        contrasts=contrast_names,
        contrast_types=contrasts,
        out_dir=out_dir,
        prefix=f"sub-{subject}_task-{task_label}"
    )

    del model

    # break

if filter_file is not None:
    [x.cleanup() for x in [bids_dir, derivs_dir]]