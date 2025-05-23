
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

input.add_argument("--include", "-i", metavar="FILE", nargs=2,
                   help="JSON file with which runs per sub to include.")

input.add_argument("--fmriprep_only", "-f", action="store_true",
                   help="Use only subs with fMRIPREP output for selected"
                        "task and space.")

# Other options

parser.add_argument("--skip", nargs="+",
                    help="Skip these participants.")

parser.add_argument("--space", default="MNI152NLin2009cAsym",
                    help="Name of output space to work in, without 'space-'.")

parser.add_argument("--overwrite", "-o", action="store_true",
                    help="Overwrite results if they exist?")

parser.add_argument("--test", "-1", action="store_true",
                    help="Only run one participant, as a test.")

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
                                 "in mm.")


args = parser.parse_args()

# Main loop =====

# Give CLI args better names
bids_dataset = args.bids_dir
derivatives_folder = args.derivatives_dir
task_label = args.task
space_label = args.space
failstop = args.stop_on_failure

inclusion_file = args.include[0]
inclusion_min = int(args.include[1])

contrasts_file = args.contrast_file
out_dir = args.out_dir

hrf_deriv = args.hrf_deriv
hrf_disp = args.hrf_disp

# Read in JSON if given

def any_runs_of_x(sub, task, n=1):

    tasks = sub["tasks"]
    task_names = [x["task"] for x in tasks]

    task_ = f"task-{task}"
    if task_ not in task_names:
        return(False)

    runs = [x["runs"] for x in tasks if x["task"] == task_]
    included = [not x["exclude"] for x in runs[0]]
    n_included = sum(included)

    if n_included >= n:
        return(True)
    else:
        return(False)

if inclusion_file is not None:

    # with open(inclusion_file, 'r') as f:
    #     include1 = json.load(f)
    #     include = include1["motion"]

    # sub_labels = [x["sub"].replace("sub-", "") for x in include if
    #                   task_label in x["task"]]

    # Get list of subj
    with open(inclusion_file, "r") as f:
        include1 = json.load(f)

    all_subs = [x["sub"] for x in include1]
    subs_with_task = [x["sub"] for x in include1 if
                      any_runs_of_x(x, task=task_label, n=inclusion_min)]

    print(f"Dropped {len(all_subs) - len(subs_with_task)} participants"
          f"without task-{task_label} (or insufficient low-motion runs).")

    include = [x for x in include1 if x["sub"] in subs_with_task]

    sub_labels = [x.replace("sub-", "") for x in subs_with_task]
    sub_labels.sort()

elif args.fmriprep_only:

    the_glob = os.path.join(derivatives_folder, "sub-*", "func",
                            f"*_task-{task_label}_*_space-{space_label}"
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
print("Contrasts:")
pp.pprint(contrasts)
print()
print("Participants:")
pp.pprint(sub_labels, compact=True)
print()

# Settings

settings = {}

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

    with open(f"{out_dir}/settings.json", "w") as f:
        json.dump(settings, f, indent=4)

else:

    with open(f"{out_dir}/settings.json", "r") as f:
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

if inclusion_file is not None:

    print()
    print("Creating temporary directories to exclude runs ...")

    bids_dir = tf.TemporaryDirectory(delete = False, suffix="_bids")
    bids_temp = bids_dir.name

    derivs_dir = tf.TemporaryDirectory(delete = False, suffix="_derivs")
    derivs_temp = derivs_dir.name

    print(f"BIDS: {bids_dataset}\nDerivatives: {derivatives_folder}\n")

    for subtaskruns in include:

        sub = subtaskruns["sub"]

        runs = [x["runs"] for x in subtaskruns["tasks"]
                if x["task"] == f"task-{task_label}"]

        runs_to_use = [x["run"] for x in runs[0] if not x["exclude"]]

        bids_dest = Path(bids_temp) / sub / "func"
        derivs_dest = Path(derivs_temp) / sub / "func"

        print(f"{sub} {task_label} {runs_to_use}")
        [os.makedirs(x) for x in [bids_dest, derivs_dest]]

        for run in runs_to_use:

            # Only need event files from BIDS
            bfiles = glob.glob(f"{bids_dataset}/{sub}/func/"
                               f"{sub}_task-{task_label}_{run}*_events.tsv")

            # Copy relevant files from derivative - this could probably
            #   be pared down to save space/copy time, but it's pretty fast
            #   as-is.
            dfiles = glob.glob(f"{derivatives_folder}/{sub}/func/"
                               f"{sub}_task-{task_label}_{run}_*")

            [copy(x, derivs_dest) for x in dfiles]
            [copy(x, bids_dest) for x in bfiles]

    bids_dataset = bids_dir.name
    derivatives_folder = derivs_dir.name

# CREATE MODELS

(

    models,
    models_run_imgs,
    models_events,
    models_confounds,

) = first_level_from_bids(

    bids_dataset,
    task_label=task_label,
    space_label=space_label,
    sub_labels=sub_labels,
    derivatives_folder=derivatives_folder,

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

# The models are not being generated in alphabetical order fix that here

models_orig = [m.subject_label for m in models]
new_order = np.argsort(models_orig).tolist()

models_ord = [models[i] for i in new_order]
models_run_imgs_ord = [models_run_imgs[i] for i in new_order]
models_events_ord = [models_events[i] for i in new_order]
models_confounds_ord = [models_confounds[i] for i in new_order]

def check_dir(dir, img, contrasts):

    n_img = len(img)
    ok = True

    print()
    print(f"Checking output dir {dir} ...")

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

    print(f"  Result of checking {dir} is: {"OK" if ok else "NOT OK"}")

    return(ok)

# Run the models

for subject, model, imgs, event, confound in zip(
    sub_labels, models_ord, models_run_imgs_ord, models_events_ord,
    models_confounds_ord
):

    out_sub = f"{out_dir}/sub-{subject}/"

    # If there is only one run, then we can't loop over a single pandas DF;
    #   so put this in a list to loop over it in the confound-checking stage.
    if type(confound) is not list:
        confound = [confound]

    if not check_dir(out_sub, imgs, contrast_names):

        os.makedirs(out_sub, exist_ok=True)

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

            if not failstop:
                print(f"Continuing past ValueError for {subject} ... ")
                continue

        save_glm_to_bids(
            m,
            contrasts=contrast_names,
            contrast_types=contrasts,
            out_dir=out_dir,
            prefix=f"sub-{subject}_task-{task_label}"
        )

        del model

    # break

if inclusion_file is not None:
    [x.cleanup() for x in [bids_dir, derivs_dir]]