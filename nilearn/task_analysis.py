
import argparse as ap
import json
import pprint as pp
import numpy as np
import os
import glob
import re

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
parser.add_argument("space",
                    help="Name of output space to work in, without 'space-'.")
parser.add_argument("contrast_file", 
                    help="JSON file with contrasts listed.")
parser.add_argument("out_dir", type=str,
                    help="Directory to save maps to. ")

parser.add_argument("--subs", "-s", nargs="+",
                    help="Subs to analyze. Without 'sub-'.")

parser.add_argument("--overwrite", "-o", action="store_true",
                    help="Overwrite results if they exist?")


args = parser.parse_args()

# Main loop =====

# Give CLI args better names
bids_dataset = args.bids_dir
derivatives_folder = args.derivatives_dir
sub_labels = args.subs
task_label = args.task
space_label = args.space
contrasts_file = args.contrast_file
out_dir = args.out_dir

## Get contrasts from file

with open(contrasts_file, 'r') as f:
    contrasts = json.load(f)

# Get contrast labels

contrast_names = list(contrasts.keys())
contrast_names2 = [x.replace(" - ", "Minus") for x in contrast_names]

pp.pprint(contrasts)

# Calculate the first-level model

if sub_labels is not None:
    print(f"Participants: {sub_labels}\n")
else:
    dirs = glob.glob(os.path.join(bids_dataset, "sub-*"))
    sub_labels = [re.search("sub-[^_]*", x).group(0).replace("sub-", "") 
                  for x in dirs]

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
    smoothing_fwhm=5.0,
    derivatives_folder=derivatives_folder,
    drift_model=None,

    # Add more things to regression strategy (including deriv, power2)
    # Set motion threshold to 0.5 mm for task 
    confounds_strategy=("motion", "compcor", "scrub", "high_pass"),

    confounds_motion="full",
    confounds_fd_threshold=0.5,

    # Running settings
    # verbose=2,
    minimize_memory=False,

    # Set these values to None to infer
    t_r=None,
    slice_time_ref=None,
    hrf_model = "spm + derivative + dispersion"
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
    
    # pp.pprint(event)
    # pp.pprint(confound[0])

    for i in confound:

        # Check that there's not more predictors than observations (i.e.
        #   frames). If there is, matrix error will occur.

        n_frames, n_confounds = i.shape
        print(f"Design matrix is {n_frames}x{n_confounds}")
        print(f"    {n_confounds} should be < {n_frames}")

    # Check if the model has already been run and saved
    out_sub = f"{out_dir}/sub-{subject}/"

    if not check_dir(out_sub, imgs, contrast_names):

        os.makedirs(out_sub, exist_ok=True)

        m = model.fit(imgs, events=event, confounds=confound)
        
        save_glm_to_bids(
            m,
            contrasts=contrast_names,
            contrast_types=contrasts,
            out_dir=out_dir,
            prefix=f"sub-{subject}_task-{task_label}"
        )

        del model

    # break