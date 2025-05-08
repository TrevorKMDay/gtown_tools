
import argparse as ap

import nilearn as nil
from nilearn.glm.first_level import first_level_from_bids

parser = ap.ArgumentParser()

parser.add_argument("bids_dir", type=str,
                    help="BIDS source directory.")
parser.add_argument("derivatives_dir", type=str,
                    help="fMRIPREP outputs directory.")
parser.add_argument("task",
                    help="Name of task to process, without 'task-'.")
parser.add_argument("space",
                    help="Name of output space to work in, without 'space-'.")

args = parser.parse_args()

bids_dataset = args.bids_dir
task_label = args.task
space_label = args.space
derivatives_folder=args.derivatives_dir

(

    models,
    models_run_imgs,
    models_events,
    models_confounds,

) = first_level_from_bids(

    bids_dataset,
    task_label=task_label,
    space_label=space_label,
    derivatives_folder=derivatives_folder,

    smoothing_fwhm=5.0,
    drift_model=None,

    # Running settings
    verbose=3,
    minimize_memory=False,

    # Set these values to None to infer
    t_r=None,
    slice_time_ref=None,

)
