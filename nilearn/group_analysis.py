import pandas as pd
import nibabel as nib
import argparse as ap
import json
import pprint as pp
import re

from pathlib import Path
from nilearn.glm.second_level import SecondLevelModel
from nilearn.interfaces.bids import save_glm_to_bids

# Get arguments ====

parser = ap.ArgumentParser()

input_args = parser.add_mutually_exclusive_group(required=True)

parser.add_argument("--list_of_files", "-f", metavar="FILE", nargs='+',
                    help="List of fixed-effect files to use "
                         "('stat-effect_statmap'). "
                         "Sets task label to 'task-all'.")

input_args.add_argument("--results_dir", "-d", type=str, metavar="DIR",
                        help="Directory with fixed-effects (first-level) "
                             "results. If set, uses all 'stat-effect_statmap' "
                             "files.")

parser.add_argument("--collapse_tasks", action="store_true",
                    default=False,
                    help="Combine all tasks into one?")


parser.add_argument("--list_of_subs", "-s", metavar="SUB", nargs="+",
                        help="List of subjects to run all contrasts for. "
                             "Format: LABEL SUB1 [SUB2 ...]")

args = parser.parse_args()

if args.list_of_subs is not None:
    group_label = args.list_of_subs[0]
    subs = args.list_of_subs[1::]
    subs.sort()
else:
    group_label = None
    subs  = None

# Check input arguments ====

if args.results_dir is not None:

    results_dir=args.results_dir
    results_path = Path(results_dir)
    all_fixed_files = list(results_path.glob("**/"
                                             "sub-*stat-effect_statmap.nii.gz"))

elif args.list_of_files is not None:

    all_fixed_files = input_args.list_of_files

# Define function =====

def run_secondlevel(fixed_files, task, con, group=None):

    n_maps = len(fixed_files)
    print(f" Found {n_maps} files for task {task}, contrast {con}")

    design_matrix = pd.DataFrame([1] * n_maps, columns=['Intercept'])

    slm = SecondLevelModel(minimize_memory=False)

    slm.fit(second_level_input=[nib.load(x) for x in fixed_files],
                                design_matrix=design_matrix)

    if group is not None:
        group_str = f"_group-{group}"
    else:
        group_str = ""

    save_glm_to_bids(
        slm,
        contrasts="Intercept",
        contrast_types={"Intercept": "t"},
        out_dir=f"{results_dir}",
        prefix=f"task-{task}_contrast-{con}{group_str}",
     )

# Main loop =====

# Identify contrasts
# (Sorry: This converts a PosixPath to string [.stem], extracts the contrast
#   and then removes the `contrast-` prefix)
contrasts = list(set([re.search('contrast-[^_]*', i.stem).group(0).replace("contrast-", "")
                     for i in all_fixed_files]))

tasks = list(set([re.search('task-[^_]*', i.stem).group(0).replace("task-", "")
                  for i in all_fixed_files]))

print(f"Tasks: {tasks}")
print(f"Contrasts: {contrasts}")

if len(tasks) > 1 and args.collapse_tasks:

    print(f"I found {len(tasks)} tasks, but collapsing them all")

    for con in contrasts:

        effect_maps_glob = f"**/sub-*_contrast-{con}_stat-effect_statmap.nii.gz"
        fixed_files = list(results_path.glob(effect_maps_glob))

        run_secondlevel(fixed_files=fixed_files, task="all", con=con)

else:

    for task in tasks:

        for con in contrasts:

            effect_maps_glob = f"**/sub-*_task-{task}_contrast-{con}_stat-effect_statmap.nii.gz"
            fixed_files = list(results_path.glob(effect_maps_glob))
            fixed_files.sort()

            if subs is not None:

                # pp.pprint(fixed_files)
                # pp.pprint(subs)
                n_orig = len(fixed_files)

                fixed_files = [x for x in fixed_files
                               if re.search("sub-[^_/]*", str(x)).group(0)
                               in subs]

                # for x in fixed_files:
                #     print(re.search("sub-[^_/]*", str(x)).group(0))

                print(f"Reduced list of files from {n_orig} to "
                      f"{len(fixed_files)}")

            run_secondlevel(fixed_files=fixed_files, task=task, con=con,
                            group=group_label)