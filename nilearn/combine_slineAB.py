import argparse as ap
import glob
import pprint as pp
import re
import pandas as pd
import os
import nibabel as nib

from nilearn.glm.second_level import SecondLevelModel

parser = ap.ArgumentParser()

parser.add_argument("directory", type=str, help="Directory to work in.")

parser.add_argument("task_regex", type=str, help="Regex to match tasks.")

args = parser.parse_args()

directory = args.directory
regex = args.task_regex

# Run combine

# from nilearn.datasets import fetch_localizer_contrasts

# n_subjects = 16
# data = fetch_localizer_contrasts(
#     ["left vs right button press"],
#     n_subjects,
#     get_tmaps=True,
#     legacy_format=False,
# )

# print(data)

t_maps = glob.glob(f"{directory}/"
                   f"sub-*_task-{regex}_contrast-*_stat-t_statmap.nii.gz")

# Extract the contrast names
contrasts = set([re.search("contrast-[^_]*", x)[0].replace("contrast-", "") \
                 for x in t_maps])

print(f"Found contrasts: {contrasts}")

for contrast in contrasts: 

    # Subset files; make sure contrast is terminated, eg. contrast-xx_ to a
    #   avoid contrast names bing detected at the start of 'xMinusY' strings
    second_level_input = [f for f in t_maps if contrast + "_" in f]

    # pp.pprint(second_level_input)

    design_matrix = pd.DataFrame(
        [1] * len(second_level_input),
        columns=["intercept"],
    )

    second_level_model = SecondLevelModel(smoothing_fwhm=8.0, n_jobs=2)
    second_level_model = second_level_model.fit(
        second_level_input,
        design_matrix=design_matrix,
    )

    # Create maps

    results = second_level_model.compute_contrast(
        second_level_contrast="intercept",
        output_type="all",
    )

    for result in results.keys():

        # The imaging data to save
        data = results[result]

        if result == "effect_size":
            label = "effect"
        elif result == "effect_variance":
            label = "variance"
        elif result == "p_value":
            label = "p"
        elif result == "stat":
            label = "t"
        elif result == "z_score":
            label = "z"

        bn = os.path.basename(second_level_input[0])
        sub = re.search("sub-[^_]*", bn)[0]

        fname = f"{sub}_contrast-{contrast}_stat-{label}_statmap.nii.gz"

        nib.save(data, f"{directory}/{fname}")

        print(f"Wrote {fname}")

        