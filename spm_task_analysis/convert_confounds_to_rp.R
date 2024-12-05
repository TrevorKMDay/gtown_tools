suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)

# Confounds filepath
confounds_fp <- args[1]
output_file <- args[2]

# trans in mm, rot in rad
# confounds_fp <- "~/Projects/LBS/ferrara_lbs/code/analysis_spm/fmriprep_input/sub-TDCh121/func/sub-TDCh121_task-easy_run-1_desc-confounds_timeseries.tsv"

# Read in confound data
confounds <- read_tsv(confounds_fp, na = "n/a", show_col_types = FALSE) %>%
  select(matches("trans_[xyz]$"), matches("rot_[xyz]$")) %>%
  mutate(
    # For some reason, fMRIPREP has these oppositely signed
    across(c(trans_x, rot_z), ~ -1 * .x)
  )

write.table(confounds, output_file, row.names = FALSE, col.names = FALSE)