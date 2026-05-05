library(tidyverse)
library(R.matlab)

# all_spmT_files <- list.files("/Volumes/thufir/kLat_SPM/",
#                               pattern = "SPM.mat", recursive = TRUE,
#                               full.names = TRUE)
#
# optcens_spmT_files <- str_subset(all_spmT_files, "optcens/01_censored/")

source("~/code/preproc_spm/get_cens_frames.R")

args <- commandArgs(trailingOnly = TRUE)
# print(args)

if (length(args) > 0) {

  for (file in args) {

    output_file <- str_replace(file, "SPM.mat", "cens_frames.csv")

    # If output file doesn't exist, or is older than the file to extract
    #   number of frames from, get the number of censored frames.
    if (!file.exists(output_file) |
        file.info(file)$mtime > file.info(output_file)$mtime) {

      result <- get_cens_frames(readMat(file)$SPM)
      write_csv(result, output_file)
      message(str_glue("Wrote {output_file}"))

    } else {

      message(str_glue("{output_file} is newer than {file}, not rerunning ",
                      "calculation."))

    }

  }

}

# get_cens_frames(readMat(optcens_spmT_files[1])$SPM)

# test <- tibble(f = optcens_spmT_files) %>%
#   mutate(
#     cens_frames = map(f, ~get_cens_frames(readMat(.x)$SPM),
#                       .progress = TRUE)
#   )
#
# test2 <- test %>%
#   unnest(cens_frames) %>%
#   mutate(
#     pct = signif(n_cens / n_frames * 100, 3)
#   )
