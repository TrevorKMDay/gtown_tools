library(tidyverse)
library(hdf5r)

# all_spmT_files <- list.files("/Volumes/thufir/kLat_SPM/",
#                               pattern = "SPM.mat", recursive = TRUE,
#                               full.names = TRUE)
#
# optcens_spmT_files <- str_subset(all_spmT_files, "optcens/01_censored/")

get_cens_frames <- function(SPM) {

  xX_names <- unlist(SPM[[8]][[7]])
  n_scan <- SPM[[3]][1, 1]

  basis_functions <- str_detect(xX_names, "[*]bf")
  rotations <- str_detect(xX_names, "R[1-6]$")

  cols_to_keep <- !(basis_functions | rotations)
  cols_to_keep[length(cols_to_keep)] <- FALSE

  result <- tribble(
    ~n_frames, ~n_cens,
    n_scan, sum(cols_to_keep)
  )

  return(result)

}

args <- commandArgs(trailingOnly = TRUE)

for (file in args) {

  result <- get_cens_frames(readMat(file)$SPM)

  output_file <- str_replace(file, "SPM.mat", "cens_frames.csv")

  write_csv(result, output_file)
  message(str_glue("Wrote {output_file}"))

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
