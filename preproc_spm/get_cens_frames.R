require(tidyverse)

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
