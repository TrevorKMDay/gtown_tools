suppressMessages(library(tidyverse))

# Yes, I know this isn't the right tool for this, but MATLAB was crashing
# every time I asked it to open an SPM.mat file.

args <- commandArgs(trailingOnly = TRUE)
# print(length(args))

if (length(args) != 3) {
  stop("Not enough arguments, supply input and output")
}

dir <- args[1]
task <- paste0("task-", args[2])
output <- args[3]

message(paste("Moving files from", dir, "to", output))

# dir <- "~/code/spm_task_analysis/foo/sub-TDCh121/func/ ...
#         sub-TDCh121_task-easy_run-1_desc-smoothedmasked/results/optcens/ ...
#         01_censored"

# Get data
spm <- suppressMessages(R.matlab::readMat(paste0(dir, "/SPM.mat"))$SPM)

xcon_loc <- which(rownames(spm) == "xCon")
vbeta_loc <- which(rownames(spm) == "Vbeta")

# Contrast names ====

xCon <- spm[[xcon_loc]] 

# Extract names from SPM object, and replace '>' with '_m_' after removing
#   pre-existing underscores.
xCon_names <- sapply(1:dim(xCon)[3], function(x) xCon[,,x]$name) %>%
  str_remove(" - .*") %>%
  str_remove_all("_") %>%
  tolower() %>%
  str_replace(">", "M")

# Beta names ====

Vbeta <- spm[[vbeta_loc]]

# Extract names from SPM object, removing extra detail
Vbeta_names <- sapply(1:dim(Vbeta)[3], function(x) Vbeta[,,x]$descrip) %>%
  str_remove("spm_spm:beta [(][0-9]*[)] - Sn[(]") %>%
  str_replace("[)] ", "b") %>%
  str_remove("[*]bf[(]1[)]")

# Rename con/spmT files ====

sub <- str_extract(dir, "sub-[^_/]*")

timestamp <- format(now(), "%y%m%d%H%M%S")
logfile <- paste0(output, "/", sub, "_date-", timestamp, "_logfile.txt")
dir.create(output, showWarnings = FALSE)
file.create(logfile)

con_files <- list.files(dir, pattern = "con_.*.nii", full.names = TRUE)
spmT_files <- list.files(dir, pattern = "spmT_.*.nii", full.names = TRUE)

stopifnot(length(con_files) > 1)
stopifnot(length(spmT_files) > 1)

dir.create(paste0(output, "/", sub, "/"), recursive = TRUE, 
           showWarnings = FALSE)

for (i in 1:length(con_files)) {
  
  con_src <- con_files[i]
  spmT_src <- spmT_files[i]
  
  con_new_name <- paste0(sub, "_", task, "_con-", xCon_names[i], 
                          "_statmap.nii")
  spmT_new_name <- paste0(sub, "_", task, "_con-", xCon_names[i], 
                          "_stat-t_statmap.nii")
  
  con_dest <- paste0(output, "/", sub, "/", con_new_name)
  spmT_dest <- paste0(output, "/", sub, "/", spmT_new_name)
  
  file.copy(con_src, con_dest, overwrite = FALSE)
  file.copy(spmT_src, spmT_dest, overwrite = FALSE)
  
  cat(paste0(con_src, " -> ", "_con-", xCon_names[i], "_statmap.nii"), 
      file = logfile, fill = TRUE, append = TRUE)
  cat(paste0(spmT_src, " -> ", "_con-", xCon_names[i], "_stat-t_statmap.nii"), 
      file = logfile, fill = TRUE, append = TRUE)
  
}

# Rename beta files ====

beta_files <- list.files(dir, pattern = "beta_.*.nii", full.names = TRUE)

stopifnot(length(beta_files) > 1)

for (i in 1:length(beta_files)) {
  
  beta_src <- beta_files[i]
  
  new_name <- paste0(sub, "_", task, "_con-", Vbeta_names[i], 
                     "_stat-beta_statmap.nii")
  
  beta_dest <- paste0(output, "/", sub, "/", new_name)
  
  file.copy(beta_src, beta_dest, overwrite = FALSE)
  
  cat(paste0(beta_src, " -> ", "_con-", Vbeta_names[i], 
             "_stat-beta_statmap.nii"), 
      file = logfile, fill = TRUE, append = TRUE)
  
}

