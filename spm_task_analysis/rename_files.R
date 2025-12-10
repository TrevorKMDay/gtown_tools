suppressMessages(library(tidyverse))
library(optparse)

# Yes, I know this isn't the right tool for this, but MATLAB was crashing
# every time I asked it to open an SPM.mat file.

option_list <- list(

    make_option(c("-a", "--all_only"), action = "store_true", default = FALSE,
                help = "Only copy files from 'All Sessions'"),

    make_option(c("-k", "--keep_con_all"), action = "store_true",
                default = FALSE,
                help = "Don't skip con-all files."),

    make_option(c("-t", "--spmT"), action = "store_true", default = TRUE,
                help = "Copy and rename spmT files"),

    make_option(c("-b", "--beta"), action = "store_true", default = FALSE,
                help = "Copy and rename beta files"),

    make_option(c("-c", "--con"), action = "store_true", default = FALSE,
                help = "Copy and rename con files")

  )

args <- parse_args2(OptionParser(option_list = option_list))
print(args)

all_only <- args$options$all_only
keep_con_all <- args$options$keep_con_all

if (length(args$args) < 3) {
  stop("Not enough arguments, supply input, task, output")
}

dir <- args$args[1]
output <- args$args[3]

# If task name supplied w/o 'task-', prepend it.
task <- args$args[2]
if (!str_detect(task, "task-"))
  task <- str_glue("task-{task}")

message(str_glue("Copying files from\n    {dir}\n    ->\n    {output}"))

# dir <- "~/code/spm_task_analysis/foo/sub-TDCh121/func/ ...
#         sub-TDCh121_task-easy_run-1_desc-smoothedmasked/results/optcens/ ...
#         01_censored"

# dir <- "/Volumes/thufir/kLat_SPM/sub-kLat002/func/task-ADDT_postOCT"

# Get data ====

spm <- suppressMessages(R.matlab::readMat(paste0(dir, "/SPM.mat"))$SPM)

xcon_loc <- which(rownames(spm) == "xCon")
vbeta_loc <- which(rownames(spm) == "Vbeta")

sub <- str_extract(dir, "sub-[^_/]*")

if (sub == "") {
  sub <- "sub-UNK"
  message("Can't find subject label, using 'sub-UNK'")
}

# Contrast names ====

xCon <- spm[[xcon_loc]]

format_con_name <- function(contrast) {

  parts <- str_split_1(contrast, ">")
  parts[1] <- tolower(parts[1])

  if (parts[2] == "Rest") {
    result <- parts[1]
  } else {
    parts[2] <- str_to_title(parts[2])
    result <- str_glue("{parts[1]}Minus{parts[2]}")
  }

  return(result)

}

# Extract names from SPM object, and replace '>' with 'Minus'
xCon_names <- tibble(str = sapply(1:dim(xCon)[3],
                                  function(x) xCon[,,x]$name)) %>%
  separate_wider_delim(str, delim = " - ",
                       names = c("contrast", "session")) %>%
  mutate(
    contrast = Vectorize(format_con_name)(contrast),
    session = tolower(str_remove(session, " ?Session[s ]?"))
  )

# Rename con/spmT files ====

timestamp <- format(now(), "%y%m%d%H%M%S")
logfile <- paste0(output, "/", sub, "_date-", timestamp, "_logfile.txt")
dir.create(output, showWarnings = FALSE)
file.create(logfile)

con_files <- list.files(dir, pattern = "con_.*.nii", full.names = TRUE)
spmT_files <- list.files(dir, pattern = "spmT_.*.nii", full.names = TRUE)

stopifnot(length(con_files) > 1)
stopifnot(length(spmT_files) > 1)

dir.create(paste0(output), recursive = TRUE, showWarnings = FALSE)

for (i in 1:length(con_files)) {

  # Skip this one if all_only is set and this is the all session
  ses <- xCon_names$session[i]
  if (all_only & ses != "all") next

  # If this is the 'all>rest' contrast, and keep_con_all isn't set, skip it
  con <- xCon_names$contrast[i]
  if (con == "all" & !keep_con_all) next

  if (args$options$con) {

    con_src <- con_files[i]

    con_new_name <- str_glue("{sub}_{task}_con-{con}_",
                             "ses-{ses}_stat-con_",
                             "statmap.nii")

    con_dest <- str_glue("{output}/{con_new_name}")

    file.copy(con_src, con_dest, overwrite = FALSE)

    cat(paste0(con_src, " -> ", con_new_name),
        file = logfile, fill = TRUE, append = TRUE)

  }

  if (args$options$spmT) {

    spmT_src <- spmT_files[i]

    spmT_new_name <- str_glue("{sub}_{task}_con-{con}_",
                             "ses-{ses}_stat-t_",
                             "statmap.nii")

    spmT_dest <- paste0(str_glue("{output}/{spmT_new_name}"))
    file.copy(spmT_src, spmT_dest, overwrite = FALSE)

    cat(paste0(spmT_src, " -> ", spmT_new_name),
        file = logfile, fill = TRUE, append = TRUE)

  }

}

# Rename beta files ====

# This needs to be refactored too

if (args$options$beta) {

  ## Beta names ====

  Vbeta <- spm[[vbeta_loc]]

  # Extract names from SPM object, removing extra detail
  Vbeta_names <- sapply(1:dim(Vbeta)[3], function(x) Vbeta[,,x]$descrip) %>%
    str_remove("spm_spm:beta [(][0-9]*[)] - Sn[(]") %>%
    str_replace("[)] ", "b") %>%
    str_remove("[*]bf[(]1[)]")

  # Do the copying

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

}

