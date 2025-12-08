suppressPackageStartupMessages(library(tidyverse))

inputs1 <- commandArgs(trailingOnly = TRUE)
inputs <- str_subset(inputs1, ".*[.]tsv")

# Check that at least one TSV file was supplied
stopifnot(length(inputs) > 0)

extract_fields <- function(filenames) {

  filenames <- sapply(filenames, function(x) str_split_1(x, " ")[1])

  fields <- str_extract_all(filenames, "[^_]+-[^_]+", simplify = TRUE) %>%
    str_remove("-[^_]+$") %>%
    as.vector() %>%
    unique()

  # Remove empty string from list
  fields <- fields[!(fields %in% c(""))]

  n_fields <- length(fields)

  message(paste("Found", n_fields, "fields:", paste(fields, collapse = " ")))

  return(fields)

}

rewrite <- function(data0, fname, fields = NA) {

  # message(input)
  # data0 <- read_delim(input, "\t", show_col_types = FALSE)

  # Extract failed image and ROI from error string
  errors <- data0 %>%
    filter(
      str_detect(`Input image`, "An error occurred")
    ) %>%
    select(`Input image`) %>%
    mutate(
      `Input image` = str_remove(`Input image`,
                                 "An error occurred when processing ") %>%
        str_remove("in combination with ") %>%
        str_remove(" - skipping..."),
      desc = "ERROR",
    ) %>%
    separate_wider_delim(`Input image`, delim = " ",
                         names = c("input_image", "inclusive_mask"))

  new_data <- data0 %>%
    select(-any_of(c("...1"))) %>%
    filter(
      !str_detect(`Input image`, "An error occurred")
    ) %>%
    select_all(~tolower(.) %>%
                str_replace_all(" ", "_") %>%
                str_remove_all("[^a-z_]")) %>%
    separate_wider_delim(input_image, delim = " ",
                        names = c("input_image", "desc"),
                        too_many = "merge", too_few = "align_start") %>%
    mutate(
      contrast_index = str_extract(desc, "contrast [0-9]*") %>%
        str_remove("contrast ") %>%
        as.numeric(),
      contrast_name = str_extract(desc, "[^ ]*>[^ ]*"),
      desc = str_remove_all(desc, "contrast [0-9*]: .[^ ]*>[^ ]* - |^[(]|[)]$"),
      across(starts_with("li"), ~replace(.x, is.nan(.x), NA))
    ) %>%
    select(source_path, input_image, starts_with("contrast"), everything())

  result <- bind_rows(new_data, errors) %>%
    arrange(input_image)

  if (length(fields) > 1 & !all(is.na(fields))) {

    message("Extracting fields from input_image")

    for (f in fields) {

      keys <- str_extract(result$input_image, paste0(f, "-[^_]*")) %>%
        str_remove(paste0(f, "-"))

      result[, f] <- keys

    }

    result <- result %>%
      select(all_of(fields), everything(), -input_image)

  }

  output_file <- str_replace(fname, ".tsv", "-reformatted.csv")
  message(paste("Writing to", output_file))

  stopifnot(fname != output_file)

  write_csv(result, output_file)

}

for (i in inputs) {

  tsv <- read_tsv(i, show_col_types = FALSE)
  fields <- extract_fields(tsv$`Input image`)

  # print(tsv$`Input image`)

  rewrite(data0 = tsv, fname = i, fields = fields)

}
