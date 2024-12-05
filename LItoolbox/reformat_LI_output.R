suppressPackageStartupMessages(library(tidyverse))

inputs1 <- commandArgs(trailingOnly = TRUE)
inputs <- str_subset(inputs1, ".*[.]tsv")

# Check that at least one TSV file was supplied
stopifnot(length(inputs) > 0)

rewrite <- function(input) {
  
  message(input)

  data0 <- read_delim(input, "\t", show_col_types = FALSE)

  new_data <- data0 %>%
    select(-any_of(c("...1"))) %>%
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

  output_file <- str_replace(input, ".tsv", "-reformatted.csv")
  message(output_file)

  stopifnot(input != output_file)

  write_csv(new_data, output_file)

}

for (i in inputs) { rewrite(i) }
