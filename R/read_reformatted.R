read_reformatted <- function(refmt, names, roi = NA) {

  x <- refmt %>%
    select(input_image, li_wm) %>%
    separate_wider_delim(input_image, delim = "_", names = names,
                         too_few = "align_start", too_many = "merge") %>%
    mutate(
      across(where(is.character), ~str_remove(., "^[^-]*-")),
      roi = roi,
      li_wm = -1 * as.numeric(li_wm),
      contrast = str_replace(contrast, "Minus", "-") %>%
        str_replace("Plus", "+")
    )

  return(x)

}
