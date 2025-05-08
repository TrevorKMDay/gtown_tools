suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)

# args <- c("/Volumes/thufir/kLat_fmriprep/", 0.9, 0.2, ".")

fmriprep <- args[1]

settings <- tibble(fd_thresh = as.numeric(args[2]),
                   fd_pct = as.numeric(args[3]),
                   datetime = now())

out <- args[4]


# Get confound data =====

cfiles <- list.files(fmriprep, pattern = "_desc-confounds_timeseries.tsv",
                     recursive = TRUE, full.names = TRUE)

confounds <- tibble(f = cfiles) %>%
  mutate(
    dat = map(f, read_tsv, show_col_types = FALSE, na = "n/a",
              .progress = TRUE)
  )

confounds2 <- confounds %>%
  mutate(
    bn = basename(f),
    n_frames = map_int(dat, nrow),
    mean_FD = map_dbl(dat, ~mean(.x$framewise_displacement, na.rm = TRUE)),
    n_FD_outlier = map_dbl(dat, ~sum(.x$framewise_displacement > settings$fd_thresh,
                               na.rm = TRUE)),
    pct_FD_outlier = n_FD_outlier / n_frames,
    exclude = pct_FD_outlier > settings$fd_pct[1]
  ) %>%
  separate_wider_delim(bn, delim = "_", names = c("sub", "task", "run", NA, NA),
                        cols_remove = FALSE)

# Save table of values
confounds2 %>%
  mutate(
    across(where(is.numeric), ~round(.x, 3))
  ) %>%
  write_csv(paste0(out, "/run_motion.csv"))

confounds2 %>%
  group_by(sub, task) %>%
  summarize(
    n = n(),
    included = sum(!exclude)
  ) %>%
  write_csv(paste0(out, "/sub_inclusion.csv"))

# Plots =====

png(paste0(out, "/mean_FD.png"), width = 400, height = 300, units = "px")

ggplot(confounds2, aes(x = mean_FD, fill = task)) +
  geom_histogram(binwidth = 0.2, boundary = 0) +
  geom_vline(xintercept = settings$fd_thresh, color = "red") +
  theme_bw() +
  facet_wrap(vars(task)) +
  labs(x = "Mean FD (mm)", fill = "# frames")

dev.off()

png(paste0(out, "/pct_FD_outlier.png"), width = 400, height = 300,
    units = "px")

ggplot(confounds2, aes(x = pct_FD_outlier, fill = task)) +
  geom_histogram(binwidth = 0.05, boundary = 0) +
  geom_vline(xintercept = settings$fd_pct, color = "red") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_bw() +
  facet_wrap(vars(task)) +
  labs(x = "% FD outlier", fill = "# frames")

dev.off()


png(paste0(out, "/FD_scatterplot.png"), width = 400, height = 300,
    units = "px")

ggplot(confounds2, aes(x = pct_FD_outlier, y = mean_FD, color = task)) +
  geom_vline(xintercept = settings$fd_pct, color = "red") +
  geom_hline(yintercept = settings$fd_thresh, color = "red") +
  geom_point(alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_bw() +
  labs(x = "% FD outlier", y = "Mean FD (mm)", fill = "# frames")

dev.off()

png(paste0(out, "/FD_boxplot.png"), width = 400, height = 300,
    units = "px")

ggplot(confounds2, aes(x = task, y = pct_FD_outlier, fill = task)) +
  geom_boxplot(notch = TRUE) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

dev.off()


message()
message("Effect of task on % outlier")
message("---------------------------")

aov(pct_FD_outlier ~ task, data = confounds2) %>%
  anova()

message()
message("Effect of task on mean FD")
message("-------------------------")

aov(mean_FD ~ task, data = confounds2) %>%
  anova()

message()
message("Summary by task")
message("---------------")

confounds2 %>%
  group_by(task) %>%
  summarize(
    n = n(),
    excluded = sum(exclude),
    mean_FD = mean(mean_FD),
    pct_FD_outlier = mean(pct_FD_outlier)
  ) %>%
  mutate(
    pct_excl = excluded / n
  )

# Number of participants with n runs of each task
confounds_summary <- confounds2 %>%
  group_by(sub, task) %>%
  filter(
    !exclude
  ) %>%
  summarize(
    n = n()
  ) %>%
  group_by(task, n) %>%
  summarize(
    n_subs = n()
  )

# Create list of runs to keep
runs_to_keep <- confounds2 %>%
  filter(
    !exclude
  ) %>%
  select(-exclude) %>%
  group_by(sub, task) %>%
  reframe(
    runs = list(run)
  ) %>%
  arrange(sub, task)

runs_to_exclude <- confounds2 %>%
  filter(
    exclude
  ) %>%
  select(-exclude) %>%
  group_by(sub, task) %>%
  reframe(
    runs = list(run)
  ) %>%
  arrange(sub, task)

# Summary statistics

confounds_fd <- confounds2 %>%
  filter(
    pct_FD_outlier <= 0.2
  )

mean_surv_fd <- weighted.mean(confounds_fd$mean_FD, confounds_fd$n_frames)
sd_surv_fd <- Hmisc::wtd.var(confounds_fd$mean_FD, confounds_fd$n_frames) %>%
  sqrt()

message(paste0("Mean surviving FD (all scans): ",
               round(mean_surv_fd, 3), " ",
               "(", round(sd_surv_fd, 3), ")"))

# Write out decisions ====

list("settings" = settings, "motion" = runs_to_keep) %>%
  toJSON(pretty = TRUE) %>%
  write_lines(paste0(out, "/included_runs.json"))

list("settings" = settings, "motion" = runs_to_exclude) %>%
  toJSON(pretty = TRUE) %>%
  write_lines(paste0(out, "/excluded_runs.json"))
