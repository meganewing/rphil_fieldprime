## =========================================================================
## Manila clam acute stress survival analysis
## Kaplan-Meier curves + Cox Proportional Hazards, following the workflow in:
## https://ahuffmyer.github.io/posts/2026-06-03-analyzing-manchester-stress-survival.html
##
## Plates are read in 24-well format (rows A-D, cols 1-6); columns 1 & 6 are
## blank (no clam). 0 = alive, 1 = dead. Plate names look like "LR-CC1-40" or
## "LR-CC1-30" -> treatment code ("CC") + replicate ("1") + temperature ("40"/"30").
## =========================================================================

library(tidyverse)
library(survival)
library(ggsurvfit)
library(gtsummary)
library(cardx)
library(cowplot)
library(emmeans)

## ========================== USER CONFIGURATION ===========================
## Everything you're likely to change between datasets lives here.

input_file    <- "data/LR_manila-40C.csv"   # path to the raw plate-reader CSV
output_prefix <- "output/LR_manila_40C"       # prefix for all saved plots/tables

# Treatment design: which treatment codes exist, and how they map onto the
# 2x2 immune-priming x thermal-priming factorial design (1 = yes, 0 = no).
treatment_design <- tribble(
  ~treatment, ~immune_priming, ~thermal_priming,
  "CC",        0,               0,
  "EE",        0,               1,
  "CP",        1,               0,
  "EP",        1,               1
)

# Human-readable labels and plot colors, keyed by treatment code
treatment_labels <- c(
  CC = "Control",
  EE = "Thermal priming",
  CP = "Immune priming",
  EP = "Thermal + Immune priming"
)

treatment_colors <- c(
  CC = "#1f77b4",  # blue
  EE = "#d62728",  # red
  CP = "#2ca02c",  # green
  EP = "#ff7f0e"   # orange
)

## Sheet layout constants (only change these if your plate-reader template changes)
plate_start_cols  <- c(1, 9, 17, 25, 33, 41, 49, 57)  # 1-indexed column of each plate's name
data_well_offsets <- 2:5                               # offsets of the 4 real (non-blank) wells
hr_pattern        <- "^([0-9.]+)HR$"                   # e.g. "4HR", "29.5HR"
plate_pattern     <- "^LR-([A-Z]{2})([12])-[0-9]+$"    # e.g. "LR-CC1-40" -> treatment "CC", replicate "1"
## ===========================================================================


# ---- 1. Read the raw file exactly as it is laid out (no header) ----------
raw <- read.csv(
  input_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  colClasses = "character",
  fill = TRUE
)

# ---- 2. Walk the sheet and pull out every well reading, at every timepoint -
# Each plate block is 8 columns wide:
#   [plate name] [well1=blank] [well2] [well3] [well4] [well5] [well6=blank] [spacer]
# 8 plates are laid out side-by-side; blocks repeat down the sheet, one per
# time check (HR row + header row + rows A-D + blank spacer row).
records <- list()
i <- 1

for (r in seq_len(nrow(raw))) {
  first_cell <- trimws(raw[r, 1])
  m <- regmatches(first_cell, regexec(hr_pattern, first_cell))[[1]]
  if (length(m) == 0) next                # not a time-block header row -> skip
  if (r + 5 > nrow(raw)) next              # guard against a truncated trailing block
  
  hour <- as.numeric(m[2])
  header_row <- raw[r + 1, ]              # row with plate names ("LR-CC1-40", 1..6, ...)
  data_rows  <- raw[(r + 2):(r + 5), ]    # rows A, B, C, D of well statuses
  
  for (pc in plate_start_cols) {
    plate_name <- trimws(header_row[[pc]])
    pm <- regmatches(plate_name, regexec(plate_pattern, plate_name))[[1]]
    if (length(pm) == 0) next             # empty/blank block -> skip
    
    treatment <- pm[2]                    # e.g. CC, EE, CP, EP
    replicate <- pm[3]                    # 1 or 2
    
    for (row_idx in seq_len(nrow(data_rows))) {
      well_row <- LETTERS[row_idx]
      for (off in data_well_offsets) {
        val <- trimws(data_rows[row_idx, pc + off])
        if (val %in% c("0", "1")) {
          records[[i]] <- data.frame(
            hour = hour,
            treatment = treatment,
            replicate = replicate,
            well_id = paste0(well_row, off),   # e.g. "A2"
            status = as.integer(val)
          )
          i <- i + 1
        }
      }
    }
  }
}

well_data <- bind_rows(records)

stopifnot(
  "No data parsed - check input_file / plate_pattern / hr_pattern against the sheet layout" =
    nrow(well_data) > 0
)

# ---- 3. Collapse each individual clam's time series into one event record -
# individual_id uniquely identifies a physical clam: treatment + replicate + well position
# (well position is consistent for the same clam across all timepoints on a given plate)
well_data <- well_data %>%
  mutate(individual_id = paste(treatment, replicate, well_id, sep = "_"))

surv_data <- well_data %>%
  arrange(individual_id, hour) %>%
  group_by(individual_id, treatment, replicate) %>%
  summarise(
    # time of first observed death; NA if never observed dead
    death_hour = if (any(status == 1)) min(hour[status == 1]) else NA_real_,
    last_hour  = max(hour),
    .groups = "drop"
  ) %>%
  mutate(
    time   = if_else(!is.na(death_hour), death_hour, last_hour),
    status = if_else(!is.na(death_hour), 1L, 0L)   # 1 = died, 0 = censored (survived to last check)
  ) %>%
  select(individual_id, treatment, replicate, time, status) %>%
  left_join(treatment_design, by = "treatment") %>%
  mutate(
    treatment = factor(treatment, levels = names(treatment_labels)),
    treatment_label = recode(treatment, !!!treatment_labels),
    immune_priming  = factor(immune_priming,  levels = c(0, 1), labels = c("No", "Yes")),
    thermal_priming = factor(thermal_priming, levels = c(0, 1), labels = c("No", "Yes"))
  )

write_csv(surv_data, paste0(output_prefix, "_survival_individuals.csv"))

# ---- 4. Kaplan-Meier survival curves, by treatment -------------------------
km_fit <- survival::survfit(Surv(time, status) ~ treatment, data = surv_data)
str(km_fit)

my_theme <- theme_classic() +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    legend.position = "right"
  )

plot_km <- ggsurvfit::survfit2(Surv(time, status) ~ treatment, data = surv_data) %>%
  ggsurvfit::ggsurvfit(linewidth = 1.1) +
  add_confidence_interval() +
  add_risktable() +
  labs(
    x = "Hours",
    y = "Survival probability",
    title = paste("Manila clam survival —", output_prefix)
  ) +
  ylim(0, 1) +
  my_theme +
  scale_colour_manual(values = treatment_colors, labels = treatment_labels, name = "Treatment") +
  scale_fill_manual(values = treatment_colors, labels = treatment_labels, name = "Treatment")

plot_km
ggsave(paste0(output_prefix, "_KM_by_treatment.png"), plot = plot_km, width = 8, height = 6, dpi = 300)

# ---- 5. Cox Proportional Hazards models -------------------------------------

# Primary model for treatment comparisons:
# Replicate is included as an adjustment term so treatment comparisons account
# for the systematic difference between replicate plates.
#
# CC is the reference treatment, and replicate 1 is the reference replicate.
cox_treatment_rep_adj <- survival::coxph(
  Surv(time, status) ~ treatment + replicate,
  data = surv_data
)

summary(cox_treatment_rep_adj)
cox_treatment_rep_adj %>% tbl_regression(exp = TRUE)

# Pairwise treatment comparisons:
# This gives all six unique treatment comparisons:
# CC vs EE, CC vs CP, CC vs EP, EE vs CP, EE vs EP, CP vs EP.
#
# emmeans computes the estimated marginal means on the Cox model's log-hazard
# scale and then forms pairwise contrasts. exp = TRUE reports hazard ratios.
# Holm adjustment controls the family-wise error rate across the six comparisons.
if (!requireNamespace("emmeans", quietly = TRUE)) {
  stop("Package 'emmeans' is required for the pairwise treatment comparisons. ",
       "Install it with install.packages('emmeans').")
}

emm_treatment <- emmeans::emmeans(
  cox_treatment_rep_adj,
  ~ treatment
)

pairwise_treatment <- emmeans::contrast(
  emm_treatment,
  method = "pairwise",
  adjust = "holm"
)

pairwise_treatment_table <- as.data.frame(
  summary(pairwise_treatment, infer = c(TRUE, TRUE), type = "response")
) %>%
  dplyr::rename(
    comparison = contrast,
    HR = ratio,
    CI_low = asymp.LCL,
    CI_high = asymp.UCL,
    p_value_adj = p.value
  ) %>%
  dplyr::select(comparison, HR, CI_low, CI_high, p_value_adj)

print(pairwise_treatment_table)

# Save the six pairwise comparisons as a CSV
write_csv(
  pairwise_treatment_table,
  paste0(output_prefix, "_Cox_pairwise_treatment_comparisons.csv")
)

# (b) 2x2 factorial: immune priming x thermal priming main effects + interaction
# This model is retained as a separate analysis because it answers a different
# question: whether immune priming, thermal priming, and their interaction
# affect mortality risk.
cox_factorial <- survival::coxph(
  Surv(time, status) ~ immune_priming * thermal_priming,
  data = surv_data
)
summary(cox_factorial)
cox_factorial %>% tbl_regression(exp = TRUE)

# Replicate-adjusted factorial model:
# This accounts for the replicate effect while retaining the factorial design.
cox_factorial_rep_adj <- survival::coxph(
  Surv(time, status) ~ immune_priming * thermal_priming + replicate,
  data = surv_data
)
summary(cox_factorial_rep_adj)
cox_factorial_rep_adj %>% tbl_regression(exp = TRUE)

# (c) Treatment x replicate interaction model:
# Retained as a diagnostic model to test whether treatment effects differ
# between replicate plates. This is NOT the primary model for the six
# pairwise treatment comparisons above.
cox_treatment_rep_interaction <- survival::coxph(
  Surv(time, status) ~ treatment * replicate,
  data = surv_data
)
summary(cox_treatment_rep_interaction)
cox_treatment_rep_interaction %>% tbl_regression(exp = TRUE)

# ---- 6. Supplementary plot: survival by immune x thermal priming ----------
plot_factorial <- ggsurvfit::survfit2(
  Surv(time, status) ~ immune_priming + thermal_priming, data = surv_data
) %>%
  ggsurvfit::ggsurvfit(linewidth = 1.1) +
  labs(
    x = "Hours",
    y = "Survival probability",
    title = paste("Survival by immune priming x thermal priming —", output_prefix),
    color = "Immune, Thermal priming"
  ) +
  ylim(0, 1) +
  my_theme

plot_factorial
ggsave(paste0(output_prefix, "_KM_by_factorial_design.png"), plot = plot_factorial, width = 9, height = 6, dpi = 300)

# ---- 7. Supplementary plot: survival by replicate within each treatment ---
plot_rep <- ggsurvfit::survfit2(Surv(time, status) ~ treatment + replicate, data = surv_data) %>%
  ggsurvfit::ggsurvfit(linewidth = 0.9) +
  labs(
    x = "Hours",
    y = "Survival probability",
    title = paste("Survival by treatment and replicate plate —", output_prefix)
  ) +
  ylim(0, 1) +
  my_theme

plot_rep
ggsave(paste0(output_prefix, "_KM_by_treatment_replicate.png"), plot = plot_rep, width = 9, height = 6, dpi = 300)