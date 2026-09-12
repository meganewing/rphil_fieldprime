# =============================================================================
# Resazurin Metabolic Assay Analysis – APC Cockle Thermal Stress Experiment
# =============================================================================
# Treatments:
#   CC – control            EE – elevated temperature priming
#   CP – immune priming     EP – combined temperature + immune priming
#
# Exposure: 11°C for 4 h → 26°C for 4 h
# Fluorescence (528/20 ex, 590/20 em) read every ~20–40 min (T0–T8).
#
# Normalisation pipeline (Huffmyer et al.):
#   1. Fold-change relative to each well's own T0 (applied to samples & blanks)
#   2. Subtract mean blank fold-change per plate × timepoint
#   3. Divide by animal area (mm²) → size-normalised metabolism
#   4. Trapezoid AUC over real elapsed hours (parsed from file headers)
#
# Inputs:  data/AP_cockle_final/APC-*-T*.txt   (plate exports, T0–T8)
#          data/AP_cockle_final/APC_PLATES.csv  (layout + area mm²)
# Outputs: output/APC_cockle/
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pracma)       # trapz()
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(cowplot)
})

# ── Directories ---------------------------------------------------------------
data_dir <- file.path("data", "AP_cockle_final")
out_dir  <- file.path("output", "APC_cockle")
fig_dir  <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Temperature to analyse ----------------------------------------------------
# Set to 11 or 26. Only plates from this exposure temperature will be included.
TARGET_TEMP <- 26

# ── Colour palette (Okabe-Ito, colourblind-safe) ------------------------------
treatment_levels <- c("CC", "EE", "CP", "EP")
okabe_ito <- c(CC = "#009E73", EE = "#E69F00", CP = "#56B4E9", EP = "#D55E00")

# =============================================================================
# Helper functions
# =============================================================================

normalize_well_id <- function(x) {
  x <- toupper(trimws(x))
  m <- str_match(x, "^([A-Z]+)0*([0-9]+)$")
  paste0(m[, 2], as.integer(m[, 3]))
}

# Parse the T-index from filename  e.g. APC-CC1-26-T3.txt -> 3
parse_t_index <- function(path) {
  hit <- str_match(basename(path), "(?i)-T([0-9]+)\\.txt$")
  as.integer(hit[, 2])
}

# Parse plate ID  e.g. APC-CC1-26-T3.txt -> full id "APC-CC1-26"
# Also extracts:
#   base_plate_id : "APC-CC1"  (treatment + replicate number, no temp)
#   temp_c        : 26         (exposure temperature)
#
# The same physical individuals appear in both the 11°C and 26°C files.
# Size normalisation uses base_plate_id + well_id so that a single area
# value (measured once per individual) applies to reads at both temperatures.
parse_plate_id <- function(path) {
  hit <- str_match(basename(path),
    "(?i)^(APC-([A-Za-z0-9]+)-([0-9]+))-T[0-9]+\\.txt$")
  # hit[,2] = full plate id (e.g. APC-CC1-26)
  # hit[,3] = treatment+rep  (e.g. CC1)
  # hit[,4] = temperature    (e.g. 26)
  list(
    plate_id      = ifelse(is.na(hit[, 2]), "unknown", hit[, 2]),
    base_plate_id = ifelse(is.na(hit[, 3]), "unknown",
                           paste0("APC-", hit[, 3])),
    temp_c        = suppressWarnings(as.integer(hit[, 4]))
  )
}

# Parse wall-clock datetime from file header  "Date\t8/20/2026\nTime\t3:36:43 PM"
parse_wall_time <- function(lines) {
  date_line <- lines[str_detect(lines, "^Date\\t")]
  time_line <- lines[str_detect(lines, "^Time\\t")]
  if (length(date_line) == 0 || length(time_line) == 0) return(NA_real_)
  date_str <- str_split(date_line[1], "\\t")[[1]][2] |> trimws()
  time_str <- str_split(time_line[1], "\\t")[[1]][2] |> trimws()
  dt <- tryCatch(
    as.POSIXct(paste(date_str, time_str),
               format = "%m/%d/%Y %I:%M:%S %p", tz = "UTC"),
    error = function(e) NA_POSIXct_
  )
  as.numeric(dt)   # seconds since epoch
}

extract_results_block <- function(lines) {
  idx <- which(trimws(lines) == "Results")
  if (length(idx) == 0) stop("No 'Results' section found")
  idx <- idx[1]
  header_tokens <- str_split(lines[idx + 1], "\\t")[[1]] |> trimws()
  col_ids <- header_tokens[header_tokens != "" &
                             str_detect(header_tokens, "^[0-9]+$")]
  j <- idx + 2
  data_lines <- character()
  while (j <= length(lines)) {
    line <- lines[j]
    if (trimws(line) == "") break
    if (!str_detect(line, "^[A-Za-z]\\t")) break
    data_lines <- c(data_lines, line)
    j <- j + 1
  }
  list(col_ids = col_ids, data_lines = data_lines)
}

parse_plate_export <- function(path) {
  lines      <- readLines(path, warn = FALSE)
  wall_sec   <- parse_wall_time(lines)
  res        <- extract_results_block(lines)

  ids <- parse_plate_id(path)

  map_dfr(res$data_lines, function(line) {
    tokens <- str_split(line, "\\t")[[1]] |> trimws()
    tokens <- tokens[tokens != ""]
    row_letter <- tokens[1]
    nums <- suppressWarnings(as.numeric(tokens[-1]))
    valid_idx <- which(!is.na(nums))
    if (length(valid_idx) == 0) return(tibble())
    vals <- nums[valid_idx]
    n    <- min(length(vals), length(res$col_ids))
    tibble(
      row_id   = toupper(row_letter),
      col_id   = as.integer(res$col_ids[seq_len(n)]),
      well_id  = normalize_well_id(
        paste0(toupper(row_letter), res$col_ids[seq_len(n)])),
      value    = vals[seq_len(n)]
    )
  }) %>%
    mutate(
      plate_id      = ids$plate_id,       # e.g. APC-CC1-26
      base_plate_id = ids$base_plate_id,  # e.g. APC-CC1 (links 11°C & 26°C)
      temp_c        = ids$temp_c,         # 11 or 26
      t_index       = parse_t_index(path),
      wall_sec      = wall_sec
    )
}

trapezoid_auc <- function(time_hr, value) {
  ok <- is.finite(time_hr) & is.finite(value)
  t  <- time_hr[ok]; v <- value[ok]
  if (length(t) < 2) return(NA_real_)
  ord <- order(t); t <- t[ord]; v <- v[ord]
  sum(diff(t) * (head(v, -1) + tail(v, -1)) / 2)
}

sig_label <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}

add_sig_brackets <- function(p, pairs_df, group_levels, y_vals) {
  sig_pairs <- pairs_df %>% mutate(label = sig_label(p.value)) %>%
    filter(label != "ns")
  if (nrow(sig_pairs) == 0) return(p)
  y_max  <- max(y_vals, na.rm = TRUE)
  y_rng  <- diff(range(y_vals, na.rm = TRUE))
  step   <- y_rng * 0.14
  for (i in seq_len(nrow(sig_pairs))) {
    parts <- str_split(as.character(sig_pairs$contrast[i]), " - ", 2)[[1]]
    x1 <- match(trimws(parts[1]), group_levels)
    x2 <- match(trimws(parts[2]), group_levels)
    if (is.na(x1) || is.na(x2)) next
    bar_y <- y_max + i * step
    p <- p +
      annotate("segment", x = x1, xend = x2, y = bar_y, yend = bar_y,
               colour = "black", linewidth = 0.6) +
      annotate("segment", x = x1, xend = x1,
               y = bar_y, yend = bar_y - step * 0.3,
               colour = "black", linewidth = 0.6) +
      annotate("segment", x = x2, xend = x2,
               y = bar_y, yend = bar_y - step * 0.3,
               colour = "black", linewidth = 0.6) +
      annotate("text", x = (x1 + x2) / 2, y = bar_y + step * 0.15,
               label = sig_pairs$label[i], size = 4.5)
  }
  p
}

# =============================================================================
# 1. Load all plate export files (T0–T8)
# =============================================================================
message("\n── 1. Loading plate files (T0–T8) ──────────────────────────────────")

plate_files <- list.files(data_dir,
  pattern = "(?i)^APC-.*-T[0-9]+\\.txt$",
  full.names = TRUE)

if (length(plate_files) == 0)
  stop("No plate files found in: ", data_dir)

plate_raw <- map_dfr(plate_files, function(path) {
  tryCatch(parse_plate_export(path),
           error = function(e) {
             message("  Parse error: ", basename(path), " — ", e$message)
             tibble()
           })
})

message("  Loaded ", nrow(plate_raw), " rows | ",
        n_distinct(plate_raw$plate_id), " plates | ",
        n_distinct(plate_raw$t_index), " timepoints (T",
        min(plate_raw$t_index), "–T", max(plate_raw$t_index), ")")

# =============================================================================
# 2. Convert wall-clock times to elapsed hours from each plate's T0
# =============================================================================
message("\n── 2. Computing elapsed hours from T0 ───────────────────────────────")

# Anchor elapsed time to each plate's own T0 read (separately for 11°C and
# 26°C plates — they are run on different days/times).
# base_plate_id links the same individuals across both temperature phases,
# but time_hr is measured from each phase's own T0.
t0_wall <- plate_raw %>%
  group_by(plate_id) %>%
  filter(t_index == min(t_index)) %>%
  summarise(t0_wall_sec = first(wall_sec), .groups = "drop")

plate_raw <- plate_raw %>%
  left_join(t0_wall, by = "plate_id") %>%
  mutate(time_hr = (wall_sec - t0_wall_sec) / 3600)

# Sanity check: print elapsed times per plate
time_check <- plate_raw %>%
  select(plate_id, t_index, time_hr) %>%
  distinct() %>%
  arrange(plate_id, t_index)

message("  Elapsed time schedule (first plate shown):")
print(time_check %>% filter(plate_id == first(plate_id)))

# =============================================================================
# 3. Load layout / size file
# =============================================================================
message("\n── 3. Loading layout / size file ───────────────────────────────────")

layout_path <- file.path(data_dir, "APC_PLATES.csv")
if (!file.exists(layout_path))
  stop("Layout file not found: ", layout_path)

layout_raw <- read_csv(layout_path, col_types = cols(.default = "c"),
                        show_col_types = FALSE)
names(layout_raw) <- names(layout_raw) |>
  str_to_lower() |> str_replace_all("[^a-z0-9]+", "_") |>
  str_replace_all("_+", "_") |> str_remove("_$")

layout_clean <- layout_raw %>%
  mutate(
    plate_id      = toupper(trimws(plate)),
    well_id       = normalize_well_id(well),
    is_blank      = is.na(size) | toupper(trimws(size)) %in% c("NA", ""),
    area_mm2      = suppressWarnings(as.numeric(size)),
    treatment     = factor(toupper(trimws(treatment)), levels = treatment_levels),
    rep_group     = trimws(replicate_group),
    # Extract base_plate_id (e.g. APC-CC1) to enable cross-temperature matching.
    # The same individuals appear in both -11 and -26 plates, so a single
    # area measurement applies to reads at both temperatures.
    base_plate_id = str_remove(plate_id, "-[0-9]+$")
  ) %>%
  select(plate_id, base_plate_id, well_id, is_blank, area_mm2,
         treatment, rep_group, imagej_label)

# Build a deduplicated area lookup keyed on base_plate_id + well_id.
# Where sizes differ between the two temperature entries for the same
# individual (they shouldn't), take the mean.
area_lookup <- layout_clean %>%
  filter(!is_blank, is.finite(area_mm2)) %>%
  group_by(base_plate_id, well_id) %>%
  summarise(area_mm2 = mean(area_mm2, na.rm = TRUE), .groups = "drop")

message("  Layout rows: ", nrow(layout_clean),
        " | Plates: ",    n_distinct(layout_clean$plate_id),
        " | Blanks: ",    sum(layout_clean$is_blank, na.rm = TRUE))

# =============================================================================
# 4. Merge plates with layout
# =============================================================================
message("\n── 4. Merging plate data with layout ───────────────────────────────")

# Step 1: join treatment/blank metadata using the full plate_id
#   (e.g. APC-CC1-26 matches APC-CC1-26 rows in layout_clean)
meta_lookup <- layout_clean %>%
  select(plate_id, well_id, is_blank, treatment, rep_group, imagej_label) %>%
  distinct()

dat <- plate_raw %>%
  left_join(meta_lookup, by = c("plate_id", "well_id")) %>%
  # Step 2: join area on base_plate_id + well_id so that both the 11°C
  #   and 26°C reads for the same individual get the same area value
  left_join(area_lookup, by = c("base_plate_id", "well_id")) %>%
  mutate(
    is_blank  = replace_na(is_blank, FALSE),
    # sample_id uses base_plate_id + well so the same animal has the same
    # ID across both temperature phases
    sample_id = if_else(!is_blank,
                        paste(base_plate_id, well_id, sep = "_"),
                        NA_character_)
  )

# Plate consistency check: every plate should have the same well count at every T
well_counts <- dat %>%
  group_by(plate_id, t_index) %>%
  summarise(n_wells = n_distinct(well_id), .groups = "drop")

expected_n <- as.integer(names(which.max(table(well_counts$n_wells))))
bad_plates  <- well_counts %>% filter(n_wells != expected_n) %>%
  pull(plate_id) %>% unique()

if (length(bad_plates) > 0) {
  warning("Dropping ", length(bad_plates),
          " plate(s) with inconsistent well counts: ",
          paste(bad_plates, collapse = ", "))
  dat <- dat %>% filter(!plate_id %in% bad_plates)
} else {
  message("  Plate consistency check passed (", expected_n,
          " wells per plate per timepoint).")
}

# Verify area linkage: report how many non-blank wells have a valid area
area_check <- dat %>%
  filter(!is_blank, t_index == min(t_index)) %>%
  summarise(
    n_total     = n(),
    n_with_area = sum(is.finite(area_mm2)),
    n_missing   = sum(!is.finite(area_mm2))
  )
message("  Area linkage: ", area_check$n_with_area, "/", area_check$n_total,
        " non-blank wells have a valid area_mm2 value (",
        area_check$n_missing, " missing).")
message("  NOTE: area values are shared between 11°C and 26°C plates via ",
        "base_plate_id + well_id (same individuals, one size measurement).")

# ── Filter to target temperature ----------------------------------------------
dat <- dat %>% filter(temp_c == TARGET_TEMP)
n_indiv <- n_distinct(dat$sample_id, na.rm = TRUE)
message("\n  Retaining only ", TARGET_TEMP, "°C plates: ",
        n_distinct(dat$plate_id), " plates, ",
        n_indiv, " unique individuals.")

# =============================================================================
# 5. Blank correction via fold-change normalisation (Huffmyer et al.)
# =============================================================================
message("\n── 5. Blank-correcting via fold-change normalisation ────────────────")

# Step 1 – fold-change relative to each well's own T0 reading
t0_vals <- dat %>%
  filter(is.finite(value)) %>%
  group_by(plate_id, well_id) %>%
  slice_min(t_index, n = 1, with_ties = FALSE) %>%
  select(plate_id, well_id, value_t0 = value) %>%
  ungroup()

dat_fc <- dat %>%
  left_join(t0_vals, by = c("plate_id", "well_id")) %>%
  mutate(fold_change = if_else(
    is.finite(value_t0) & value_t0 > 0,
    value / value_t0,
    NA_real_))

# Step 2 – mean blank fold-change per plate × timepoint
blank_fc_ref <- dat_fc %>%
  filter(is_blank) %>%
  group_by(plate_id, t_index) %>%
  summarise(mean_blank_fc = mean(fold_change, na.rm = TRUE), .groups = "drop")

# Step 3 – subtract blank fold-change from sample fold-change
samples <- dat_fc %>%
  filter(!is_blank, !is.na(treatment)) %>%
  left_join(blank_fc_ref, by = c("plate_id", "t_index")) %>%
  mutate(corrected_fc = fold_change - mean_blank_fc)

message("  Sample rows after blank correction: ", nrow(samples),
        " | Individuals: ", n_distinct(samples$sample_id))

# =============================================================================
# 6. Raw fluorescence plots
# =============================================================================
message("\n── 6. Plotting raw fluorescence ─────────────────────────────────────")

p_raw_plates <- dat %>%
  filter(is.finite(value)) %>%
  mutate(grp = if_else(is_blank, "blank", as.character(treatment))) %>%
  ggplot(aes(x = time_hr, y = value,
             group = paste(plate_id, well_id), colour = grp)) +
  geom_line(alpha = 0.5) + geom_point(size = 1, alpha = 0.6) +
  facet_wrap(~ plate_id, ncol = 4) +
  scale_colour_manual(values = c(okabe_ito, blank = "black"),
                      na.value = "grey70", name = "Treatment") +
  labs(x = "Elapsed time (h)", y = "Raw fluorescence (RFU)",
       title = "Raw fluorescence by plate (all timepoints)") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7))

ggsave(file.path(fig_dir, "01_raw_fluor_by_plate.png"),
       p_raw_plates, width = 14, height = 8)

# Mean raw by treatment
raw_trt_sum <- dat %>%
  filter(!is_blank, is.finite(value), !is.na(treatment)) %>%
  group_by(treatment, time_hr) %>%
  summarise(mean_val = mean(value), se_val = sd(value)/sqrt(n()), .groups="drop")

p_raw_trt <- ggplot(raw_trt_sum,
    aes(x = time_hr, y = mean_val, colour = treatment,
        fill = treatment, group = treatment)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = okabe_ito, name = "Treatment") +
  scale_fill_manual(values = okabe_ito, name = "Treatment") +
  labs(x = "Elapsed time (h)", y = "Mean raw fluorescence (RFU ± SE)",
       title = "Mean raw fluorescence by treatment") +
  theme_classic(base_size = 13)

ggsave(file.path(fig_dir, "02_raw_fluor_mean_by_treatment.png"),
       p_raw_trt, width = 7, height = 5)

# =============================================================================
# 7. Blank-corrected fold-change plots
# =============================================================================
message("\n── 7. Plotting blank-corrected fold-change ───────────────────────────")

p_bc_ind <- samples %>%
  ggplot(aes(x = time_hr, y = corrected_fc,
             group = sample_id, colour = treatment)) +
  geom_line(alpha = 0.45) + geom_point(size = 1, alpha = 0.6) +
  facet_wrap(~ treatment) +
  scale_colour_manual(values = okabe_ito, name = "Treatment") +
  labs(x = "Elapsed time (h)", y = "Blank-corrected fold-change",
       title = "Individual blank-corrected metabolic traces by treatment") +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "03_bc_fc_individual_by_treatment.png"),
       p_bc_ind, width = 10, height = 5)

bc_sum <- samples %>%
  group_by(treatment, time_hr) %>%
  summarise(mean_val = mean(corrected_fc, na.rm = TRUE),
            se_val   = sd(corrected_fc, na.rm = TRUE) /
                       sqrt(sum(!is.na(corrected_fc))),
            .groups  = "drop")

p_bc_mean <- ggplot(bc_sum,
    aes(x = time_hr, y = mean_val, colour = treatment,
        fill = treatment, group = treatment)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = okabe_ito, name = "Treatment") +
  scale_fill_manual(values = okabe_ito, name = "Treatment") +
  labs(x = "Elapsed time (h)",
       y = "Mean blank-corrected fold-change (± SE)",
       title = "Mean blank-corrected metabolic activity by treatment") +
  theme_classic(base_size = 13)

ggsave(file.path(fig_dir, "04_bc_fc_mean_by_treatment.png"),
       p_bc_mean, width = 7, height = 5)

# =============================================================================
# 8. Size-normalised metabolism  (corrected_fc / area_mm²)
# =============================================================================
message("\n── 8. Computing size-normalised metabolism ───────────────────────────")

metabolism_df <- samples %>%
  mutate(
    metabolism = if_else(
      is.finite(area_mm2) & area_mm2 > 0 & is.finite(corrected_fc),
      corrected_fc / area_mm2,
      NA_real_))

n_valid <- metabolism_df %>%
  filter(t_index == 0, is.finite(metabolism)) %>% nrow()
message("  Animals with valid area data at T0: ", n_valid)

p_metab_ind <- metabolism_df %>%
  filter(is.finite(metabolism)) %>%
  ggplot(aes(x = time_hr, y = metabolism,
             group = sample_id, colour = treatment)) +
  geom_line(alpha = 0.45) + geom_point(size = 1, alpha = 0.6) +
  facet_wrap(~ treatment) +
  scale_colour_manual(values = okabe_ito, name = "Treatment") +
  labs(x = "Elapsed time (h)", y = "Metabolism (fold-change / mm²)",
       title = "Individual size-normalised metabolism by treatment") +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "05_metabolism_individual_by_treatment.png"),
       p_metab_ind, width = 10, height = 5)

metab_sum <- metabolism_df %>%
  group_by(treatment, time_hr) %>%
  summarise(mean_val = mean(metabolism, na.rm = TRUE),
            se_val   = sd(metabolism, na.rm = TRUE) /
                       sqrt(sum(!is.na(metabolism))),
            n        = sum(!is.na(metabolism)),
            .groups  = "drop")

p_metab_mean <- ggplot(metab_sum,
    aes(x = time_hr, y = mean_val, colour = treatment,
        fill = treatment, group = treatment)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = okabe_ito, name = "Treatment") +
  scale_fill_manual(values = okabe_ito, name = "Treatment") +
  labs(x = "Elapsed time (h)",
       y = "Metabolism (fold-change / mm² ± SE)",
       title = "Mean size-normalised metabolic activity by treatment") +
  theme_classic(base_size = 13)

ggsave(file.path(fig_dir, "06_metabolism_mean_by_treatment.png"),
       p_metab_mean, width = 7, height = 5)

# =============================================================================
# 9. AUC – trapezoid rule over real elapsed hours
# =============================================================================
message("\n── 9. Computing AUC (trapezoid rule, real elapsed hours) ────────────")

auc_df <- metabolism_df %>%
  filter(is.finite(time_hr), is.finite(metabolism)) %>%
  group_by(sample_id, treatment, rep_group, plate_id, well_id) %>%
  summarise(
    AUC          = trapezoid_auc(time_hr, metabolism),
    n_timepoints = n(),
    .groups      = "drop"
  ) %>%
  filter(is.finite(AUC))

message("  AUC computed for ", nrow(auc_df), " individuals")

auc_summary <- auc_df %>%
  group_by(treatment) %>%
  summarise(n      = n(),
            mean   = mean(AUC),
            sd     = sd(AUC),
            se     = sd / sqrt(n),
            median = median(AUC),
            .groups = "drop")

message("\n  AUC summary table:")
print(auc_summary)

write_csv(auc_df,      file.path(out_dir, "auc_all_individuals.csv"))
write_csv(auc_summary, file.path(out_dir, "auc_summary.csv"))
write_csv(metabolism_df, file.path(out_dir, "metabolism_full.csv"))

# =============================================================================
# 10. Statistical analysis – 2×2 factorial design
# =============================================================================
# Treatments map onto two binary priming factors:
#   heat_priming   : EE, EP = TRUE  |  CC, CP = FALSE
#   immune_priming : CP, EP = TRUE  |  CC, EE = FALSE
#
# This allows testing:
#   (1) Main effect of heat priming
#   (2) Main effect of immune priming
#   (3) Heat × immune interaction (synergy / antagonism)
#   (4) Tukey pairwise comparisons among all 4 groups
# =============================================================================
message("\n── 10. Statistical analysis (2×2 factorial) ─────────────────────────")

auc_df <- auc_df %>%
  mutate(
    heat_priming   = factor(if_else(treatment %in% c("EE", "EP"), "yes", "no"),
                            levels = c("no", "yes")),
    immune_priming = factor(if_else(treatment %in% c("CP", "EP"), "yes", "no"),
                            levels = c("no", "yes"))
  )

# ── AUC: 2×2 factorial linear model ─────────────────────────────────────────
model_auc <- lm(AUC ~ heat_priming * immune_priming, data = auc_df)
anova_res  <- anova(model_auc)

message("\n  2×2 factorial ANOVA on AUC:")
print(anova_res)

# Main effects (marginal means collapsed across the other factor)
emm_heat   <- emmeans(model_auc, ~ heat_priming)
emm_immune <- emmeans(model_auc, ~ immune_priming)
emm_inter  <- emmeans(model_auc, ~ heat_priming * immune_priming)

pairs_heat   <- as.data.frame(pairs(emm_heat,   adjust = "none"))
pairs_immune <- as.data.frame(pairs(emm_immune, adjust = "none"))
# All pairwise (Tukey) among the 4 treatment groups
pairs_trt    <- as.data.frame(pairs(emm_inter,  adjust = "tukey")) %>%
  # relabel contrasts to treatment names for readability
  mutate(contrast = str_replace_all(contrast,
    c("no no" = "CC", "yes no" = "EE", "no yes" = "CP", "yes yes" = "EP")))

message("\n  Main effect – heat priming:")
print(pairs_heat)

message("\n  Main effect – immune priming:")
print(pairs_immune)

message("\n  Tukey pairwise comparisons (all 4 groups):")
print(pairs_trt)

write_csv(as.data.frame(anova_res), file.path(out_dir, "anova_2x2_results.csv"))
write_csv(pairs_heat,   file.path(out_dir, "contrast_heat_priming.csv"))
write_csv(pairs_immune, file.path(out_dir, "contrast_immune_priming.csv"))
write_csv(pairs_trt,    file.path(out_dir, "pairwise_tukey_treatment.csv"))

# Interaction contrast: is the EP effect additive or synergistic?
# (EP - CP) vs (EE - CC)  — i.e. does adding heat priming to immune-primed
# animals produce the same gain as adding it to naive animals?
interaction_contrast <- contrast(emm_inter,
  interaction = c(heat_priming = "pairwise", immune_priming = "pairwise"))
message("\n  Interaction contrast (synergy / antagonism test):")
print(as.data.frame(interaction_contrast))
write_csv(as.data.frame(interaction_contrast),
          file.path(out_dir, "interaction_contrast.csv"))

# ── Time-series LMM: 2×2 factorial × time ────────────────────────────────────
message("\n  Fitting time-series linear mixed model (2×2 × time) …")

# ── Time-series LMM: 2×2 factorial × time ────────────────────────────────────
# t_index (0, 1, 2 … 8) is used as the time axis rather than the parsed
# wall-clock hours.  Wall-clock timestamps vary by a few minutes between
# plates (each plate is read sequentially), so using rounded hours produces
# many near-duplicate factor levels (e.g. "0.419" vs "0.454") that are all
# nominally the same timepoint but appear as separate levels in the model.
# This creates hundreds of empty treatment × time cells and makes the
# interaction terms unestimable.  t_index is already a clean shared integer
# (all plates have T0, T1, … T8) and avoids this entirely.
# real elapsed hours (time_hr) are kept for the AUC trapezoid calculation
# and axis labelling, but t_index drives the statistical model.
message("\n  Fitting time-series linear mixed model (2×2 × timepoint) …")

ts_df <- metabolism_df %>%
  filter(is.finite(metabolism)) %>%
  mutate(
    timepoint      = factor(t_index),          # clean 0–8 integer index
    individual     = factor(sample_id),
    heat_priming   = factor(if_else(treatment %in% c("EE", "EP"), "yes", "no"),
                            levels = c("no", "yes")),
    immune_priming = factor(if_else(treatment %in% c("CP", "EP"), "yes", "no"),
                            levels = c("no", "yes"))
  )

ts_model <- tryCatch(
  lmer(metabolism ~ timepoint * heat_priming * immune_priming + (1 | individual),
       data = ts_df),
  error = function(e) { message("  LMM failed: ", e$message); NULL })

if (!is.null(ts_model)) {
  ts_anova <- anova(ts_model, type = 3, ddf = "Satterthwaite")
  message("\n  Time-series 2×2 factorial ANOVA (LMM):")
  print(ts_anova)
  write_csv(as.data.frame(ts_anova),
            file.path(out_dir, "timeseries_lmm_anova.csv"))

  # Marginal means across timepoints for each treatment combination
  emm_ts <- emmeans(ts_model, ~ heat_priming * immune_priming | timepoint)
  ts_pairs <- as.data.frame(pairs(emm_ts, adjust = "tukey")) %>%
    mutate(contrast = str_replace_all(contrast,
      c("no no" = "CC", "yes no" = "EE", "no yes" = "CP", "yes yes" = "EP")))
  write_csv(ts_pairs, file.path(out_dir, "timeseries_pairwise_by_timepoint.csv"))
}

# =============================================================================
# 11. AUC box plots
# =============================================================================
message("\n── 11. AUC box plots ─────────────────────────────────────────────────")

df_plot <- auc_df %>%
  mutate(treatment = factor(treatment, levels = treatment_levels))

# ── Fig 7: All 4 groups with Tukey brackets ───────────────────────────────────
p_auc <- ggplot(df_plot, aes(x = treatment, y = AUC, fill = treatment)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.8) +
  scale_fill_manual(values = okabe_ito, guide = "none") +
  labs(x = "Treatment",
       y = "Metabolism (AUC; fold-change·h / mm²)",
       title = "Size-normalised metabolic AUC by treatment") +
  theme_classic(base_size = 13)

p_auc <- add_sig_brackets(p_auc, pairs_trt, treatment_levels, df_plot$AUC)
ggsave(file.path(fig_dir, "07_auc_boxplot_treatment.png"),
       p_auc, width = 6, height = 5)
print(p_auc)

# ── Fig 8: Mean ± SE bar chart ────────────────────────────────────────────────
p_auc_bar <- ggplot(auc_summary,
    aes(x = factor(treatment, levels = treatment_levels),
        y = mean, fill = treatment)) +
  geom_col(alpha = 0.75, colour = "white") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25) +
  scale_fill_manual(values = okabe_ito, guide = "none") +
  labs(x = "Treatment", y = "Mean AUC (fold-change·h / mm² ± SE)",
       title = "Mean metabolic AUC by treatment") +
  theme_classic(base_size = 13)

ggsave(file.path(fig_dir, "08_auc_bar_treatment.png"),
       p_auc_bar, width = 6, height = 5)

# ── Fig 9: 2×2 factorial interaction plot ─────────────────────────────────────
# Shows means ± SE with immune_priming on x-axis, heat_priming as line colour.
# This is the canonical way to visualise a 2×2 interaction.
factorial_sum <- auc_df %>%
  group_by(heat_priming, immune_priming) %>%
  summarise(mean = mean(AUC), se = sd(AUC)/sqrt(n()), n = n(),
            .groups = "drop") %>%
  mutate(
    heat_label   = if_else(heat_priming   == "yes",
                           "Heat primed (EE/EP)", "Not heat primed (CC/CP)"),
    immune_label = if_else(immune_priming == "yes",
                           "Immune primed\n(CP/EP)", "Not immune primed\n(CC/EE)")
  )

p_factorial <- ggplot(factorial_sum,
    aes(x = immune_label, y = mean,
        colour = heat_label, group = heat_label)) +
  geom_line(linewidth = 1) +
  geom_pointrange(aes(ymin = mean - se, ymax = mean + se),
                  size = 0.8, linewidth = 1) +
  scale_colour_manual(values = c("Heat primed (EE/EP)"       = "#E69F00",
                                 "Not heat primed (CC/CP)"   = "#009E73"),
                      name = NULL) +
  labs(x = "Immune priming",
       y = "Mean AUC (fold-change·h / mm² ± SE)",
       title = "2×2 factorial interaction plot",
       subtitle = "Parallel lines = additive; crossing/diverging = interaction") +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "09_auc_factorial_interaction.png"),
       p_factorial, width = 6, height = 5)
print(p_factorial)

# ── Fig 10: Main effects box plots ───────────────────────────────────────────
p_heat <- ggplot(auc_df,
    aes(x = heat_priming, y = AUC, fill = heat_priming)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.45, size = 1.5) +
  scale_fill_manual(values = c(no = "#009E73", yes = "#E69F00"),
                    labels = c(no = "No (CC/CP)", yes = "Yes (EE/EP)"),
                    name = "Heat priming") +
  labs(x = "Heat priming", y = "AUC (fold-change·h / mm²)",
       title = "Main effect: heat priming") +
  theme_classic(base_size = 13)

p_immune <- ggplot(auc_df,
    aes(x = immune_priming, y = AUC, fill = immune_priming)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.45, size = 1.5) +
  scale_fill_manual(values = c(no = "#009E73", yes = "#56B4E9"),
                    labels = c(no = "No (CC/EE)", yes = "Yes (CP/EP)"),
                    name = "Immune priming") +
  labs(x = "Immune priming", y = "AUC (fold-change·h / mm²)",
       title = "Main effect: immune priming") +
  theme_classic(base_size = 13)

p_main_effects <- plot_grid(p_heat, p_immune, ncol = 2)
ggsave(file.path(fig_dir, "10_auc_main_effects.png"),
       p_main_effects, width = 10, height = 5)
print(p_main_effects)

# ── Fig 11: AUC by replicate group within treatment ───────────────────────────
p_rep <- ggplot(auc_df,
    aes(x = rep_group, y = AUC, colour = treatment)) +
  geom_point(size = 2.5, alpha = 0.75) +
  facet_wrap(~ treatment, scales = "free_x") +
  scale_colour_manual(values = okabe_ito, guide = "none") +
  labs(x = "Replicate group", y = "AUC (fold-change·h / mm²)",
       title = "Individual AUC by replicate group within treatment") +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "11_auc_by_rep_group.png"),
       p_rep, width = 9, height = 5)

# =============================================================================
# Done
# =============================================================================
message("\n══════════════════════════════════════════════════════════════════════")
message("Analysis complete.")
message("Figures  → ", fig_dir)
message("CSV data → ", out_dir)
message("══════════════════════════════════════════════════════════════════════\n")
