# =============================================================================
# Power Analysis – APC Cockle Resazurin Experiment (2×2 Factorial Design)
# =============================================================================
# Factors:
#   heat_priming   (A): no = CC/CP,  yes = EE/EP
#   immune_priming (B): no = CC/EE,  yes = CP/EP
#
# Primary metric: size-normalised AUC (fold-change·h / mm²)
#
# Analyses:
#   1. 2×2 factorial ANOVA power (main effects + interaction)
#   2. Pairwise t-test power (Tukey-corrected, all 6 pairs)
#   3. Focused contrasts: main effect of A, main effect of B, A×B interaction
#   4. Sensitivity analysis: power across a range of effect sizes
#   5. Simulation-based validation (non-parametric bootstrap)
#   6. Power curves and minimum-n recommendation table
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pwr)
  library(emmeans)
  library(broom)
  library(cowplot)
})

set.seed(2026)

# ── Directories ---------------------------------------------------------------
out_dir <- file.path("output", "APC_cockle")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

okabe_ito <- c(CC = "#009E73", EE = "#E69F00", CP = "#56B4E9", EP = "#D55E00")
treatment_levels <- c("CC", "EE", "CP", "EP")

# =============================================================================
# 1. Load AUC data (produced by 01_resazurin_analysis.R)
# =============================================================================
auc_path <- file.path(out_dir, "auc_all_individuals.csv")

if (!file.exists(auc_path))
  stop("Run 01_resazurin_analysis.R first to generate: ", auc_path)

auc_df <- read_csv(auc_path, col_types = cols(.default = "c")) %>%
  mutate(
    AUC            = as.numeric(AUC),
    treatment      = factor(treatment, levels = treatment_levels),
    heat_priming   = factor(if_else(treatment %in% c("EE", "EP"), "yes", "no"),
                            levels = c("no", "yes")),
    immune_priming = factor(if_else(treatment %in% c("CP", "EP"), "yes", "no"),
                            levels = c("no", "yes"))
  ) %>%
  filter(is.finite(AUC))

message("AUC data loaded: ", nrow(auc_df), " individuals | ",
        n_distinct(auc_df$treatment), " treatments")

# =============================================================================
# 2. Observed descriptive statistics
# =============================================================================
message("\n── 2. Observed descriptive statistics ──────────────────────────────")

obs_trt <- auc_df %>%
  group_by(treatment, heat_priming, immune_priming) %>%
  summarise(n = n(), mean = mean(AUC), sd = sd(AUC),
            se = sd/sqrt(n), cv_pct = 100*sd/mean,
            .groups = "drop")

message("\n  Per-treatment summary:")
print(obs_trt)

# Pooled within-group SD (root mean square of group SDs)
pooled_sd <- sqrt(mean(obs_trt$sd^2, na.rm = TRUE))
message("\n  Pooled within-group SD: ", round(pooled_sd, 6))

# Bootstrap 95% CI on pooled SD
boot_sd <- replicate(2000, {
  df_b <- auc_df %>% group_by(treatment) %>%
    slice_sample(prop = 1, replace = TRUE) %>% ungroup()
  sds <- df_b %>% group_by(treatment) %>%
    summarise(s = sd(AUC, na.rm = TRUE), .groups = "drop") %>% pull(s)
  sqrt(mean(sds^2, na.rm = TRUE))
})
sd_ci <- quantile(boot_sd, c(0.025, 0.975))
message("  Bootstrap 95% CI for pooled SD: [",
        round(sd_ci[1], 6), ", ", round(sd_ci[2], 6), "]")

# Marginal means for main effects
means_heat   <- auc_df %>% group_by(heat_priming) %>%
  summarise(mean = mean(AUC), .groups = "drop")
means_immune <- auc_df %>% group_by(immune_priming) %>%
  summarise(mean = mean(AUC), .groups = "drop")

delta_heat   <- diff(means_heat$mean)          # yes - no
delta_immune <- diff(means_immune$mean)         # yes - no

# Interaction: (EP - CP) - (EE - CC)  [synergy above additive expectation]
m <- setNames(obs_trt$mean, obs_trt$treatment)
delta_interaction <- (m["EP"] - m["CP"]) - (m["EE"] - m["CC"])

message("\n  Marginal mean difference (heat priming):   ", round(delta_heat, 6))
message("  Marginal mean difference (immune priming): ", round(delta_immune, 6))
message("  Interaction contrast (synergy term):       ", round(delta_interaction, 6))

# =============================================================================
# 3. Effect sizes
# =============================================================================
message("\n── 3. Effect sizes ─────────────────────────────────────────────────")

# Cohen's d for main effects (marginal means / pooled SD)
d_heat        <- abs(delta_heat)   / pooled_sd
d_immune      <- abs(delta_immune) / pooled_sd
d_interaction <- abs(delta_interaction) / pooled_sd

# Cohen's f for the 2×2 ANOVA (overall)
# f = sigma_means / sigma_within, where sigma_means = SD of the 4 cell means
f_anova <- sd(obs_trt$mean) / pooled_sd

# Pairwise Cohen's d (all 6 pairs)
pairs_grid <- combn(treatment_levels, 2, simplify = FALSE)
cohens_d_table <- map_dfr(pairs_grid, function(pair) {
  g1 <- obs_trt %>% filter(treatment == pair[1])
  g2 <- obs_trt %>% filter(treatment == pair[2])
  d  <- abs(g1$mean - g2$mean) / pooled_sd
  tibble(group1 = pair[1], group2 = pair[2],
         delta_AUC = abs(g1$mean - g2$mean), cohens_d = d,
         magnitude = case_when(d < 0.2 ~ "negligible", d < 0.5 ~ "small",
                               d < 0.8 ~ "medium", TRUE ~ "large"))
})

message("\n  Summary of effect sizes:")
message("    Cohen's d – heat priming main effect:   ", round(d_heat, 3))
message("    Cohen's d – immune priming main effect: ", round(d_immune, 3))
message("    Cohen's d – interaction contrast:       ", round(d_interaction, 3))
message("    Cohen's f – overall 2×2 ANOVA:          ", round(f_anova, 3))
message("\n  Pairwise Cohen's d:")
print(cohens_d_table %>% mutate(across(where(is.numeric), \(x) round(x, 4))))

write_csv(cohens_d_table, file.path(out_dir, "power_cohens_d.csv"))

# =============================================================================
# 4. Parametric power analysis
# =============================================================================
message("\n── 4. Parametric power analysis ────────────────────────────────────")

alpha      <- 0.05
target_pw  <- 0.80
n_seq      <- 3:80   # n per treatment cell (= n per group in 2×2)
n_groups   <- 4

# ── 4a. Overall 2×2 ANOVA (any effect at all) ────────────────────────────────
anova_curve <- tibble(
  n_per_group = n_seq,
  power = map_dbl(n_seq, function(n) {
    pwr.anova.test(k = n_groups, n = n, f = f_anova,
                   sig.level = alpha)$power
  })
)

min_n_anova <- anova_curve %>% filter(power >= target_pw) %>%
  slice(1) %>% pull(n_per_group)
if (length(min_n_anova) == 0) min_n_anova <- NA_integer_

message("  Overall 2×2 ANOVA (f = ", round(f_anova, 3),
        "): min n = ", min_n_anova)

# ── 4b. Focused contrasts: main effects and interaction ───────────────────────
# For a 2×2 design each main effect compares two groups of size 2n (marginal)
# so the effective n per contrast arm = 2 × n_per_cell
focused_contrasts <- tibble(
  contrast     = c("Heat priming main effect",
                   "Immune priming main effect",
                   "Interaction (synergy)"),
  cohens_d     = c(d_heat, d_immune, d_interaction),
  # Bonferroni correction for 3 planned contrasts
  alpha_adj    = alpha / 3
) %>%
  mutate(
    power_curve = map2(cohens_d, alpha_adj, function(d, a) {
      tibble(
        n_per_cell = n_seq,
        # effective n per arm = 2 * n_per_cell (marginal means pool two cells)
        power = map_dbl(n_seq, function(n) {
          pwr.t.test(n = 2 * n, d = d, sig.level = a,
                     type = "two.sample", alternative = "two.sided")$power
        })
      )
    }),
    min_n = map2_int(power_curve, cohens_d, function(pc, d) {
      hit <- pc %>% filter(power >= target_pw) %>% slice(1) %>%
        pull(n_per_cell)
      if (length(hit) == 0) NA_integer_ else as.integer(hit)
    })
  )

message("\n  Focused contrasts (Bonferroni α = ", round(alpha/3, 4), "):")
focused_contrasts %>%
  select(contrast, cohens_d, min_n) %>%
  mutate(cohens_d = round(cohens_d, 3)) %>%
  print()

# ── 4c. All pairwise t-tests (Tukey-corrected) ───────────────────────────────
n_pairs <- choose(n_groups, 2)   # = 6

pairwise_power <- map_dfr(pairs_grid, function(pair) {
  d <- cohens_d_table %>%
    filter(group1 == pair[1], group2 == pair[2]) %>% pull(cohens_d)

  tibble(
    comparison  = paste(pair, collapse = " vs "),
    n_per_group = n_seq,
    cohens_d    = d,
    power = map_dbl(n_seq, function(n) {
      pwr.t.test(n = n, d = d, sig.level = alpha / n_pairs,
                 type = "two.sample", alternative = "two.sided")$power
    })
  )
}) %>%
  group_by(comparison, cohens_d) %>%
  mutate(min_n = {
    hit <- n_per_group[power >= target_pw]
    if (length(hit) == 0) NA_integer_ else as.integer(hit[1])
  }) %>%
  ungroup()

pairwise_min <- pairwise_power %>%
  select(comparison, cohens_d, min_n) %>% distinct() %>%
  arrange(cohens_d)

message("\n  Pairwise min n (Tukey-corrected α = ",
        round(alpha / n_pairs, 4), "):")
print(pairwise_min %>% mutate(cohens_d = round(cohens_d, 3)))

# =============================================================================
# 5. Sensitivity analysis (power over a range of effect sizes)
# =============================================================================
message("\n── 5. Sensitivity analysis ─────────────────────────────────────────")

d_range <- seq(0.1, max(c(cohens_d_table$cohens_d, 1.5), na.rm = TRUE),
               by = 0.05)

sensitivity_df <- crossing(cohens_d = d_range, n_per_group = 3:60) %>%
  mutate(
    power = map2_dbl(cohens_d, n_per_group, function(d, n) {
      pwr.t.test(n = n, d = d, sig.level = alpha / n_pairs,
                 type = "two.sample", alternative = "two.sided")$power
    })
  )

sensitivity_min <- sensitivity_df %>%
  filter(power >= target_pw) %>%
  group_by(cohens_d) %>%
  slice_min(n_per_group, n = 1) %>% ungroup()

# =============================================================================
# 6. Simulation-based power (2×2 ANOVA, non-parametric validation)
# =============================================================================
message("\n── 6. Simulation-based power (2×2 ANOVA) ───────────────────────────")

sim_means <- setNames(obs_trt$mean, obs_trt$treatment)

simulate_2x2_power <- function(n_per_cell, means, sd_within,
                                n_sim = 1000, alpha = 0.05) {
  sig_count <- 0
  for (i in seq_len(n_sim)) {
    dat_sim <- map_dfr(names(means), function(trt) {
      tibble(
        AUC       = rnorm(n_per_cell, mean = means[trt], sd = sd_within),
        treatment = trt
      )
    }) %>%
      mutate(
        heat_priming   = if_else(treatment %in% c("EE", "EP"), "yes", "no"),
        immune_priming = if_else(treatment %in% c("CP", "EP"), "yes", "no")
      )
    p_val <- anova(lm(AUC ~ heat_priming * immune_priming,
                      data = dat_sim))$`Pr(>F)`
    # Significant if any term (main effects or interaction) is significant
    if (any(p_val[1:3] < alpha, na.rm = TRUE)) sig_count <- sig_count + 1
  }
  sig_count / n_sim
}

sim_n_seq <- seq(4, 50, by = 2)
message("  Running ", length(sim_n_seq),
        " sample sizes × 1 000 simulations …")

sim_curve <- tibble(
  n_per_cell  = sim_n_seq,
  power_sim   = map_dbl(sim_n_seq, function(n) {
    simulate_2x2_power(n, sim_means, pooled_sd, n_sim = 1000)
  }),
  power_param = map_dbl(sim_n_seq, function(n) {
    pwr.anova.test(k = n_groups, n = n, f = f_anova,
                   sig.level = alpha)$power
  })
)

min_n_sim <- sim_curve %>% filter(power_sim >= target_pw) %>%
  slice(1) %>% pull(n_per_cell)
if (length(min_n_sim) == 0) min_n_sim <- NA_integer_

# =============================================================================
# 7. Minimum-n recommendation table
# =============================================================================
message("\n── 7. Minimum-n recommendation table ───────────────────────────────")

focused_min <- focused_contrasts %>%
  transmute(
    test          = contrast,
    method        = paste0("Parametric t-test (Bonferroni α = ",
                           round(alpha_adj, 4), ")"),
    effect_size   = paste0("d = ", round(cohens_d, 3)),
    # min_n is per cell; for main effects each arm pools 2 cells so actual
    # animals needed per cell = min_n
    min_n_per_cell = min_n,
    note          = "n per treatment group (4 groups total)"
  )

final_table <- bind_rows(
  tibble(test = "Overall 2×2 ANOVA",
         method = "Parametric (pwr.anova.test)",
         effect_size = paste0("f = ", round(f_anova, 3)),
         min_n_per_cell = min_n_anova,
         note = "n per treatment group; detects any effect"),
  tibble(test = "Overall 2×2 ANOVA",
         method = "Simulation (1 000 iterations)",
         effect_size = paste0("f = ", round(f_anova, 3)),
         min_n_per_cell = min_n_sim,
         note = "n per treatment group; detects any effect"),
  focused_min,
  pairwise_min %>%
    transmute(test = paste("Pairwise:", comparison),
              method = paste0("Parametric t-test (Tukey α = ",
                              round(alpha / n_pairs, 4), ")"),
              effect_size = paste0("d = ", round(cohens_d, 3)),
              min_n_per_cell = min_n,
              note = "n per treatment group")
)

message("\n  Final minimum-n table (", target_pw * 100, "% power, α = ", alpha, "):")
print(final_table)
write_csv(final_table, file.path(out_dir, "power_minimum_n_table.csv"))

# =============================================================================
# 8. Figures
# =============================================================================
message("\n── 8. Generating power figures ─────────────────────────────────────")

# ── Fig A: Overall ANOVA power curve (parametric + simulation) ───────────────
p_anova <- ggplot() +
  geom_line(data = anova_curve,
            aes(x = n_per_group, y = power, linetype = "Parametric"),
            colour = "#0072B2", linewidth = 1) +
  geom_point(data = sim_curve,
             aes(x = n_per_cell, y = power_sim, shape = "Simulation"),
             colour = "#D55E00", size = 2.2, alpha = 0.85) +
  geom_hline(yintercept = target_pw, linetype = "dashed", colour = "grey40") +
  { if (!is.na(min_n_anova))
      geom_vline(xintercept = min_n_anova, linetype = "dotted",
                 colour = "#0072B2") } +
  { if (!is.na(min_n_anova))
      annotate("text", x = min_n_anova + 1, y = 0.08,
               label = paste0("n = ", min_n_anova), hjust = 0,
               colour = "#0072B2", size = 3.5) } +
  scale_linetype_manual(name = NULL, values = c(Parametric = "solid")) +
  scale_shape_manual(name = NULL, values = c(Simulation = 16)) +
  labs(x = "Sample size per treatment group (n per cell)",
       y = "Statistical power",
       title = "Power to detect any treatment effect (2×2 factorial ANOVA)",
       subtitle = paste0("f = ", round(f_anova, 3),
                         "; pooled SD = ", round(pooled_sd, 5),
                         "; α = ", alpha)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(0, 80, 10)) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "PA_power_anova_curve.png"),
       p_anova, width = 7, height = 5)
print(p_anova)

# ── Fig B: Focused contrast power curves (main effects + interaction) ─────────
focused_curves_long <- focused_contrasts %>%
  select(contrast, power_curve) %>%
  unnest(power_curve)

p_focused <- ggplot(focused_curves_long,
    aes(x = n_per_cell, y = power, colour = contrast)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = target_pw, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(
    values = c("Heat priming main effect"   = "#E69F00",
               "Immune priming main effect" = "#56B4E9",
               "Interaction (synergy)"      = "#D55E00"),
    name = "Contrast") +
  labs(x = "Sample size per treatment group (n per cell)",
       y = "Statistical power",
       title = "Power for main effects and interaction contrast",
       subtitle = paste0("Bonferroni α = ", round(alpha / 3, 4),
                         " (3 planned contrasts); effective n = 2 × n per cell")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(0, 80, 10)) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10))

ggsave(file.path(fig_dir, "PB_power_focused_contrasts.png"),
       p_focused, width = 8, height = 5)
print(p_focused)

# ── Fig C: Pairwise power curves ─────────────────────────────────────────────
pair_colours <- c(
  "CC vs EE" = "#E69F00", "CC vs CP" = "#56B4E9",
  "CC vs EP" = "#D55E00", "EE vs CP" = "#009E73",
  "EE vs EP" = "#CC79A7", "CP vs EP" = "#0072B2"
)

p_pairwise <- ggplot(pairwise_power,
    aes(x = n_per_group, y = power, colour = comparison)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = target_pw, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = pair_colours, name = "Comparison") +
  labs(x = "Sample size per treatment group",
       y = "Statistical power",
       title = "Power for pairwise treatment comparisons",
       subtitle = paste0("Two-sample t-test; Tukey-corrected α = ",
                         round(alpha / n_pairs, 4))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(0, 80, 10)) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "PC_power_pairwise_curves.png"),
       p_pairwise, width = 8, height = 5)
print(p_pairwise)

# ── Fig D: Sensitivity heatmap ────────────────────────────────────────────────
p_heat_map <- ggplot(sensitivity_df %>% filter(n_per_group <= 60),
    aes(x = n_per_group, y = cohens_d, fill = power)) +
  geom_tile() +
  geom_contour(aes(z = power), breaks = c(0.5, 0.8, 0.9),
               colour = "white", linewidth = 0.4) +
  geom_label(data = sensitivity_min %>% filter(n_per_group <= 60),
             aes(label = n_per_group), fill = "white", size = 2.3,
             label.padding = unit(0.1, "lines")) +
  # Mark observed pairwise d values
  geom_hline(data = cohens_d_table,
             aes(yintercept = cohens_d,
                 colour = paste(group1, "vs", group2)),
             linetype = "dashed", linewidth = 0.5, show.legend = TRUE) +
  scale_fill_viridis_c(name = "Power", limits = c(0, 1),
                        breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = unname(pair_colours), name = "Observed d") +
  scale_x_continuous(breaks = seq(5, 60, 5)) +
  labs(x = "Sample size per treatment group",
       y = "Cohen's d (pairwise effect size)",
       title = "Power sensitivity: n × effect size",
       subtitle = "Labels = min n for 80% power; dashed lines = observed pairwise d") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "PD_power_sensitivity_heatmap.png"),
       p_heat_map, width = 9, height = 6)
print(p_heat_map)

# ── Fig E: Observed group means ± SD ─────────────────────────────────────────
p_obs <- ggplot(obs_trt,
    aes(x = treatment, y = mean, colour = treatment)) +
  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd),
                  linewidth = 0.9, size = 0.9) +
  geom_text(aes(y = mean - sd - 0.0001 * max(mean),
                label = paste0("n=", n)), vjust = 1.5, size = 3.2) +
  scale_colour_manual(values = okabe_ito, guide = "none") +
  labs(x = "Treatment", y = "Mean AUC ± SD (fold-change·h / mm²)",
       title = "Observed group means and SDs used for power estimation") +
  theme_classic(base_size = 13)

ggsave(file.path(fig_dir, "PE_observed_means_sd.png"),
       p_obs, width = 6, height = 5)

# =============================================================================
# 9. Recommendations printout
# =============================================================================
message("\n══════════════════════════════════════════════════════════════════════")
message("POWER ANALYSIS SUMMARY  (2×2 Factorial: Heat × Immune Priming)")
message("══════════════════════════════════════════════════════════════════════")
message(sprintf("  Target power:                  %.0f%%",   target_pw * 100))
message(sprintf("  Significance level (α):        %.3f",    alpha))
message(sprintf("  Number of treatment groups:    %d",       n_groups))
message(sprintf("  Pooled within-group SD:        %.5f",    pooled_sd))
message(sprintf("  Bootstrap SD 95%% CI:          [%.5f, %.5f]",
                sd_ci[1], sd_ci[2]))
message("")
message("  Effect sizes:")
message(sprintf("    Overall ANOVA (Cohen's f):           %.3f  (%s)",
                f_anova,
                case_when(f_anova < 0.1 ~ "negligible", f_anova < 0.25 ~ "small",
                          f_anova < 0.4 ~ "medium", TRUE ~ "large")))
message(sprintf("    Heat priming main effect (d):        %.3f", d_heat))
message(sprintf("    Immune priming main effect (d):      %.3f", d_immune))
message(sprintf("    Interaction / synergy contrast (d):  %.3f", d_interaction))
message("")
message("  Minimum n per treatment group for 80% power (α = 0.05):")
message(sprintf("    Overall 2×2 ANOVA (parametric):  n = %s",
                ifelse(is.na(min_n_anova), ">80", min_n_anova)))
message(sprintf("    Overall 2×2 ANOVA (simulation):  n = %s",
                ifelse(is.na(min_n_sim), ">50 simulated", min_n_sim)))
for (i in seq_len(nrow(focused_contrasts))) {
  message(sprintf("    %-38s n = %s",
                  paste0(focused_contrasts$contrast[i], ":"),
                  ifelse(is.na(focused_contrasts$min_n[i]), ">80",
                         focused_contrasts$min_n[i])))
}
message("")
message("  IMPORTANT CAVEATS:")
message("  • Effect sizes are estimated from a small pilot dataset.")
message("    Treat minimum-n values as lower bounds — plan for ~20% extra")
message("    animals to account for mortality and exclusions.")
message("  • Plate = treatment in this design; plate and treatment effects")
message("    are confounded. Future designs should distribute animals")
message("    across multiple plates per treatment.")
message("  • The interaction contrast (synergy term) has the largest")
message(sprintf("    required n because its effect size is small (d = %.3f).", d_interaction))
message("    If detecting synergy is a key goal, power this contrast")
message("    explicitly rather than relying on the overall ANOVA.")
message("══════════════════════════════════════════════════════════════════════\n")

message("Figures  → ", fig_dir)
message("CSV data → ", out_dir)
