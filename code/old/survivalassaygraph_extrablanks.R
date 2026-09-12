library(tidyverse)
library(scales)

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------

dat <- read.csv(
  "data/AP_manila-35C.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Find time-point rows
# ------------------------------------------------------------

time_rows <- which(
  grepl(
    "^[0-9]+\\.?[0-9]*HR$",
    trimws(dat[[1]]),
    ignore.case = TRUE
  )
)

all_data <- list()

# ------------------------------------------------------------
# Parse each timepoint and plate
# ------------------------------------------------------------

for (i in seq_along(time_rows)) {
  
  time_row <- time_rows[i]
  plate_row <- time_row + 1
  well_rows <- (time_row + 2):(time_row + 5)
  
  # Time in hours
  time_hr <- as.numeric(
    str_remove(
      trimws(dat[time_row, 1]),
      regex("HR", ignore_case = TRUE)
    )
  )
  
  # ----------------------------------------------------------
  # Find plate labels on this row
  # ----------------------------------------------------------
  
  plate_cols <- which(
    grepl(
      "^AP-(EE|CP|EP)[12]-35$|^AP-CC1-35$",
      trimws(as.character(dat[plate_row, ])),
      ignore.case = TRUE
    )
  )
  
  for (p in plate_cols) {
    
    plate_name <- trimws(
      as.character(dat[plate_row, p])
    )
    
    # Treatment
    treatment <- str_extract(
      plate_name,
      "(?<=AP-)[A-Z]{2}"
    )
    
    # Replicate
    replicate <- as.numeric(
      str_extract(
        plate_name,
        "[12](?=-35)"
      )
    )
    
    # --------------------------------------------------------
    # Six columns belonging to this plate
    # --------------------------------------------------------
    
    six_cols <- (p + 1):(p + 6)
    
    plate_values <- dat[
      well_rows,
      six_cols,
      drop = FALSE
    ]
    
    # --------------------------------------------------------
    # Columns 1 and 6 are blank.
    # Actual clam wells are columns 2-5.
    # --------------------------------------------------------
    
    clam_values <- plate_values[, 2:5, drop = FALSE]
    
    # --------------------------------------------------------
    # Convert ONLY exact 0 and 1 values to numeric.
    # Anything else = NA.
    # --------------------------------------------------------
    
    clam_values <- as.matrix(clam_values)
    
    clam_values <- trimws(
      as.character(clam_values)
    )
    
    clam_values <- ifelse(
      clam_values == "0",
      0,
      ifelse(
        clam_values == "1",
        1,
        NA_real_
      )
    )
    
    # --------------------------------------------------------
    # Create well-level data
    # --------------------------------------------------------
    
    wells <- expand.grid(
      row = LETTERS[1:4],
      column = 2:5
    )
    
    wells$status <- as.vector(clam_values)
    
    wells$time_hr <- time_hr
    wells$treatment <- treatment
    wells$replicate <- replicate
    wells$plate <- plate_name
    
    all_data[[length(all_data) + 1]] <- wells
  }
}

# Combine all wells
mortality_data <- bind_rows(all_data)


# ------------------------------------------------------------
# IMPORTANT CHECK:
# Number of actual wells and dead wells for every plate
# at every timepoint
# ------------------------------------------------------------

plate_summary <- mortality_data %>%
  group_by(
    time_hr,
    treatment,
    replicate,
    plate
  ) %>%
  summarise(
    
    # Every 1 = one dead clam
    dead = sum(status == 1, na.rm = TRUE),
    
    # Every 0 OR 1 = one actual clam
    total = sum(
      !is.na(status)
    ),
    
    # Number alive
    alive = sum(status == 0, na.rm = TRUE),
    
    # Mortality
    mortality = dead / total,
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# Combine plates within each treatment
#
# CC has one plate.
# EE, CP, and EP have two plates.
#
# Therefore CC gets its actual denominator,
# while the other treatments get the sum of both plates.
# ------------------------------------------------------------

mortality_summary <- plate_summary %>%
  group_by(
    time_hr,
    treatment
  ) %>%
  summarise(
    
    dead = sum(dead),
    total = sum(total),
    alive = sum(alive),
    
    mortality = dead / total,
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# Treatment order
# ------------------------------------------------------------

mortality_summary <- mortality_summary %>%
  mutate(
    treatment = factor(
      treatment,
      levels = c("CC", "EE", "CP", "EP")
    )
  )


# ------------------------------------------------------------
# Diagnostic table
#
# This lets you see exactly what denominator is being used.
# ------------------------------------------------------------

print(
  plate_summary %>%
    arrange(
      time_hr,
      treatment,
      replicate
    )
)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p <- ggplot(
  mortality_summary,
  aes(
    x = time_hr,
    y = mortality,
    color = treatment,
    group = treatment
  )
) +
  
  geom_line(
    linewidth = 1.2
  ) +
  
  geom_point(
    size = 3
  ) +
  
  scale_color_manual(
    values = c(
      CC = "blue",
      EE = "red",
      CP = "green",
      EP = "orange"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    labels = percent_format(accuracy = 1)
  ) +
  
  scale_x_continuous(
    breaks = sort(
      unique(
        mortality_summary$time_hr
      )
    )
  ) +
  
  labs(
    x = "Time at 35 °C (hours)",
    y = "Cumulative mortality",
    color = "Treatment"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    legend.position = "right"
  )

print(p)


# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------

ggsave(
  "figures/AP_survival/clam_mortality_35C.png",
  plot = p,
  width = 12,
  height = 6,
  dpi = 300
)