library(tidyverse)
library(scales)

dat <- read.csv(
  "data/LR_manila-40C.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Find time-point rows
# ------------------------------------------------------------

time_rows <- which(
  grepl("^[0-9]+\\.?[0-9]*HR$", dat[[1]], ignore.case = TRUE)
)

all_data <- list()

for (i in seq_along(time_rows)) {
  
  time_row <- time_rows[i]
  plate_row <- time_row + 1
  well_rows <- (time_row + 2):(time_row + 5)
  
  time_hr <- as.numeric(
    str_remove(dat[time_row, 1], "HR")
  )
  
  # Find all plate labels
  plate_cols <- which(
    grepl(
      "^LR-[A-Z]{2}[12]-40$",
      dat[plate_row, ],
      ignore.case = TRUE
    )
  )
  
  for (p in plate_cols) {
    
    plate_name <- dat[plate_row, p]
    
    # Treatment: CC, EE, CP, or EP
    treatment <- str_extract(
      plate_name,
      "(?<=LR-)[A-Z]{2}"
    )
    
    # Replicate: 1 or 2
    replicate <- as.numeric(
      str_extract(
        plate_name,
        "[12](?=-40)"
      )
    )
    
    # --------------------------------------------------------
    # Get the SIX plate columns
    # --------------------------------------------------------
    
    six_cols <- (p + 1):(p + 6)
    
    plate_values <- dat[
      well_rows,
      six_cols
    ]
    
    # Columns 1 and 6 are blanks.
    # Keep ONLY columns 2, 3, 4, and 5.
    clam_values <- plate_values[, 2:5]
    
    # Convert to numeric
    clam_values <- apply(
      clam_values,
      2,
      as.numeric
    )
    
    # Create well identifiers
    wells <- expand.grid(
      row = LETTERS[1:4],
      column = 2:5
    )
    
    wells$dead <- as.vector(clam_values)
    
    wells$time_hr <- time_hr
    wells$treatment <- treatment
    wells$replicate <- replicate
    wells$plate <- plate_name
    
    all_data[[length(all_data) + 1]] <- wells
  }
}

mortality_data <- bind_rows(all_data)


# ------------------------------------------------------------
# Calculate mortality
# ------------------------------------------------------------

mortality_summary <- mortality_data %>%
  group_by(time_hr, treatment) %>%
  summarise(
    dead = sum(dead, na.rm = TRUE),
    total = 32,
    mortality = dead / total,
    .groups = "drop"
  )


# ------------------------------------------------------------
# Treatment order
# ------------------------------------------------------------

mortality_summary$treatment <- factor(
  mortality_summary$treatment,
  levels = c("CC", "EE", "CP", "EP")
)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

ggplot(
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
    breaks = sort(unique(mortality_summary$time_hr))
  ) +
  
  labs(
    x = "Time at 40 °C (hours)",
    y = "Cumulative mortality",
    color = "Treatment"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    legend.position = "right"
  )


ggplot(
  mortality_summary,
  aes(
    x = time_hr,
    y = mortality,
    color = treatment,
    group = treatment
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
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
    breaks = sort(unique(mortality_summary$time_hr))
  ) +
  labs(
    x = "Time at 40 °C (hours)",
    y = "Cumulative mortality",
    color = "Treatment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right"
  )

# Save figure as a high-resolution PNG
ggsave(
  "figures/LR_survival/clam_mortality_40C.png",
  width = 12,
  height = 6,
  dpi = 300
)