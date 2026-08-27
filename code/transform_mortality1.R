# =============================================================================
# transform_mortality.R
# Transforms plate-based clam mortality CSV data into tidy output.
# All timepoints present in the file are automatically detected and each
# gets its own column in the output (named exactly as in the raw data,
# e.g. "4HR", "29.5HR").
#
# Usage:
#   source("transform_mortality.R")
#   OR from terminal: Rscript transform_mortality.R
#
# Configure the variables in the CONFIG section below before running.
# =============================================================================

# =============================================================================
# CONFIG — edit these variables for each new input file
# =============================================================================

INPUT_FILE <- "data/LR_manila-40C.csv"

# Output file path. Set to NULL to auto-generate (data/<input_name>_mortality.csv)
OUTPUT_FILE <- NULL  # e.g. "/path/to/output.csv"

# Treatment lookup: treatment code -> c(Immune_Priming, Heat_Priming)
TREATMENT_INFO <- list(
  "CC" = c("No",  "No"),
  "EE" = c("No",  "Yes"),
  "CP" = c("Yes", "No"),
  "EP" = c("Yes", "Yes")
)

# =============================================================================
# PARSE — no edits needed below this line
# =============================================================================

parse_plate_id <- function(plate_id) {
  # Parse e.g. "LR-CC1-40" -> list(site, treatment_code, replicate, assay_temp)
  m <- regmatches(plate_id, regexec("^([^-]+)-([A-Z]+)(\\d+)-(\\d+)$", plate_id))[[1]]
  if (length(m) == 0) stop(paste("Cannot parse plate ID:", plate_id))
  list(
    site           = m[2],
    treatment_code = m[3],
    replicate      = as.integer(m[4]),
    assay_temp     = paste0(m[5], "C")
  )
}

transform_mortality <- function(input_file, output_file = NULL) {

  # Read entire file as character matrix (no header, all strings)
  raw <- read.csv(input_file, header = FALSE, colClasses = "character",
                  stringsAsFactors = FALSE, fill = TRUE, check.names = FALSE)
  mat <- as.matrix(raw)
  mat[is.na(mat)] <- ""
  mat <- trimws(mat)

  nrows <- nrow(mat)
  ncols <- ncol(mat)

  # ------------------------------------------------------------------
  # 1. Find the first plate-header row (col 1 matches plate ID pattern)
  # ------------------------------------------------------------------
  plate_pattern <- "^[A-Za-z]+-[A-Z]+[0-9]+-[0-9]+$"
  first_plate_row <- NA
  for (i in seq_len(nrows)) {
    if (grepl(plate_pattern, mat[i, 1])) {
      first_plate_row <- i
      break
    }
  }
  if (is.na(first_plate_row)) stop("Could not find any plate ID row in the file.")

  # ------------------------------------------------------------------
  # 2. Discover all plates and their column offsets
  # ------------------------------------------------------------------
  header_row <- mat[first_plate_row, ]
  plates <- list()
  for (j in seq_len(ncols)) {
    if (grepl(plate_pattern, header_row[j])) {
      plates <- c(plates, list(list(plate_id = header_row[j], col = j)))
    }
  }
  cat(sprintf("Found %d plates: %s\n", length(plates),
              paste(sapply(plates, `[[`, "plate_id"), collapse = ", ")))

  # ------------------------------------------------------------------
  # 3. Discover ALL timepoint blocks in file order
  # ------------------------------------------------------------------
  hr_pattern <- "^[0-9]+\\.?[0-9]*HR$"
  timepoint_blocks <- list()
  for (i in seq_len(nrows)) {
    val <- mat[i, 1]
    if (grepl(hr_pattern, val, ignore.case = TRUE)) {
      plate_hdr_row <- i + 1
      if (plate_hdr_row <= nrows) {
        timepoint_blocks <- c(timepoint_blocks,
                              list(list(hr = val, plate_hdr_row = plate_hdr_row)))
      }
    }
  }

  # Deduplicate: keep only the first occurrence of each HR label
  seen_hrs      <- c()
  unique_blocks <- list()
  for (blk in timepoint_blocks) {
    hr_upper <- toupper(blk$hr)
    if (!hr_upper %in% seen_hrs) {
      seen_hrs      <- c(seen_hrs, hr_upper)
      unique_blocks <- c(unique_blocks, list(blk))
    }
  }

  tp_cols <- sapply(unique_blocks, `[[`, "hr")
  cat(sprintf("Found %d timepoints: %s\n", length(tp_cols), paste(tp_cols, collapse = ", ")))

  # ------------------------------------------------------------------
  # 4. Discover non-blank well column positions for each plate
  # ------------------------------------------------------------------
  get_well_cols <- function(plate_col) {
    well_info <- list()
    for (offset in 1:7) {
      abs_col <- plate_col + offset
      if (abs_col > ncols) break
      cell <- mat[first_plate_row, abs_col]
      if (cell == "" || tolower(cell) == "blank") next
      well_label <- sub("\\.0$", "", cell)
      well_info <- c(well_info, list(list(label = well_label, col = abs_col)))
    }
    well_info
  }

  # ------------------------------------------------------------------
  # 5. Build output records
  #
  # Replicate structure:
  #   Silo_Replicate  = plate number (1 or 2); one silo per plate
  #                     plate 1 holds silo 1 clams (rows A & B)
  #                     plate 2 holds silo 2 clams (rows C & D)
  #   Field_Replicate = each row letter is its own field box replicate,
  #                     labelled as silo+row (e.g. "1A", "1B", "2C", "2D")
  #   Assay_Replicate = full plate ID (the physical 24-well assay plate)
  # ------------------------------------------------------------------
  ROW_LABELS <- c("A", "B", "C", "D")
  meta_cols  <- c("PlateID_well", "Silo_Replicate", "Field_Replicate", "Assay_Replicate",
                  "Immune_Priming", "Heat_Priming", "Assay_Temp")

  records <- list()

  for (pl in plates) {
    plate_id  <- pl$plate_id
    plate_col <- pl$col

    parsed  <- parse_plate_id(plate_id)
    tx_info <- TREATMENT_INFO[[parsed$treatment_code]]
    if (is.null(tx_info)) tx_info <- c("Unknown", "Unknown")

    silo_rep  <- parsed$replicate   # plate number = silo number
    well_cols <- get_well_cols(plate_col)

    for (row_label in ROW_LABELS) {
      row_offset <- which(ROW_LABELS == row_label)  # A=1, B=2, C=3, D=4

      # Field replicates A/B/C/D map to plate rows as follows:
      # Plate 1: rows A&B = field rep A, rows C&D = field rep B
      # Plate 2: rows A&B = field rep C, rows C&D = field rep D
      field_rep_labels <- list("1" = c("A"="A","B"="A","C"="B","D"="B"),
                               "2" = c("A"="C","B"="C","C"="D","D"="D"))
      field_rep <- field_rep_labels[[as.character(silo_rep)]][[row_label]]

      for (wc in well_cols) {
        cup_id <- paste0(plate_id, "_", row_label, wc$label)

        rec <- list(
          "PlateID_well"    = cup_id,
          "Silo_Replicate"  = silo_rep,
          "Field_Replicate" = field_rep,
          "Assay_Replicate" = plate_id,
          "Immune_Priming"  = tx_info[1],
          "Heat_Priming"    = tx_info[2],
          "Assay_Temp"      = parsed$assay_temp
        )

        for (blk in unique_blocks) {
          data_row <- blk$plate_hdr_row + row_offset
          if (data_row > nrows) {
            rec[[blk$hr]] <- NA
            next
          }
          cell_val <- mat[data_row, wc$col]
          if (cell_val == "" || tolower(cell_val) == "blank") {
            rec[[blk$hr]] <- NA
          } else {
            num <- suppressWarnings(as.numeric(cell_val))
            rec[[blk$hr]] <- if (!is.na(num)) as.integer(num) else cell_val
          }
        }

        records <- c(records, list(rec))
      }
    }
  }

  # ------------------------------------------------------------------
  # 6. Assemble data frame
  # ------------------------------------------------------------------
  out_df <- do.call(rbind, lapply(records, as.data.frame,
                                  stringsAsFactors = FALSE, check.names = FALSE))
  out_df <- out_df[, c(meta_cols, tp_cols)]

  # Drop rows where ALL timepoint columns are NA
  all_na <- apply(out_df[, tp_cols, drop = FALSE], 1, function(r) all(is.na(r)))
  out_df <- out_df[!all_na, ]

  # ------------------------------------------------------------------
  # 7. Write output
  # ------------------------------------------------------------------
  if (is.null(output_file)) {
    base        <- tools::file_path_sans_ext(basename(input_file))
    output_file <- file.path(getwd(), "data", paste0(base, "_mortality.csv"))
  }

  write.csv(out_df, output_file, row.names = FALSE, quote = TRUE)
  cat(sprintf("\nSaved %d rows to: %s\n", nrow(out_df), output_file))
  invisible(out_df)
}

# =============================================================================
# RUN
# =============================================================================

result <- transform_mortality(INPUT_FILE, OUTPUT_FILE)
cat("\nPreview of output (first 10 rows):\n")
print(head(result, 10))
