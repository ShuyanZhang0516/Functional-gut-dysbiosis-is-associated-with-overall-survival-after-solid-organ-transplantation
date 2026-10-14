
# =============================================================================
# KO-BASED FUNCTIONAL DYSBIOSIS ANALYSIS
# =============================================================================
#
# Purpose:
#   1. Load KO abundance data and sample metadata
#   2. Filter KOs by prevalence
#   3. Perform CLR transformation
#   4. Calculate functional dysbiosis scores
#   5. Correlate SOTR KO profiles with dysbiosis scores
#   6. Perform PCA and generate a visualization
#
# Required input files:
#   - KO abundance CSV: Sample column + KO columns (e.g. K00001)
#   - Metadata CSV: Sample and Group columns
#
# Optional input:
#   - KO annotation CSV: KO column + functional annotation columns
#
# Configure file paths using environment variables before running.
#
# =============================================================================


# =============================================================================
# STEP 0. CONFIGURATION AND PACKAGES
# =============================================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tibble",
  "ggplot2",
  "ggrepel",
  "ggnewscale"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install missing packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(ggnewscale)
})

# Configure paths with environment variables.
KO_ABUNDANCE_FILE <- Sys.getenv("KO_ABUNDANCE_FILE", "")
METADATA_FILE <- Sys.getenv("METADATA_FILE", "")
KO_ANNOTATION_FILE <- Sys.getenv("KO_ANNOTATION_FILE", "")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", "results")

# Group labels.
SOTR_LABEL <- Sys.getenv("SOTR_LABEL", "SOTR")
HC_LABEL <- Sys.getenv("HC_LABEL", "HC")

# Analysis thresholds.
PREVALENCE_THRESHOLD_INITIAL <- 0.05
PREVALENCE_THRESHOLD_FINAL <- 0.20
FDR_THRESHOLD <- 0.05

# Check required input paths.
if (!nzchar(KO_ABUNDANCE_FILE)) {
  stop("Set KO_ABUNDANCE_FILE to your KO abundance CSV path.")
}

if (!nzchar(METADATA_FILE)) {
  stop("Set METADATA_FILE to your metadata CSV path.")
}

if (!file.exists(KO_ABUNDANCE_FILE)) {
  stop("KO abundance file not found: ", KO_ABUNDANCE_FILE)
}

if (!file.exists(METADATA_FILE)) {
  stop("Metadata file not found: ", METADATA_FILE)
}

# Create output folders.
dir.create(
  file.path(RESULTS_DIR, "intermediate"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(RESULTS_DIR, "tables"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(RESULTS_DIR, "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

# Helper function for writing output files.
write_output <- function(data, filename, folder = "tables") {
  output_path <- file.path(RESULTS_DIR, folder, filename)
  write_csv(data, output_path, na = "")
  message("Saved: ", output_path)
}


# =============================================================================
# STEP 1. LOAD AND VALIDATE INPUT DATA
# =============================================================================

message("\nSTEP 1: Loading input data")

ko_data <- read_csv(
  KO_ABUNDANCE_FILE,
  show_col_types = FALSE
)

metadata <- read_csv(
  METADATA_FILE,
  show_col_types = FALSE
)

if (!"Sample" %in% names(ko_data)) {
  stop("KO abundance file must contain a Sample column.")
}

if (!all(c("Sample", "Group") %in% names(metadata))) {
  stop("Metadata file must contain Sample and Group columns.")
}

if (anyDuplicated(ko_data$Sample)) {
  stop("Duplicate sample IDs found in KO abundance file.")
}

if (anyDuplicated(metadata$Sample)) {
  stop("Duplicate sample IDs found in metadata.")
}

if (anyNA(ko_data$Sample) || anyNA(metadata$Sample)) {
  stop("Missing sample IDs found in input data.")
}

# Select KO columns with standard KEGG Orthology identifiers.
ko_cols <- names(ko_data)[grepl("^K[0-9]{5}$", names(ko_data))]

if (length(ko_cols) == 0) {
  stop(
    "No KO columns found. Expected column names such as K00001."
  )
}

# Keep only the two study groups.
metadata <- metadata %>%
  filter(Group %in% c(SOTR_LABEL, HC_LABEL))

# Match sample IDs across the two input tables.
common_samples <- ko_data$Sample[
  ko_data$Sample %in% metadata$Sample
]

if (length(common_samples) == 0) {
  stop("No matching samples found between the input files.")
}

ko_data <- ko_data %>%
  filter(Sample %in% common_samples) %>%
  arrange(match(Sample, common_samples))

metadata <- metadata %>%
  filter(Sample %in% common_samples) %>%
  arrange(match(Sample, common_samples))

if (!identical(ko_data$Sample, metadata$Sample)) {
  stop("Sample alignment failed.")
}

if (!all(c(SOTR_LABEL, HC_LABEL) %in% metadata$Group)) {
  stop("Both SOTR and HC groups must be represented.")
}

message("Matched samples: ", nrow(ko_data))
message("KO features found: ", length(ko_cols))


# =============================================================================
# STEP 2. CLEAN KO ABUNDANCE DATA
# =============================================================================

message("\nSTEP 2: Cleaning KO abundance data")

for (ko in ko_cols) {
  ko_data[[ko]] <- suppressWarnings(as.numeric(ko_data[[ko]]))
}

ko_matrix <- as.matrix(ko_data[, ko_cols, drop = FALSE])
storage.mode(ko_matrix) <- "numeric"

if (any(is.infinite(ko_matrix), na.rm = TRUE)) {
  stop("Infinite values detected in KO abundance data.")
}

if (any(ko_matrix < 0, na.rm = TRUE)) {
  stop("Negative KO abundances are not allowed.")
}

# Preserve the original workflow: interpret missing abundance values as zero.
ko_matrix[is.na(ko_matrix)] <- 0

# Remove samples with zero total KO abundance.
valid_samples <- rowSums(ko_matrix) > 0

message(
  "Zero-total samples removed: ",
  sum(!valid_samples)
)

ko_matrix <- ko_matrix[valid_samples, , drop = FALSE]
metadata <- metadata[valid_samples, , drop = FALSE]
sample_ids <- ko_data$Sample[valid_samples]

if (nrow(ko_matrix) == 0) {
  stop("No samples remain after removing zero-total samples.")
}

if (!all(c(SOTR_LABEL, HC_LABEL) %in% metadata$Group)) {
  stop("Both groups must remain after sample filtering.")
}


# =============================================================================
# STEP 3. KO PREVALENCE FILTERING
# =============================================================================

message("\nSTEP 3: Filtering KOs by prevalence")

# Prevalence = proportion of samples with abundance > 0.
sotr_index <- metadata$Group == SOTR_LABEL
hc_index <- metadata$Group == HC_LABEL

prevalence_sotr <- colMeans(
  ko_matrix[sotr_index, , drop = FALSE] > 0
)

prevalence_hc <- colMeans(
  ko_matrix[hc_index, , drop = FALSE] > 0
)

# First filter:
# Retain a KO if prevalence is at least 5% in either group.
keep_initial <- (
  prevalence_sotr >= PREVALENCE_THRESHOLD_INITIAL |
  prevalence_hc >= PREVALENCE_THRESHOLD_INITIAL
)

ko_matrix_initial <- ko_matrix[, keep_initial, drop = FALSE]

if (ncol(ko_matrix_initial) == 0) {
  stop("No KOs passed the initial prevalence threshold.")
}

# Second filter:
# Retain KOs with combined prevalence strictly greater than 20%.
prevalence_combined <- colMeans(ko_matrix_initial > 0)

keep_final <- prevalence_combined > PREVALENCE_THRESHOLD_FINAL

ko_matrix_final <- ko_matrix_initial[, keep_final, drop = FALSE]

if (ncol(ko_matrix_final) == 0) {
  stop("No KOs passed the final prevalence threshold.")
}

# Save the filtered abundance table.
filtered_abundance <- data.frame(
  Sample = sample_ids,
  ko_matrix_final,
  check.names = FALSE
)

write_output(
  filtered_abundance,
  "ko_abundance_filtered.csv",
  "intermediate"
)

write_output(
  metadata,
  "metadata_matched.csv",
  "intermediate"
)

# Save prevalence statistics for all original KO features.
prevalence_table <- tibble(
  KO = colnames(ko_matrix),
  prevalence_SOTR = as.numeric(prevalence_sotr),
  prevalence_HC = as.numeric(prevalence_hc),
  retained_initial = keep_initial
) %>%
  mutate(
    prevalence_combined_after_initial_filter = NA_real_,
    retained_final = FALSE
  )

prevalence_table$prevalence_combined_after_initial_filter[
  match(names(prevalence_combined), prevalence_table$KO)
] <- as.numeric(prevalence_combined)

prevalence_table$retained_final[
  match(names(keep_final), prevalence_table$KO)
] <- as.logical(keep_final)

write_output(
  prevalence_table,
  "ko_prevalence.csv",
  "intermediate"
)

message("Samples retained: ", nrow(ko_matrix_final))
message("KOs retained: ", ncol(ko_matrix_final))


# =============================================================================
# STEP 4. RELATIVE ABUNDANCE AND CLR TRANSFORMATION
# =============================================================================

message("\nSTEP 4: Performing CLR transformation")

# Convert abundances to relative abundance per sample.
relative_abundance <- ko_matrix_final /
  rowSums(ko_matrix_final)

# Replace zeros using one-fifth of the smallest positive relative abundance.
positive_values <- relative_abundance[relative_abundance > 0]

if (length(positive_values) == 0) {
  stop("No positive relative abundance values available.")
}

pseudocount <- min(positive_values) / 5

relative_abundance[relative_abundance == 0] <- pseudocount

# Re-normalize each sample after zero replacement.
relative_abundance <- relative_abundance /
  rowSums(relative_abundance)

# CLR = log(value) - mean(log(values)) within each sample.
log_relative_abundance <- log(relative_abundance)

clr_matrix <- sweep(
  log_relative_abundance,
  1,
  rowMeans(log_relative_abundance),
  FUN = "-"
)

clr_data <- data.frame(
  Sample = sample_ids,
  clr_matrix,
  check.names = FALSE
)

write_output(
  clr_data,
  "ko_clr.csv",
  "intermediate"
)

message("CLR transformation complete.")
message("Zero-replacement pseudocount: ", signif(pseudocount, 5))


# =============================================================================
# STEP 5. FUNCTIONAL DYSBIOSIS SCORE
# =============================================================================

message("\nSTEP 5: Calculating functional dysbiosis scores")

# CLR samples and metadata must be in identical order.
metadata <- metadata[
  match(clr_data$Sample, metadata$Sample),
  ,
  drop = FALSE
]

if (anyNA(metadata$Sample) ||
    !identical(clr_data$Sample, metadata$Sample)) {
  stop("Sample mismatch between CLR data and metadata.")
}

clr_features <- setdiff(names(clr_data), "Sample")

clr_matrix <- as.matrix(
  clr_data[, clr_features, drop = FALSE]
)

storage.mode(clr_matrix) <- "numeric"

# Calculate group centroids in CLR space.
centroid_sotr <- colMeans(
  clr_matrix[metadata$Group == SOTR_LABEL, , drop = FALSE]
)

centroid_hc <- colMeans(
  clr_matrix[metadata$Group == HC_LABEL, , drop = FALSE]
)

# Euclidean distance from each sample to a centroid.
distance_to_centroid <- function(x, centroid) {
  difference <- sweep(x, 2, centroid, FUN = "-")
  sqrt(rowSums(difference^2))
}

distance_sotr <- distance_to_centroid(
  clr_matrix,
  centroid_sotr
)

distance_hc <- distance_to_centroid(
  clr_matrix,
  centroid_hc
)

# Positive scores indicate closer proximity to the SOTR centroid.
# Negative scores indicate closer proximity to the HC centroid.
dysbiosis_scores <- data.frame(
  Sample = clr_data$Sample,
  Group = metadata$Group,
  DistanceToSOTR = distance_sotr,
  DistanceToHC = distance_hc,
  FunctionalDysbiosis = distance_hc - distance_sotr
)

# Standardize scores within SOTR samples, as in the original workflow.
sotr_score_index <- dysbiosis_scores$Group == SOTR_LABEL

dysbiosis_scores$FunctionalDysbiosis_scaled <- NA_real_

if (sum(sotr_score_index) < 2) {
  stop("At least two SOTR samples are required to standardize scores.")
}

sotr_score_sd <- sd(
  dysbiosis_scores$FunctionalDysbiosis[sotr_score_index]
)

if (!is.finite(sotr_score_sd) || sotr_score_sd == 0) {
  stop("SOTR dysbiosis scores have zero or invalid standard deviation.")
}

dysbiosis_scores$FunctionalDysbiosis_scaled[sotr_score_index] <-
  as.numeric(
    scale(
      dysbiosis_scores$FunctionalDysbiosis[sotr_score_index]
    )
  )

write_output(
  dysbiosis_scores,
  "functional_dysbiosis_scores.csv"
)


# =============================================================================
# STEP 6. KO CORRELATION WITH FUNCTIONAL DYSBIOSIS
# =============================================================================

message("\nSTEP 6: Correlating KO profiles with dysbiosis scores")

sotr_scores <- dysbiosis_scores %>%
  filter(Group == SOTR_LABEL) %>%
  select(Sample, FunctionalDysbiosis_scaled)

# Explicitly match by Sample ID instead of relying on row order.
sotr_clr <- clr_data %>%
  inner_join(sotr_scores, by = "Sample")

correlation_results <- lapply(clr_features, function(ko) {

  x <- sotr_clr[[ko]]
  y <- sotr_clr$FunctionalDysbiosis_scaled

  valid <- is.finite(x) & is.finite(y)

  if (
    sum(valid) < 3 ||
    length(unique(x[valid])) < 2 ||
    length(unique(y[valid])) < 2
  ) {
    return(data.frame(
      KO = ko,
      correlation = NA_real_,
      p_value = NA_real_,
      n = sum(valid)
    ))
  }

  test <- cor.test(
    x[valid],
    y[valid],
    method = "pearson"
  )

  data.frame(
    KO = ko,
    correlation = unname(test$estimate),
    p_value = test$p.value,
    n = sum(valid)
  )
})

correlation_table <- bind_rows(correlation_results) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(FDR, p_value)

write_output(
  correlation_table,
  "ko_dysbiosis_correlations.csv"
)

significant_correlations <- correlation_table %>%
  filter(!is.na(FDR), FDR < FDR_THRESHOLD)

write_output(
  significant_correlations,
  "ko_dysbiosis_correlations_significant.csv"
)

message(
  "KOs with FDR < ",
  FDR_THRESHOLD,
  ": ",
  nrow(significant_correlations)
)

# Optional: join functional annotation.
if (
  nzchar(KO_ANNOTATION_FILE) &&
  file.exists(KO_ANNOTATION_FILE)
) {

  annotations <- read_csv(
    KO_ANNOTATION_FILE,
    show_col_types = FALSE
  )

  if (!"KO" %in% names(annotations)) {
    stop("KO annotation file must contain a KO column.")
  }

  annotated_correlations <- correlation_table %>%
    left_join(annotations, by = "KO")

  write_output(
    annotated_correlations,
    "ko_dysbiosis_correlations_annotated.csv"
  )

  write_output(
    filter(
      annotated_correlations,
      !is.na(FDR),
      FDR < FDR_THRESHOLD
    ),
    "ko_dysbiosis_correlations_significant_annotated.csv"
  )
}


# =============================================================================
# STEP 7. PRINCIPAL COMPONENT ANALYSIS (PCA)
# =============================================================================

message("\nSTEP 7: Performing PCA")

# Center CLR features, without scaling each feature to unit variance.
pca_result <- prcomp(
  clr_matrix,
  center = TRUE,
  scale. = FALSE
)

if (ncol(pca_result$x) < 2) {
  stop("PCA did not produce two principal components.")
}

variance_explained <- (
  pca_result$sdev^2 / sum(pca_result$sdev^2)
) * 100

# Sample scores.
pca_scores <- data.frame(
  Sample = clr_data$Sample,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = metadata$Group
) %>%
  left_join(
    dysbiosis_scores %>%
      select(Sample, FunctionalDysbiosis_scaled),
    by = "Sample"
  )

# KO loadings.
pca_loadings <- data.frame(
  KO = rownames(pca_result$rotation),
  PC1 = pca_result$rotation[, 1],
  PC2 = pca_result$rotation[, 2]
)

write_output(pca_scores, "pca_scores.csv")
write_output(pca_loadings, "pca_loadings.csv")

# Calculate group centroids in PCA score space.
pca_centroids <- pca_scores %>%
  group_by(Group) %>%
  summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  )

# Select up to 10 most positive and 10 most negative KO correlations.
ranked_correlations <- correlation_table %>%
  filter(is.finite(correlation)) %>%
  arrange(desc(correlation))

top_positive <- head(ranked_correlations$KO, 10)

top_negative <- head(
  ranked_correlations$KO[
    order(ranked_correlations$correlation)
  ],
  10
)

selected_kos <- unique(c(top_positive, top_negative))

loading_subset <- pca_loadings %>%
  filter(KO %in% selected_kos) %>%
  left_join(
    correlation_table %>% select(KO, correlation),
    by = "KO"
  ) %>%
  mutate(
    KO_group = ifelse(correlation >= 0, "Positive", "Negative")
  )


# Build PCA plot.
p <- ggplot() +

  # Healthy controls.
  geom_point(
    data = filter(pca_scores, Group == HC_LABEL),
    aes(x = PC1, y = PC2),
    color = "darkorange",
    size = 2.5,
    alpha = 0.85
  ) +

  # SOTR samples colored by scaled dysbiosis score.
  geom_point(
    data = filter(pca_scores, Group == SOTR_LABEL),
    aes(
      x = PC1,
      y = PC2,
      color = FunctionalDysbiosis_scaled
    ),
    size = 2.5,
    na.rm = TRUE
  ) +

  scale_color_gradientn(
    colors = c("darkblue", "yellow"),
    name = "Functional dysbiosis"
  ) +

  # Group centroids.
  geom_point(
    data = pca_centroids,
    aes(x = PC1, y = PC2, fill = Group),
    shape = 21,
    size = 4.5,
    color = "black"
  ) +

  scale_fill_manual(
    values = c(
      "HC" = "darkorange",
      "SOTR" = "steelblue"
    )
  ) +

  geom_text_repel(
    data = pca_centroids,
    aes(x = PC1, y = PC2, label = Group),
    inherit.aes = FALSE
  ) +

  labs(
    title = "PCA of KO functional profiles",
    x = sprintf("PC1 (%.1f%% variance)", variance_explained[1]),
    y = sprintf("PC2 (%.1f%% variance)", variance_explained[2])
  ) +

  theme_bw()

ggsave(
  filename = file.path(
    RESULTS_DIR,
    "figures",
    "functional_pca.png"
  ),
  plot = p,
  width = 10,
  height = 7,
  dpi = 300
)

message("PCA plot saved.")


# =============================================================================
# STEP 8. ANALYSIS SUMMARY
# =============================================================================

message("\nSTEP 8: Analysis summary")

summary_table <- tibble(
  Metric = c(
    "Samples analysed",
    "SOTR samples",
    "HC samples",
    "KOs after prevalence filtering",
    "KOs with FDR below threshold",
    "Initial prevalence threshold",
    "Final combined prevalence threshold",
    "FDR threshold"
  ),
  Value = as.character(c(
    nrow(clr_data),
    sum(metadata$Group == SOTR_LABEL),
    sum(metadata$Group == HC_LABEL),
    ncol(clr_matrix),
    nrow(significant_correlations),
    PREVALENCE_THRESHOLD_INITIAL,
    PREVALENCE_THRESHOLD_FINAL,
    FDR_THRESHOLD
  ))
)

write_output(
  summary_table,
  "analysis_summary.csv"
)

print(summary_table)

message("\nAnalysis pipeline completed.")
message("Results directory: ", RESULTS_DIR)

# =============================================================================
# END OF SCRIPT
# =============================================================================

