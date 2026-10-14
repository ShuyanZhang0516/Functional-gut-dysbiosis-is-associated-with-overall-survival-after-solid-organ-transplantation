# Combined survival analyses:
#   A. Overall functional dysbiosis and mortality
#   B. Individual KO associations with mortality
#
# Shared cohort:
#   - SOTR samples matched to mortality and clinical metadata
#   - months_since_tx > 11
#   - KO analysis additionally excludes samples with all-zero KO profiles
#
# Main outputs:
#   results/tables/
#   results/figures/
#   results/intermediate/

# ============================================================
# 0. Setup
# ============================================================

required_packages <- c(
  "readr", "dplyr", "tibble", "tidyr",
  "ggplot2", "survival", "survminer", "cowplot"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install missing packages before running: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)
library(cowplot)

# Required environment variables:
# MORTALITY_FILE
# SAMPLE_ID_FILE
# METADATA_FILE
# MEDICATION_FILE
# KO_ABUNDANCE_FILE
# DYSBIOSIS_FILE
#
# Optional:
# RESULTS_DIR (default: results)
#
# Example:
# Sys.setenv(
#   MORTALITY_FILE = "data/mortality_transplantationtime.csv",
#   SAMPLE_ID_FILE = "data/Sample_ID.csv",
#   METADATA_FILE = "data/Meta_cross_organ_v4.rds",
#   MEDICATION_FILE = "data/medication_regimen_ID.csv",
#   KO_ABUNDANCE_FILE = "data/merge_ko_coverage_all_normalised.csv",
#   DYSBIOSIS_FILE = "data/functional_dysbiosis_scores.csv",
#   RESULTS_DIR = "results"
# )

get_required_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (value == "") {
    stop("Please set environment variable: ", name)
  }
  if (!file.exists(value)) {
    stop("Input file does not exist: ", value)
  }
  value
}

mortality_file <- get_required_env("MORTALITY_FILE")
sample_id_file <- get_required_env("SAMPLE_ID_FILE")
metadata_file <- get_required_env("METADATA_FILE")
medication_file <- get_required_env("MEDICATION_FILE")
ko_abundance_file <- get_required_env("KO_ABUNDANCE_FILE")
dysbiosis_file <- get_required_env("DYSBIOSIS_FILE")

results_dir <- Sys.getenv("RESULTS_DIR", unset = "results")
tables_dir <- file.path(results_dir, "tables")
figures_dir <- file.path(results_dir, "figures")
intermediate_dir <- file.path(results_dir, "intermediate")

for (directory in c(tables_dir, figures_dir, intermediate_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

# Analysis settings
MIN_MONTHS_SINCE_TX <- 11
MIN_KO_PREVALENCE <- 0.50
ZERO_REPLACEMENT_DIVISOR <- 5

# ============================================================
# 1. Read and validate inputs
# ============================================================

message("Reading input files...")

mortality <- read_csv(mortality_file, show_col_types = FALSE)
sample_id <- read_csv(sample_id_file, show_col_types = FALSE)
metadata <- as_tibble(readRDS(metadata_file))
medication <- read_csv(medication_file, show_col_types = FALSE)
ko_abundance <- read_csv(ko_abundance_file, show_col_types = FALSE)
dysbiosis <- read_csv(dysbiosis_file, show_col_types = FALSE)

check_columns <- function(data, required, data_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      data_name, " is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

check_columns(
  mortality,
  c("ID", "months_since_tx", "event", "survival_time_months"),
  "Mortality data"
)
check_columns(sample_id, c("ID", "Sample"), "Sample ID data")
check_columns(metadata, c("ID", "Age", "GESLACHT", "BMI"), "Metadata")
check_columns(
  medication,
  c("ID", "medication_regimen"),
  "Medication data"
)
check_columns(ko_abundance, "Sample", "KO abundance data")
check_columns(
  dysbiosis,
  c("Sample", "Group", "FunctionalDysbiosis"),
  "Functional dysbiosis data"
)

# ============================================================
# 2. Prepare shared clinical cohort
# ============================================================

# Retain the first mortality record per participant, as in the original
# analysis. Confirm that this is appropriate for your study design.
mortality <- mortality %>%
  distinct(ID, .keep_all = TRUE)

sample_id <- sample_id %>%
  distinct(ID, Sample)

metadata <- metadata %>%
  select(ID, Age, GESLACHT, BMI) %>%
  distinct(ID, .keep_all = TRUE)

medication <- medication %>%
  select(ID, medication_regimen) %>%
  distinct(ID, .keep_all = TRUE)

clinical_all <- mortality %>%
  left_join(sample_id, by = "ID") %>%
  filter(!is.na(Sample), Sample != "") %>%
  left_join(metadata, by = "ID") %>%
  left_join(medication, by = "ID") %>%
  filter(months_since_tx > MIN_MONTHS_SINCE_TX) %>%
  distinct(ID, .keep_all = TRUE)

if (nrow(clinical_all) == 0) {
  stop("No participants remain after clinical cohort filtering.")
}

if (anyDuplicated(clinical_all$Sample)) {
  stop(
    "Duplicate Sample IDs occur in the clinical cohort. ",
    "Resolve the sample-to-participant mapping before analysis."
  )
}

clinical_all <- clinical_all %>%
  mutate(
    event = as.integer(event),
    survival_time_months = as.numeric(survival_time_months),
    months_since_tx = as.numeric(months_since_tx),
    Age = as.numeric(Age),
    BMI = as.numeric(BMI),
    GESLACHT = factor(GESLACHT),
    medication_regimen = factor(medication_regimen)
  )

if (any(!clinical_all$event %in% c(0L, 1L), na.rm = TRUE)) {
  stop("event must be coded as 0 = censored and 1 = event.")
}

if (any(
  clinical_all$survival_time_months <= 0,
  na.rm = TRUE
)) {
  stop("survival_time_months must be positive.")
}

write_csv(
  clinical_all,
  file.path(intermediate_dir, "clinical_cohort_over_11_months.csv")
)

message("Clinical cohort after filtering: ", nrow(clinical_all))

# ============================================================
# 3. Analysis A: functional dysbiosis and mortality
# ============================================================

message("Starting functional dysbiosis survival analysis...")

# Restrict dysbiosis data to SOTR samples and match by sample ID.
sotr_dysbiosis <- dysbiosis %>%
  filter(toupper(as.character(Group)) == "SOTR") %>%
  select(Sample, FunctionalDysbiosis) %>%
  distinct(Sample, .keep_all = TRUE) %>%
  mutate(FunctionalDysbiosis = as.numeric(FunctionalDysbiosis))

dysbiosis_data <- clinical_all %>%
  inner_join(sotr_dysbiosis, by = "Sample") %>%
  filter(
    !is.na(FunctionalDysbiosis),
    is.finite(FunctionalDysbiosis)
  ) %>%
  mutate(
    age_c = Age - mean(Age, na.rm = TRUE),
    bmi_c = BMI - mean(BMI, na.rm = TRUE)
  )

if (nrow(dysbiosis_data) == 0) {
  stop(
    "No SOTR samples matched between the clinical cohort and ",
    "functional dysbiosis file."
  )
}

if (sd(dysbiosis_data$FunctionalDysbiosis) == 0) {
  stop("FunctionalDysbiosis has zero variance in the matched cohort.")
}

dysbiosis_data <- dysbiosis_data %>%
  mutate(
    FunctionalDysbiosis_scaled = as.numeric(
      scale(FunctionalDysbiosis)
    ),
    high_low = if_else(
      FunctionalDysbiosis >= median(FunctionalDysbiosis, na.rm = TRUE),
      "High", "Low"
    ),
    high_low = factor(high_low, levels = c("Low", "High"))
  )

write_csv(
  dysbiosis_data,
  file.path(intermediate_dir, "functional_dysbiosis_survival_cohort.csv")
)

# ------------------------------------------------------------
# 3A. Event percentage by dysbiosis group
# ------------------------------------------------------------

event_summary <- dysbiosis_data %>%
  group_by(high_low) %>%
  summarise(
    N = n(),
    Events = sum(event == 1, na.rm = TRUE),
    Event_Percentage = 100 * mean(event == 1, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  event_summary,
  file.path(tables_dir, "dysbiosis_event_summary.csv")
)

event_plot <- ggplot(
  event_summary,
  aes(x = high_low, y = Event_Percentage, fill = high_low)
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = sprintf("%.1f%%", Event_Percentage)),
    vjust = -0.4,
    size = 4
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Mortality by functional dysbiosis group",
    x = "Functional dysbiosis group",
    y = "Participants with event (%)",
    fill = "Group"
  ) +
  theme_classic()

ggsave(
  file.path(figures_dir, "dysbiosis_event_percentage.png"),
  event_plot,
  width = 6,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 3B. Kaplan-Meier analysis
# ------------------------------------------------------------

km_data <- dysbiosis_data %>%
  filter(
    !is.na(survival_time_months),
    survival_time_months > 0,
    !is.na(event),
    !is.na(high_low)
  )

km_fit <- survfit(
  Surv(survival_time_months, event) ~ high_low,
  data = km_data
)

km_test <- survdiff(
  Surv(survival_time_months, event) ~ high_low,
  data = km_data
)

km_p <- pchisq(
  km_test$chisq,
  df = length(km_test$n) - 1,
  lower.tail = FALSE
)

write_csv(
  tibble(
    Test = "Log-rank test",
    Chi_square = unname(km_test$chisq),
    P_value = km_p,
    N = nrow(km_data),
    Events = sum(km_data$event == 1)
  ),
  file.path(tables_dir, "dysbiosis_kaplan_meier_logrank.csv")
)

km_plot <- ggsurvplot(
  km_fit,
  data = km_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Follow-up time (months)",
  ylab = "Survival probability",
  legend.title = "Dysbiosis group",
  legend.labs = c("Low", "High"),
  palette = c("#377EB8", "#E41A1C"),
  ggtheme = theme_classic()
)

ggsave(
  file.path(figures_dir, "dysbiosis_kaplan_meier.png"),
  km_plot$plot,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(figures_dir, "dysbiosis_kaplan_meier_risk_table.png"),
  arrange_ggsurvplots(km_plot, print = FALSE),
  width = 8,
  height = 8,
  dpi = 300
)

combined_plot <- plot_grid(
  km_plot$plot,
  event_plot,
  ncol = 2,
  rel_widths = c(1.3, 1)
)

ggsave(
  file.path(figures_dir, "dysbiosis_survival_combined.png"),
  combined_plot,
  width = 12,
  height = 5.5,
  dpi = 300
)

# ------------------------------------------------------------
# 3C. Standard Cox model and time-dependent Cox model
# ------------------------------------------------------------

dysbiosis_model_data <- dysbiosis_data %>%
  filter(
    !is.na(survival_time_months),
    survival_time_months > 0,
    !is.na(event),
    !is.na(FunctionalDysbiosis_scaled),
    !is.na(months_since_tx),
    !is.na(age_c),
    !is.na(bmi_c),
    !is.na(GESLACHT),
    !is.na(medication_regimen)
  )

if (nrow(dysbiosis_model_data) < 20 ||
    dplyr::n_distinct(dysbiosis_model_data$event) < 2) {
  stop("Insufficient complete cases or events for dysbiosis Cox models.")
}

cox_dysbiosis_standard <- coxph(
  Surv(survival_time_months, event) ~
    FunctionalDysbiosis_scaled +
    months_since_tx +
    age_c + bmi_c + GESLACHT + medication_regimen,
  data = dysbiosis_model_data,
  ties = "efron",
  x = TRUE
)

cox_dysbiosis_td <- coxph(
  Surv(survival_time_months, event) ~
    FunctionalDysbiosis_scaled +
    tt(months_since_tx) +
    age_c + bmi_c + GESLACHT + medication_regimen,
  data = dysbiosis_model_data,
  tt = function(x, t, ...) x * log(t + 1),
  ties = "efron",
  x = TRUE
)

extract_cox_results <- function(model, model_name) {
  coef_table <- summary(model)$coefficients
  ci <- confint(model)

  tibble(
    Term = rownames(coef_table),
    Model = model_name,
    Estimate = coef_table[, "coef"],
    Hazard_Ratio = exp(coef_table[, "coef"]),
    StdError = coef_table[, "se(coef)"],
    Z = coef_table[, "z"],
    P_value = coef_table[, "Pr(>|z|)"],
    CI_Lower = exp(ci[, 1]),
    CI_Upper = exp(ci[, 2])
  )
}

dysbiosis_cox_results <- bind_rows(
  extract_cox_results(cox_dysbiosis_standard, "Standard Cox"),
  extract_cox_results(cox_dysbiosis_td, "Time-dependent covariate Cox")
)

write_csv(
  dysbiosis_cox_results,
  file.path(tables_dir, "dysbiosis_cox_results.csv")
)

# Proportional hazards diagnostics for the standard Cox model.
ph_test <- cox.zph(cox_dysbiosis_standard)

ph_table <- as.data.frame(ph_test$table) %>%
  rownames_to_column("Term")

write_csv(
  ph_table,
  file.path(tables_dir, "dysbiosis_standard_cox_ph_test.csv")
)

pdf(
  file.path(figures_dir, "dysbiosis_standard_cox_ph_diagnostics.pdf"),
  width = 8,
  height = 7
)
plot(ph_test)
dev.off()

# Dysbiosis-specific effect estimates for a compact comparison.
dysbiosis_hr_comparison <- dysbiosis_cox_results %>%
  filter(Term == "FunctionalDysbiosis_scaled") %>%
  select(
    Model, Hazard_Ratio, CI_Lower, CI_Upper, P_value
  )

write_csv(
  dysbiosis_hr_comparison,
  file.path(tables_dir, "dysbiosis_hazard_ratio_comparison.csv")
)

# ============================================================
# 4. Analysis B: individual KO associations with mortality
# ============================================================

message("Starting KO-by-KO survival analysis...")

# Keep clinical samples with KO abundance profiles.
if (anyDuplicated(ko_abundance$Sample)) {
  stop(
    "Duplicate Sample IDs in KO abundance table. ",
    "Resolve duplicates before analysis."
  )
}

# Adjust this expression if the KO columns use a different naming scheme.
ko_cols <- names(ko_abundance)[grepl("^K[0-9]{5}$", names(ko_abundance))]

if (length(ko_cols) == 0) {
  stop(
    "No KO columns found using '^K[0-9]{5}$'. ",
    "Update the KO column selection rule to match your file."
  )
}

ko_input <- ko_abundance %>%
  filter(Sample %in% clinical_all$Sample)

if (nrow(ko_input) == 0) {
  stop("No overlap between KO abundance samples and clinical samples.")
}

ko_matrix <- as.matrix(ko_input[, ko_cols, drop = FALSE])
storage.mode(ko_matrix) <- "numeric"
rownames(ko_matrix) <- ko_input$Sample

if (any(is.infinite(ko_matrix), na.rm = TRUE)) {
  stop("KO abundance matrix contains infinite values.")
}

ko_matrix[is.na(ko_matrix)] <- 0

if (any(ko_matrix < 0)) {
  stop("KO abundances must be non-negative.")
}

# Remove all-zero KO profiles.
nonzero_samples <- rowSums(ko_matrix) > 0
message("KO samples before all-zero filtering: ", nrow(ko_matrix))
message("All-zero samples removed: ", sum(!nonzero_samples))

ko_matrix <- ko_matrix[nonzero_samples, , drop = FALSE]

ko_clinical <- clinical_all %>%
  filter(Sample %in% rownames(ko_matrix))

ko_clinical <- ko_clinical[
  match(rownames(ko_matrix), ko_clinical$Sample),
  ,
  drop = FALSE
]

stopifnot(identical(ko_clinical$Sample, rownames(ko_matrix)))

write_csv(
  ko_clinical,
  file.path(intermediate_dir, "ko_survival_included_samples.csv")
)

# ------------------------------------------------------------
# 4A. KO prevalence filtering
# ------------------------------------------------------------

ko_prevalence <- colMeans(ko_matrix > 0, na.rm = TRUE)

prevalence_table <- tibble(
  KO = names(ko_prevalence),
  Prevalence = as.numeric(ko_prevalence)
) %>%
  arrange(desc(Prevalence))

write_csv(
  prevalence_table,
  file.path(tables_dir, "ko_survival_prevalence.csv")
)

keep_kos <- names(ko_prevalence)[
  ko_prevalence > MIN_KO_PREVALENCE
]

if (length(keep_kos) == 0) {
  stop("No KOs passed the prevalence threshold.")
}

ko_matrix <- ko_matrix[, keep_kos, drop = FALSE]

message("KOs passing prevalence filter: ", ncol(ko_matrix))

# ------------------------------------------------------------
# 4B. Relative abundance and CLR transformation
# ------------------------------------------------------------

relative_abundance <- sweep(
  ko_matrix,
  1,
  rowSums(ko_matrix),
  FUN = "/"
)

positive_values <- relative_abundance[relative_abundance > 0]

if (length(positive_values) == 0) {
  stop("No positive relative abundance values are available.")
}

zero_replacement <- min(positive_values) / ZERO_REPLACEMENT_DIVISOR

relative_abundance[relative_abundance == 0] <- zero_replacement

# Re-normalize after replacing zeros.
relative_abundance <- sweep(
  relative_abundance,
  1,
  rowSums(relative_abundance),
  FUN = "/"
)

log_abundance <- log(relative_abundance)

ko_clr <- sweep(
  log_abundance,
  1,
  rowMeans(log_abundance),
  FUN = "-"
)

stopifnot(all(is.finite(ko_clr)))

# ------------------------------------------------------------
# 4C. Fit one Cox model per KO
# ------------------------------------------------------------

ko_results <- list()

failed_kos <- tibble(
  KO = character(),
  ErrorMessage = character()
)

n_kos <- ncol(ko_clr)

for (i in seq_len(n_kos)) {

  ko <- colnames(ko_clr)[i]

  if (i %% 100 == 0 || i == 1 || i == n_kos) {
    message("Fitting KO ", i, "/", n_kos, ": ", ko)
  }

  dat <- ko_clinical
  dat$KO_value <- as.numeric(ko_clr[, ko])
  dat$KO_scaled <- as.numeric(scale(dat$KO_value))

  required_model_cols <- c(
    "survival_time_months", "event", "KO_scaled",
    "months_since_tx", "Age", "BMI",
    "GESLACHT", "medication_regimen", "ID"
  )

  dat <- dat %>%
    filter(if_all(all_of(required_model_cols), ~ !is.na(.x))) %>%
    mutate(
      age_c = Age - mean(Age, na.rm = TRUE),
      bmi_c = BMI - mean(BMI, na.rm = TRUE)
    ) %>%
    filter(
      is.finite(KO_scaled),
      is.finite(survival_time_months),
      is.finite(months_since_tx),
      survival_time_months > 0,
      event %in% c(0L, 1L)
    )

  if (nrow(dat) < 20 || dplyr::n_distinct(dat$event) < 2) {
    failed_kos <- bind_rows(
      failed_kos,
      tibble(
        KO = ko,
        ErrorMessage = "Insufficient complete cases or only one event class"
      )
    )
    next
  }

  if (!is.finite(sd(dat$KO_scaled)) || sd(dat$KO_scaled) == 0) {
    failed_kos <- bind_rows(
      failed_kos,
      tibble(KO = ko, ErrorMessage = "KO has zero variance")
    )
    next
  }

  fit_result <- tryCatch({

    model <- coxph(
      Surv(survival_time_months, event) ~
        KO_scaled +
        tt(months_since_tx) +
        age_c + bmi_c + GESLACHT + medication_regimen,
      data = dat,
      tt = function(x, t, ...) x * log(t + 1),
      ties = "efron",
      x = TRUE
    )

    coef_table <- summary(model)$coefficients
    ci <- confint(model)

    if (!"KO_scaled" %in% rownames(coef_table)) {
      stop("KO coefficient was not estimable.")
    }

    tibble(
      KO = ko,
      N = nrow(dat),
      Events = sum(dat$event == 1),
      Estimate = unname(coef_table["KO_scaled", "coef"]),
      Hazard_Ratio = exp(unname(coef_table["KO_scaled", "coef"])),
      StdError = unname(coef_table["KO_scaled", "se(coef)"]),
      Z = unname(coef_table["KO_scaled", "z"]),
      P_value = unname(coef_table["KO_scaled", "Pr(>|z|)"]),
      CI_Lower = exp(unname(ci["KO_scaled", 1])),
      CI_Upper = exp(unname(ci["KO_scaled", 2]))
    )

  }, error = function(e) {

    failed_kos <<- bind_rows(
      failed_kos,
      tibble(KO = ko, ErrorMessage = conditionMessage(e))
    )

    NULL
  })

  if (!is.null(fit_result)) {
    ko_results[[ko]] <- fit_result
  }
}

# ------------------------------------------------------------
# 4D. Multiple-testing correction and result tables
# ------------------------------------------------------------

if (length(ko_results) == 0) {
  write_csv(
    failed_kos,
    file.path(tables_dir, "ko_survival_failed_models.csv")
  )
  stop("No KO models completed successfully.")
}

ko_survival_results <- bind_rows(ko_results) %>%
  mutate(
    p_adj_BH = p.adjust(P_value, method = "BH"),
    p_adj_bonferroni = p.adjust(P_value, method = "bonferroni")
  ) %>%
  arrange(P_value)

write_csv(
  ko_survival_results,
  file.path(tables_dir, "ko_survival_cox_results.csv")
)

write_csv(
  filter(ko_survival_results, p_adj_BH < 0.05),
  file.path(tables_dir, "ko_survival_significant_BH.csv")
)

write_csv(
  filter(ko_survival_results, p_adj_bonferroni < 0.05),
  file.path(tables_dir, "ko_survival_significant_bonferroni.csv")
)

write_csv(
  failed_kos,
  file.path(tables_dir, "ko_survival_failed_models.csv")
)

# ------------------------------------------------------------
# 4E. KO result summary
# ------------------------------------------------------------

ko_summary <- tibble(
  Metric = c(
    "Samples analysed in KO analysis",
    "KOs passing prevalence filter",
    "KO models successfully fitted",
    "KO models failed",
    "KOs with nominal P < 0.05",
    "KOs with BH-adjusted P < 0.05",
    "KOs with Bonferroni-adjusted P < 0.05"
  ),
  Value = c(
    nrow(ko_matrix),
    ncol(ko_matrix),
    nrow(ko_survival_results),
    nrow(failed_kos),
    sum(ko_survival_results$P_value < 0.05, na.rm = TRUE),
    sum(ko_survival_results$p_adj_BH < 0.05, na.rm = TRUE),
    sum(ko_survival_results$p_adj_bonferroni < 0.05, na.rm = TRUE)
  )
)

write_csv(
  ko_summary,
  file.path(tables_dir, "ko_survival_summary.csv")
)

# ============================================================
# 5. Combined run summary
# ============================================================

combined_summary <- tibble(
  Analysis = c(
    "Functional dysbiosis",
    "Individual KO"
  ),
  Samples_analysed = c(
    nrow(dysbiosis_model_data),
    nrow(ko_matrix)
  ),
  Events = c(
    sum(dysbiosis_model_data$event == 1),
    sum(ko_clinical$event[match(rownames(ko_matrix), ko_clinical$Sample)] == 1,
        na.rm = TRUE)
  )
)

write_csv(
  combined_summary,
  file.path(tables_dir, "survival_analysis_summary.csv")
)

print(combined_summary)
print(ko_summary)

message("Both survival analyses completed.")
message("Tables: ", normalizePath(tables_dir))
message("Figures: ", normalizePath(figures_dir))
message("Intermediate files: ", normalizePath(intermediate_dir))

