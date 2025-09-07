# This script evaluates the performance of a previously calculated Polygenic Risk Score (PRS)
# for the SWE cohort and creates new sample metadata with PRS-based subgroup labels.

# --- 1. Setup: Load Libraries and Define File Paths ---
library(dplyr)
library(ggplot2)
library(pROC)

# --- Input Files ---
swe_fam_raw_file <- "/path/to/raw_cohort_data/SWE.fam"
swe_fam_qcd_file <- "/path/to/swe_cohort/SWE_eur_afterQC.fam"
pcs_swe_file <- "/path/to/swe_cohort/PC6628.bigstatsr.txt"
prs_files_prefix_swe <- "/path/to/prs_results/SWE_prs_auto_p" # Note: Prefix only

# --- Output Files ---
auc_results_file <- "/path/to/output/SWE_AUC_results.txt"
swe_meta_with_prs_file <- "/path/to/output/SWE_metadata_with_PRS.txt"
swe_subgroups_file <- "/path/to/output/subgroups_SWE.txt"
swe_subgroups_summary_file <- "/path/to/output/subgroup_summary_SWE.txt"

# --- 2. Load and Prepare Metadata for the SWE Cohort ---
# (This section remains unchanged)
metaSWE <- read.delim(swe_fam_raw_file, sep = ' ', header = FALSE) %>%
  mutate(family.ID = V1, sample.ID = V2, vcfid = paste0(paste0(family.ID, "_", sample.ID), "_", paste0(family.ID, "_", sample.ID)), sex = V5, case_control = V6 - 1) %>%
  select(family.ID, sample.ID, vcfid, sex, case_control)

pcs_swe <- read.delim(pcs_swe_file, sep = ' ')
colnames(pcs_swe) <- c('family.ID', 'vcfid', paste0("V", 1:10))

newmetaSWE <- read.delim(swe_fam_qcd_file, sep = ' ', header = FALSE) %>%
  mutate(family.ID = V1, sample.ID = V2, vcfid = V2) %>%
  select(family.ID, sample.ID, vcfid) %>%
  left_join(metaSWE[, 3:5], by = "vcfid") %>%
  left_join(pcs_swe, by = c("family.ID", "vcfid")) %>%
  mutate(sex = if_else(sex == 2, 0, sex))


# --- 3. PRS Performance Evaluation (AUC Calculation) ---

# Define the p-value threshold(s) for loading PRS results.
# The commented-out line can be used to evaluate multiple PRS versions as a sensitivity analysis.
prs_p_thresholds <- "1.0" # Standard analysis: Use the p=1.0 PRS
# prs_p_thresholds <- c("1.0", "0.5", "0.3", "0.2", "0.1", "0.05") # Sensitivity analysis: Uncomment to evaluate multiple thresholds

# Generate file paths based on the thresholds
prs_files_swe <- paste0(prs_files_prefix_swe, prs_p_thresholds, ".txt")

# Load all specified PRS files
prs_list_swe <- lapply(prs_files_swe, read.delim, header = FALSE)

# Calculate AUC for each loaded PRS
swe_auc_scores <- sapply(prs_list_swe, function(prs) auc(newmetaSWE$case_control, prs$V1))

# Combine AUC results into a data frame and save
swe_auc_results <- data.frame(
  p_thresh = as.numeric(prs_p_thresholds),
  SWE_AUC = swe_auc_scores
)
write.table(swe_auc_results, auc_results_file, quote = FALSE, row.names = FALSE)

# Plot the AUC results
ggplot(swe_auc_results, aes(x = as.character(p_thresh), y = SWE_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_cartesian(ylim = c(0.75, 0.87)) +
  ggtitle("SWE Cohort PRS Performance (AUC)") +
  xlab("GWAS P-value Threshold for SNP Inclusion") +
  ylab("AUC")


# --- 4. Append PRS and Define Subgroups ---

# Append the PRS to the metadata for subgrouping.
# If multiple thresholds were tested, this now correctly uses the FIRST one in the list (p=1.0 by default).
newmetaSWE$PRS <- prs_list_swe[[1]]$V1

# Save the metadata with the appended PRS scores
write.table(newmetaSWE, swe_meta_with_prs_file, quote = FALSE, row.names = FALSE)

# --- (The rest of the script for subgroup definition remains the same) ---

# --- Function to Define PRS-based Subgroups ---
define_prs_subgroups <- function(metadata) {
  # ... (function content is unchanged) ...
  subgroup_meta <- metadata %>%
    rename(prs = PRS,
           X1 = V1, X2 = V2, X3 = V3, X4 = V4, X5 = V5,
           X6 = V6, X7 = V7, X8 = V8, X9 = V9, X10 = V10)
  prsSCZ_sorted <- sort(subgroup_meta$prs[subgroup_meta$case_control == 1])
  prsCtr_sorted <- sort(subgroup_meta$prs[subgroup_meta$case_control == 0])
  subgroup_meta <- subgroup_meta %>%
    mutate(
      prsSCZ = case_when(case_control == 1 ~ prs),
      prsCtr = case_when(case_control == 0 ~ prs),
      TB10all = case_when(prs > quantile(prs, 0.875) & prs < quantile(prs, 0.975) ~ 1,
                          prs < quantile(prs, 0.125) & prs > quantile(prs, 0.025) ~ 0),
      TB15all = case_when(prs > quantile(prs, 0.825) & prs < quantile(prs, 0.975) ~ 1,
                          prs < quantile(prs, 0.175) & prs > quantile(prs, 0.025) ~ 0),
      TB20all = case_when(prs > quantile(prs, 0.775) & prs < quantile(prs, 0.975) ~ 1,
                          prs < quantile(prs, 0.225) & prs > quantile(prs, 0.025) ~ 0)
    )
  sprs <- sort(subgroup_meta$prs)
  BL <- quantile(prsSCZ_sorted, 0.025, na.rm = TRUE)
  TU <- quantile(prsCtr_sorted, 0.975, na.rm = TRUE)
  BU10 <- sprs[which.min(abs(sprs - BL)) + round(length(sprs) * 0.10)]
  BU15 <- sprs[which.min(abs(sprs - BL)) + round(length(sprs) * 0.15)]
  BU20 <- sprs[which.min(abs(sprs - BL)) + round(length(sprs) * 0.20)]
  TL10 <- sprs[which.min(abs(sprs - TU)) - round(length(sprs) * 0.10)]
  TL15 <- sprs[which.min(abs(sprs - TU)) - round(length(sprs) * 0.15)]
  TL20 <- sprs[which.min(abs(sprs - TU)) - round(length(sprs) * 0.20)]
  subgroup_meta <- subgroup_meta %>%
    mutate(
      T10 = case_when(prs >= TL10 & prs <= TU ~ case_control),
      B10 = case_when(prs >= BL & prs <= BU10 ~ case_control),
      T15 = case_when(prs >= TL15 & prs <= TU ~ case_control),
      B15 = case_when(prs >= BL & prs <= BU15 ~ case_control),
      T20 = case_when(prs >= TL20 & prs <= TU ~ case_control),
      B20 = case_when(prs >= BL & prs <= BU20 ~ case_control)
    )
  return(subgroup_meta)
}

# --- 5. Generate and Save Final Subgroup Files for SWE Cohort ---
subgroups_swe <- define_prs_subgroups(newmetaSWE)
write.table(subgroups_swe, file = swe_subgroups_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

subgroup_cols_swe <- colnames(subgroups_swe[, c("case_control", "TB10all", "TB15all", "TB20all", "T10", "B10", "T15", "B15", "T20", "B20")])
summary_swe <- data.frame(
  subgroup = subgroup_cols_swe,
  n = sapply(subgroups_swe[subgroup_cols_swe], function(x) sum(!is.na(x))),
  n_case = sapply(subgroups_swe[subgroup_cols_swe], function(x) sum(x == 1, na.rm = TRUE)),
  n_control = sapply(subgroups_swe[subgroup_cols_swe], function(x) sum(x == 0, na.rm = TRUE))
)
write.table(summary_swe, file = swe_subgroups_summary_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

print("Script finished. SWE metadata with PRS and subgroup labels have been saved.")