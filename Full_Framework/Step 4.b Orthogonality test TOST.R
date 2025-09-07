# This script performs an orthogonality test to determine if a Polygenic Risk Score (PRS)
# is statistically independent of (i.e., not correlated with) specific biological pathway
# activities, as measured by Gene Set Variation Analysis (GSVA) scores.
#
# The statistical method used is the Two One-Sided T-tests (TOST) for equivalence,
# which tests if a correlation is statistically equivalent to zero within a defined margin (delta).

# --- 1. Setup: Load Libraries and Define Paths ---

# Load necessary libraries
library(tidyverse)
library(stringr)
library(TOSTER) # For equivalence testing

# --- Input Files ---
# Metadata files for different cohorts
meta_2722_file <- "/path/to/metadata/subgrp_ADcmeta0217.txt"
meta_2854_file <- "/path/to/metadata/subgrp_repmeta0217.txt"
meta_swe_file <- "/path/to/metadata/subgrp_SWEmeta0217.txt"

# GSVA score matrices (genesets x individuals)
gsva_ad_go_file <- "/path/to/gsva_matrices/go_2722_inAPOE_TWAS_in_PRS_X_matrix.txt"
gsva_ad_canonical_file <- "/path/to/gsva_matrices/canonical_2854_inAPOE_TWAS_in_PRS_X_matrix.txt"
gsva_ad_brain1_file <- "/path/to/gsva_matrices/brain_cell_type_function_AD_pathways_2722_exAPOE_TWAS_in_PRS_X_matrix.txt"
gsva_ad_brain2_file <- "/path/to/gsva_matrices/brain_cell_type_function_AD_pathways_2854_inAPOE_TWAS_in_PRS_X_matrix.txt"
gsva_scz_go_file <- "/path/to/gsva_matrices/go_SWE_TWAS_in_PRS_X_matrix.txt"

# --- Output Files ---
output_dir <- "/path/to/orthogonality_test_results/"


# --- 2. Define Reusable Analysis Function ---

#' Performs a TOST equivalence test for the correlation between PRS and GSVA scores.
#'
#' @param pathway_name The exact name of the pathway in the GSVA matrix.
#' @param gsva_data The full GSVA data frame (pathways x samples).
#' @param gsva_path_info A data frame mapping GSVA rownames to tissue and pathway names.
#' @param metadata The metadata data frame for the cohort.
#' @param prs_column_name The name of the column containing the PRS scores.
#' @param delta The equivalence margin for the correlation test.
#' @return A numeric vector of p-values for the 13 tissues.
perform_tost_analysis <- function(pathway_name, gsva_data, gsva_path_info, metadata, prs_column_name, delta) {
  
  # Extract the GSVA scores for the specified pathway, transposed to samples x tissues
  gsva_scores <- t(gsva_data[gsva_path_info$pathway == pathway_name, ])
  colnames(gsva_scores) <- gsva_path_info$tissue[1:13]
  
  # Get the PRS scores from the metadata
  prs_scores <- metadata[[prs_column_name]]
  
  # Loop through each of the 13 tissues and perform the equivalence test
  p_values <- numeric(13)
  for (i in 1:13) {
    tissue_gsva_scores <- gsva_scores[, i]
    res <- z_cor_test(
      x = prs_scores,
      y = tissue_gsva_scores,
      alternative = "e", # "e" for equivalence
      null = c(delta, -delta) # The equivalence bounds
    )
    p_values[i] <- res$p.value
  }
  
  return(p_values)
}


# --- 3. Load Data ---
meta2722 <- read.delim(meta_2722_file, header = TRUE)
meta2854 <- read.delim(meta_2854_file, header = TRUE)
metaSWE <- read.delim(meta_swe_file, header = TRUE)


# --- 4. Power Analysis ---
# Perform a power analysis to determine the sample size needed to detect
# equivalence within a given margin (delta), or vice versa.
power_res_d050 <- power_z_cor(alpha = 0.05, power = 0.70, rho = 0, null = c(0.05, -0.05))
power_res_d033 <- power_z_cor(alpha = 0.05, power = 0.70, rho = 0, null = c(0.033, -0.033))
print("Power analysis result for delta = 0.05:")
print(power_res_d050)
print("Power analysis result for delta = 0.033:")
print(power_res_d033)


# --- 5. Run Analyses for Each Dataset ---

## Analysis 1: AD Gene Ontology Pathways
print("Running TOST for AD GO pathways...")
GSVA_GO <- read.delim(gsva_ad_go_file)
GSVA_GO_path <- data.frame(str_split_fixed(rownames(GSVA_GO), ":", n = 2), stringsAsFactors = FALSE)
colnames(GSVA_GO_path) <- c("tissue", "pathway")
tissues <- unique(GSVA_GO_path$tissue)

ad_go_results <- data.frame(
  tissue = tissues,
  AB_eq_p = perform_tost_analysis("BP_amyloid_beta_clearance", GSVA_GO, GSVA_GO_path, meta2722, "Kprs", 0.05),
  BP_eq_p = perform_tost_analysis("BP_divalent_inorganic_cation_homeostasis", GSVA_GO, GSVA_GO_path, meta2722, "Kprs", 0.05),
  MF_eq_p = perform_tost_analysis("MF_lipoprotein_particle_receptor_binding", GSVA_GO, GSVA_GO_path, meta2722, "Kprs", 0.05)
)
write.table(ad_go_results, paste0(output_dir, "AD_GO_delta_0.05_2722kprs.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


## Analysis 2: AD Canonical Pathway
print("Running TOST for AD Canonical pathway...")
GSVA_canonical <- read.delim(gsva_ad_canonical_file)
GSVA_canonical_path <- data.frame(str_split_fixed(rownames(GSVA_canonical), ":", n = 2), stringsAsFactors = FALSE)
colnames(GSVA_canonical_path) <- c("tissue", "pathway")

ad_kegg_results <- data.frame(
  tissue = unique(GSVA_canonical_path$tissue),
  KEGG_eq_p = perform_tost_analysis("KEGG_CALCIUM_SIGNALING_PATHWAY", GSVA_canonical, GSVA_canonical_path, meta2854, "Sprs", 0.05)
)
write.table(ad_kegg_results, paste0(output_dir, "AD_KEGG_delta_0.05_2854sprs.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


## Analysis 3 & 4: AD Brain-specific Pathways
print("Running TOST for AD Brain pathways...")
GSVA_brain1 <- read.delim(gsva_ad_brain1_file)
GSVA_brain1_path <- data.frame(str_split_fixed(rownames(GSVA_brain1), ":", n = 2), stringsAsFactors = FALSE)
colnames(GSVA_brain1_path) <- c("tissue", "pathway")

ad_brain_results1 <- data.frame(
  tissue = unique(GSVA_brain1_path$tissue),
  microglia_eq_p = perform_tost_analysis("Microglia_Cell_death_&_apoptosis", GSVA_brain1, GSVA_brain1_path, meta2722, "Sprs", 0.05)
)
write.table(ad_brain_results1, paste0(output_dir, "AD_brain_delta_0.05_2722sprs.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

GSVA_brain2 <- read.delim(gsva_ad_brain2_file)
GSVA_brain2_path <- data.frame(str_split_fixed(rownames(GSVA_brain2), ":", n = 2), stringsAsFactors = FALSE)
colnames(GSVA_brain2_path) <- c("tissue", "pathway")

ad_brain_results2 <- data.frame(
  tissue = unique(GSVA_brain2_path$tissue),
  MYELIN_eq_p = perform_tost_analysis("GOBP_MYELIN_MAINTENANCE", GSVA_brain2, GSVA_brain2_path, meta2854, "Sprs", 0.05)
)
write.table(ad_brain_results2, paste0(output_dir, "AD_brain_delta_0.05_2854sprs.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


## Analysis 5: SCZ Gene Ontology Pathway
print("Running TOST for SCZ GO pathway...")
GSVA_GO_SWE <- read.delim(gsva_scz_go_file)
GSVA_GO_path_SWE <- data.frame(str_split_fixed(rownames(GSVA_GO_SWE), ":", n = 2), stringsAsFactors = FALSE)
colnames(GSVA_GO_path_SWE) <- c("tissue", "pathway")

scz_go_results <- data.frame(
  tissue = unique(GSVA_GO_path_SWE$tissue),
  BP_SWE_eq_p = perform_tost_analysis("BP_muscle_tissue_development", GSVA_GO_SWE, GSVA_GO_path_SWE, metaSWE, "prs", 0.033)
)
write.table(scz_go_results, paste0(output_dir, "SCZ_delta_0.033_SWEprs.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

print("All analyses complete.")