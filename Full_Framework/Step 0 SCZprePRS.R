# This script prepares GWAS summary statistics and target genotype datasets
# for a Polygenic Risk Score (PRS) calculation using the bigsnpr package.
# The main steps include:
# 1. Formatting and filtering the GWAS summary statistics.
# 2. Processing two separate target genotype cohorts (e.g., discovery and replication).
#    - Reading genotype data from .bed files.
#    - Performing mean imputation for missing values.
#    - Matching the genotype data with the summary statistics.
#    - Saving the prepared data objects for downstream PRS analysis.

# --- Load Libraries ---
library(bigsnpr)
library(bigstatsr)
library(bigparallelr)
library(data.table)
library(tidyverse)

# --- Summary Statistics Preparation ---

# Download and load a reference map of HapMap3 SNPs, which is often used for LD calculations.
print("Loading reference SNP map...")
info <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/25503788",
  fname = "map_hm3_ldpred2.rds"
))

# Load the primary GWAS summary statistics (e.g., from PGC3 SCZ GWAS).
print("Loading and formatting GWAS summary statistics...")
sumstats <- bigreadr::fread2("/path/to/gwas/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv")
sumstats <- sumstats %>%
  mutate(
    chr = CHROM,
    pos = POS,
    a1 = A1,
    a0 = A2,
    n_eff = NEFF,
    beta = BETA,
    beta_se = SE,
    rsid = ID,
    pvalue = PVAL
  ) %>%
  select(chr, pos, rsid, a1, a0, n_eff, beta, beta_se, pvalue)

# Load an adjusted set of summary statistics (e.g., from EraSOR) and merge with the primary set.
adj_sumstats <- bigreadr::fread2("/path/to/adjusted/entire.EraSOR.adjusted.SCZ.assoc")
adj_sumstats <- adj_sumstats %>%
  rename(rsid = SNP) %>%
  left_join(sumstats, by = "rsid")

# Filter the summary statistics:
# 1. Keep only SNPs present in the HapMap3 reference panel.
# 2. Join with the reference panel to get allele frequencies (e.g., af_UKBB).
# 3. Recalculate beta and beta standard error from Z-scores, which ensures consistency.
print("Filtering and recalculating beta values for summary statistics...")
adj_sumstats_hapmap <- adj_sumstats[adj_sumstats$rsid %in% info$rsid, ] %>%
  left_join(info[, c("rsid", "af_UKBB")]) %>%
  mutate(
    beta = Z / sqrt(2 * af_UKBB * (1 - af_UKBB) * (n_eff + Z^2)),
    beta_se = 1 / sqrt(2 * af_UKBB * (1 - af_UKBB) * (n_eff + Z^2))
  ) %>%
  select(-A1, -A2)

# Save the fully processed summary statistics to a file.
print("Saving processed summary statistics...")
write.table(
  adj_sumstats_hapmap,
  file = "/path/to/output/processed.EraSOR.adjusted.SCZ.assoc.txt",
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)


# --- Function to Process a Target Genotype Cohort ---

#' Prepares a target genotype dataset for PRS calculation.
#'
#' @param bed_file Path to the target cohort's .bed file.
#' @param summary_stats The processed GWAS summary statistics data frame.
#' @param output_rdata_file Path to save the final RData object.
#' @param genetic_map_dir Directory containing genetic map files for `snp_asGeneticPos`.
#' @param n_cores Number of cores to use for parallel processing.
#' @return Does not return a value; saves the prepared data to the specified file.
process_target_cohort <- function(bed_file, summary_stats, output_rdata_file, genetic_map_dir, n_cores) {
  
  print(paste("Processing target cohort from:", bed_file))
  
  # Attach the bigSNP object from the .bed, .bim, and .fam files.
  rds_file <- sub("\\.bed$", ".rds", bed_file)
  snp_readBed(bed_file)
  obj_bigSNP <- snp_attach(rds_file)
  
  # Extract genotype data and map information.
  G <- obj_bigSNP$genotypes
  CHR <- obj_bigSNP$map$chromosome
  POS <- obj_bigSNP$map$physical.pos
  
  # Perform simple and fast imputation for missing genotypes (replaces NA with mean).
  print("Imputing missing genotypes...")
  Gipt <- snp_fastImputeSimple(G, method = "mean2")
  
  # Get genetic positions from physical positions.
  print("Calculating genetic positions...")
  POS2 <- snp_asGeneticPos(CHR, POS, dir = genetic_map_dir, ncores = n_cores)
  
  # Match the summary statistics with the genotype data's map.
  # This ensures that alleles are aligned between the two datasets.
  print("Matching summary statistics with genotype map...")
  map <- setNames(obj_bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
  df_beta <- snp_match(summary_stats, map, join_by_pos = FALSE)
  
  # Save all necessary prepared objects into a single RData file.
  print(paste("Saving prepared data to:", output_rdata_file))
  save(df_beta, G, Gipt, POS2, file = output_rdata_file)
}


# --- Process Discovery and Replication Cohorts ---

# Define the number of cores for parallel computation.
NCORES <- 35

# Process the discovery cohort.
process_target_cohort(
  bed_file = "/path/to/discovery_cohort/genotypes.bed",
  summary_stats = adj_sumstats_hapmap,
  output_rdata_file = "/path/to/output/discovery_cohort_prepared.Rdata",
  genetic_map_dir = "/path/to/genetic_map_directory/",
  n_cores = NCORES
)

# Process the replication cohort.
process_target_cohort(
  bed_file = "/path/to/replication_cohort/genotypes.bed",
  summary_stats = adj_sumstats_hapmap,
  output_rdata_file = "/path/to/output/replication_cohort_prepared.Rdata",
  genetic_map_dir = "/path/to/genetic_map_directory/",
  n_cores = NCORES
)

print("Data preparation complete.")