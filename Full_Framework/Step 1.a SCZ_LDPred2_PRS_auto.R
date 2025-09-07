# This script calculates Polygenic Risk Scores (PRS) using the LDpred2-auto method
# from the bigsnpr package. The process is repeated for several different
# p-value thresholds to filter the input summary statistics.

# --- 1. Setup and Data Loading ---

# Load the required library
library(bigsnpr)

# Define key parameters and file paths
# Input file containing the prepared data from the previous step
input_rdata <- "/path/to/prepared_data/cohort_prepared.Rdata"
# Directory for storing large temporary files created during correlation calculation
tmp_dir <- "/path/to/tmp-data/"
# Base path for saving the final PRS results
output_prs_prefix <- "/path/to/output/mycohort_prs_auto"
# Number of cores to use for parallel computations
NCORES <- 40

# Load the prepared data objects:
# - df_beta: A data frame of GWAS summary statistics matched to the target cohort.
# - G: The original genotype matrix (a bigSNP object).
# - Gipt: The imputed genotype matrix.
# - POS2: A vector of genetic positions for the SNPs.
load(input_rdata)

# Store the original, unfiltered summary statistics
df_beta_original <- df_beta


# --- 2. Iterative PRS Calculation for Different P-value Thresholds ---

# Define the set of p-value thresholds to iterate through.
# A value of 1.0 means no p-value filtering is applied initially. Normallyl p_thresholds was set for 1.0 for LDPred2.
#p_thresholds <- c(1.0, 0.5, 0.3, 0.2, 0.1, 0.05)
p_thresholds <- 1.0

# Loop over each p-value threshold to calculate a corresponding PRS
for (p_val_threshold in p_thresholds) {

  cat("\n--- Starting PRS calculation for p-value threshold:", p_val_threshold, "---\n")

  # --- 2a. Filter Summary Statistics ---
  # Create a subset of the summary statistics based on the current p-value threshold
  df_beta <- df_beta_original[df_beta_original$pvalue <= p_val_threshold, ]
  cat("Number of SNPs after filtering:", nrow(df_beta), "\n")

  # --- 2b. Calculate LD (Linkage Disequilibrium) Matrices ---
  # This step calculates the correlation between SNPs, which is crucial for LDpred2.
  # NOTE: This is a computationally intensive step that is repeated for each p-value subset.

  # Create a temporary file for the on-disk correlation matrix
  tmp <- tempfile(tmpdir = tmp_dir)
  
  # Calculate LD matrix from imputed genotypes (used for LDpred2)
  print("Calculating LD matrix from imputed genotypes...")
  for (chr in 1:22) {
    cat(chr, ".. ", sep = "")
    ind_chr_gwas <- which(df_beta$chr == chr)
    ind_chr_geno <- df_beta$`_NUM_ID_`[ind_chr_gwas]
    
    corr0 <- snp_cor(Gipt, ind.col = ind_chr_geno, size = 3 / 1000,
                     infos.pos = POS2[ind_chr_geno], ncores = NCORES)
    
    if (chr == 1) {
      ld <- Matrix::colSums(corr0^2)
      corr <- as_SFBM(corr0, tmp, compact = TRUE)
    } else {
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
    }
  }

  # Calculate LD matrix from original genotypes (used for LD Score Regression)
  print("\nCalculating LD matrix from original genotypes for LDSC...")
  tmp_o <- tempfile(tmpdir = tmp_dir)
  for (chr in 1:22) {
    cat(chr, ".. ", sep = "")
    ind_chr_gwas <- which(df_beta$chr == chr)
    ind_chr_geno <- df_beta$`_NUM_ID_`[ind_chr_gwas]

    corr0_o <- snp_cor(G, ind.col = ind_chr_geno, size = 3 / 1000,
                       infos.pos = POS2[ind_chr_geno], ncores = NCORES)

    if (chr == 1) {
      ld_o <- Matrix::colSums(corr0_o^2)
      corr_o <- as_SFBM(corr0_o, tmp_o, compact = TRUE)
    } else {
      ld_o <- c(ld_o, Matrix::colSums(corr0_o^2))
      corr_o$add_columns(corr0_o, nrow(corr_o))
    }
  }

  # --- 2c. Estimate Heritability using LD Score Regression ---
  print("\nRunning LD Score regression to estimate h2...")
  ldsc <- snp_ldsc(ld_score = ld_o,
                   ld_size = length(ld_o),
                   chi2 = (df_beta$beta / df_beta$beta_se)^2,
                   sample_size = df_beta$n_eff)
  h2_est <- ldsc[["h2"]]
  cat("Heritability estimate (h2):", h2_est, "\n")

  # --- 2d. Run LDpred2-auto ---
  # This is the main algorithm that automatically infers SNP effect sizes (betas)
  # by modeling LD and learning from the data.
  print("Running LDpred2-auto...")
  set.seed(116) # for reproducibility
  
  # A shrinkage coefficient to adjust for potential mismatch between the LD reference
  # and the GWAS summary statistics LD structure.
  coef_shrink <- 0.95
  
  multi_auto <- snp_ldpred2_auto(
    corr, df_beta, h2_init = h2_est,
    vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
    ncores = NCORES,
    allow_jump_sign = FALSE,
    shrink_corr = coef_shrink
  )

  # --- 2e. Post-process Betas and Calculate PRS ---
  # The auto method returns multiple chains; we select the robust ones and average them.
  print("Finalizing SNP effects and calculating PRS...")
  range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
  keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))
  
  # Average the beta estimates from the selected chains
  beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
  
  # Calculate the final PRS for each individual by multiplying genotypes by the estimated betas
  pred_auto <- big_prodVec(Gipt, beta_auto, ind.col = df_beta[["_NUM_ID_"]])

  # --- 2f. Save the Results ---
  output_file <- sprintf("%s_p%s.txt", output_prs_prefix, p_val_threshold)
  print(paste("Saving PRS scores to:", output_file))
  
  write.table(pred_auto, file = output_file, sep = "\t",
              col.names = FALSE, row.names = FALSE, quote = FALSE)
}

print("\nAll PRS calculations are complete.")