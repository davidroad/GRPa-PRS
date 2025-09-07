# This script performs a combined association analysis of GSVA scores across multiple tissues
# against various phenotypic subgroups. It uses a likelihood ratio test to determine the
# significance of pathway associations and then generates heatmaps to visualize the results.

# --- Load Libraries ---
library(tidyverse)
library(pheatmap)

# --- Argument Parsing ---
# The script expects the following command-line arguments:
# 1. PC_num: The number of principal components to use as covariates (e.g., "3", "5", "10").
# 2. GP: A group parameter, likely related to subgroup definitions (e.g., "10", "15", "20").
# 3. apoe: A string indicating APOE adjustment status (e.g., "_no_adjAPOE", "_adjAPOE").
# 4. cov_file: Path to the covariate and phenotype file (e.g., "/path/to/covariate_file.txt").
# 5. status: A status or dataset name for file paths (e.g., "2722inAPOE").
# 6. APOE_status: A more specific dataset identifier for input files (e.g., "disc2722_MASHR").
# 7. GMT_name: The name of the gene set collection used (e.g., "go", "canonical", "brain").
# 8. output_path: The base directory for saving combined results (e.g., "/path/to/combined_results/").
# 9. input_path: The directory where per-tissue GSVA and DEG results are stored (e.g., "/path/to/gsva_results/").
args <- commandArgs(trailingOnly = TRUE)

PC_num <- args[1]
GP <- args[2]
apoe <- args[3]
cov_file <- args[4]
status <- args[5]
APOE_status <- args[6]
GMT_name <- args[7]
output_path <- args[8]
input_path <- args[9]

# --- Directory and Parameter Setup ---
setwd(output_path)
dir.create(file.path(status, GMT_name), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(status, GMT_name))

# Define tissues of interest
tissues <- c("Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia",
             "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9",
             "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia",
             "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra")

# Define traits (outcomes) for association testing
traits <- c("case_control", paste0("TB", GP, "all"), paste0("TB", GP, "SCZ"), paste0("TB", GP, "Ctr"), paste0("T", GP), paste0("B", GP))

# --- Function Definition ---

#' Performs a logistic regression likelihood ratio test.
#'
#' @param meta A dataframe with covariates and the outcome variable.
#' @param X A dataframe of predictor variables (e.g., GSVA scores across tissues).
#' @param PC_num The number of PCs to include in the model.
#' @param outcome The name of the outcome column in the meta dataframe.
#' @return A numeric vector containing the degrees of freedom and the p-value.
logisticLRT <- function(meta, X, PC_num, outcome) {
    y <- meta %>% dplyr::select(all_of(outcome)) %>% pull()

    if (apoe == "_no_adjAPOE") {
        null_data <- meta %>% dplyr::select(starts_with("X") & ends_with(as.character(1:PC_num)), sex)
    } else {
        null_data <- meta %>% dplyr::select(starts_with("X") & ends_with(as.character(1:PC_num)), sex, apoe2, apoe4)
    }
    
    reg_data <- cbind(null_data, X)
    
    # Fit null model (covariates only) and full model (covariates + predictors)
    null_model <- glm(y ~ ., data = null_data, family = binomial(link = "logit"))
    reg_model <- glm(y ~ ., data = reg_data, family = binomial(link = "logit"))
    
    # Perform likelihood ratio test
    chisq <- null_model$deviance - reg_model$deviance
    df <- null_model$df.residual - reg_model$df.residual
    pvalue <- pchisq(chisq, df, lower.tail = FALSE)
    
    return(c(df, pvalue))
}


# --- Main Association Analysis ---

# Load covariate data
cov <- read.delim(cov_file, as.is = TRUE, check.names = FALSE)

# Load GSVA scores for all tissues
exp <- list()
pathway <- list()
for (i in 1:length(tissues)) {
    file_path <- file.path(input_path, paste0("gsvaOut_", GMT_name, "_mixdiff_T_", APOE_status, "_", tissues[i]))
    exp[[i]] <- read.delim(file_path, check.names = FALSE)
    pathway[[i]] <- exp[[i]][, 1]
}

# Get a list of all unique pathways across all tissues
pathway_uni <- unique(unlist(pathway))
print(paste("Found", length(pathway_uni), "unique pathways."))

# Loop through each trait to perform association analysis
for (outcome in traits) {
    print(paste("Running analysis for outcome:", outcome))
    
    # Initialize a dataframe to store results
    result <- data.frame()
    
    # Test each pathway for association with the outcome
    for (j in 1:length(pathway_uni)) {
        pw <- pathway_uni[j]
        
        # Consolidate scores for the current pathway from all tissues
        x <- data.frame()
        for (i in 1:length(tissues)) {
            x <- rbind(x, exp[[i]] %>% filter(id == pw))
        }
        
        X <- t(x[, -1])
        colnames(X) <- paste0("tissue", 1:ncol(X))
        
        # Run the likelihood ratio test and append results
        lrt_result <- logisticLRT(cov, X, PC_num = PC_num, outcome = outcome)
        result <- rbind(result, c(pw, lrt_result))
    }
    
    colnames(result) <- c("pathway", "n_tissue", "raw_p")
    
    # Format results and calculate FDR adjusted p-value
    result <- result %>%
        mutate(
            n_tissue = as.numeric(n_tissue),
            raw_p = as.numeric(raw_p)
        ) %>%
        mutate(
            adj_p = p.adjust(raw_p, method = "fdr")
        ) %>%
        arrange(adj_p)
        
    # --- Add Directionality Information from DEG Results ---
    
    # Load corresponding DEG analysis results
    GSVA_deg <- list()
    for (i in 1:length(tissues)) {
        deg_path <- file.path(input_path, status, paste0(GMT_name, "_DEG"), paste0("PC", PC_num), outcome, paste0(APOE_status, "_", tissues[i], ".txt"))
        GSVA_deg[[i]] <- read.delim(deg_path, check.names = FALSE, as.is = TRUE)
    }

    # Extract effect direction and logFC for each pathway
    direction_out <- c()
    logFC <- c()
    direction_out_combined <- c()
    logFC_combined <- c()
    
    for (q in 1:length(pathway_uni)) {
        tmp_direction <- c()
        tmp_FC <- c()
        for (p in 1:length(tissues)) {
            match_row <- match(pathway_uni[q], rownames(GSVA_deg[[p]]))
            tmp_direction <- c(tmp_direction, GSVA_deg[[p]][match_row, 2]) # Assuming column 2 is effect direction
            tmp_FC <- c(tmp_FC, GSVA_deg[[p]][match_row, 1])      # Assuming column 1 is logFC
        }
        
        # Get direction and logFC from the tissue with the max absolute effect
        direction_out <- c(direction_out, tmp_direction[which.max(abs(tmp_direction))])
        logFC <- c(logFC, tmp_FC[which.max(abs(tmp_FC))])
        
        # Calculate combined direction and logFC by summing non-NA values
        direction_out_combined <- c(direction_out_combined, sum(tmp_direction[!is.na(tmp_direction)]))
        logFC_combined <- c(logFC_combined, sum(tmp_FC[!is.na(tmp_FC)]))
    }

    # Combine association results with directionality
    result <- data.frame(result, direction_out, logFC, direction_out_combined, logFC_combined, check.names = FALSE)
    
    # Write final results to file
    output_filename <- paste0(status, "_", GMT_name, "_PC", PC_num, "_", outcome, "_", GP, "_combined.txt")
    write.table(result, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}


# --- Heatmap Generation ---
# This section generates heatmaps summarizing the association results across all tested traits.
# It creates plots based on different significance levels (FDR or raw p-value).

print("Generating summary heatmaps...")

for (use_fdr in c(TRUE, FALSE)) {

    # Identify significant pathways across all traits
    Term_list <- c()
    for (outcome in traits) {
        file_name <- paste0(status, "_", GMT_name, "_PC", PC_num, "_", outcome, "_", GP, "_combined.txt")
        file <- read.delim(file_name, check.names = FALSE, as.is = TRUE)
        
        if (use_fdr) {
            threshold <- 0.2
            significant_terms <- file[file$adj_p < threshold, "pathway"]
        } else {
            threshold <- 0.001
            significant_terms <- file[file$raw_p < threshold, "pathway"]
        }
        Term_list <- c(Term_list, significant_terms)
    }
    Term_list <- unique(Term_list)

    if (length(Term_list) == 0) {
        print(paste("No significant terms found for FDR =", use_fdr, ". Skipping heatmap generation."))
        next
    }

    # --- Create Heatmap for Filtered (Significant) Pathways Only ---
    
    # Create a summary matrix to hold p-values
    summary_matrix <- matrix(1, ncol = length(traits), nrow = length(Term_list))
    colnames(summary_matrix) <- traits
    rownames(summary_matrix) <- Term_list
    
    # Populate the matrix with signed p-values
    for (outcome in traits) {
        file_name <- paste0(status, "_", GMT_name, "_PC", PC_num, "_", outcome, "_", GP, "_combined.txt")
        file <- read.delim(file_name, check.names = FALSE, as.is = TRUE)
        
        matched_rows <- match(Term_list, file$pathway)
        
        # Use raw p-value and combined direction for the heatmap value
        p_values <- file[matched_rows, "raw_p"]
        signs <- sign(file[matched_rows, "logFC_combined"])
        summary_matrix[, outcome] <- p_values * signs
    }

    # Transform p-values to -log10 scale for visualization
    summary_matrix_log10 <- -log10(abs(summary_matrix)) * sign(summary_matrix)
    summary_matrix_log10[is.na(summary_matrix_log10)] <- 0 # Handle NAs if any pathway was missing
    
    # Filter out columns that have no significant results
    summary_matrix_log10_filtered <- summary_matrix_log10[, colSums(abs(summary_matrix_log10)) != 0, drop = FALSE]
    summary_matrix_log10_filtered <- abs(summary_matrix_log10_filtered) # Use absolute values for heatmap color

    if (ncol(summary_matrix_log10_filtered) > 0) {
        # Define file names based on whether FDR or raw p-value was used
        fdr_tag <- if (use_fdr) paste0("FDR_", threshold) else paste0("rawP_", threshold)
        
        # Save the data table
        write.table(summary_matrix_log10_filtered, file = paste0("summary_heatmap_filtered_", fdr_tag, ".txt"), quote = FALSE, sep = "\t")

        # Save the heatmap as a PDF
        pdf(paste0("summary_heatmap_filtered_", fdr_tag, ".pdf"), 10, 15)
        tryCatch({
            p <- pheatmap(summary_matrix_log10_filtered,
                          cluster_rows = FALSE,
                          angle_col = 315,
                          cellwidth = 18,
                          cellheight = 18,
                          display_numbers = matrix(ifelse(summary_matrix_log10_filtered > -log10(threshold), "*", ""), nrow(summary_matrix_log10_filtered)))
            print(p)
        }, error = function(e) { print(paste("Error plotting heatmap:", e)) })
        dev.off()
        
        # Save a one-color version of the heatmap
        pdf(paste0("summary_heatmap_filtered_one_color_", fdr_tag, ".pdf"), 10, 15)
        tryCatch({
            p <- pheatmap(summary_matrix_log10_filtered,
                          color = colorRampPalette(c("white", "firebrick1"))(50),
                          cluster_rows = FALSE,
                          angle_col = 315,
                          cellwidth = 18,
                          cellheight = 18,
                          display_numbers = matrix(ifelse(summary_matrix_log10_filtered > -log10(threshold), "*", ""), nrow(summary_matrix_log10_filtered)))
            print(p)
        }, error = function(e) { print(paste("Error plotting heatmap:", e)) })
        dev.off()
    }
}