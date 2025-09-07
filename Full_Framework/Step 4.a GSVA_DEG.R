# This script performs differential analysis on Gene Set Variation Analysis (GSVA) scores
# to identify significantly altered pathways between cases and controls. It uses the limma package
# for differential analysis and generates volcano plots for visualization.

# --- Load Libraries ---
library(limma)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

# --- Argument Parsing ---
# The script expects the following command-line arguments:
# 1. The number of principal components to use as covariates (e.g., "PC3", "PC5", "PC10").
# 2. The working directory containing the GSVA results (e.g., "/path/to/gsva/results/").
# 3. The path to the covariate and phenotype file (e.g., "/path/to/covariate_file.txt").
# 4. The enrichment option name, used to filter input files (e.g., "brain_cell_type_function_AD_pathways").
# 5. A status/dataset identifier, also for filtering input files (e.g., "my_dataset_ID").
# 6. The name of the main output folder to be created within the working directory (e.g., "output_folder_name").
args <- commandArgs(trailingOnly = TRUE)

pc_num <- args[1]
work_directory <- args[2]
cov_file_path <- args[3]
enrich_option <- args[4]
dataset_status <- args[5]
output_folder <- args[6]

# --- Directory Setup ---
# Set the primary working directory and create output subdirectories.
setwd(work_directory)
system(paste0("mkdir -p ", file.path(work_directory, output_folder)))
system(paste0("mkdir -p ", file.path(work_directory, output_folder, paste0(enrich_option, "_DEG"))))

# --- File and Data Preparation ---
# Find GSVA result files matching the specified enrichment option and dataset status.
all_files <- dir()
files_to_process <- all_files[grep(enrich_option, all_files)]
files_to_process <- files_to_process[grep(dataset_status, files_to_process)]
print("Files to be processed:")
print(files_to_process)

# Define significance cutoffs
logFC_cutoff <- 0.1
adjPvalue_cutoff <- 0.05

# Define the groups/subgroups for comparison from the covariate file columns.
groups <- c("case_control", "TB10all", "TB10SCZ", "TB10Ctr", "TB15all", "TB15SCZ", "TB15Ctr", "TB20all", "TB20SCZ", "TB20Ctr", "T10", "B10", "T15", "B15", "T20", "B20")

# Load the covariate file once before the loop.
cov <- read.delim(cov_file_path, as.is = TRUE, check.names = FALSE)

# --- Main Processing Loop ---
# Loop through each GSVA result file.
for (i in files_to_process) {
    
    # Load the GSVA scores matrix.
    gsva_scores <- read.delim(i, check.names = FALSE, row.names = 1)
    
    # Clean up the filename for use in output files.
    i_name <- gsub(paste0("gsvaOut_", enrich_option, "_mixdiff_T_"), "", i)

    # Loop through each predefined group for differential analysis.
    for (j in 1:length(groups)) {
        subgroup <- groups[j]
        print(paste("Processing file:", i_name, "for subgroup:", subgroup))

        # --- Subgroup-specific Directory Setup ---
        setwd(file.path(work_directory, output_folder, paste0(enrich_option, "_DEG")))
        system(paste0("mkdir -p ", pc_num))
        system(paste0("mkdir -p ", file.path(pc_num, subgroup)))

        # Match and filter samples present in both the GSVA scores and the covariate file for the current subgroup.
        # Include only samples where the group status is 0 (control) or 1 (case).
        match_idx <- which(cov[, colnames(cov) == subgroup] %in% c(0, 1))
        
        rt_filter <- gsva_scores[, match_idx]
        cov_filter <- cov[match_idx, ]

        # Prepare the design matrix for limma, including covariates.
        type1 <- factor(cov_filter[, colnames(cov_filter) == subgroup])
        sex <- factor(cov_filter$sex)

        design <- NULL
        if (pc_num == "PC10") {
            design <- model.matrix(~ 0 + type1 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex, data = cov_filter)
        } else if (pc_num == "PC5") {
            design <- model.matrix(~ 0 + type1 + X1 + X2 + X3 + X4 + X5 + sex, data = cov_filter)
        } else if (pc_num == "PC3") {
            design <- model.matrix(~ 0 + type1 + X1 + X2 + X3 + sex, data = cov_filter)
        } else {
            stop("Invalid number of PCs specified. Please use 'PC3', 'PC5', or 'PC10'.")
        }

        colnames(design)[1:2] <- c("Control", "Case")
        row.names(design) <- colnames(rt_filter)
        
        # --- Differential Analysis with limma ---
        contrast.matrix <- makeContrasts("Case-Control", levels = design)
        fit <- lmFit(rt_filter, design)
        fit <- contrasts.fit(fit, contrast.matrix)
        fit <- eBayes(fit)
        
        # Get the full results table.
        all_results <- topTable(fit, coef = 1, n = Inf, sort.by = "none")

        # Add standard error and calculate Cohen's d for effect size.
        all_results$beta_se <- fit$stdev.unscaled * fit$sigma
        all_results$effect_size <- all_results$logFC / fit$sigma
        
        # Sort results by absolute effect size.
        all_results <- all_results[order(abs(all_results$effect_size), decreasing = TRUE), ]
        
        # Write the full results table to a file.
        write.table(all_results, file = file.path(pc_num, subgroup, paste0(i_name, ".txt")), sep = "\t", quote = FALSE, col.names = TRUE)

        # Filter for significantly differentiated pathways.
        df <- all_results[all_results$adj.P.Val < adjPvalue_cutoff & abs(all_results$logFC) > logFC_cutoff, ]
        write.table(df, file = file.path(pc_num, subgroup, paste0("GSVA_DEG_", enrich_option, "_", i_name, ".txt")), sep = "\t", quote = FALSE, row.names = TRUE)

        # --- Volcano Plot Generation ---
        volcano_data <- all_results
        volcano_data$pathway <- rownames(all_results)
        volcano_data$FDR <- -log10(volcano_data$adj.P.Val)

        volcano_sig_data <- df
        volcano_sig_data$pathway <- rownames(df)
        volcano_sig_data$FDR <- -log10(volcano_sig_data$adj.P.Val)

        mycol <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))

        p <- ggplot() +
            geom_point(data = volcano_data, aes(x = logFC, y = FDR), color = mycol[2], size = 1) +
            geom_point(data = volcano_sig_data, aes(x = logFC, y = FDR), color = mycol[1], size = 1) +
            theme_bw() +
            labs(x = "Log (Fold Change)", y = "-Log10 (FDR)") +
            geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), col = mycol[3], linetype = "dotted", size = 1) +
            geom_hline(yintercept = -log10(adjPvalue_cutoff), col = mycol[3], linetype = "dotted", size = 1) +
            theme(
                axis.title = element_text(size = 15, color = "black", face = "bold"),
                axis.text = element_text(size = 10, color = "black", face = "bold")
            )

        # Create a second plot with labels for significant points.
        p1 <- p + geom_text(data = volcano_sig_data, aes(x = logFC, y = FDR, label = pathway), hjust = -0.1, vjust = 0, size = 3)

        pdf(file.path(pc_num, subgroup, paste0(enrich_option, "_", i_name, "_volcano.pdf")), width = 7, height = 7)
        print(p)
        print(p1)
        dev.off()
    }
}