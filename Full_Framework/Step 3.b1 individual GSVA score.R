# This script performs Gene Set Variation Analysis (GSVA) on predicted gene expression data.
# It requires several command-line arguments to specify input files, the enrichment analysis type, and paths.

# Load required libraries
library(GSVA)
library(GSEABase)

# --- Argument Parsing ---
# The script expects the following command-line arguments:
# 1. An index (e.g., "1") corresponding to the input file to be processed.
# 2. An enrichment option: "GO", "canonical", "brain_cell_type_synapse_immune", or "brain_cell_type_function_AD_pathways".
# 3. The base path for output files (e.g., "/path/to/output/directory/").
# 4. The path to the directory containing the input TWAS (Transcriptome-Wide Association Study) data (e.g., "/path/to/twas/results/").
# 5. The working directory for the script (e.g., "/path/to/script/location/").
# 6. A status or dataset name for creating a unique output subfolder (e.g., "dataset_name").
args <- commandArgs(trailingOnly = TRUE)

input_index <- as.numeric(args[1])
enrich_option <- args[2]
output_path <- args[3]
twas_path <- args[4]
work_path <- args[5]
status <- args[6]

# Set the working directory and create the output directory if it doesn't exist
setwd(work_path)
system(paste0("mkdir -p ", output_path, status))

# --- Gene Set Preparation ---
# Load gene sets from GMT files.

# 1. Canonical Pathways (KEGG, REACTOME, BIOCARTA)
gmt_pathway_file <- "/path/to/c2.all.v7.4.symbols.gmt"
all_pathway_sets <- getGmt(gmt_pathway_file,
                           collectionType = BroadCollection(category = "c2"),
                           geneIdType = SymbolIdentifier())
# Filter for canonical pathway gene sets
canonical_geneset_pathway <- all_pathway_sets[c(grep("^KEGG", names(all_pathway_sets)),
                                                grep("^REACTOME", names(all_pathway_sets)),
                                                grep("^BIOCARTA", names(all_pathway_sets)))]

# 2. Brain Cell Type and Function Gene Sets
gmt_cell_function_file <- "/path/to/brain_cell_type_function_synapse_immune.gmt"
geneset_cell_type_function <- getGmt(gmt_cell_function_file,
                                     collectionType = BroadCollection(category = "c5"),
                                     geneIdType = SymbolIdentifier())

gmt_ad_pathways_file <- "/path/to/brain_cell_type_function_AD_pathways.gmt"
geneset_brain_cell_type_function_AD_pathways <- getGmt(gmt_ad_pathways_file,
                                                       collectionType = BroadCollection(category = "c5"),
                                                       geneIdType = SymbolIdentifier())

# 3. Gene Ontology (GO) Gene Sets (Non-redundant)
gmt_go_file <- "/path/to/GO_merged_non_redundant_full_term_name_gene_symbol_rename.gmt"
geneset_go <- getGmt(gmt_go_file,
                     collectionType = BroadCollection(category = "c5"),
                     geneIdType = SymbolIdentifier())

# --- Reference Gene Expression Data Preparation ---
# Load and format a reference gene expression matrix (e.g., GTEx median TPM)
# to map gene symbols.
gtex_tpm_file <- "/path/to/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
h19 <- read.delim(gtex_tpm_file, as.is = TRUE, head = FALSE)
h19 <- h19[-c(1, 2), ] # Remove header lines
colnames(h19) <- h19[1, ]
h19 <- h19[-1, ]

# --- Main Processing Loop ---
# List all files in the TWAS input directory and process the one specified by the input_index.
files <- dir(twas_path)
for (i in files[input_index]) {
    
    file_name <- gsub("_predict.txt", "", i)
    print(paste("Processing file:", file_name))

    # Read and format the predicted expression matrix
    mashr_result <- read.delim(paste0(twas_path, i), as.is = TRUE, check.names = FALSE)
    mashr_result <- mashr_result[, -1]
    rownames(mashr_result) <- mashr_result[, 1]
    mashr_result <- mashr_result[, -1]

    # Transpose the matrix for GSVA (genes as rows, samples as columns)
    new_file_filter_t <- t(mashr_result)
    new_file_filter_t <- as.matrix(new_file_filter_t)

    # Match ENSEMBL IDs (from rows) to Gene Symbols using the GTEx reference
    idx_19 <- match(rownames(new_file_filter_t), h19[, 1])
    rownames(new_file_filter_t) <- h19[idx_19, 2]

    # --- GSVA Execution ---
    # Based on the 'enrich_option', run GSVA with the corresponding gene set.
    setwd(output_path)

    if (enrich_option == "GO") {
        gsvaOut_go <- gsva(new_file_filter_t, geneset_go,
                           min.sz = 5,       # Minimum size of the gene sets.
                           max.sz = 500,     # Maximum size of the gene sets.
                           verbose = TRUE)   # Provides information about each calculation step.
        gsvaOut_go <- rbind(id = colnames(gsvaOut_go), gsvaOut_go)
        write.table(gsvaOut_go, file = paste0(status, "/gsvaOut_go_mixdiff_T_", file_name, ".txt"), sep = "\t", quote = FALSE, col.names = FALSE)
    }

    if (enrich_option == "canonical") {
        gsvaOut_canonical <- gsva(new_file_filter_t, canonical_geneset_pathway,
                                  min.sz = 5,
                                  max.sz = 500,
                                  verbose = TRUE)
        gsvaOut_canonical <- rbind(id = colnames(gsvaOut_canonical), gsvaOut_canonical)
        write.table(gsvaOut_canonical, file = paste0(status, "/gsvaOut_canonical_mixdiff_T_", file_name, ".txt"), sep = "\t", quote = FALSE, col.names = FALSE)
    }

    if (enrich_option == "brain_cell_type_synapse_immune") {
        gsvaOut_brain_cell_type_function_synapse_immune <- gsva(new_file_filter_t, geneset_cell_type_function,
                                                                min.sz = 5,
                                                                max.sz = 500,
                                                                verbose = TRUE)
        gsvaOut_brain_cell_type_function_synapse_immune <- rbind(id = colnames(gsvaOut_brain_cell_type_function_synapse_immune), gsvaOut_brain_cell_type_function_synapse_immune)
        write.table(gsvaOut_brain_cell_type_function_synapse_immune, file = paste0(status, "/gsvaOut_brain_cell_type_function_synapse_immune_mixdiff_T_", file_name, ".txt"), sep = "\t", quote = FALSE, col.names = FALSE)
    }

    if (enrich_option == "brain_cell_type_function_AD_pathways") {
        gsvaOut_brain_cell_type_function_AD_pathways <- gsva(new_file_filter_t, geneset_brain_cell_type_function_AD_pathways,
                                                             min.sz = 5,
                                                             max.sz = 500,
                                                             verbose = TRUE)
        gsvaOut_brain_cell_type_function_AD_pathways <- rbind(id = colnames(gsvaOut_brain_cell_type_function_AD_pathways), gsvaOut_brain_cell_type_function_AD_pathways)
        write.table(gsvaOut_brain_cell_type_function_AD_pathways, file = paste0(status, "/gsvaOut_brain_cell_type_function_AD_pathways_mixdiff_T_", file_name, ".txt"), sep = "\t", quote = FALSE, col.names = FALSE)
    }
}