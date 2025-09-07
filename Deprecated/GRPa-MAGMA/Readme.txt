## GRPa-MAGMA Readme
### Step 1: run the related MetaXcan code to prepare for further analysis
Firstly, PrediXcan imputed gene expression data in 13 brain tissues based on the genotype (.vcf). Then, MultiXcan was first used to integrate TWAS across brain regions and identify the associations between genes and traits or subgroup indicators we generated based on the PRS and diagnosis. 

### Step 2: run MAGMA for the GRPa to do GSEA in different gene sets
After the annotation, we regenerate MAGMA correlation matrix for each subgroup with their multixcan p-value. Then, run the magma to conduct GSEA.

### Step 3: summarize the GSEA result 
GSEA results in different conditions were stored to plot the heatmap for visualization.
