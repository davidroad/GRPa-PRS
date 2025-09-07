# GRPa-PRS
## GRPa-PRS: A risk stratification method to identify genetically-regulated pathways in polygenic diseases
https://pmc.ncbi.nlm.nih.gov/articles/PMC10327215/ 

<img width="1098" height="716" alt="image" src="https://github.com/user-attachments/assets/384b211f-35c5-468e-a48d-d32d03adbbf2" />



**GRPa-PRS workflow and study design**

![image](https://github.com/davidroad/GRPa-PRS/assets/4857356/09ef6ac1-d922-41d6-a879-52ac36c3b91c)

**Illustration of six strata comparisons and genetic factor distribution in GRPas**

Polygenic risk scores (PRS) are tools used to evaluate an individual’s susceptibility to polygenic diseases based on their genetic profile. A considerable proportion of people carry a high genetic risk but evade the disease. On the other hand, some individuals with a low risk of eventually developing the disease. We hypothesized that unknown counterfactors might be involved in reversing the PRS prediction, which might provide new insights into the pathogenesis, prevention, and early intervention of diseases.


<img width="920" height="774" alt="image" src="https://github.com/user-attachments/assets/ca72b473-3b44-4a88-b05c-1cb1cd7e91cc" />


# GRPa-PRS Workflow using SCZ Dataset as an Example

### **Step 0: Data Preparation for PRS**
* Load and format GWAS summary statistics (e.g., an EraSOR adjusted version) and process target cohort genotype data from `.bed` files to prepare materials for PRS calculation.
* Prepare `.gmt` files to define the gene sets for pathway analyses.

### **Step 1.a: PRS Calculation**
* Calculate Polygenic Risk Scores (PRS) using the `LDpred2` method. This is often run for multiple GWAS p-value thresholds to find the best-performing model. The input is the output from Step 0.

### **Step 1.b: PRS Evaluation and Subgroup Definition**
* Evaluate the different PRS models from Step 1a using AUC analysis to determine predictive performance.
* Select the best-performing PRS and use its quantiles to define high-risk and low-risk subgroups.
* Compile a final metadata (`META info`) file that includes covariates and subgroup labels for use in all subsequent analyses.

### **Step 2: Predict Individual-Level Gene Expression (GReX)**
* Use PrediXcan to predict individual-level gene expression from genotype data for multiple tissues (e.g., 13 brain tissues). The output is a matrix of GReX per individual per tissue.

### **Step 3.a1: Multi-Tissue TWAS with MultiXcan**
* Perform a multi-tissue gene-based association test (TWAS) using the predicted GReX (from Step 2) and the binary subgroups (from Step 1.b). The output is a table of gene-level p-values for each trait association.

### **Step 3.a2: MAGMA Gene Set Enrichment Analysis**
* Use MAGMA to test for enrichment of biological pathways within the TWAS results.
* Input is the gene-level p-value list from the MultiXcan analysis (Step 3.a1) and a gene set (`.gmt`) file.
* The analysis can be run as a competitive test (outputting a `.gsa.out` file) or a self-contained test (outputting a `.gsa.self.out` file).

### **Step 3.b1: GSVA Score Calculation**
* Calculate pathway-level scores using Gene Set Variation Analysis (GSVA).
* The inputs are the individual GReX matrices (from Step 2) and a gene set (`.gmt`) file.
* The output is a GSVA score for each gene set, individual, and tissue.

### **Step 3.b2: GSVA Likelihood Ratio Test**
* Use a likelihood ratio test to assess the combined association of pathway GSVA scores from multiple tissues with the clinical subgroups defined in Step 1.b.
* The results are often visualized with heatmaps to show significant association patterns.

### **Step 4.a: Differential Pathway Analysis**
* Identify which pathways are differentially active between the subgroups from Step 1.b. This test is performed on the GSVA scores (from Step 3.b1), not on individual genes.

### **Step 4.b: Orthogonality Test**
* Conduct a TOST equivalence test to confirm that the PRS (from Step 1.b) and significant pathway GSVA scores (from Step 3.b1) are independent (orthogonal).
