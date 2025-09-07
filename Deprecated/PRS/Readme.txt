## PRS Readme

### Step 1: prepare for your PRS calculation
The code includes preparing your summary statistics and target genotype (.bed) in a way that can be directly used in the PRS calculation.

### Step 2: PRS calculation by ldpred2_auto
The code includes loading what you prepared and completing the PRS calculation. In this example, the P-value threshold was set to 0.5.

### Step 3: subgroup based on the PRS and diagnosis
The code based on a data frame (newmetaSWE) includes at least the diagnosis label and PRS you have calculated and other covariates you are interested in. You can use this code to generate a data frame that includes the subgroup label mentioned in our manuscript, which can be used in the following analysis.
