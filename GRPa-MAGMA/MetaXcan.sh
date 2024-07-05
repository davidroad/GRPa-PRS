conda activate imlabtools

tissue=(Amygdala Anterior_cingulate_cortex_BA24 Caudate_basal_ganglia Cerebellar_Hemisphere Cerebellum Cortex Frontal_Cortex_BA9 Hippocampus Hypothalamus Nucleus_accumbens_basal_ganglia Putamen_basal_ganglia Spinal_cord_cervical_c-1 Substantia_nigra)

# use PrediXcan to impute gene expression data in 13 brain tissues based on the genotype (.vcf)

for i in $(seq 0 12)
do
python3 ./software/Predict.py \
--model_db_path ./eqtl/mashr/mashr_Brain_${tissue[i]}.db \
--model_db_snp_key varID \
--vcf_genotypes /data2/xli37/data/SCZ_SWE/SWE_eur_afterQC.vcf \
--vcf_mode genotyped \
--liftover ../summary-gwas-imputation/data/liftover/hg19ToHg38.over.chain.gz \
--on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
--prediction_output /data2/xli37/data/mashr/SCZ/SWE6628_MASHR_Brain_${tissue[i]}_predict.txt \
--prediction_summary_output /data2/xli37/data/mashr/SCZ/SWE6628_MASHR_Brain_${tissue[i]}_summary.txt \
--verbosity 9 \
--throw
done

# use MultiXcan to generate gene p-value 
subgroup=(case_control TB10all T10 B10 TB10SCZ TB10Ctr TB15all TB15SCZ TB15Ctr T15 B15 TB20all TB20SCZ TB20Ctr T20 B20)

for i in $(seq 0 15)
do
python3 ./software/MulTiXcan.py \
--expression_folder /data2/xli37/data/mashr/SCZ/SWE/ \
--expression_pattern "SWE6628_MASHR_(.*)_predict.txt" \
--input_phenos_file /data2/xli37/data/TWAS/subgrp_SWEmeta0217.txt \
--input_phenos_column ${subgroup[i]} \
--covariates_file /data2/xli37/data/TWAS/subgrp_SWEmeta0217.txt \
--covariates sex X1 X2 X3 X4 X5 \
--mode logistic \
--verbosity 1 \
--output /data2/xli37/data/TWAS/MultiXcan/SWE/pc5/MASHR_mt_${subgroup[i]}.txt \
--throw
done
