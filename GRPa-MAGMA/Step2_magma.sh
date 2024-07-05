# annotate
./magma --annotate --snp-loc ../SCZ_SWE/SWE_eur_afterQC.bim --gene-loc geneloc/Gencode26_GRCh37_52648.gene.loc --out SWE_newloc
./magma --bfile   ../SCZ_SWE/SWE_eur_afterQC --gene-annot SWE_newloc.genes.annot --out SWE_newloc

subgroup=(case_control TB10all TB10SCZ TB10Ctr T10 B10 TB15all TB15SCZ TB15Ctr T15 B15 TB20all TB20SCZ TB20Ctr T20 B20)

# MAGMA
for i in $(seq 0 15)
do
## Regenerate corresponding correlation matrix for each subgroup with their multixcan p-value
Rscript replace_TWAS_correlation_matrix_Gene_symbol.R \
SWE_newloc.genes.raw  \
SWE_newloc.genes.out  \
./geneloc/Gencode26_GRCh37_52648.gene.loc  \
/data2/xli37/data/TWAS/MultiXcan/SWE/pc5/MASHR_mt_${subgroup[i]}.txt  \
/data2/xli37/data/MAGMA/  \
/data2/xli37/data/MAGMA/GSEA/SWE/SWE_new_loc.${subgroup[i]}

## run GSEA for 3 different gene sets
./magma \
--model  self-contained \
--gene-results /data2/xli37/data/MAGMA/GSEA/SWE/SWE_new_loc.${subgroup[i]}.genes.raw \
--set-annot /data2/ydai2/data/TP/brain_cell_type_function_synapse_immune.gmt \
--out  /data2/xli37/data/MAGMA/GSEA/SWE/${subgroup[i]}.brain

./magma \
--model  self-contained \
--gene-results /data2/xli37/data/MAGMA/GSEA/SWE/SWE_new_loc.${subgroup[i]}.genes.raw \
--set-annot /home/xli37/AD_PRS/GO_merged_non_redundant_full_term_name_gene_symbol_rename.gmt \
--out  /data2/xli37/data/MAGMA/GSEA/SWE/${subgroup[i]}.GO

./magma \
--model  self-contained \
--gene-results /data2/xli37/data/MAGMA/GSEA/SWE/SWE_new_loc.${subgroup[i]}.genes.raw \
--set-annot /home/xli37/AD_PRS/c2.canonical.v7.4.symbols.gmt \
--out  /data2/xli37/data/MAGMA/GSEA/SWE/${subgroup[i]}.canonical
done
