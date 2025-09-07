library(bigsnpr)
library(bigstatsr)
library(bigparallelr)
library(data.table)
library(tidyverse)


info <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/25503788",
  fname = "map_hm3_ldpred2.rds"))

sumstats <- bigreadr::fread2("/data2/zhaolab-shared/data/GWAS_SCZ_Trubetskoy/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv")
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
adj.sumstats <- bigreadr::fread2("/data2/xli37/data/overlap_sample/entire.EraSOR.adjusted.SCZ.assoc") 
adj.sumstats <- adj.sumstats %>% rename(rsid = SNP) %>% left_join(sumstats) 
adj.sumstats_hapmap <- adj.sumstats[adj.sumstats$rsid %in% info$rsid,] %>% left_join(info[,c(5,6)]) %>% 
  mutate(beta =  Z / sqrt(2*af_UKBB*(1 - af_UKBB)*(n_eff +Z^2)),
         beta_se = 1 / sqrt(2*af_UKBB*(1 - af_UKBB)*(n_eff +Z^2))) %>% select(-A1, -A2)

write.table(adj.sumstats_hapmap, file = "/data2/xli37/data/overlap_sample/entire.EraSOR.adjusted.SCZ.assoc.wposbeta.txt", sep = "\t",col.names = TRUE, row.names = FALSE, quote = FALSE)


snp_readBed("/data2/xli37/data/SCZ_SWE/SWE_eur_afterQC.bed")
obj.bigSNP = snp_attach("/data2/xli37/data/SCZ_SWE/SWE_eur_afterQC.rds")

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
rsID <- obj.bigSNP$map$marker.ID

NCORES <- 35

Gipt <- snp_fastImputeSimple(G, method = "mean2") 

map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
POS2 <- snp_asGeneticPos(CHR, POS, dir = "/data2/xli37/data/AD_PRS/position", ncores = NCORES)

df_beta <- snp_match(adj.sumstats_hapmap, map, join_by_pos = FALSE)

save(df_beta, G, Gipt, POS2, file = "/data2/xli37/data/AD_PRS/review/SWEforPadj.Rdata")
