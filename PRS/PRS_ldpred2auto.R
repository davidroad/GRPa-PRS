library(bigsnpr)

# require: df_beta, G, Gipt, POS2
load("/data2/xli37/data/AD_PRS/review/SWEforPadj.Rdata")
NCORES = 40

# include only p-value < 0.5
df_beta = df_beta[df_beta$pvalue <=0.5,]

df_beta = df_beta

tmp <- tempfile(tmpdir = "/data2/xli37/data/tmp-data")

for (chr in 1:22) {
  
  # Compute the correlation for each chromosome
  # print(chr)
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'G'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  
  corr0 <- snp_cor(Gipt, ind.col = ind.chr2, size = 3 / 1000, # window size of 3 cM
                   infos.pos = POS2[ind.chr2], ncores = NCORES)
  
  # create the on-disk sparse genome-wide correlation matrix
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

tmp <- tempfile(tmpdir = "/data2/xli37/data/tmp-data")

for (chr in 1:22) {
  
  # Compute the correlation for each chromosome
  # print(chr)
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'G'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  
  corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000, # window size of 3 cM
                   infos.pos = POS2[ind.chr2], ncores = NCORES)
  
  # create the on-disk sparse genome-wide correlation matrix
  if (chr == 1) {
    ld_O <- Matrix::colSums(corr0^2)
    corr_O <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld_O <- c(ld_O, Matrix::colSums(corr0^2))
    corr_O$add_columns(corr0, nrow(corr_O))
  }
}

#LD score regression
ldsc <- snp_ldsc(   ld_O,
                    length(ld_O),
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff)
h2_est <- ldsc[["h2"]]
h2_se <- ldsc[["h2_se"]]

# auto

coef_shrink <- 0.95  # reduce this up to 0.4 if you have some (large) mismatch with the LD ref

set.seed(116)  # to get the same result every time
# takes less than 2 min with 4 cores
multi_auto <- snp_ldpred2_auto(
  corr, df_beta, h2_init = h2_est,
  vec_p_init = seq_log(1e-4, 0.2, length.out = 30), ncores = NCORES,
  # use_MLE = FALSE,  # uncomment if you have convergence issues or when power is low (need v1.11.9)
  allow_jump_sign = FALSE, shrink_corr = coef_shrink)


(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
pred_auto <- big_prodVec(Gipt, beta_auto,  ind.col = df_beta[["_NUM_ID_"]])

write.table(pred_auto, file = "/data2/xli37/data/AD_PRS/review/SWE_prs_auto_p0.5.txt", sep = "\t",col.names = FALSE, row.names = FALSE, quote = FALSE)
