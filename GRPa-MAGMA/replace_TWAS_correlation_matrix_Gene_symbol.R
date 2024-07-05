##Replace TWAS genes z score to correlation matrix of MAGMA
library(readr)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
input.genes.raw = args[1] ## "Wightman_new_loc.genes.raw"
input.genes.out = args[2] ## "Wightman_new_loc.genes.out"
gene.loc = args[3] ## "Gencode26_GRCh37_52648.gene.loc"
TWAS_file = args[4] ## "./code/MASHR_mt_case_control.txt"
work_path = args[5] ##./
output_name = args[6] ## "Wightman"

setwd(work_path)
#dat_raw <- read_lines("Wightman_new_loc.genes.raw")
dat_raw <-  read_lines(input.genes.raw)
dat_top <- dat_raw[1:2]
dat = dat_raw[-c(1:2)]
annotation4magma <- read.delim(gene.loc, header = F)[,c(1,6)]
colnames(annotation4magma) <- c("GENE", "ID")
outADc <- read.table(input.genes.out, header = T) %>% left_join(annotation4magma)

  ADc_TWAS <- read.delim(TWAS_file) %>% rename(TWAS.P = pvalue,GENE =gene) %>% select(TWAS.P, GENE)
  ADc_TWAS_filtered <- ADc_TWAS[!is.na(ADc_TWAS$TWAS.P),]
  ADc_update <- outADc %>% left_join(ADc_TWAS_filtered)
  ADc_update$TWAS.P[ADc_update$TWAS.P == 1] <- 0.9999999
  ADc_update <- ADc_update%>%
    mutate(TWAS.Z = qnorm(1 - TWAS.P) #  qf(1 - TWAS.P, df1 = 13, df2 = 2722 - 13 )
	) %>% rename(P0 = START, P1 = STOP)
  ADc_update$TWAS.Z[which(ADc_update$TWAS.Z == Inf)] <- 10 # max(ADc_update$TWAS.Z[ADc_update$TWAS.Z != Inf & !is.na(ADc_update$TWAS.Z)])
	

ADc_update_gene = ADc_update$GENE[!is.na(ADc_update$TWAS.P)]
print(paste0("Overlap TWAS and correlation gene number:", length(ADc_update_gene),"!") )
## regenerate the gene.raw file using the existing TWAS gene list and existing correlation
dat_tmp <- dat
genes <- vector("list", length(dat))
correlations <- vector("list", length(dat))
correlation_genes <- vector("list", length(dat))
row_to_remove <- {}
replaced_row <- 0
for (i in 1:length(dat)) {
  row_data <- strsplit(dat[i], " ")[[1]]
  genes[[i]] <- row_data[1]
  tmp <- {}
  tmp_gene <- {}
  tmp_correlations <- {}
  if (genes[[i]] %in% ADc_update_gene) {
	if (length(row_data) > 9) {
		correlations[[i]] <- as.numeric(row_data[10:length(row_data)])
		for (p in (i-length(correlations[[i]])): (i-1)) {
			tmp <- c(tmp,genes[[p]])
		}
		correlation_genes[[i]] <- tmp
		if (sum(correlation_genes[[i]] %in% ADc_update_gene)>0) {
			replaced_row = replaced_row +1
			tmp_gene = correlation_genes[[i]][correlation_genes[[i]] %in% ADc_update_gene]
			tmp_correlations = correlations[[i]][correlation_genes[[i]] %in% ADc_update_gene]
			dat_tmp[i] = paste(c(ADc_update$ID[match(genes[[i]], ADc_update$GENE)],row_data[2:8], ADc_update$TWAS.Z[match(genes[[i]], ADc_update$GENE)],tmp_correlations),collapse=" ")
		}else{
			replaced_row = replaced_row +1
			dat_tmp[i] = paste(c(ADc_update$ID[match(genes[[i]], ADc_update$GENE)],row_data[2:8], ADc_update$TWAS.Z[match(genes[[i]], ADc_update$GENE)]),collapse=" ")
		}
	} else {
		replaced_row = replaced_row + 1
		dat_tmp[i] = paste(c(ADc_update$ID[match(genes[[i]], ADc_update$GENE)],row_data[2:8], ADc_update$TWAS.Z[match(genes[[i]], ADc_update$GENE)]),collapse=" ")
	}
  }else {
    row_to_remove = c(row_to_remove,i)
  }
}
print(paste0("Replaced row number:", replaced_row,"!"))
length(row_to_remove)

##replaced_row + row_to_remove should equals to length(dat)
dat_out = dat_tmp[-row_to_remove]

writeLines(paste(dat_top,collapse="\n"), paste0(output_name,".genes.raw"))
write.table(dat_out, file = paste0(output_name,".genes.raw"),append = TRUE, row.names=F,col.names=F,quote=F)
