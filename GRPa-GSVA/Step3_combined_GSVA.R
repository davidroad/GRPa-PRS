library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
PC_num = args[1] ##3, 5, 10
GP = args[2] ##10,15,20
apoe = args[3] ## _no_adjAPOE, adjAPOE
cov_file = args[4] ## "/collab2/ydai2/24_GRPa/code/subgrp_ADcmeta0217.txt"
status = args[5] ## "2722inAPOE" or "2722exAPOE" or "2854inAPOE" or "2854exAPOE"
APOE_status = args[6] ## "disc2722_MASHR"
GMT_name = args[7] ## brain, GO, canonical pathway
output_path = args[8] ## "/collab2/ydai2/24_GRPa/GSVA/combined_result/AD/"
input_path = args[9] ## "/collab2/ydai2/24_GRPa/GSVA/result/AD"
##go,canonical,brain,cell_type_function
#numPC ="3"
#apoe = "APOE"
#all_tissue = "Brain_only"
#model = "brain"

setwd(output_path)
system(paste0("mkdir -p ",status))
system(paste0("mkdir -p ",status,"/",GMT_name,"/"))

setwd(paste0(output_path,status,"/",GMT_name,"/"))
##example_code##
##Rscript combined_GSVA_1854.R 3 _no_APOE Brain_only go
apoe_out = ""
 if(apoe=="_no_adjAPOE"){
  apoe_out = "_no_adjAPOE"
  }else{
  apoe_out = "_adjAPOE"
  }

logisticLRT = function(meta, X , PC_num, outcome){
  y = meta %>% dplyr::select(outcome) %>% as_vector()
  if(apoe=="_no_adjAPOE"){
  null_data = meta %>% dplyr::select(X1:paste0("X", PC_num), sex)
  }					#,apoeGenotype
  else{
  null_data = meta %>% dplyr::select(X1:paste0("X", PC_num), sex, apoe2, apoe4)
  }
  reg_data = cbind(null_data, X)
  
  null_model = glm(data = null_data, y ~., family =  binomial(link = logit))
  reg_model = glm(data = reg_data, y ~., family =  binomial(link = logit))
  chisq = null_model$deviance - reg_model$deviance
  df = null_model$df.residual - reg_model$df.residual
  pvalue = pchisq(chisq, df , lower.tail=FALSE)
  return(c(df,pvalue))
}

# tissues you want to include
element <- c("Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hypothalamus","Brain_Hippocampus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra")

cov = read.delim( cov_file,as.is=T,check.names=F)
traits <- c("case_control",paste0("TB",GP,"all"),paste0("TB",GP,"SCZ"),paste0("TB",GP,"Ctr"),paste0("T",GP),paste0("B",GP))

##treat PRS as the covariates
exp = list()
pathway = list()
for (i in 1:length(element)) {
  exp[[i]] = read.delim(paste0(input_path,"/gsvaOut_",GMT_name,"_mixdiff_T_",APOE_status,"_",element[i]),check.names=F)
  pathway[[i]] = exp[[i]][,1]
}

pathway_uni = unique(unlist(pathway))
print(length(pathway_uni))
for (outcome in traits){


result = data.frame()

for (j in 1:length(pathway_uni)) {
  pw = pathway_uni[j]
  x = data.frame()
  for (i in 1:length(element)) {
    x = rbind(x,exp[[i]] %>% filter(id == pw))
  }
  X = t(x[,-1])
  colnames(X) = c(paste0("tissue",1:ncol(X)))
  result = rbind (result, c(pw, logisticLRT(cov, X, PC_num = PC_num, outcome = outcome))) # specify the outcome you want to explore
}

colnames(result) = c("pathway", "n_tissue", "raw_p")

result = result %>% mutate(
  n_tissue = as.numeric(n_tissue),
  raw_p = as.numeric(raw_p)
) %>% 
  mutate(
    adj_p =  p.adjust(raw_p, method = "fdr")
  ) %>% arrange(adj_p)

##add direction
GSVA = list()
for (i in 1:length(element)) {
  GSVA[[i]] = read.delim(paste0(input_path,"/",status,"/",GMT_name,"_DEG", "/PC",PC_num,"/",outcome,"/",APOE_status,"_",element[i],".txt"),check.names=F,as.is=T)
  #pathway[[i]] = exp[[i]][,1]
}
logFC <- {}
direction_out <- {}
direction_out_combined <- {}
logFC_combined <- {}
for (q in 1:length(pathway_uni)){
	tmp <- {}
	tmp_FC <- {}
	for (p in 1:length(element)){
		tmp <- c(tmp,GSVA[[p]][match(pathway_uni[q],rownames(GSVA[[p]])),2])
		tmp_FC <- c(tmp_FC,GSVA[[p]][match(pathway_uni[q],rownames(GSVA[[p]])),1])
	}
direction_out <- c(direction_out,tmp[which.max(abs(tmp))])
direction_out_combined <- c(direction_out_combined,sum(tmp[!is.na(tmp)]))
logFC <- c(logFC,tmp_FC[which.max(abs(tmp_FC))])
logFC_combined <- c(logFC_combined,sum(tmp_FC[!is.na(tmp_FC)]))
}

result <- data.frame(result,direction_out,logFC,direction_out_combined,logFC_combined,check.names=F)

write.table(result, file= paste0(status,"_",GMT_name,"_PC",PC_num,"_",outcome,"_",GP,"_combined.txt"),sep="\t",quote=F,row.names=F)

}

##add figures##
###1/31
##example for cell_type_function##
library(pheatmap)
Term = GMT_name
##cell_type_function
##"SelfT10","SelfB10" removed
##canonical,go,brain_cell_type_function_AD_pathways together
for (FDR in c(TRUE, FALSE)) {

##without WB
###without LogFC filtration
#FDR = FALSE
###use the absolute value###
##canonical,go,brain_cell_type_function_AD_pathways together
for (GP in c(GP)) {
	for (PC in c(paste0("PC",PC_num))) {
#		for (Term in c("go","canonical","brain_cell_type_function_AD_pathways")){
		for (Term in c(GMT_name)){
print(c(GP,PC,Term))
#GP = 10
library(pheatmap)
##cell_type_function
#,"SelfT10","SelfB10"

#FDR = TRUE 
#traits <- c("case_control","STB10all","STB10AD","STB10Ctr","KTB10all","KTB10AD","KTB10Ctr","WTB10all","WTB10AD","WTB10Ctr","SelfTB10all","SelfTB10AD","SelfTB10Ctr","ST10","SB10","KT10","KB10","WT10","WB10")
#traits <- c("case_control","STB10all","STB10AD","STB10Ctr","KTB10all","KTB10AD","KTB10Ctr","WTB10all","WTB10AD","WTB10Ctr","SelfTB10all","SelfTB10AD","SelfTB10Ctr","ST10","SB10","KT10","KB10","WT10","WB10")
traits <- c("case_control",paste0("TB",GP,"all"),paste0("TB",GP,"SCZ"),paste0("TB",GP,"Ctr"),paste0("T",GP),paste0("B",GP))
#model = c("_Brain_whole_blood","_Brain_only")
##get qualified term 
#threshold = 0.05/(167*16) # 0.001,0.0005
#FC_threshold = 0.1478
FC_threshold = 0.0
#FC_threshold = 0.196
exp_threshold = 0.0
#FC_combine_threshold = 0.416
#FC_combine_threshold = 0.838
#FC_combine_threshold = 0.15
#FC_combine_threshold = 0.05

#PC = "PC5" # by default, option,"PC3","PC5","PC10"
#Term="brain_cell_type_function_AD_pathways" # "go","brain","canonical","cell_type_function"
#Term = "canonical"
Term_list <-{}
all_file <- {}
    for (outcome in traits) {
##/data2/ydai2/projects/22_AD_TWAS_PRS/GSVA/2492_inAPOE_result/brain_cell_type_function_AD_pathways_no_adjAPOE_DEG/PC5/diag
      file_name=paste0(status,"_",GMT_name,"_PC",PC_num,"_",outcome,"_",GP,"_combined.txt")
      file = read.delim(file_name, check.names=F,as.is=T)
#	  FC_combine_threshold = quantile(x=abs(file$logFC_combined),probs=c(0.5))
	  FC_combine_threshold = 0
	  all_file = rbind(all_file,file)
#	  threshold = 0.05/(nrow(file))
	  if (FDR) {
		threshold = 0.2
#		Term_list <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		Term_list <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold)  &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
	  }else{
		threshold = 0.001
#		Term_list <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		Term_list <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold)  &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])

	  }
#      Term_list <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs(file[,5])> exp_threshold) & (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold) & ((sign(file[,6])* sign(file[,7]))>0) ),1])
#     Term_list <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
#      Term_list <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs( file[,5] + sign(file[,5]) * abs(file[,6])/2 ) > exp_threshold) & (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold)  ),1])

    }


##first use raw pvalue to cut a threshold
quantile(x=abs(all_file$direction_out[all_file$raw_p<0.05/nrow(file)]),probs=c(0.05,0.25,0.5,0.95,0.975,1))
quantile(x=abs(all_file$logFC[all_file$raw_p<0.05/nrow(file)] ),probs=c(0.05,0.25,0.4,0.5,0.75,0.8,0.95,0.975,1))
quantile(x=abs(all_file$logFC_combined[all_file$raw_p<0.05/nrow(file)] ),probs=c(0.05,0.25,0.4,0.5,0.75,0.8,0.95,0.975,1))
Term_list <- unique(Term_list)
#construct one matrix to hold#
sum_table_col <- c(paste0(traits))
sum_table_cell_type_function <- matrix(1,ncol=length(sum_table_col ),nrow=length(Term_list))
colnames(sum_table_cell_type_function) <- sum_table_col
rownames(sum_table_cell_type_function) <- Term_list
##fill in the table
    for (outcome in traits) {
      file_name=paste0(status,"_",GMT_name,"_PC",PC_num,"_",outcome,"_",GP,"_combined.txt")
      file = read.delim(file_name, check.names=F,as.is=T)
#	  FC_combine_threshold = quantile(x=abs(file$logFC_combined),probs=c(0.5))
	  FC_combine_threshold = 0 
#      tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs(file[,5])> exp_threshold)& (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold) & ((sign(file[,6])* sign(file[,7]))>0)),1]
#      tmp <- file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold) & (abs(file[,5])> exp_threshold)& (sign(file[,6])== sign(file[,8])) & ((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold)) ),1]
#       tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs( file[,5] + sign(file[,5]) * abs(file[,6])/2 ) > exp_threshold) & (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold)  ),1]
	
	  if (FDR) {
		threshold = 0.2
#		tmp <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		tmp <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
	  }else{
		threshold = 0.001
#		tmp <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		tmp <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])

	  }
        if (length(tmp)>=1){
            idx_row = match(tmp,Term_list)
            idx_col = match(paste0(outcome),sum_table_col)
			if (FDR) {
				sum_table_cell_type_function[idx_row,idx_col] <- file[match(tmp,file[,1]),3] *sign(file[match(tmp,file[,1]),8])
			}else{
				sum_table_cell_type_function[idx_row,idx_col] <- file[match(tmp,file[,1]),3] *sign(file[match(tmp,file[,1]),8])
			}
        }
	}

colnames(sum_table_cell_type_function) <- sum_table_col
rownames(sum_table_cell_type_function) <- Term_list
sum_table_cell_type_function_log10 <- -log10(abs(sum_table_cell_type_function))*sign(sum_table_cell_type_function)
sum_table_cell_type_function_log10_filter <- rbind(sum_table_cell_type_function_log10[,!colSums(sum_table_cell_type_function_log10)==0])
sum_table_cell_type_function_log10_filter <- abs(sum_table_cell_type_function_log10_filter)
#sum_table_cell_type_function_log10_filter <- sum_table_cell_type_function_log10_filter[!rowSums(sum_table_cell_type_function_log10)==0,]
#print(pheatmap(sum_table_cell_type_function_log10_filter,cellwidth = 18, cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter))))
##tsea.plot(sum_table_cell_type_function_log10_filter)

colnames(sum_table_cell_type_function_log10_filter) <- sapply(colnames(sum_table_cell_type_function_log10_filter),function(x){gsub("_no_adjAPOE","",x)})
rownames(sum_table_cell_type_function_log10_filter) <- Term_list

if (FDR){
write.table(sum_table_cell_type_function_log10_filter,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_without_filtration_noWB_0.2.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_without_filtration_noWB_0.2.pdf"),10,15)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_filter,cluster_rows =F,cellwidth = 18, angle_col = 315,cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter)))), error = function(e)NULL)
dev.off()
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_without_filtration_noWB_one_color_0.2.pdf"),10,15)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_filter,color = colorRampPalette(c("white", "firebrick1"))(50),cluster_rows =F,cellwidth = 18, angle_col = 315,cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter)))), error = function(e)NULL)
dev.off()

}else{
write.table(sum_table_cell_type_function_log10_filter,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_without_filtration_noWB_0.001.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_without_filtration_noWB_0.001.pdf"),10,15)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_filter,cluster_rows =F,cellwidth = 18, angle_col = 315,cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter)))), error = function(e)NULL)
dev.off()
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_without_filtration_noWB_one_color_0.001.pdf"),10,15)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_filter,color = colorRampPalette(c("white", "firebrick1"))(50),cluster_rows =F,cellwidth = 18, angle_col = 315,cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter)))), error = function(e)NULL)
dev.off()
}

##fill in the full box
Term_list <- unique(Term_list)
#construct one matrix to hold#
sum_table_col <- c(paste0(traits))
sum_table_cell_type_function_full <- matrix(1,ncol=length(sum_table_col ),nrow=length(Term_list))
colnames(sum_table_cell_type_function_full) <- sum_table_col
rownames(sum_table_cell_type_function_full) <- Term_list
##fill in the table
    for (outcome in traits) {
      file_name=paste0(status,"_",GMT_name,"_PC",PC_num,"_",outcome,"_",GP,"_combined.txt")
      file = read.delim(file_name, check.names=F,as.is=T)

#      tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs(file[,5])> exp_threshold)& (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold) & ((sign(file[,6])* sign(file[,7]))>0)),1]
#       tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs( file[,5] + sign(file[,5]) * abs(file[,6])/2 ) > exp_threshold) & (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold)  ),1]
		tmp  <- file[match(Term_list,file[,1]),1]
        if (length(tmp)>=1){
			for (m in tmp){
				idx_row = m
				idx_col = match(paste0(outcome),sum_table_col)
				if (FDR) {
					sum_table_cell_type_function_full[idx_row,idx_col] <- file[match(m,file[,1]),3] *sign(file[match(m,file[,1]),8])
				}else{
					sum_table_cell_type_function_full[idx_row,idx_col] <- file[match(m,file[,1]),3] *sign(file[match(m,file[,1]),8])
				}
			}
		}
	}

colnames(sum_table_cell_type_function_full) <- sum_table_col
rownames(sum_table_cell_type_function_full) <- Term_list
sum_table_cell_type_function_log10_full <- -log10(abs(sum_table_cell_type_function_full))*sign(sum_table_cell_type_function_full)
#sum_table_cell_type_function_log10_full <- sum_table_cell_type_function_log10_full[,match(colnames(sum_table_cell_type_function_log10_filter),colnames(sum_table_cell_type_function_log10_full))]
#sum_table_cell_type_function_log10_filter <- sum_table_cell_type_function_log10_filter[!rowSums(sum_table_cell_type_function_log10)==0,]
sum_table_cell_type_function_log10_full <- rbind(sum_table_cell_type_function_log10_full[,c(1,2,3,4,5,6)])
rownames(sum_table_cell_type_function_log10_full) <- Term_list
#pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,gaps_col = cumsum(c(4,4,4,3,3)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(0.05/(nrow(file))), "*", ""), nrow(sum_table_cell_type_function_log10_full)))
if (length(nrow(sum_table_cell_type_function_log10_full))>0){
colnames(sum_table_cell_type_function_log10_full) <- sapply(colnames(sum_table_cell_type_function_log10_full),function(x){gsub("_no_adjAPOE","",x)})
#print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,gaps_col = cumsum(c(1,4,4,4,3,3)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14))

sum_table_cell_type_function_log10_full <- abs(sum_table_cell_type_function_log10_full)


if (FDR){
write.table(sum_table_cell_type_function_log10_full,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full_without_filtration_0.2.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full_without_filtration_0.2.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,2)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()
}else{
write.table(sum_table_cell_type_function_log10_full,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full_without_filtration_0.001.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full_without_filtration_0.001.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,2)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()

}

}else{
next
}
##filter WB
sum_table_cell_type_function_log10_full <- rbind(sum_table_cell_type_function_log10_full[,c(1,2,3,4,5,6)])
rownames(sum_table_cell_type_function_log10_full) <- Term_list


#pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,gaps_col = cumsum(c(4,4,4,3,3)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(0.05/(nrow(file))), "*", ""), nrow(sum_table_cell_type_function_log10_full)))
if (length(nrow(sum_table_cell_type_function_log10_full))>0){
colnames(sum_table_cell_type_function_log10_full) <- sapply(colnames(sum_table_cell_type_function_log10_full),function(x){gsub("_no_adjAPOE","",x)})
#print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,gaps_col = cumsum(c(1,4,4,4,3,3)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14))

sum_table_cell_type_function_log10_full <- abs(sum_table_cell_type_function_log10_full)

sum_table_cell_type_function_log10_full_out <- sum_table_cell_type_function_log10_full
sum_table_cell_type_function_log10_full_out = matrix(as.integer(sum_table_cell_type_function_log10_full_out>-log10(threshold)),ncol=ncol(sum_table_cell_type_function_log10_full))
colnames(sum_table_cell_type_function_log10_full_out) <- colnames(sum_table_cell_type_function_log10_full)
rownames(sum_table_cell_type_function_log10_full_out) <- rownames(sum_table_cell_type_function_log10_full)



if (FDR){
write.table(sum_table_cell_type_function_log10_full_out,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full_without_filtration_noWB_binary_0.2.txt"),quote=F,sep="\t")
write.table(sum_table_cell_type_function_log10_full,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full_without_filtration_noWB_0.2.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full_without_filtration_noWB_0.2.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,2)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full_without_filtration_noWB_one_color_0.2.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,color = colorRampPalette(c("white", "firebrick1"))(50),angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,3,3)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()
}else{
write.table(sum_table_cell_type_function_log10_full_out,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full_without_filtration_noWB_binary_0.001.txt"),quote=F,sep="\t")
write.table(sum_table_cell_type_function_log10_full,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full_without_filtration_noWB_0.001.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full_without_filtration_noWB_0.001.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,2)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full_without_filtration_noWB_one_color_0.001.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,color = colorRampPalette(c("white", "firebrick1"))(50),angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,3,3)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()



}
}else{
next
}
}
}

}



for (GP in c(GP)) {
	for (PC in c(paste0("PC",PC_num))) {
#		for (Term in c("go","canonical","brain_cell_type_function_AD_pathways")){
		for (Term in c(GMT_name)){
#GP = 10
###without LogFC filtration
#FDR = FALSE
print(c(GP,PC,Term))
#GP = 10
library(pheatmap)
traits <- c("case_control",paste0("TB",GP,"all"),paste0("TB",GP,"SCZ"),paste0("TB",GP,"Ctr"),paste0("T",GP),paste0("B",GP))

FC_threshold = 0.0
#FC_threshold = 0.196
exp_threshold = 0.0


Term_list <-{}
all_file <- {}
    for (outcome in traits) {
##/data2/ydai2/projects/22_AD_TWAS_PRS/GSVA/2722_inAPOE_result/brain_cell_type_function_AD_pathways_no_adjAPOE_DEG/PC5/diag
      file_name=paste0(status,"_",GMT_name,"_PC",PC_num,"_",outcome,"_",GP,"_combined.txt")
	 
      file = read.delim(file_name, check.names=F,as.is=T)
#	  FC_combine_threshold = quantile(x=abs(file$logFC_combined),probs=c(0.5))
	  FC_combine_threshold = 0
	  all_file = rbind(all_file,file)
	  if (FDR) {
		threshold = 0.05
#		Term_list <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		Term_list <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold)  &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
	  }else{
		threshold = 0.05/(nrow(file))
#		Term_list <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		Term_list <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold)  &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])

	  }
#      Term_list <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs(file[,5])> exp_threshold) & (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold) & ((sign(file[,6])* sign(file[,7]))>0) ),1])

    }


##first use raw pvalue to cut a threshold
quantile(x=abs(all_file$direction_out[all_file$raw_p<0.05/nrow(file)]),probs=c(0.05,0.25,0.5,0.95,0.975,1))
quantile(x=abs(all_file$logFC[all_file$raw_p<0.05/nrow(file)] ),probs=c(0.05,0.25,0.4,0.5,0.75,0.8,0.95,0.975,1))
quantile(x=abs(all_file$logFC_combined[all_file$raw_p<0.05/nrow(file)] ),probs=c(0.05,0.25,0.4,0.5,0.75,0.8,0.95,0.975,1))
Term_list <- unique(Term_list)
#construct one matrix to hold#
sum_table_col <- c(paste0(traits))
sum_table_cell_type_function <- matrix(1,ncol=length(sum_table_col ),nrow=length(Term_list))
colnames(sum_table_cell_type_function) <- sum_table_col
rownames(sum_table_cell_type_function) <- Term_list
##fill in the table
    for (outcome in traits) {
      file_name=paste0(status,"_",GMT_name,"_PC",PC_num,"_",outcome,"_",GP,"_combined.txt")
      file = read.delim(file_name, check.names=F,as.is=T)
#	  FC_combine_threshold = quantile(x=abs(file$logFC_combined),probs=c(0.5))
	  FC_combine_threshold = 0 
#      tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs(file[,5])> exp_threshold)& (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold) & ((sign(file[,6])* sign(file[,7]))>0)),1]
#      tmp <- file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold) & (abs(file[,5])> exp_threshold)& (sign(file[,6])== sign(file[,8])) & ((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold)) ),1]
#       tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs( file[,5] + sign(file[,5]) * abs(file[,6])/2 ) > exp_threshold) & (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold)  ),1]
	
	  if (FDR) {
		threshold = 0.05
#		tmp <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		tmp <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
	  }else{
		threshold = 0.05/(nrow(file))
#		tmp <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		tmp <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])

	  }
        if (length(tmp)>=1){
            idx_row = match(tmp,Term_list)
            idx_col = match(paste0(outcome),sum_table_col)
			if (FDR) {
				sum_table_cell_type_function[idx_row,idx_col] <- file[match(tmp,file[,1]),3] *sign(file[match(tmp,file[,1]),8])
			}else{
				sum_table_cell_type_function[idx_row,idx_col] <- file[match(tmp,file[,1]),3] *sign(file[match(tmp,file[,1]),8])
			}
        }
	}


colnames(sum_table_cell_type_function) <- sum_table_col
rownames(sum_table_cell_type_function) <- Term_list
sum_table_cell_type_function_log10 <- -log10(abs(sum_table_cell_type_function))*sign(sum_table_cell_type_function)
sum_table_cell_type_function_log10_filter <- rbind(sum_table_cell_type_function_log10[,!colSums(sum_table_cell_type_function_log10)==0])
sum_table_cell_type_function_log10_filter <- abs(sum_table_cell_type_function_log10_filter)
#sum_table_cell_type_function_log10_filter <- sum_table_cell_type_function_log10_filter[!rowSums(sum_table_cell_type_function_log10)==0,]
#print(pheatmap(sum_table_cell_type_function_log10_filter,cellwidth = 18, cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter))))
##tsea.plot(sum_table_cell_type_function_log10_filter)

colnames(sum_table_cell_type_function_log10_filter) <- sapply(colnames(sum_table_cell_type_function_log10_filter),function(x){gsub("_no_adjAPOE","",x)})
rownames(sum_table_cell_type_function_log10_filter) <- Term_list

if (FDR){
write.table(sum_table_cell_type_function_log10_filter,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_without_filtration.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_without_filtration.pdf"),10,15)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_filter,cluster_rows =F,cellwidth = 18, angle_col = 315,cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter)))), error = function(e)NULL)
dev.off()
}else{
write.table(sum_table_cell_type_function_log10_filter,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_without_filtration.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_without_filtration.pdf"),10,15)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_filter,cluster_rows =F,cellwidth = 18, angle_col = 315,cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter)))), error = function(e)NULL)
dev.off()
}

##fill in the full box
Term_list <- unique(Term_list)
#construct one matrix to hold#
sum_table_col <- c(paste0(traits))
sum_table_cell_type_function_full <- matrix(1,ncol=length(sum_table_col ),nrow=length(Term_list))
colnames(sum_table_cell_type_function_full) <- sum_table_col
rownames(sum_table_cell_type_function_full) <- Term_list
##fill in the table
    for (outcome in traits) {
      file_name=paste0(status,"_",GMT_name,"_PC",PC_num,"_",outcome,"_",GP,"_combined.txt")
					  
      file = read.delim(file_name, check.names=F,as.is=T)

#      tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs(file[,5])> exp_threshold)& (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold) & ((sign(file[,6])* sign(file[,7]))>0)),1]
#       tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs( file[,5] + sign(file[,5]) * abs(file[,6])/2 ) > exp_threshold) & (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold)  ),1]
		tmp  <- file[match(Term_list,file[,1]),1]
        if (length(tmp)>=1){
			for (m in tmp){
				idx_row = m
				idx_col = match(paste0(outcome),sum_table_col)
				if (FDR) {
					sum_table_cell_type_function_full[idx_row,idx_col] <- file[match(m,file[,1]),3] *sign(file[match(m,file[,1]),8])
				}else{
					sum_table_cell_type_function_full[idx_row,idx_col] <- file[match(m,file[,1]),3] *sign(file[match(m,file[,1]),8])
				}
			}
		}
	}


colnames(sum_table_cell_type_function_full) <- sum_table_col
rownames(sum_table_cell_type_function_full) <- Term_list
sum_table_cell_type_function_log10_full <- -log10(abs(sum_table_cell_type_function_full))*sign(sum_table_cell_type_function_full)
#sum_table_cell_type_function_log10_full <- sum_table_cell_type_function_log10_full[,match(colnames(sum_table_cell_type_function_log10_filter),colnames(sum_table_cell_type_function_log10_full))]
#sum_table_cell_type_function_log10_filter <- sum_table_cell_type_function_log10_filter[!rowSums(sum_table_cell_type_function_log10)==0,]
sum_table_cell_type_function_log10_full <- rbind(sum_table_cell_type_function_log10_full[,c(1,2,3,4,5,6)])
rownames(sum_table_cell_type_function_log10_full) <- Term_list
#pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,gaps_col = cumsum(c(4,4,4,3,3)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(0.05/(nrow(file))), "*", ""), nrow(sum_table_cell_type_function_log10_full)))
if (length(nrow(sum_table_cell_type_function_log10_full))>0){
colnames(sum_table_cell_type_function_log10_full) <- sapply(colnames(sum_table_cell_type_function_log10_full),function(x){gsub("_no_adjAPOE","",x)})
#print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,gaps_col = cumsum(c(1,4,4,4,3,3)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14))

sum_table_cell_type_function_log10_full <- abs(sum_table_cell_type_function_log10_full)




if (FDR){
write.table(sum_table_cell_type_function_log10_full,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full_without_filtration.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full_without_filtration.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,2)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()
}else{
write.table(sum_table_cell_type_function_log10_full,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full_without_filtration.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full_without_filtration.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,2)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()



}
}else{
next
}
}
}
}



##without WB
###without LogFC filtration
#FDR = FALSE
###use the absolute value###
##canonical,go,brain_cell_type_function_AD_pathways together
for (GP in c(GP)) {
	for (PC in c(paste0("PC",PC_num))) {
#		for (Term in c("go","canonical","brain_cell_type_function_AD_pathways")){
		for (Term in c(GMT_name)){
print(c(GP,PC,Term))
#GP = 10
library(pheatmap)
##cell_type_function
#,"SelfT10","SelfB10"

#FDR = TRUE 
#traits <- c("case_control","STB10all","STB10AD","STB10Ctr","KTB10all","KTB10AD","KTB10Ctr","WTB10all","WTB10AD","WTB10Ctr","SelfTB10all","SelfTB10AD","SelfTB10Ctr","ST10","SB10","KT10","KB10","WT10","WB10")
#traits <- c("case_control","STB10all","STB10AD","STB10Ctr","KTB10all","KTB10AD","KTB10Ctr","WTB10all","WTB10AD","WTB10Ctr","SelfTB10all","SelfTB10AD","SelfTB10Ctr","ST10","SB10","KT10","KB10","WT10","WB10")
traits <- c("case_control",paste0("TB",GP,"all"),paste0("TB",GP,"SCZ"),paste0("TB",GP,"Ctr"),paste0("T",GP),paste0("B",GP))
#model = c("_Brain_whole_blood","_Brain_only")
##get qualified term 
#threshold = 0.05/(167*16) # 0.001,0.0005
#FC_threshold = 0.1478
FC_threshold = 0.0
#FC_threshold = 0.196
exp_threshold = 0.0
#FC_combine_threshold = 0.416
#FC_combine_threshold = 0.838
#FC_combine_threshold = 0.15
#FC_combine_threshold = 0.05

#PC = "PC5" # by default, option,"PC3","PC5","PC10"
#Term="brain_cell_type_function_AD_pathways" # "go","brain","canonical","cell_type_function"
#Term = "canonical"
Term_list <-{}
all_file <- {}
    for (outcome in traits) {
##/data2/ydai2/projects/22_AD_TWAS_PRS/GSVA/2492_inAPOE_result/brain_cell_type_function_AD_pathways_no_adjAPOE_DEG/PC5/diag
      file_name=paste0(status,"_",GMT_name,"_PC",PC_num,"_",outcome,"_",GP,"_combined.txt")
      file = read.delim(file_name, check.names=F,as.is=T)
#	  FC_combine_threshold = quantile(x=abs(file$logFC_combined),probs=c(0.5))
	  FC_combine_threshold = 0
	  all_file = rbind(all_file,file)
#	  threshold = 0.05/(nrow(file))
	  if (FDR) {
		threshold = 0.05
#		Term_list <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		Term_list <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold)  &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
	  }else{
		threshold = 0.05/(nrow(file))
#		Term_list <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		Term_list <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold)  &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])

	  }
#      Term_list <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs(file[,5])> exp_threshold) & (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold) & ((sign(file[,6])* sign(file[,7]))>0) ),1])
#     Term_list <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
#      Term_list <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs( file[,5] + sign(file[,5]) * abs(file[,6])/2 ) > exp_threshold) & (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold)  ),1])

    }


##first use raw pvalue to cut a threshold
quantile(x=abs(all_file$direction_out[all_file$raw_p<0.05/nrow(file)]),probs=c(0.05,0.25,0.5,0.95,0.975,1))
quantile(x=abs(all_file$logFC[all_file$raw_p<0.05/nrow(file)] ),probs=c(0.05,0.25,0.4,0.5,0.75,0.8,0.95,0.975,1))
quantile(x=abs(all_file$logFC_combined[all_file$raw_p<0.05/nrow(file)] ),probs=c(0.05,0.25,0.4,0.5,0.75,0.8,0.95,0.975,1))
Term_list <- unique(Term_list)
#construct one matrix to hold#
sum_table_col <- c(paste0(traits))
sum_table_cell_type_function <- matrix(1,ncol=length(sum_table_col ),nrow=length(Term_list))
colnames(sum_table_cell_type_function) <- sum_table_col
rownames(sum_table_cell_type_function) <- Term_list
##fill in the table
    for (outcome in traits) {
      file_name=paste0(status,"_",GMT_name,"_PC",PC_num,"_",outcome,"_",GP,"_combined.txt")
      file = read.delim(file_name, check.names=F,as.is=T)
#	  FC_combine_threshold = quantile(x=abs(file$logFC_combined),probs=c(0.5))
	  FC_combine_threshold = 0 
#      tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs(file[,5])> exp_threshold)& (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold) & ((sign(file[,6])* sign(file[,7]))>0)),1]
#      tmp <- file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold) & (abs(file[,5])> exp_threshold)& (sign(file[,6])== sign(file[,8])) & ((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold)) ),1]
#       tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs( file[,5] + sign(file[,5]) * abs(file[,6])/2 ) > exp_threshold) & (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold)  ),1]
	
	  if (FDR) {
		threshold = 0.05
#		tmp <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		tmp <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
	  }else{
		threshold = 0.05/(nrow(file))
#		tmp <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		tmp <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])

	  }
        if (length(tmp)>=1){
            idx_row = match(tmp,Term_list)
            idx_col = match(paste0(outcome),sum_table_col)
			if (FDR) {
				sum_table_cell_type_function[idx_row,idx_col] <- file[match(tmp,file[,1]),3] *sign(file[match(tmp,file[,1]),8])
			}else{
				sum_table_cell_type_function[idx_row,idx_col] <- file[match(tmp,file[,1]),3] *sign(file[match(tmp,file[,1]),8])
			}
        }
	}

colnames(sum_table_cell_type_function) <- sum_table_col
rownames(sum_table_cell_type_function) <- Term_list
sum_table_cell_type_function_log10 <- -log10(abs(sum_table_cell_type_function))*sign(sum_table_cell_type_function)
sum_table_cell_type_function_log10_filter <- rbind(sum_table_cell_type_function_log10[,!colSums(sum_table_cell_type_function_log10)==0])
sum_table_cell_type_function_log10_filter <- abs(sum_table_cell_type_function_log10_filter)
#sum_table_cell_type_function_log10_filter <- sum_table_cell_type_function_log10_filter[!rowSums(sum_table_cell_type_function_log10)==0,]
#print(pheatmap(sum_table_cell_type_function_log10_filter,cellwidth = 18, cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter))))
##tsea.plot(sum_table_cell_type_function_log10_filter)

colnames(sum_table_cell_type_function_log10_filter) <- sapply(colnames(sum_table_cell_type_function_log10_filter),function(x){gsub("_no_adjAPOE","",x)})
rownames(sum_table_cell_type_function_log10_filter) <- Term_list

if (FDR){
write.table(sum_table_cell_type_function_log10_filter,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_without_filtration_noWB.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_without_filtration_noWB.pdf"),10,15)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_filter,cluster_rows =F,cellwidth = 18, angle_col = 315,cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter)))), error = function(e)NULL)
dev.off()
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_without_filtration_noWB_one_color.pdf"),10,15)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_filter,color = colorRampPalette(c("white", "firebrick1"))(50),cluster_rows =F,cellwidth = 18, angle_col = 315,cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter)))), error = function(e)NULL)
dev.off()

}else{
write.table(sum_table_cell_type_function_log10_filter,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_without_filtration_noWB.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_without_filtration_noWB.pdf"),10,15)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_filter,cluster_rows =F,cellwidth = 18, angle_col = 315,cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter)))), error = function(e)NULL)
dev.off()
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_without_filtration_noWB_one_color.pdf"),10,15)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_filter,color = colorRampPalette(c("white", "firebrick1"))(50),cluster_rows =F,cellwidth = 18, angle_col = 315,cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter)))), error = function(e)NULL)
dev.off()
}

##fill in the full box
Term_list <- unique(Term_list)
#construct one matrix to hold#
sum_table_col <- c(paste0(traits))
sum_table_cell_type_function_full <- matrix(1,ncol=length(sum_table_col ),nrow=length(Term_list))
colnames(sum_table_cell_type_function_full) <- sum_table_col
rownames(sum_table_cell_type_function_full) <- Term_list
##fill in the table
    for (outcome in traits) {
      file_name=paste0(status,"_",GMT_name,"_PC",PC_num,"_",outcome,"_",GP,"_combined.txt")
      file = read.delim(file_name, check.names=F,as.is=T)

#      tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs(file[,5])> exp_threshold)& (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold) & ((sign(file[,6])* sign(file[,7]))>0)),1]
#       tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs( file[,5] + sign(file[,5]) * abs(file[,6])/2 ) > exp_threshold) & (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold)  ),1]
		tmp  <- file[match(Term_list,file[,1]),1]
        if (length(tmp)>=1){
			for (m in tmp){
				idx_row = m
				idx_col = match(paste0(outcome),sum_table_col)
				if (FDR) {
					sum_table_cell_type_function_full[idx_row,idx_col] <- file[match(m,file[,1]),3] *sign(file[match(m,file[,1]),8])
				}else{
					sum_table_cell_type_function_full[idx_row,idx_col] <- file[match(m,file[,1]),3] *sign(file[match(m,file[,1]),8])
				}
			}
		}
	}

colnames(sum_table_cell_type_function_full) <- sum_table_col
rownames(sum_table_cell_type_function_full) <- Term_list
sum_table_cell_type_function_log10_full <- -log10(abs(sum_table_cell_type_function_full))*sign(sum_table_cell_type_function_full)
#sum_table_cell_type_function_log10_full <- sum_table_cell_type_function_log10_full[,match(colnames(sum_table_cell_type_function_log10_filter),colnames(sum_table_cell_type_function_log10_full))]
#sum_table_cell_type_function_log10_filter <- sum_table_cell_type_function_log10_filter[!rowSums(sum_table_cell_type_function_log10)==0,]
sum_table_cell_type_function_log10_full <- rbind(sum_table_cell_type_function_log10_full[,c(1,2,3,4,5,6)])
rownames(sum_table_cell_type_function_log10_full) <- Term_list


#pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,gaps_col = cumsum(c(4,4,4,3,3)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(0.05/(nrow(file))), "*", ""), nrow(sum_table_cell_type_function_log10_full)))
if (length(nrow(sum_table_cell_type_function_log10_full))>0){
colnames(sum_table_cell_type_function_log10_full) <- sapply(colnames(sum_table_cell_type_function_log10_full),function(x){gsub("_no_adjAPOE","",x)})
#print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,gaps_col = cumsum(c(1,4,4,4,3,3)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14))

sum_table_cell_type_function_log10_full <- abs(sum_table_cell_type_function_log10_full)

sum_table_cell_type_function_log10_full_out <- sum_table_cell_type_function_log10_full
sum_table_cell_type_function_log10_full_out = matrix(as.integer(sum_table_cell_type_function_log10_full_out>-log10(threshold)),ncol=ncol(sum_table_cell_type_function_log10_full))
colnames(sum_table_cell_type_function_log10_full_out) <- colnames(sum_table_cell_type_function_log10_full)
rownames(sum_table_cell_type_function_log10_full_out) <- rownames(sum_table_cell_type_function_log10_full)



if (FDR){
write.table(sum_table_cell_type_function_log10_full_out,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full_without_filtration_noWB_binary.txt"),quote=F,sep="\t")
write.table(sum_table_cell_type_function_log10_full,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full_without_filtration_noWB.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full_without_filtration_noWB.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,2)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full_without_filtration_noWB_one_color.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,color = colorRampPalette(c("white", "firebrick1"))(50),angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,2)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()
}else{
write.table(sum_table_cell_type_function_log10_full_out,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full_without_filtration_noWB_binary.txt"),quote=F,sep="\t")
write.table(sum_table_cell_type_function_log10_full,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full_without_filtration_noWB.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full_without_filtration_noWB.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,2)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full_without_filtration_noWB_one_color.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,color = colorRampPalette(c("white", "firebrick1"))(50),angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,2)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()



}
}else{
next
}
}
}
}

###
#FDR = FALSE
###use the absolute value###
##canonical,go,brain_cell_type_function_AD_pathways together
for (GP in c(GP)) {
	for (PC in c(paste0("PC",PC_num))) {
#		for (Term in c("go","canonical","brain_cell_type_function_AD_pathways")){
		for (Term in c(GMT_name)){
print(c(GP,PC,Term))
#GP = 10
library(pheatmap)
traits <- c("case_control",paste0("TB",GP,"all"),paste0("TB",GP,"SCZ"),paste0("TB",GP,"Ctr"),paste0("T",GP),paste0("B",GP))
#GP = 10
#model = c("_Brain_whole_blood","_Brain_only")
##get qualified term 
#threshold = 0.05/(167*16) # 0.001,0.0005
#FC_threshold = 0.1478
FC_threshold = 0.0
#FC_threshold = 0.196
exp_threshold = 0.0

Term_list <-{}
all_file <- {}
    for (outcome in traits) {
##/data2/ydai2/projects/22_AD_TWAS_PRS/GSVA/2722_inAPOE_result/brain_cell_type_function_AD_pathways_no_adjAPOE_DEG/PC5/diag
      file_name=paste0(status,"_",GMT_name,"_PC",PC_num,"_",outcome,"_",GP,"_combined.txt")
      file = read.delim(file_name, check.names=F,as.is=T)
	  FC_combine_threshold = quantile(x=abs(file$logFC_combined),probs=c(0.5))
	  all_file = rbind(all_file,file)
#	  threshold = 0.05/(nrow(file))
	  if (FDR) {
		threshold = 0.05
#		Term_list <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		Term_list <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold)  &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
	  }else{
		threshold = 0.05/(nrow(file))
#		Term_list <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		Term_list <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold)  &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])

	  }

    }


##first use raw pvalue to cut a threshold
quantile(x=abs(all_file$direction_out[all_file$raw_p<0.05/nrow(file)]),probs=c(0.05,0.25,0.5,0.95,0.975,1))
quantile(x=abs(all_file$logFC[all_file$raw_p<0.05/nrow(file)] ),probs=c(0.05,0.25,0.4,0.5,0.75,0.8,0.95,0.975,1))
quantile(x=abs(all_file$logFC_combined[all_file$raw_p<0.05/nrow(file)] ),probs=c(0.05,0.25,0.4,0.5,0.75,0.8,0.95,0.975,1))
Term_list <- unique(Term_list)
#construct one matrix to hold#
sum_table_col <- c(paste0(traits))
sum_table_cell_type_function <- matrix(1,ncol=length(sum_table_col ),nrow=length(Term_list))
colnames(sum_table_cell_type_function) <- sum_table_col
rownames(sum_table_cell_type_function) <- Term_list
##fill in the table
    for (outcome in traits) {
      file_name=paste0(status,"_",GMT_name,"_PC",PC_num,"_",outcome,"_",GP,"_combined.txt")
      file = read.delim(file_name, check.names=F,as.is=T)
	  FC_combine_threshold = quantile(x=abs(file$logFC_combined),probs=c(0.5))
#      tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs(file[,5])> exp_threshold)& (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold) & ((sign(file[,6])* sign(file[,7]))>0)),1]
#      tmp <- file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold) & (abs(file[,5])> exp_threshold)& (sign(file[,6])== sign(file[,8])) & ((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold)) ),1]
#       tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs( file[,5] + sign(file[,5]) * abs(file[,6])/2 ) > exp_threshold) & (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold)  ),1]
	
	  if (FDR) {
		threshold = 0.05
#		tmp <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		tmp <- c(Term_list,file[((file[,4]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
	  }else{
		threshold = 0.05/(nrow(file))
#		tmp <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &  (sign(file[,6])== sign(file[,8])) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])
		tmp <- c(Term_list,file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold)  & (abs(file[,5])> exp_threshold) &((abs(file[,6])> FC_combine_threshold)|(abs(file[,8])> FC_combine_threshold))  ),1])

	  }
        if (length(tmp)>=1){
            idx_row = match(tmp,Term_list)
            idx_col = match(paste0(outcome),sum_table_col)
			if (FDR) {
				sum_table_cell_type_function[idx_row,idx_col] <- file[match(tmp,file[,1]),3] *sign(file[match(tmp,file[,1]),8])
			}else{
				sum_table_cell_type_function[idx_row,idx_col] <- file[match(tmp,file[,1]),3] *sign(file[match(tmp,file[,1]),8])
			}
        }
	}
}

colnames(sum_table_cell_type_function) <- sum_table_col
rownames(sum_table_cell_type_function) <- Term_list
sum_table_cell_type_function_log10 <- -log10(abs(sum_table_cell_type_function))*sign(sum_table_cell_type_function)
sum_table_cell_type_function_log10_filter <- rbind(sum_table_cell_type_function_log10[,!colSums(sum_table_cell_type_function_log10)==0])
sum_table_cell_type_function_log10_filter <- abs(sum_table_cell_type_function_log10_filter)
#sum_table_cell_type_function_log10_filter <- sum_table_cell_type_function_log10_filter[!rowSums(sum_table_cell_type_function_log10)==0,]
#print(pheatmap(sum_table_cell_type_function_log10_filter,cellwidth = 18, cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter))))
##tsea.plot(sum_table_cell_type_function_log10_filter)

colnames(sum_table_cell_type_function_log10_filter) <- sapply(colnames(sum_table_cell_type_function_log10_filter),function(x){gsub("_no_adjAPOE","",x)})
rownames(sum_table_cell_type_function_log10_filter) <- Term_list

if (FDR){
write.table(sum_table_cell_type_function_log10_filter,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR.pdf"),10,15)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_filter,cluster_rows =F,cellwidth = 18, angle_col = 315,cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter)))), error = function(e)NULL)
dev.off()
}else{
write.table(sum_table_cell_type_function_log10_filter,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_abs.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_abs.pdf"),10,15)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_filter,cluster_rows =F,cellwidth = 18, angle_col = 315,cellheight = 18,display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_filter) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_filter)))), error = function(e)NULL)
dev.off()
}

##fill in the full box
Term_list <- unique(Term_list)
#construct one matrix to hold#
sum_table_col <- c(paste0(traits))
sum_table_cell_type_function_full <- matrix(1,ncol=length(sum_table_col ),nrow=length(Term_list))
colnames(sum_table_cell_type_function_full) <- sum_table_col
rownames(sum_table_cell_type_function_full) <- Term_list
##fill in the table
    for (outcome in traits) {
      file_name=paste0(status,"_",GMT_name,"_PC",PC_num,"_",outcome,"_",GP,"_combined.txt")
      file = read.delim(file_name, check.names=F,as.is=T)

#      tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs(file[,5])> exp_threshold)& (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold) & ((sign(file[,6])* sign(file[,7]))>0)),1]
#       tmp <- file[((file[,3]<threshold) & (abs(file[,6]) > FC_threshold) & (sign(file[,6])== sign(file[,8])) & (abs( file[,5] + sign(file[,5]) * abs(file[,6])/2 ) > exp_threshold) & (sign(file[,5])== sign(file[,7])) & (abs(file[,8])> FC_combine_threshold)  ),1]
		tmp  <- file[match(Term_list,file[,1]),1]
        if (length(tmp)>=1){
			for (m in tmp){
				idx_row = m
				idx_col = match(paste0(outcome),sum_table_col)
				if (FDR) {
					sum_table_cell_type_function_full[idx_row,idx_col] <- file[match(m,file[,1]),3] *sign(file[match(m,file[,1]),8])
				}else{
					sum_table_cell_type_function_full[idx_row,idx_col] <- file[match(m,file[,1]),3] *sign(file[match(m,file[,1]),8])
				}
			}
		}
	}


colnames(sum_table_cell_type_function_full) <- sum_table_col
rownames(sum_table_cell_type_function_full) <- Term_list
sum_table_cell_type_function_log10_full <- -log10(abs(sum_table_cell_type_function_full))*sign(sum_table_cell_type_function_full)
#sum_table_cell_type_function_log10_full <- sum_table_cell_type_function_log10_full[,match(colnames(sum_table_cell_type_function_log10_filter),colnames(sum_table_cell_type_function_log10_full))]
#sum_table_cell_type_function_log10_filter <- sum_table_cell_type_function_log10_filter[!rowSums(sum_table_cell_type_function_log10)==0,]
sum_table_cell_type_function_log10_full <- rbind(sum_table_cell_type_function_log10_full[,c(1,2,3,4,5,6)])
rownames(sum_table_cell_type_function_log10_full) <- Term_list
#pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,gaps_col = cumsum(c(4,4,4,3,3)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(0.05/(nrow(file))), "*", ""), nrow(sum_table_cell_type_function_log10_full)))
if (length(nrow(sum_table_cell_type_function_log10_full))>0){
colnames(sum_table_cell_type_function_log10_full) <- sapply(colnames(sum_table_cell_type_function_log10_full),function(x){gsub("_no_adjAPOE","",x)})
#print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,gaps_col = cumsum(c(1,4,4,4,3,3)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14))

sum_table_cell_type_function_log10_full <- abs(sum_table_cell_type_function_log10_full)

if (FDR){
write.table(sum_table_cell_type_function_log10_full,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_","abs_FDR_full.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,2)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()
}else{
write.table(sum_table_cell_type_function_log10_full,file=paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full.txt"),quote=F,sep="\t")
pdf(paste0("sum_table_",PC,"_",Term,"_",GP,"_abs_full.pdf"),15,10)
tryCatch(print(pheatmap(sum_table_cell_type_function_log10_full,angle_col = 315,cluster_cols = F,cluster_rows =F,gaps_col = cumsum(c(1,3,2)),display_numbers = matrix(ifelse(abs(sum_table_cell_type_function_log10_full) > -log10(threshold), "*", ""), nrow(sum_table_cell_type_function_log10_full)),cellwidth = 14, cellheight = 14)), error = function(e)NULL)
dev.off()

}
}else{
next
}
}
}

}
