#heatmap
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(limma)
#Differential analysis for GSVA results
args = commandArgs(trailingOnly=TRUE)
PC_num = args[1] # "PC10"
work_directory = args[2] # "/collab2/ydai2/24_GRPa/GSVA/result/AD/"
cov_file = args[3] ## "/collab2/ydai2/24_GRPa/code/subgrp_ADcmeta0217.txt"
enrich_option = args[4] ## "brain_cell_type_function_AD_pathways"
APOE_status = args[5] ## "disc2722_MASHR"
folder = args[6] ## "2722inAPOE"
system(paste0("mkdir -p ",work_directory,folder))
system(paste0("mkdir -p ",work_directory,folder,"/",enrich_option,"_DEG/"))
setwd(work_directory)
files <- dir()
files_canonical <- files[grep(enrich_option,files)]
print(files_canonical)
files_canonical <- files_canonical[grep(APOE_status,files_canonical)]
print(files_canonical)
for (i in files_canonical) {

logFCcutoff<-0.1
adjPvalueCutoff<-0.05 #0.01, 0.001
##for each canonical file, do DEG analysis##
##process phenotype data##
##cov label 15,16,74
#PC_num = "PC10"
#PC_num = "PC5"

##cov_16 <- read.delim("/collab2/ydai2/21_AD_TWAS/PrediXcan_imputed_expression/Phenotype/58516/covs.txt",as.is=T,check.names=F)
#groups <- c("clinic","patho","subtype","Sgroup4_5","Sgroup4_10","Sgroup4_15","Sgroup4_20","Sgroup5","Sgroup6","Sgroup7","Sgroup8","Kgroup4_5","Kgroup4_10","Kgroup4_15","Kgroup4_20","Kgroup5","Kgroup6","Kgroup7","Kgroup8")
groups <- c("case_control","TB10all","TB10SCZ","TB10Ctr","TB15all","TB15SCZ","TB15Ctr","TB20all","TB20SCZ","TB20Ctr","T10","B10","T15","B15","T20","B20")
for (j in 1:length(groups)){
print(groups[j])
cov <- read.delim(cov_file,as.is=T,check.names=F)
#cov <- cov[!is.na(cov$apoeGenotype),]
setwd(work_directory)
rt<-read.delim(i,check.names = F,row.names = 1)
i_name <- gsub(paste0("gsvaOut_", enrich_option,"_mixdiff_T_"),"",i)

#rt <- rt[,!is.na(cov$apoeGenotype)]
#cov$age_new <- cov$age_at_visit_max
#cov$age_new[!is.na(cov$age_first_ad_dx)] <- cov$age_first_ad_dx[!is.na(cov$age_first_ad_dx)]
#cov$age_new[which(cov$age_new=="90+")] = 90
#cov$age_new <- as.numeric(cov$age_new)
subgroup = groups[j]
print(subgroup)
setwd(paste0(work_directory,folder,"/",enrich_option,"_DEG/"))
system(paste0("mkdir -p ",PC_num))
system(paste0("mkdir -p ",PC_num,"/",subgroup))

match_idx <- c(which(cov[,which(colnames(cov)==groups[j])]==0),which(cov[,which(colnames(cov)==groups[j])]==1))
#for (m in 1:nrow(cov)) {match_idx <- c(match_idx,grep(cov[m,1],colnames(rt)))}
rt_filter <- rt[,match_idx]
cov_filter<- cov[match_idx,]
type1<-factor(cov_filter[,which(colnames(cov)==groups[j])])
PC1 <- cov_filter$X1
PC2 <- cov_filter$X2
PC3 <- cov_filter$X3
PC4 <- cov_filter$X4
PC5 <- cov_filter$X5
PC6 <- cov_filter$X6
PC7 <- cov_filter$X7
PC8 <- cov_filter$X8
PC9 <- cov_filter$X9
PC10 <- cov_filter$X10
sex <- factor(cov_filter$sex)
#age_new <- cov_filter$age_new

if(PC_num == "PC10") {
design=model.matrix(~0+type1+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+sex)
}
if (PC_num == "PC5") {
design=model.matrix(~0+type1+PC1+PC2+PC3+PC4+PC5+sex)
}
if (PC_num == "PC3") {
design=model.matrix(~0+type1+PC1+PC2+PC3+sex)
}

print(head(design))
colnames(design)[1:2]<-c("Control","Case")
row.names(design)<-colnames(rt_filter) #
contrast.matrix<-makeContrasts("Case-Control",levels=design)#
fit=lmFit(rt_filter,design)
fit =contrasts.fit(fit,contrast.matrix) #
fit<-eBayes(fit)
all<-topTable(fit,coef = 1 ,n=Inf,sort.by="none" )
all<-data.frame(all, "residual_se"=fit$sigma,"se"=fit$stdev.unscaled,check.names=F)
colnames(all)[ncol(all)] <- "beta_se"
##get cohen's d
effect_size= (all$logFC) / fit$sigma
print(effect_size)
all<-data.frame(all, "effect_size"= effect_size, check.names=F)
all<-all[order(abs(all$effect_size),decreasing=T),]
write.table(all,file=paste0(PC_num,"/",subgroup,"/",i_name,".txt"),sep="\t",quote=F,col.names = T)

df<-all[all$adj.P.Val<adjPvalueCutoff & abs(all$logFC)>logFCcutoff,]
dim(df)
#diffname<-c()
diffname<-row.names(df) #25+43+9+17+9
length(diffname)
diffname<-diffname[!duplicated(diffname)]
write.table(df,file=paste0(PC_num,"/",subgroup,"/","GSVA_DEG_",enrich_option,"_",i_name,".txt"),sep="\t",quote=F,row.names=T)

#anno$Case_Control_order <-factor(anno$Case_Control,levels=c(0,1))
#anno<-anno[order(anno$Case_Control_order),]
#rt_filter_order<-rt_filter[,anno$ID]

volcano = all
volcano$pathway = rownames(all)
volcano$FDR = -log(volcano$adj.P.Val, 10)
volcano.sig = df
volcano.sig$pathway = rownames(df)
volcano.sig$FDR = -log(volcano.sig$adj.P.Val, 10)
mycol <- c(brewer.pal(9,"Set1"),brewer.pal(8, "Set2"),brewer.pal(12,"Set3")[c(-2,-12)],brewer.pal(12,"Paired"))
p <- ggplot(volcano) +
            	           geom_point(data = volcano, aes(x = logFC, y = FDR), color = mycol[2], cex = 1) +
                          geom_point(data = volcano.sig, aes(x = logFC, y = FDR), color = mycol[1], cex = 1) +
                          theme_bw() + xlab("Log (fold change)") + ylab("-Log10 (FDR)") +
                          geom_vline(xintercept = 0.1, col = mycol[3], linetype = "dotted", size = 1) +
                          geom_vline(xintercept = -0.1, col = mycol[3], linetype = "dotted", size = 1) +
                          geom_hline(yintercept = -log(0.05, 10), col = mycol[3], linetype = "dotted", size = 1)

p <- p + theme(axis.title.x = element_text(size = 15, color = "black", face = "bold", angle = 0))
p <- p + theme(axis.title.y = element_text(size = 15, color = "black", face = "bold", angle = 90))
p <- p + theme(axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))
p <- p + theme(axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))

p1 <- p + geom_text(data = volcano.sig, aes(x = logFC, y = FDR, label = pathway), hjust = -0.5, vjust = 0)
#p2 <- p + ylim(0, 15)

pdf(paste0(PC_num,"/",subgroup,"/",enrich_option,"_",i_name,"_volcano.pdf"), 7, 7)
print (p)
print (p1)
dev.off()
#names(Type)<-colnames(rt_filter_order)
#Type<-as.data.frame(Type)
}
}

