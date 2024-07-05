##input two arguments, first: samples from "1-13" ; second: enrich options from c("GO","canonical","brain")
library(GSVA)
library(GSEABase)
library(pheatmap)
library(GSVAdata)
data(c2BroadSets)
args = commandArgs(trailingOnly=TRUE)
input <- as.numeric(args[1])
enrich_option <- args[2]
output_path <- args[3] #"/collab2/ydai2/24_GRPa/GSVA/result/AD/"
TWAS_path <- args[4] # "/collab2/ydai2/24_GRPa/TWAS/AD/disc_inAPOE/"
work_path <- args[5] # "/collab2/ydai2/24_GRPa/GSVA/code"
status <- args[6] # "2722inAPOE"
setwd(work_path)
system(paste0("mkdir -p ",output_path,status))

###RUN GSVA
##get canonical_pathway
gmtfiles_pathway <-"c2.all.v7.4.symbols.gmt"
geneset_pathway <- getGmt(gmtfiles_pathway,collectionType = BroadCollection(category = "c2"),
                          geneIdType = SymbolIdentifier())
canonical_geneset_pathway <- geneset_pathway[c(grep("^KEGG", names(geneset_pathway)), grep("^REACTOME", names(geneset_pathway)), grep("^BIOCARTA", names(geneset_pathway)))]
#dmGWAS <- "dm_GWAS_genes_symbol.gmt"
#geneset_dmGWAS <- getGmt(dmGWAS,collectionType=BroadCollection(category="c5"), geneIdType = SymbolIdentifier())
cell_type_function <- "brain_cell_type_function_synapse_immune.gmt"
geneset_cell_type_function <- getGmt(cell_type_function,collectionType=BroadCollection(category="c5"), geneIdType = SymbolIdentifier())
brain_cell_type_function_AD_pathways <- "brain_cell_type_function_AD_pathways.gmt"
geneset_brain_cell_type_function_AD_pathways <- getGmt(brain_cell_type_function_AD_pathways,collectionType=BroadCollection(category="c5"),geneIdType = SymbolIdentifier())
##get 
#gmtfiles_go <-"c5.go.v7.4.symbols.gmt"
#non redundant
gmtfiles_go <- "GO_merged_non_redundant_full_term_name_gene_symbol_rename.gmt"
geneset_go <- getGmt(gmtfiles_go,collectionType=BroadCollection(category="c5"),
                     geneIdType = SymbolIdentifier())

h19 <- read.delim("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",as.is=T,head=F)
h19 <- h19[-c(1,2),]
colnames(h19) <- h19[1,]
h19 <- h19[-1,]
files <- dir(TWAS_path)
for (i in files[input]) {

  file_name = gsub("_predict.txt","",i)
  print(file_name)

  MASHR_result <- read.delim(paste0(TWAS_path,i),as.is=T,check.names=F)
  MASHR_result <- MASHR_result[,-1]
  new_file_filter <- MASHR_result
  rownames(new_file_filter) <- new_file_filter[,1]
  new_file_filter <- new_file_filter[,-1]
  new_file_filter_t <- t(new_file_filter)
  new_file_filter_t <- as.matrix(new_file_filter_t,ncol = nrow(new_file_filter))
  idx_19 <- match(rownames(new_file_filter_t),h19[,1])
  rownames(new_file_filter_t) <- h19[idx_19,2]
  setwd(output_path)
  if (enrich_option == "GO") { 
    gsvaOut_go <- gsva(new_file_filter_t,geneset_go,
                       min.sz=5,     #size of the resulting gene sets.
                       max.sz=500,
                       #                mx.diff=FALSE,#two approaches to calculate the enrichment statistic (ES) (default):TRUE
                       verbose=TRUE) #Gives information about each calculation step. Default: FALSE.
    gsvaOut_go <- rbind(id=colnames(gsvaOut_go),gsvaOut_go)
    write.table(gsvaOut_go,file=paste0("gsvaOut_go_mixdiff_T_",file_name),sep="\t",quote=F,col.names = F)
  }
  ##gsva_canonical
  if (enrich_option == "canonical") { 
    gsvaOut_canonical <- gsva(new_file_filter_t,canonical_geneset_pathway,
                              min.sz=5,     #size of the resulting gene sets.
                              max.sz=500,
                              #                mx.diff=FALSE,#two approaches to calculate the enrichment statistic (ES) (default):TRUE
                              verbose=TRUE) #Gives information about each calculation step. Default: FALSE.
    gsvaOut_canonical <- rbind(id=colnames(gsvaOut_canonical),gsvaOut_canonical)
    write.table(gsvaOut_canonical,file=paste0("gsvaOut_canonical_mixdiff_T_",file_name),sep="\t",quote=F,col.names = F)	
  }
  ##brain_cell_type_function_synapse_immune
    if (enrich_option == "brain_cell_type_synapse_immune") { 
    gsvaOut_brain_cell_type_function_synapse_immune <- gsva(new_file_filter_t,geneset_cell_type_function,
                              min.sz=5,     #size of the resulting gene sets.
                              max.sz=500,
                              #                mx.diff=FALSE,#two approaches to calculate the enrichment statistic (ES) (default):TRUE
                              verbose=TRUE) #Gives information about each calculation step. Default: FALSE.
    gsvaOut_brain_cell_type_function_synapse_immune <- rbind(id=colnames(gsvaOut_brain_cell_type_function_synapse_immune),gsvaOut_brain_cell_type_function_synapse_immune)
    write.table(gsvaOut_brain_cell_type_function_synapse_immune,file=paste0("gsvaOut_brain_cell_type_function_synapse_immune_mixdiff_T_",file_name),sep="\t",quote=F,col.names = F)	
  }
  ##brain_cell_type_function_AD_pathways
  if (enrich_option == "brain_cell_type_function_AD_pathways") { 
    gsvaOut_brain_cell_type_function_AD_pathways <- gsva(new_file_filter_t,geneset_brain_cell_type_function_AD_pathways,
                                       min.sz=5,     #size of the resulting gene sets.
                                       max.sz=500,
                                       #two approaches to calculate the enrichment statistic (ES) (default):TRUE
                                       verbose=TRUE) #Gives information about each calculation step. Default: FALSE.
    gsvaOut_brain_cell_type_function_AD_pathways <- rbind(id=colnames(gsvaOut_brain_cell_type_function_AD_pathways),gsvaOut_brain_cell_type_function_AD_pathways)
    write.table(gsvaOut_brain_cell_type_function_AD_pathways,file=paste0("gsvaOut_brain_cell_type_function_AD_pathways_mixdiff_T_",file_name),sep="\t",quote=F,col.names = F)	
  }
}
