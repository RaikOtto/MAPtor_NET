library("ggplot2")
library("stringr")
library("umap")
library("dplyr")
library("grid")
library("tibble")

expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Fröhling.S34.HGNC.DESeq2.VOOM.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X\\.)", "")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X)", "")
expr_raw[1:5,1:5]

nomatch_index = which(! (colnames(expr_raw) %in% rownames(meta_info)))
colnames(expr_raw)[nomatch_index] = str_c("X",colnames(expr_raw)[nomatch_index])

meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info$Name = str_replace(meta_info$Sample, pattern = "^(X\\.)", "")
rownames(meta_info) = meta_info$Sample

##
source("~/MAPTor_NET/Misc//Visualization_colors.R")

vis_genes = c("HOXB8","HOXB9","DLX6","HOXC12","HOXA11","HOXD13","HMGA2","HOXA6","MYBL2","IGF2BP2","SH2D2A","SLC38A5","HOXD10","HOXD11","HOXC10","HOXA5","HOXA7","HOXA9","HOXA10","HOXB5","HOXB7","HAUS7","SOCS2","SLC35F3")

### A ####

meta_data = meta_info[colnames(expr_raw),]
expr_raw_pan = expr_raw[,meta_data$Histology_Primary == "Pancreatic"]
expr = expr_raw_pan[vis_genes,]
meta_data = meta_info[colnames(expr),]

meta_data$MEN1_MT = rep("Homozygous",ncol(expr))
meta_data$MEN1_MT[meta_data$MEN1_Mut_AF < .9] = "Heterozygous"
meta_data$MEN1_MT[meta_data$MEN1_Mut_AF < .1] = "WT"
meta_data[,"MEN1_exp"] = as.double(expr_raw_pan["MEN1",])

p = pheatmap::pheatmap(
  cor(expr),
  annotation_col = meta_data[c("MEN1_exp","MEN1_MT","Study","NET_NEC_UMAP","NEC_NET","Histology_Primary","Grading")],
  annotation_colors = aka3,
  show_rownames = FALSE,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = FALSE,
  annotation_legend = TRUE,
  #fontsize_col = 7,
  clustering_method = "ward.D",
  cluster_cols = TRUE
)

#svg("~/Downloads/Figure_3_Plot_C.svg", width = 10, height = 10)
print(p)
dev.off()


### B ####

expr = expr_raw[vis_genes,]
meta_data = meta_info[colnames(expr),]

meta_data$MEN1_MT = rep("Homozygous",ncol(expr))
meta_data$MEN1_MT[meta_data$MEN1_Mut_AF < .9] = "Heterozygous"
meta_data$MEN1_MT[meta_data$MEN1_Mut_AF < .1] = "WT"
meta_data[,"MEN1_exp"] = as.double(expr_raw["MEN1",])

p = pheatmap::pheatmap(
  cor(expr),
  annotation_col = meta_data[c("MEN1_exp","MEN1_MT","Study","NET_NEC_UMAP","NEC_NET","Histology_Primary","Grading")],
  annotation_colors = aka3,
  show_rownames = FALSE,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = FALSE,
  annotation_legend = TRUE,
  #fontsize_col = 7,
  clustering_method = "average",
  cluster_cols = TRUE
)

#svg("~/Downloads/Figure_3_Plot_D.svg", width = 10, height = 10)
print(p)
dev.off()

