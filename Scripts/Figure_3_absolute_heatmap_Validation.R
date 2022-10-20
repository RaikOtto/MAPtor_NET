library("ggplot2")
library("stringr")
library("umap")
library("dplyr")
library("grid")
library("tibble")

expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Fröhling.S34.HGNC.DESeq2.VOOM.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X\\.)", "")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X)", "")
nomatch = match(colnames(expr_raw),meta_info$Sample, nomatch = 0) == 0
colnames(expr_raw)[nomatch] = str_c("X",colnames(expr_raw)[nomatch])

meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info$Name = str_replace(meta_info$Sample, pattern = "^(X\\.)", "")
rownames(meta_info) = meta_info$Sample

##
source("~/MAPTor_NET/Misc//Visualization_colors.R")

genes_dif_Exp = read.table("~/MAPTor_NET/Misc/Differentially_expressed.tsv",sep ="\t", stringsAsFactors = F, header = TRUE)
genes_dif_Exp = genes_dif_Exp[,1]


i = 76
genes_of_interest_hgnc_t = read.table("~/MAPTor_NET//Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
genes_of_interest_hgnc_t$V1[i]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
sad_genes = sad_genes[sad_genes != ""]

liver_genes = genes_of_interest_hgnc_t[70,3:ncol(genes_of_interest_hgnc_t)]
sad_genes = sad_genes[!(sad_genes %in% liver_genes)]
sad_genes = sad_genes[sad_genes %in% rownames(expr_raw)]
sad_genes = sad_genes[sad_genes %in% genes_dif_Exp]

table(sad_genes %in% genes_of_interest_hgnc_t[82,3:nrow(genes_of_interest_hgnc_t)])

vis_genes = sad_genes #c("HOXB8","HOXB9","DLX6","HOXC12","HOXA11","HOXD13","HMGA2","HOXA6","MYBL2","IGF2BP2","SH2D2A","SLC38A5","HOXD10","HOXD11","HOXC10","HOXA5","HOXA7","HOXA9","HOXA10","HOXB5","HOXB7","HAUS7","SOCS2","SLC35F3")
men1_size_vec = as.double(expr_raw["MEN1",])

### A

expr = expr_raw[vis_genes,]
meta_data = meta_info[colnames(expr),]

meta_data$MEN1_MT = rep("Homozygous")
meta_data$MEN1_MT[meta_data$MEN1_Mut_AF < .9] = "Heterozygous"
meta_data$MEN1_MT[meta_data$MEN1_Mut_AF < .1] = "WT"
meta_data[,"MEN1_exp"] = as.double(expr_raw["MEN1",meta_data$Sample])
cormat = (expr)
#cormat = cormat[rownames(cormat) != "132502",colnames(cormat) != "132502"]

p = pheatmap::pheatmap(
  cormat,
  annotation_col = meta_data[c("MEN1_exp","MEN1_MT","Study","NET_NEC_UMAP","NEC_NET","Histology_Primary","Grading")],
  annotation_colors = aka3,
  show_rownames = TRUE,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = FALSE,
  annotation_legend = TRUE,
  clustering_method = "average"
)

#svg("~/Downloads/Figure_4_Plot_A.svg", width = 10, height = 10)
print(p)
#dev.off()

# 132502
