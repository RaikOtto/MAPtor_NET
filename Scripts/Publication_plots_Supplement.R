library("stringr")
library("dplyr")
library("umap")
library("ggplot2")

# Supplementary Figure 1 Exclusion Discovery

expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Discovery_Cohort.S81.HGNC.DESeq2.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
#expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Discovery_Cohort.S64.DESeq2.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
#expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Discovery_Cohort.S64.DESeq2.vsd.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X\\.)", "")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X)", "")
expr_raw[1:5,1:5]
dim(expr_raw)

library("grid")

meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)

meta_info$Name = str_replace(meta_info$Sample, pattern = "^(X\\.)", "")
rownames(meta_info) = meta_info$Sample

no_match = match(colnames(expr_raw), meta_info$Sample, nomatch = F) == F
colnames(expr_raw)[no_match  ]
colnames(expr_raw)[no_match  ] = paste("X",colnames(expr_raw)[no_match  ],sep = "")
no_match = match(colnames(expr_raw), meta_info$Sample, nomatch = F) == F
meta_data = meta_info[ colnames(expr_raw), ]
table(no_match)
table(meta_data$Primary_Metastasis)

rownames(meta_data) = meta_data$Sample
expr_raw = expr_raw[,meta_data$Sample[meta_data$Grading != ""]]

dim(expr_raw)

table(meta_data$Grading)

##
source("~/MAPTor_NET/Misc/Visualization_colors.R")

###

i = 13
genes_of_interest_hgnc_t = read.table("~/MAPTor_NET//Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
sad_genes = sad_genes[sad_genes != ""]

liver_genes = genes_of_interest_hgnc_t[70,3:ncol(genes_of_interest_hgnc_t)]
sad_genes = sad_genes[!(sad_genes %in% liver_genes)]

#hox_genes = c("HOXA1","HOXA2","HOXA3","HOXA4","HOXA5","HOXA6","HOXA7","HOXA9","HOXA10","HOXA11","HOXA13","HOXB1","HOXB2","HOXB3","HOXB4","HOXB5","HOXB6","HOXB7","HOXB8","HOXB9","HOXB13","HOXC4","HOXC5","HOXC6","HOXC8","HOXC9","HOXC10","HOXC11","HOXC12","HOXC13","HOXD1","HOXD3","HOXD4","HOXD8","HOXD9","HOXD10","HOXD11","HOXD12","HOXD13")
#sad_genes = hox_genes

genes_of_interest_hgnc_t[i,1]

sad_genes[which(!(sad_genes %in% rownames(expr_raw)))]
table(sad_genes %in% rownames(expr_raw) )

row_var = as.double(apply(expr_raw , FUN = var, MARGIN = 1))
threshold = quantile(row_var,seq(0,1,0.01))[2]
expr_raw_vis = expr_raw[row_var >= threshold,]

expr = expr_raw_vis[ rownames(expr_raw_vis) %in% sad_genes,]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))
meta_data = meta_info[ colnames(expr), ]

## Supplementary D Heatmap

p = pheatmap::pheatmap(
  #(expr),
  cor_mat,
  #annotation_col = meta_data[c("Included","Histology_Primary","Primary_Metastasis","NET_NEC_UMAP","NEC_NET","Grading","Study")],
  annotation_col = meta_data[c("Included","Histology_Primary","Primary_Metastasis","Grading","Study")],
  annotation_colors = aka3,
  show_rownames = FALSE,
  show_colnames = TRUE,
  treeheight_row = 0,
  legend = TRUE,
  annotation_legend = TRUE,
  #fontsize_col = 7,
  clustering_method = "ward.D2"
)

svg("~/Downloads/MAPTor-NET_plots_20_10_2022/1_SM_Figure_exclusion_Discovery.svg", width = 10, height = 10)
print(p)
dev.off()

# Supplementary Figure 2 Exclusion Validation

#expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Discovery_Cohort.S81.HGNC.DESeq2.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Fröhling.S34.HGNC.DESeq2.VOOM.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
#expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Discovery_Cohort.S64.DESeq2.vsd.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X\\.)", "")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X)", "")
expr_raw[1:5,1:5]
dim(expr_raw)

library("grid")

meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)

meta_info$Name = str_replace(meta_info$Sample, pattern = "^(X\\.)", "")
rownames(meta_info) = meta_info$Sample

no_match = match(colnames(expr_raw), meta_info$Sample, nomatch = F) == F
colnames(expr_raw)[no_match  ]
colnames(expr_raw)[no_match  ] = paste("X",colnames(expr_raw)[no_match  ],sep = "")
no_match = match(colnames(expr_raw), meta_info$Sample, nomatch = F) == F
meta_data = meta_info[ colnames(expr_raw), ]
table(no_match)
table(meta_data$Primary_Metastasis)

rownames(meta_data) = meta_data$Sample
expr_raw = expr_raw[,meta_data$Sample[meta_data$Grading != ""]]

dim(expr_raw)

table(meta_data$Grading)

##
source("~/MAPTor_NET/Misc/Visualization_colors.R")

###

i = 13
genes_of_interest_hgnc_t = read.table("~/MAPTor_NET//Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
sad_genes = sad_genes[sad_genes != ""]

liver_genes = genes_of_interest_hgnc_t[70,3:ncol(genes_of_interest_hgnc_t)]
sad_genes = sad_genes[!(sad_genes %in% liver_genes)]

#hox_genes = c("HOXA1","HOXA2","HOXA3","HOXA4","HOXA5","HOXA6","HOXA7","HOXA9","HOXA10","HOXA11","HOXA13","HOXB1","HOXB2","HOXB3","HOXB4","HOXB5","HOXB6","HOXB7","HOXB8","HOXB9","HOXB13","HOXC4","HOXC5","HOXC6","HOXC8","HOXC9","HOXC10","HOXC11","HOXC12","HOXC13","HOXD1","HOXD3","HOXD4","HOXD8","HOXD9","HOXD10","HOXD11","HOXD12","HOXD13")
#sad_genes = hox_genes

genes_of_interest_hgnc_t[i,1]

sad_genes[which(!(sad_genes %in% rownames(expr_raw)))]
table(sad_genes %in% rownames(expr_raw) )

row_var = as.double(apply(expr_raw , FUN = var, MARGIN = 1))
threshold = quantile(row_var,seq(0,1,0.01))[2]
expr_raw_vis = expr_raw[row_var >= threshold,]

expr = expr_raw_vis[ rownames(expr_raw_vis) %in% sad_genes,]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))
meta_data = meta_info[ colnames(expr), ]

## Supplementary D Heatmap

p = pheatmap::pheatmap(
  #(expr),
  cor_mat,
  #annotation_col = meta_data[c("Included","Histology_Primary","Primary_Metastasis","NET_NEC_UMAP","NEC_NET","Grading","Study")],
  annotation_col = meta_data[c("Included","Histology_Primary","Primary_Metastasis","Grading","Study")],
  annotation_colors = aka3,
  show_rownames = FALSE,
  show_colnames = TRUE,
  treeheight_row = 0,
  legend = TRUE,
  annotation_legend = TRUE,
  #fontsize_col = 7,
  clustering_method = "ward.D2"
)

#svg("~/Downloads/MAPTor-NET_plots_20_10_2022/2_SM_Figure_exclusion_Validation.svg", width = 10, height = 10)
print(p)
dev.off()

