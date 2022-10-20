library("ggpubr")
library("stringr")
library("reshape2")
library("ggplot2")
library("dplyr")
library("grid")

#443 407 444 579 450 409 452 PNET08

#draw_colnames_45 <- function (coln, gaps, ...) {
#  coord = pheatmap:::find_coordinates(length(coln), gaps)
#  x = coord$coord - 0.5 * coord$size
#  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
#  return(res)}
#assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

#meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
#meta_info$NEC_NET = meta_info$NEC_NET_PCA

#expr_raw = read.table("~/Deko_Projekt/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Muraro.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = T)
expr_raw = read.table("~/Deko_Projekt/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Segerstolpe.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = T)
#expr_raw = read.table("~/Deko_Projekt/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = T)
dim(expr_raw)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "\\.", "_")
expr_raw[1:5,1:5]
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match] = paste("X",colnames(expr_raw)[no_match],sep ="")
no_match = colnames(expr_raw) %in% meta_info$Sample == F

table(no_match)
meta_data = meta_info[colnames(expr_raw),]

### downsampling

cell_types = names(table(meta_data$Subtype))
cell_type_selection_vector = meta_data$Subtype %in% c("Alpha","Beta","Gamma","Delta","Acinar","Ductal","endothelial","mesenchymal")
expr_raw = expr_raw[cell_type_selection_vector]
dim(expr_raw)
meta_data = meta_info[colnames(expr_raw),]
table(meta_data$Subtype)

for ( cell_type in cell_types){
  
  print(cell_type)
  nr_samples = as.integer(table(meta_data$Subtype)[cell_type])
  selection = which(meta_data$Subtype == names(table(meta_data$Subtype)[cell_type]))
  
  if (nr_samples > 100){
    down_sample_selection = sample(selection, size = nr_samples - 100)
    expr_raw = expr_raw[,-down_sample_selection]
  }
  
  #if (nr_samples < 100){
  #  down_sample_selection = selection
  #  expr_raw = expr_raw[,-down_sample_selection]
  #}
  meta_data = meta_info[colnames(expr_raw),]
  
}
table(meta_data$Subtype)

source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
#genes_of_interest_hgnc_t = read.table("~/Deko_Projekt/Misc//Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t = read.table("~/MAPTor_NET//Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1

liver_genes = genes_of_interest_hgnc_t[70,3:ncol(genes_of_interest_hgnc_t)]
i = 63
genes_of_interest_hgnc_t[i,1]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
#sad_genes = sad_genes[!(sad_genes %in% liver_genes)]
length(sad_genes)

expr_mat = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw));colnames(expr_mat) = colnames(expr_raw);rownames(expr_mat) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
#expr_mat = expr_mat[,meta_data[meta_data$NEC_NET %in% "NEC","Sample"]]
expr = expr_mat
row_var = apply(expr, MARGIN = 1, FUN = var)
summary(row_var)
#expr = expr[row_var > 1,]

expr[1:5,1:5]
dim(expr)

###

correlation_matrix = cor(expr)

library(umap)

tissue_vec = meta_data$Subtype
tissue_vec[tissue_vec %in% c("Alpha","Beta","Gamma","Delta")] ="Endocrine"
tissue_vec[tissue_vec %in% c("Ductal","Acinar")] ="Exocrine"
tissue_vec[tissue_vec %in% c("endothelial","mesenchymal")] ="Others"

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
custom.config$random_state
#custom.config$random_state = 350
#custom.config$n_components=2

umap_result = umap::umap(
  correlation_matrix,
  colvec = meta_data$Subtype,
  preserve.seed = FALSE,
  config=custom.config
)

colnames(umap_result$layout) = c("x","y")
umap_result$layout = as.data.frame(umap_result$layout)

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point( aes( size = 4, color = as.factor(meta_data$Subtype) ))

umap_p = umap_p + xlab("") + ylab("")
umap_p = umap_p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
umap_p

umap_p = umap_p + stat_ellipse( linetype = 1,aes( color = tissue_vec),level=.66, type ="t",size=1.5)
#umap_p = umap_p + scale_color_manual( values = c("cyan","blue","yellow","purple","brown","black","brown","green")) ##33ACFF ##FF4C33
#umap_p = umap_p + scale_color_manual( values = c("cyan","blue","yellow","purple","brown","black","white","brown","green","gray","gray")) ##33ACFF ##FF4C33
umap_p = umap_p + scale_color_manual( values = c("brown","black","black","black","brown","black","white","brown","black","gray","gray")) ##33ACFF ##FF4C33
#umap_p = umap_p + scale_color_manual( values = c("Purple","black","black","black","purple","black","brown","purple","black","brown")) ##33ACFF ##FF4C33
umap_p

umap_p = umap_p + annotate("text", x = 2.5, y = 2.5, label = "NEC",col = "darkred",size =8)
umap_p = umap_p + annotate("text", x = -4.5, y = -2.5, label = "NET",col = "#33ACFF",size =8)
umap_p


####

