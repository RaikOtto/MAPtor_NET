library("stringr")
library("umap")
library("dplyr")
library("viridis")
library("grid")
library("ggplot2")

expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Fröhling.S34.HGNC.DESeq2.VOOM.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X\\.)", "")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X)", "")
expr_raw[1:5,1:5]
dim(expr_raw)

expr_raw_mki67 = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Fröhling.S34.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
colnames(expr_raw_mki67) = str_replace(colnames(expr_raw), pattern = "^(X\\.)", "")
colnames(expr_raw_mki67) = str_replace(colnames(expr_raw), pattern = "^(X)", "")
mki67_vec = as.double(expr_raw["MKI67",])

meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)

meta_info$Name = str_replace(meta_info$Sample, pattern = "^(X\\.)", "")
rownames(meta_info) = meta_info$Sample

no_match = match(colnames(expr_raw), meta_info$Sample, nomatch = F) == F
colnames(expr_raw)[no_match  ]
colnames(expr_raw)[no_match  ] = paste("X",colnames(expr_raw)[no_match  ],sep = "")
no_match = match(colnames(expr_raw), meta_info$Sample, nomatch = F) == F
meta_data = meta_info[ colnames(expr_raw), ]
no_match
dim(meta_data)
table(meta_data$Primary_Metastasis)

rownames(meta_data) = meta_data$Sample
expr_raw = expr_raw[,meta_data$Sample[meta_data$Grading != ""]]

dim(expr_raw)
table(meta_data$Grading)

##
source("~/MAPTor_NET/Misc//Visualization_colors.R")

###

i = 13
#genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/SeneSys_gene_sets.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t = read.table("~/MAPTor_NET//Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
genes_of_interest_hgnc_t$V1[i]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
sad_genes = sad_genes[sad_genes != ""]

liver_genes = genes_of_interest_hgnc_t[70,3:ncol(genes_of_interest_hgnc_t)]
sad_genes = sad_genes[!(sad_genes %in% liver_genes)]

hox_genes = c("HOXA1","HOXA2","HOXA3","HOXA4","HOXA5","HOXA6","HOXA7","HOXA9","HOXA10","HOXA11","HOXA13","HOXB1","HOXB2","HOXB3","HOXB4","HOXB5","HOXB6","HOXB7","HOXB8","HOXB9","HOXB13","HOXC4","HOXC5","HOXC6","HOXC8","HOXC9","HOXC10","HOXC11","HOXC12","HOXC13","HOXD1","HOXD3","HOXD4","HOXD8","HOXD9","HOXD10","HOXD11","HOXD12","HOXD13")
#sad_genes = hox_genes

genes_of_interest_hgnc_t[i,1]

sad_genes[which(!(sad_genes %in% rownames(expr_raw)))]
table(sad_genes %in% rownames(expr_raw) )

expr = expr_raw[ rownames(expr_raw) %in% sad_genes,]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))
#pcr = prcomp(t(expr))
meta_data = meta_info[ colnames(expr), ]

## initiate umap

custom.config = umap.defaults
#custom.config$random_state = sample(1:1000,size = 1)
custom.config$random_state = 748
custom.config$n_components=2

umap_result = umap::umap(
  cor_mat,
  colvec = meta_data["Study"],
  preserve.seed = TRUE,
  config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")
custom.config$random_state # 748

# plot 2 gradings

mki67_size_vec = mki67_vec + -1*min(mki67_vec)
mki67_size_vec[32] = 3
#mki67_size_vec = mki67_size_vec-min(mki67_size_vec) + 1
#mki67_size_vec = mki67_size_vec / max(mki67_size_vec)
#mki67_size_vec = mki67_size_vec*1.4

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point( aes( size = mki67_size_vec, color = as.character(meta_data$Grading) ))
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$Grading), level=.5, type ="t", size=1.5)
umap_p = umap_p + scale_color_manual( values = c("#CCCCCC","#999999","#333333","gray")) 

umap_p = umap_p + theme(legend.position = "top") + xlab("") + ylab("")
#umap_p = umap_p + geom_vline( xintercept=-1, size = 2, linetype = 2) + geom_hline( yintercept = -1.25, size = 2, linetype = 2)  
umap_p = umap_p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())

#svg("~/Downloads/MAPTor-NET_plots_20_10_2022/Figure_2C_UMAP_grading_validation.svg", width = 10, height = 10)
umap_p + theme_minimal()
dev.off()

# plot 3 nec net men1


men1_colors = meta_data$MEN1
men1_colors[ (meta_data$NET_NEC_UMAP == "NET") & (meta_data$MEN1 == 0) ] = "blue"
men1_colors[ (meta_data$NET_NEC_UMAP == "NEC") & (meta_data$MEN1 == 0) ] = "red"
men1_colors[ !(men1_colors %in% c("red","blue")) ] = "white"

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point( aes( size =4, color = as.character(meta_data$NET_NEC_UMAP) ))
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$NET_NEC_UMAP), level=.5, type ="t", size=1.5)
umap_p = umap_p + geom_point( aes( size = 1, color = men1_colors )) # men 1
umap_p = umap_p + scale_color_manual( values = c("blue","red","blue","red","white"))

umap_p = umap_p + theme(legend.position = "top") + xlab("") + ylab("")
umap_p = umap_p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())

#svg("~/Downloads/MAPTor-NET_plots_20_10_2022/Figure_2_D_UMAP_net_nec_men1_validation_umap.svg", width = 10, height = 10)
umap_p + theme_minimal()
dev.off()

# plot 5 tissue of origin

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point( aes( size = 4, color = as.character(meta_data$Histology_Primary) ))
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$Histology_Primary), level=.5, type ="t", size=1.5)
umap_p = umap_p + scale_color_manual( values = c("#440154FF","#73D055FF","#FDE725FF","#1F968BFF","#39568CFF")) 

umap_p = umap_p + theme(legend.position = "top") + xlab("") + ylab("")
umap_p = umap_p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())

svg("~/Downloads/MAPTor-Net_plots_22.02.2022/3_UMAP_histology_primary.svg", width = 10, height = 10)
umap_p + theme_minimal()
dev.off()

# plot 5 primary metastasis

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point( aes( size =4, color = as.character(meta_data$Primary_Metastasis) ))
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$Primary_Metastasis), level=.5, type ="t", size=1.5)
umap_p = umap_p + scale_color_manual( values = c("orange","#333333","#CCCCCC")) 

umap_p = umap_p + theme(legend.position = "top") + xlab("") + ylab("")
umap_p = umap_p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())

#svg("~/Downloads/MAPTor-Net_plots_22.02.2022//3_UMAP_primary_metastasis.svg", width = 10, height = 10)
umap_p + theme_minimal()
dev.off()

### label plot supplement validation

umap_p = ggplot(
  umap_result$layout,
  aes(x = x, y = y, color = meta_data$Study))
umap_p = umap_p + geom_point( aes( size = 4, color = as.character(meta_data$Study) ))
#umap_p = umap_p + geom_label(aes(label = meta_data$Sample, size = NULL), nudge_y = 0.0)
umap_p = umap_p + scale_color_manual( values = c("#66FF00","#006600")) ##33ACFF ##FF4C33
umap_p = umap_p + theme(legend.position = "top") + xlab("") + ylab("")
umap_p = umap_p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
umap_p = umap_p  + geom_label_repel(
  aes(label = meta_data$Sample, size = NULL, color = meta_data$Study),
  arrow = arrow(length = unit(0.03, "npc"),
                type = "closed", ends = "last"),
  nudge_y = 1,
  segment.size  = 0.3)

svg("~/Downloads/MAPTor-NET_plots_20_10_2022/3_SM_Figure_Umap_Labels_Discovey.svg", width = 10, height = 10)
umap_p + theme_minimal()
dev.off()