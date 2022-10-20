library("stringr")
library("umap")
library("dplyr")
library("ggplot2")

#expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Discovery_Cohort.S64.DESeq2.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Fröhling.S34.HGNC.DESeq2.VOOM.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
#expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Fröhling.S34.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)

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
no_match
dim(meta_data)
table(meta_data$Primary_Metastasis)

rownames(meta_data) = meta_data$Sample

#expr_raw = expr_raw[ ,meta_data$Histology_Primary == "Pancreatic"]
#meta_data = meta_info[colnames(expr_raw),]
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
sad_genes = hox_genes

genes_of_interest_hgnc_t[i,1]

sad_genes[which(!(sad_genes %in% rownames(expr_raw)))]
table(sad_genes %in% rownames(expr_raw) )

expr = expr_raw[ rownames(expr_raw) %in% sad_genes,]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))
#pcr = prcomp(t(expr))
meta_data = meta_info[ colnames(expr), ]
meta_data$MEN1_mt = meta_data$MEN1
meta_data$MEN1_expr = log(as.double(expr_raw["MEN1",]))

## Supplementary D Heatmap

p = pheatmap::pheatmap(
  #(expr),
  cor_mat,
  annotation_col = meta_data[c("MEN1_expr","MEN1_mt","NET_NEC_PCA","NEC_NET","Grading")],
  annotation_colors = aka3,
  show_rownames = FALSE,
  show_colnames = F,
  treeheight_row = 0,
  legend = FALSE,
  annotation_legend = TRUE,
  #fontsize_col = 7,
  clustering_method = "ward.D"
)

#svg("~/Downloads/MAPTor-NET_plots_23_02_2022/4_Heatmap_discovery_39hox.svg", width = 10, height = 10)
#svg("~/Downloads/MAPTor-NET_plots_23_02_2022/4_Heatmap_discovery_39hox_pancreas_only.svg", width = 10, height = 10)
#svg("~/Downloads/MAPTor-NET_plots_23_02_2022/5_Heatmap_validation_39hox.svg", width = 10, height = 10)
print(p)
dev.off()

### labels

p = ggbiplot::ggbiplot(
  pcr,
  labels = meta_data$Sample,
  var.axes = F,
  ellipse = TRUE,
  obs.scale = 1,
  var.scale = 1,
  labels.size = 2
)
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")
#p = p + xlim(-1.9,1.3) + ylim(-.85,1.2) # Fig 5

#svg("~/Downloads/MAPtor_NET_plots_19_10_2021/Figure_0_discovery.svg", width = 10, height = 10)
p
dev.off()

## Figure 1 A

Grading_vec  = as.character(meta_data$Grading)
Grading_vec[Grading_vec %in% "G3"] = 3
Grading_vec[Grading_vec %in% "G2"] = 2
Grading_vec[Grading_vec %in% "G1"] = 1

hist_vec = meta_data$Primary_Metastasis %>% recode( "Primary" = "4") %>% recode( "Metastasis" = "6") %>% recode( "Unknown" = "5") %>% as.double() # Fig 5

MKi67_vec = as.double(expr_raw["MKI67",colnames(expr)])

MKi67_vec = ((MKi67_vec + -1*min(MKi67_vec) ) +.5)

p = ggbiplot::ggbiplot(
    pcr,
    groups = as.character(meta_data$Study),
    var.axes = F,
    ellipse = FALSE,
    obs.scale = 1,
    var.scale = 1,
    size = 0.2
    )
p = p + geom_point( aes( shape = meta_data$Grading, color = meta_data$Study ), size = MKi67_vec ) # Fig 4
p = p + scale_color_manual( values = c("#35A047","#F28500") ) #Fig 4 Master
p = p + stat_ellipse(aes(color = meta_data$Study),level = .72, size = 1.25, linetype = 2)
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")
p = p + xlim(-1.9,1.3) + ylim(-.85,1.2) # Fig 5
p

#svg("~/Downloads/MAPtor_NET_plots_19_10_2021/Figure_1_A_Discovery.svg", width = 10, height = 10)
#svg("~/Downloads/MAPtor_NET_plots_19_10_2021/Figure_1_A_Discovery_legend.svg", width = 10, height = 10)
p
dev.off()

## Figure 1 B

metastasis_vec = meta_data$Primary_Metastasis
metastasis_vec[metastasis_vec == "Primary"] = "1"
metastasis_vec[metastasis_vec != "1"] = "4"

p = ggbiplot::ggbiplot(
  pcr,
  groups = as.character(meta_data$Histology_Primary),
  ellipse = FALSE,
  var.axes = F,
  obs.scale = 1,
  var.scale = 1,)

p = p + geom_point( aes(colour= meta_data$Histology ), size = as.integer(metastasis_vec) )
#p = p + scale_y_reverse()
p = p + scale_color_manual( values = c("darkred","green","cyan","blue","black","black","black") ) 
p = p + geom_point( aes(colour= meta_data$Histology_Primary ), size = as.integer(metastasis_vec) )
#p = p + scale_y_reverse()
p = p + scale_color_manual( values = c("brown","purple","orange","black") ) 
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")
p = p + xlim(-1.9,1.3) + ylim(-.85,1.2) # Fig 5
p

#svg("~/Downloads/MAPtor_NET_plots_19_10_2021/Figure_1_B_Discovery.svg", width = 10, height = 10)
#svg("~/Downloads/MAPtor_NET_plots_19_10_2021/Figure_1_B_Discovery_legend.svg", width = 10, height = 10)
p
dev.off()

# Figure 1 C andere meta data datei

nec_net_col_vec = meta_data$NEC_NET;
nec_net_col_vec[nec_net_col_vec == "NEC"] = "red";
nec_net_col_vec[nec_net_col_vec == "NET"] = "blue"
nec_net_col_vec[nec_net_col_vec == "Ambiguous"] = "purple"

MEN1_status_col = meta_data$MEN1
MEN1_status_col[MEN1_status_col <= 0] = "MEN1_wt"
MEN1_status_col[MEN1_status_col != "MEN1_wt"] = "MEN1_mt"
MEN1_status_col[ ( nec_net_col_vec  == "blue") & (MEN1_status_col == "MEN1_wt") ] = "MEN1_wt_net"
MEN1_status_col[ ( nec_net_col_vec  == "blue") & (MEN1_status_col == "MEN1_mt") ] = "MEN1_mt_net"
MEN1_status_col[ ( nec_net_col_vec  == "red") & (MEN1_status_col == "MEN1_wt") ] = "MEN1_wt_nec"
MEN1_status_col[ ( nec_net_col_vec  == "red") & (MEN1_status_col == "MEN1_mt") ] = "MEN1_mt_nec"
MEN1_status_col[ ( nec_net_col_vec  == "purple") & (MEN1_status_col == "MEN1_wt") ] = "MEN1_wt_ambiguous"
MEN1_status_col[ ( nec_net_col_vec  == "purple") & (MEN1_status_col == "MEN1_mt") ] = "MEN1_mt_ambiguous"

size_vec = nec_net_col_vec
size_vec[size_vec != "blue"] = "4"
size_vec[size_vec != "4" ] = "1"

shape_vec = meta_data$Grading
shape_vec[shape_vec == "G1"] = 18
shape_vec[shape_vec == "G2"] = 17
shape_vec[shape_vec == "G3"] = 15
shape_vec = as.integer(shape_vec)

p = ggbiplot::ggbiplot(
  pcr,
  groups = as.character(meta_data$NET_NEC_PCA),
  ellipse = FALSE,
  var.axes = F,
  obs.scale = 1,
  var.scale = 1,)

p = p + geom_point( shape = 15 , size = as.integer(size_vec),  color = "black", stroke = .9)
p = p + geom_point( aes( color = as.factor( meta_data$NEC_NET )), shape = shape_vec, size = 3 )
#p = p + scale_y_reverse()
p = p + geom_point( aes ( color = as.factor( MEN1_status_col ) ), size = 1 )
p = p + stat_ellipse(aes(color = meta_data$NET_NEC_PCA),level = .71, size = 1.25, linetype = 2)
p = p + scale_color_manual( values = c("purple","white","white","white","purple","red","blue","red","blue"))
p = p + xlim(-1.9,1.3) + ylim(-.85,1.2) # Fig 5
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")
p

#svg("~/Downloads/MAPtor_NET_plots_19_10_2021/Figure_1_C_Discovery.svg", width = 10, height = 10)
#svg("~/Downloads/MAPtor_NET_plots_19_10_2021/Figure_1_C_Discovery_legend.svg", width = 10, height = 10)
p
dev.off()


### HOX GENE PLOT

p = ggbiplot::ggbiplot(
  pcr,
  groups = as.character(meta_data$NET_NEC_PCA),
  ellipse = FALSE,
  var.axes = F,
  obs.scale = 1,
  var.scale = 1,)

p = p + geom_point( shape = 15 , size = as.integer(size_vec),  color = "black", stroke = .9)
p = p + geom_point( aes( color = as.factor( meta_data$NEC_NET )), shape = shape_vec, size = 3 )
#p = p + scale_y_reverse()
p = p + geom_point( aes ( color = as.factor( MEN1_status_col ) ), size = 1 )
p = p + stat_ellipse(aes(color = meta_data$NET_NEC_PCA),level = .7, size = 1.25, linetype = 2)
p = p + scale_color_manual( values = c("purple","white","white","white","purple","red","blue","red","blue"))
#p = p + xlim(-1.9,1.3) + ylim(-0.75,1.1) # Fig 5
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")
p

#svg("~/Downloads/MAPtor_NET_plots_19_10_2021/Figure_1_C_Discovery_Hox_genes_legend.svg", width = 10, height = 10)
#svg("~/Downloads/MAPtor_NET_plots_19_10_2021/Figure_1_C_Discovery_Hox_genes.svg", width = 10, height = 10)
p
dev.off()

###

# Supplementary Figure 1
meta_data$Location[meta_data$Location == "Unknown"] = "missing"
p = pheatmap::pheatmap(
  expr,
  annotation_col = meta_data[c("NET_NEC_PCA","NEC_NET","Grading","Location","Histology","Study")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  treeheight_col = 0,
  treeheight_row = 0,
  legend = F,
  fontsize_col = 5,
  clustering_method = "ward.D"
)
p

#svg("~/MAPTor_NET/Results/Publication/Figures/Plot_D_Discovery.svg", width = 5, height = 5)
    print(p)
dev.off()



rank_nec = which(label_vec == "C")
rank_net = which(label_vec == "T")
wilcox.test(rank_nec,rank_net)

##### UMAP

### Figure 3 Plot E UMAP PLOT

sad_genes = hox_genes = c("HOXA1","HOXA2","HOXA3","HOXA4","HOXA5","HOXA6","HOXA7","HOXA9","HOXA10","HOXA11","HOXA13","HOXB1","HOXB2","HOXB3","HOXB4","HOXB5","HOXB6","HOXB7","HOXB8","HOXB9","HOXB13","HOXC4","HOXC5","HOXC6","HOXC8","HOXC9","HOXC10","HOXC11","HOXC12","HOXC13","HOXD1","HOXD3","HOXD4","HOXD8","HOXD9","HOXD10","HOXD11","HOXD12","HOXD13")
props = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/RepSet.S57.DESeq2.tsv",sep ="\t", header = T, as.is=TRUE,row.names = 1)
props[1:5,1:5]
colnames(props) = colnames(props) %>% str_replace("^X","")

genes_of_interest_hgnc_t = read.table("~/MAPTor_NET//Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[13,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
sad_genes = sad_genes[sad_genes != ""]
liver_genes = genes_of_interest_hgnc_t[70,3:ncol(genes_of_interest_hgnc_t)]
sad_genes = sad_genes[!(sad_genes %in% liver_genes)]
sad_genes = hox_genes
expr = expr_raw[match(sad_genes, rownames(expr_raw), nomatch = 0),]
correlation_matrix = cor(expr)
meta_data= meta_info[colnames(correlation_matrix),]

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
#custom.config$random_state = 995
custom.config$n_components=2

umap_result = umap::umap(
  correlation_matrix,
  colvec = meta_data["NET_NEC_PCA"],
  preserve.seed = TRUE,
  config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point( aes( size = 4, color = as.character(meta_data$NET_NEC_PCA) ))
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$NET_NEC_PCA), level=.5, type ="t", size=1.5)
umap_p = umap_p + scale_color_manual( values = c("darkred","blue")) ##33ACFF ##FF4C33

umap_p = umap_p + theme(legend.position = "top") + xlab("") + ylab("")
#umap_p = umap_p + geom_vline( xintercept=-1, size = 2, linetype = 2) + geom_hline( yintercept = -1.25, size = 2, linetype = 2)  
umap_p = umap_p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())

#umap_p = umap_p + annotate("text", x = 3.5, y = -3.5, label = "NEC",col = "darkred",size =8)
#umap_p = umap_p + annotate("text", x = -1, y = 2.5, label = "NET",col = "#33ACFF",size =8)
umap_p
custom.config$random_state # 523
