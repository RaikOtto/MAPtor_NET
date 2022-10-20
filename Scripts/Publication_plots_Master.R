library("stringr")
library("dplyr")
library("viridis")
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

library("grid")
meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F,comment.char = "!")

#meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F,comment.char = "!")
meta_info$Sample = str_replace(meta_info$Sample, pattern = "^(X\\.)", "")

rownames(meta_info) = meta_info$Sample

###


#expr_raw = read.table("~/MAPTor_NET/BAMs_new/Master/Master_new.S34.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1, comment.char = "!")
expr_raw = read.table("~/MAPTor_NET/BAMs_new/Master/Master_new.S34.HGNC.DESeq2.VOOM.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1, comment.char = "!")

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X\\.", "")
expr_raw[1:5,1:5]
no_match = (colnames(expr_raw) %in% meta_info$Sample) == F
colnames(expr_raw)[no_match] = paste("X",colnames(expr_raw)[no_match],sep ="")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "XX", "")
no_match = (colnames(expr_raw) %in% meta_info$Sample) == F
no_match

expr_raw = expr_raw[,str_detect(colnames(expr_raw),pattern = "_", negate = T)]
meta_data = meta_info[colnames(expr_raw),]
dim(expr_raw)

grep(rownames(expr_raw), pattern = "HOXB", value = T)

### Prep

meta_data = meta_info[ colnames(expr_raw), ]
table(meta_data$Grading)

##
source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")

###

i = 13#i = 13
genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
sad_genes = sad_genes[sad_genes != ""]

hox_genes = c("HOXA1","HOXA2","HOXA3","HOXA4","HOXA5","HOXA6","HOXA7","HOXA9","HOXA10","HOXA11","HOXA13","HOXB1","HOXB2","HOXB3","HOXB4","HOXB5","HOXB6","HOXB7","HOXB8","HOXB9","HOXB13","HOXC4","HOXC5","HOXC6","HOXC8","HOXC9","HOXC10","HOXC11","HOXC12","HOXC13","HOXD1","HOXD3","HOXD4","HOXD8","HOXD9","HOXD10","HOXD11","HOXD12","HOXD13")

sad_genes = hox_genes
#hox_genes = c("HOXB7","HOXB8","HOXB2","HOXB5","HOXD3","HOXD10","HOXD11","HOXC10","HOXD8","HOXD13","HOXD4","HOXB6","HOXD9","HOXA5","HOXB9","HOXB4","HOXB3","HOXA10")
#hox_genes = c("HOXA2","HOXA3","HOXA4","HOXA5","HOXA9","HOXA10","HOXB3","HOXB4","HOXB5","HOXB6")

liver_genes = genes_of_interest_hgnc_t[70,3:ncol(genes_of_interest_hgnc_t)]
sad_genes = sad_genes[!(sad_genes %in% liver_genes)]

genes_of_interest_hgnc_t[i,1]

sad_genes[which(!(sad_genes %in% rownames(expr_raw)))]
table(sad_genes %in% rownames(expr_raw) )

expr = expr_raw[ rownames(expr_raw) %in% sad_genes,]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))

p = ggbiplot::ggbiplot(
  pcr,
  labels = meta_data$Sample,
  var.axes = F,
  ellipse = FALSE,
  obs.scale = 1,
  var.scale = 1,
  labels.size = 1.5
)
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")
#p = p + xlim(-.75,2.25) + ylim(-.9,.8) # Fig 5
p

#svg("~/MAPTor_NET/Results/Publication/Figures/Figure_0_Validation_HOX.svg", width = 5, height = 5)
p
#dev.off()

## Figure 1 A

Grading_vec  = as.character(meta_data$Grading)
Grading_vec[Grading_vec %in% "G3"] = 3
Grading_vec[Grading_vec %in% "G2"] = 2
Grading_vec[Grading_vec %in% "G1"] = 1

hist_vec = meta_data$Primary_Metastasis %>% recode( "Primary" = "4") %>% recode( "Metastasis" = "6") %>% recode( "Local_recurrence" = "5") %>% as.double() # Fig 5

MKi67_vec = as.double(expr_raw["MKI67",colnames(expr)])
MKi67_vec = scale(MKi67_vec)

p = ggbiplot::ggbiplot(
  pcr,
  groups = as.character(meta_data$Study),
  var.axes = F,
  ellipse = TRUE,
  obs.scale = 1,
  var.scale = 1,
  size = 0.2
)

MKi67_vec = ((MKi67_vec + -1*min(MKi67_vec) ) +.5)

p = p + geom_point( aes( shape = meta_data$Grading, color = meta_data$Study ), size = MKi67_vec ) # Fig 4
p = p + scale_color_manual( values = c("#03444d") ) 
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")
p = p + xlim(-.92,2.22) + ylim(-1,.75) # Fig 5
p

svg("~/Downloads/Figure_1_A_Validation.svg", width = 10, height = 10)
p
dev.off()


# Figure 1 C

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
p = p + stat_ellipse(aes(color = meta_data$NET_NEC_PCA),level = .8, size = 1.25, linetype = 2)
p = p + scale_color_manual( values = c("white","red","blue","red","blue"))
#p = p + xlim(-.92,2.22) + ylim(-1,.75) # Fig 5
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")
p

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
p = p + stat_ellipse(aes(color = meta_data$NET_NEC_PCA),level = .8, size = 1.25, linetype = 2)
p = p + scale_color_manual( values = c("white","red","blue","red","blue"))
#p = p + xlim(-.92,2.22) + ylim(-1,.75) # Fig 5
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")
p

#svg("~/Downloads/Figure_1_C_Validation_PCA_hox_genes.svg", width = 5, height = 5)
p
dev.off()

## Supplementary D Heatmap

p = pheatmap::pheatmap(
  expr,
  #cor_mat,
  annotation_col = meta_data[c("NET_NEC_PCA","NEC_NET","MEN1","Study","Primary_Metastasis","Histology_Primary","Grading")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = F,
  treeheight_row = 0,
  legend = FALSE,
  annotation_legend = FALSE,
  fontsize_col = 7,
  clustering_method = "complete"
)

#svg("~/Downloads/Supplementary_D_Validation_Legend.svg", width = 10, height = 10)
#svg("~/Downloads/Supplementary_D_Validation_hox_genes.svg", width = 10, height = 10)
print(p)
dev.off()

### waterfall plot


library(ggplot2)

expr_raw_selection = expr_raw[,]
meta_data_selection = meta_data[,]

meta_data_selection$MEN1 = as.double(expr_raw_selection[ "MEN1",])
vis_mat = meta_data_selection[,c("Sample","MEN1","MEN1_Mut_AF","NEC_NET_PCA","Histology")]
dim(vis_mat)
colnames(vis_mat) = c("Sample","MEN1","MEN1_Mut_AF","NEC_NET","Histology")
vis_mat$Sample = factor(vis_mat$Sample, levels = vis_mat$Sample[order(as.double(vis_mat$MEN1))] )

label_vec = meta_data_selection$NEC_NET_PCA
label_vec = label_vec[order(meta_data_selection$MEN1)]

label_vec[label_vec == "NET"] = "T"
label_vec[label_vec != "T"] = "C"
col_vec = label_vec
col_vec[col_vec == "T"] = "blue"
col_vec[col_vec != "blue"] = "red"

vis_mat$MEN1 = vis_mat$MEN1 + -1*min(vis_mat$MEN1)

p = ggplot( data = vis_mat,aes( x = Sample, y = MEN1, fill = MEN1_Mut_AF ))
p = p + geom_col(position = position_dodge2( preserve = "single"),color = "black") 

p = p + geom_bar(stat="identity", position=position_dodge(), color = "black")
p = p + scale_fill_gradientn(colours = c("white","gray","black"), breaks = c(0.0,.5,1.0))

p = p + annotate("text", x=1:34,y = 8.5,parse=TRUE, label = label_vec, color = col_vec, size = 4.5 )
p = p + annotate("text", x=13,y = 7.5, label= "Wilcoxon-Smith test p-value = 0.047", size = 4.5 )

p = p+ theme_minimal() + xlab("") + ylab("") + theme(legend.position="top")
p = p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())# + ylim(0,6)
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + labs(fill = "Allel frequency mutated MEN1 transcripts")
p = p + expand_limits(x = c(0,57)) 
p

svg("~/MAPTor_NET/Results/Publication/Figures/Figure_K_Discovery.svg", width = 8, height = 5)
p 
dev.off()

rank_nec = which(label_vec == "C")
rank_net = which(label_vec == "T")
wilcox.test(rank_nec,rank_net)
