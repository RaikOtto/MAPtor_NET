library("ggplot2")
library("stringr")
library("umap")
library("dplyr")

expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Discovery_Cohort.S64.DESeq2.vsd.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)

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
expr_raw = expr_raw[,meta_data$Sample[meta_data$Grading != ""]]

#expr_raw = expr_raw[ ,meta_data$Histology == "Pancreatic"]
#meta_data = meta_info[colnames(expr_raw),]
dim(expr_raw)

table(meta_data$Grading)

##
source("~/MAPTor_NET/Misc//Visualization_colors.R")

### waterfall plot

expr_raw_selection = expr_raw[,]
meta_data_selection = meta_data[,]

expr_raw_selection = expr_raw[,meta_data$Histology == "Pancreatic"]
meta_data_selection = meta_data[meta_data$Histology == "Pancreatic",]
dim(expr_raw_selection)

meta_data_selection$MEN1 = as.double(expr_raw_selection[ "MEN1",])
vis_mat = meta_data_selection[,c("Sample","MEN1","MEN1_Mut_AF","NET_NEC_UMAP","Histology_Primary")]
dim(vis_mat)
vis_mat$Sample = factor(vis_mat$Sample, levels = vis_mat$Sample[order(as.double(vis_mat$MEN1))] )

label_vec = meta_data_selection$NET_NEC_UMAP
label_vec = label_vec[order(meta_data_selection$MEN1)]

label_vec[label_vec == "NET"] = "T"
label_vec[label_vec != "T"] = "C"
col_vec = label_vec
col_vec[col_vec == "T"] = "blue"
col_vec[col_vec != "blue"] = "red"

#vis_mat$MEN1 = vis_mat$MEN1 + -1*min(vis_mat$MEN1)

p = ggplot( data = vis_mat,aes( x = Sample, y = MEN1, fill = MEN1_Mut_AF ))
p = p + geom_col(position = position_dodge2( preserve = "single"),color = "black") 

p = p + geom_bar(stat="identity", position=position_dodge(), color = "black")
p = p + scale_fill_gradientn(colours = c("white","gray","black"), breaks = c(0.0,.5,1.0))

p = p + annotate("text", x=1:ncol(expr_raw_selection),y = 11.25,parse=TRUE, label = label_vec, color = col_vec, size = 4.5 )
p = p + annotate("text", x=15,y = 9.75, label= "Wilcoxon-Smith test p-value = 1E-2", size = 4.5 )

#p = p + annotate("text", x=1:ncol(expr_raw_selection),y = 11.25,parse=TRUE, label = label_vec, color = col_vec, size = 4.5 )
#p = p + annotate("text", x=15,y = 9.75, label= "Wilcoxon-Smith test p-value = 3E-6", size = 4.5 )

p = p+ theme_minimal() + xlab("") + ylab("") + theme(legend.position="top")
p = p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())# + ylim(0,6)
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + labs(fill = "Allel frequency mutated MEN1 transcripts")
p = p + expand_limits(x = c(0,ncol(expr_raw_selection))) 

#svg("~/Downloads//Figure_5_Discovery.svg", width = 8, height = 5)
p 
dev.off()

###

NETs = which(label_vec == "T")
NECs = which(label_vec == "C")
wilcox.test(
  x = NETs,
  y = NECs
)

####

expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Fröhling.S34.HGNC.DESeq2.VOOM.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)

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
expr_raw = expr_raw[,meta_data$Sample[meta_data$Grading != ""]]

#expr_raw = expr_raw[ ,meta_data$Histology == "Pancreatic"]
#meta_data = meta_info[colnames(expr_raw),]
dim(expr_raw)

table(meta_data$Grading)

##
source("~/MAPTor_NET/Misc//Visualization_colors.R")

### waterfall plot

expr_raw_selection = expr_raw[,]
meta_data_selection = meta_data[,]

#expr_raw_selection = expr_raw[,meta_data$Histology == "Pancreatic"]
#meta_data_selection = meta_data[meta_data$Histology == "Pancreatic",]
dim(expr_raw_selection)

meta_data_selection$MEN1 = (as.double(expr_raw_selection[ "MEN1",]))
meta_data_selection$MEN1 = meta_data_selection$MEN1 + -1*min(meta_data_selection$MEN1) + .5
vis_mat = meta_data_selection[,c("Sample","MEN1","MEN1_Mut_AF","NET_NEC_UMAP","Histology_Primary")]
dim(vis_mat)
vis_mat$Sample = factor(vis_mat$Sample, levels = vis_mat$Sample[order(as.double(vis_mat$MEN1))] )

label_vec = meta_data_selection$NET_NEC_UMAP
label_vec = label_vec[order(meta_data_selection$MEN1)]

label_vec[label_vec == "NET"] = "T"
label_vec[label_vec != "T"] = "C"
col_vec = label_vec
col_vec[col_vec == "T"] = "blue"
col_vec[col_vec != "blue"] = "red"

#vis_mat$MEN1 = vis_mat$MEN1 + -1*min(vis_mat$MEN1)

p = ggplot( data = vis_mat,aes( x = Sample, y = MEN1, fill = MEN1_Mut_AF ))
p = p + geom_col(position = position_dodge2( preserve = "single"),color = "black") 

p = p + geom_bar(stat="identity", position=position_dodge(), color = "black")
p = p + scale_fill_gradientn(colours = c("white","gray","black"), breaks = c(0.0,.5,1.0))

p = p + annotate("text", x=1:ncol(expr_raw_selection),y = 11.25,parse=TRUE, label = label_vec, color = col_vec, size = 4.5 )
p = p + annotate("text", x=15,y = 9.75, label= "Wilcoxon-Smith test p-value = 0.047", size = 4.5 )

p = p+ theme_minimal() + xlab("") + ylab("") + theme(legend.position="top")
p = p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())# + ylim(0,6)
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + labs(fill = "Allel frequency mutated MEN1 transcripts")
p = p + expand_limits(x = c(0,ncol(expr_raw_selection))) 

svg("~/Downloads/Figure_5_Validation.svg", width = 8, height = 5)
p 
dev.off()

###

NETs = which(label_vec == "T")
NECs = which(label_vec == "C")
wilcox.test(
  x = NETs,
  y = NECs
)

