library("umap")
library("ggpubr")
library("stringr")
library("reshape2")
library("ggplot2")
library("dplyr")
library("grid")

meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

expr_raw = read.table("~/MAPTor_NET/BAMs_new/RepSet_S52.PanNEN.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
#expr_raw = read.table("~/Downloads/GSE119262_non-normalized.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

expr_raw[1:5,1:5]
dim(expr_raw)
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match == TRUE] = str_replace(colnames(expr_raw)[no_match == TRUE], pattern = "^X","")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match == TRUE] = paste("X",colnames(expr_raw)[no_match == TRUE],sep ="")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[which(no_match == TRUE)]
meta_data = meta_info[colnames(expr_raw),]
table(meta_data$Study)

expr_raw = expr_raw[,meta_data$Histology_Primary == "Pancreatic" ]
meta_data = meta_info[colnames(expr_raw),]

expr_raw = expr_raw[,meta_data$NET_NEC_PCA == "NET" ]
meta_data = meta_info[colnames(expr_raw),]
dim(expr_raw)

genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t[,1]
liver_genes = genes_of_interest_hgnc_t[70,3:ncol(genes_of_interest_hgnc_t)]
i = 64
genes_of_interest_hgnc_t[i,1]

sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
sad_genes = sad_genes[!(sad_genes %in% liver_genes)]
length(sad_genes)

expr = as.data.frame(expr_raw)[rownames(expr_raw) %in% sad_genes,]
expr[1:5,1:5]
dim(expr)

#WriteXLS::WriteXLS(expr,"~/Downloads/Everolimus_metabolomic_methylomic_analyses/2_NET_Left_vs_NET_Right/S52.panNETs.NETs_only.no_study_normalization.absolute_counts.EVO_UPREGULATED.xlsx",row.names = TRUE)

#expr = expr[,!(colnames(expr) %in% c("GSM3362429","GSM3362430","GSM3362461","GSM3362462","GSM3362446") )]
meta_data = meta_info[colnames(expr),]
###

correlation_matrix = cor((expr))
pcr = prcomp(t(correlation_matrix))

dim(expr)

#svg(filename = "~/Downloads/Heatmap.svg", width = 10, height = 10)
source("~/MAPTor_NET//Scripts/Archive/Visualization_colors_thedieck.R")

CDH5_vec = CDH5_vec_ori = (as.double(expr_raw["CDH5",colnames(expr_raw) %in% colnames(expr)]))
low_threshold = quantile(CDH5_vec_ori, seq(0,1,0.01))[26]
medium_threshold = quantile(CDH5_vec_ori, seq(0,1,0.01))[76]
CDH5_vec[CDH5_vec_ori > medium_threshold ] = "High"
CDH5_vec[CDH5_vec_ori <= medium_threshold ] = "Medium"
CDH5_vec[CDH5_vec_ori <= low_threshold ] = "Low"
meta_data$CDH5 = CDH5_vec

#pdf("~/Downloads/Reproductions_21_12_2021/Figure_4_Labels.pancreas_only.pdf")
p  = pheatmap::pheatmap(
  correlation_matrix,
  #expr_raw[c("H2BP4","H3P10","SLC22A6"),],
  annotation_col = meta_data[,c("pannet_cluster","NET_NEC_PCA","Grading","Study")],
  annotation_colors = aka3,
  show_rownames = FALSE,
  show_colnames = TRUE,
  treeheight_row = 0,
  legend = T,
  fontsize_col = 5,
  clustering_method = "ward.D"
)
dev.off()

pannet_cluster = p$tree_col$labels[ p$tree_col$order]
pannet_cluster[1:which(pannet_cluster == "440")] = "Left"
pannet_cluster[pannet_cluster != "Left"] = "Right"
meta_data[p$tree_col$labels[ p$tree_col$order],"pannet_cluster"]=pannet_cluster

expr_export = expr_raw[,p$tree_col$labels[ p$tree_col$order]]
expr_export = expr_raw[,order_vec]
colnames(expr_export) = pannet_cluster

#write.table(expr_export,"~/Downloads/2_NET_Left_vs_NET_Right//S52.panNETs.NETs_only.no_study_normalization.absolute_counts.tsv",quote =F , sep ="\t", row.names = TRUE)

####

library("umap")
umap_result = umap::umap(
  correlation_matrix,
  #colvec = col_vec_nec_net,
  preserve.seed = TRUE#,
  #config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point( aes( size = 4, color = as.character(meta_data$PAKT_response) ))
umap_p
#write.table(meta_info,"~/MAPTor_NET/Misc/Meta_information.tsv",sep ="\t", quote =F ,row.names = F)

#######

vis_mat = reshape2::melt(as.matrix(expr))
new_meta_mat = meta_data[as.character(vis_mat[,2]),]
vis_mat$Cohort = new_meta_mat$Cohort
colnames(vis_mat) = c("Gene","Sample","Expression","Left_Right")

gene_differences = c()
for (gene in unique(vis_mat$Gene)){
  exp_mat = vis_mat[vis_mat$Gene == gene,"Expression"]
  exp_groups = aggregate(exp_mat, by = list(meta_data[,"VM_signature_new"]),FUN = mean)
  difference = abs(exp_groups$x[1] - exp_groups$x[2])
  gene_differences = c(gene_differences, difference)
}
names(gene_differences) = unique(vis_mat$Gene)
order_list = order(gene_differences,decreasing = T)



vis_mat$Gene = factor(
  as.character(vis_mat$Gene),
  levels= unique(as.character(vis_mat$Gene)[  order_list ]))

genes = names(table(vis_mat$Gene)); length(genes)


if (min(vis_mat$Expression) < 0)
  vis_mat$Expression = vis_mat$Expression + -1 * min(vis_mat$Expression)

p = ggplot( data = vis_mat,aes( x = Gene, y = Expression, fill = as.factor(Left_Right) ))
p = p +  geom_boxplot( )
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + theme(legend.position="top")+ xlab("") + ylab("")

p = p + scale_fill_manual(values=c("red","blue"))

#pdf("~/Downloads/Data_analysis_everolimus_04_10_2021/C_GSEA_pancreatic_NET_only_marker_gene_cohorts/13_Expression_Marker_Genes_Cohorts.pdf")
p
dev.off()

p = p + ylim(0,y_lim_val)
p = p + geom_hline(yintercept = 3.3, linetype= "dashed")

###

row_var = apply(expr, FUN = var, MARGIN = 1)
vis_mat = reshape2::melt(row_var)
vis_mat$Gene = rownames(vis_mat)
colnames(vis_mat) = c("Variance","Gene") 
vis_mat$Gene = factor(vis_mat$Gene, levels = c(vis_mat$Gene[order(vis_mat$Variance, decreasing = T)]))

vis_mat = reshape2::melt(t(expr))
colnames(vis_mat) = c("Sample","Gene","Expression")
order_vec = aggregate(vis_mat$Expression, by = list(vis_mat$Gene), FUN = var)
vis_mat$Cohort = meta_data[vis_mat$Sample,"Cohort"]

vis_mat$Gene = factor(vis_mat$Gene, levels = c(as.character(unique(vis_mat$Gene))[order(order_vec$x, decreasing = T)]))
selected_genes = order(order_vec,decreasing = TRUE)[1:20]

p = ggplot( data = vis_mat[vis_mat$Gene %in% order_vec$Group.1[1:20],],aes( x = Gene, y = Expression, fill = Cohort ))
p = p + geom_boxplot()
#p = p +  geom_bar( stat="identity")
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + theme(legend.position="top")+ xlab("") + ylab("")
p

####

row_var = rownames(expr)[order(apply(expr, MARGIN = 1, FUN = var), decreasing = TRUE)]
