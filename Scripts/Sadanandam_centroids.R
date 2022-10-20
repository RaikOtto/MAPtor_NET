### create centroids

#t_a = read.table("~/MAPTor_NET/Misc/Sadanandam_training_samples_A.tsv", sep ="\t", header =T, stringsAsFactors = F)
t_b = read.table("~/MAPTor_NET/Misc/Sadanandam_training_samples_B.tsv", sep ="\t", header =T, stringsAsFactors = F)

#int_genes = intersect( t_a[,1], t_b[,1])

#m_a = match(int_genes, t_a[,1])
#m_b = match(int_genes, t_b[,1])
#new_mat = cbind( t_a[m_a, 2:ncol(t_a)], t_b[m_b, 2:ncol(t_b)] )
new_mat = t_b[,2:ncol(t_b)]
#rownames(new_mat) = int_genes
rownames(new_mat) = t_b[,1]
new_mat[1:5,1:5]
dim(new_mat)

subtypes = colnames(new_mat)[2:ncol(new_mat)]
subtypes[str_detect(subtypes, pattern = "MLP\\.1")] = "MLP_1"
subtypes[str_detect(subtypes, pattern = "MLP\\.2")] = "MLP_2"
subtypes[str_detect(subtypes, pattern = "Intermediate")] = "Intermediate"
subtypes[str_detect(subtypes, pattern = "Insulinoma")] = "Insulinoma_like"
subtypes[str_detect(subtypes, pattern = "Normal")] = "Normal_islet_like"
subtypes_uni = unique(subtypes)

meta_data = data.frame(
    "Subtype" = subtypes
)
rownames(meta_data) = 1:(ncol(new_mat)-1)
colnames(new_mat)[2:ncol(new_mat)] = 1:(ncol(new_mat)-1)
new_mat = new_mat[,-1]

cor_mat = cor(new_mat)
pheatmap::pheatmap(
  cor_mat,
  show_rownames = F,
  show_colnames = F,
  annotation_col = meta_data
)

pca = prcomp(cor_mat)
ggbiplot::ggbiplot(
  pca,
  groups = as.character(meta_data$Subtype),
  color = meta_data$Subtype,
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F
)

exp_centroid <<- matrix(as.double(), nrow = nrow(new_mat))
for (subtype in subtypes_uni){
    data = new_mat[, subtypes == subtype ]
    exp_centroid <<- cbind( exp_centroid, rowMeans(data))
}

rownames(exp_centroid) = rownames(new_mat)
colnames(exp_centroid) = subtypes_uni

#write.table(exp_centroid,"~/MAPTor_NET/Misc/Sadanandam_gene_expr_subtypes.csv", sep ="\t", quote =F)

library("stringr")

balanced.centroid = read.table( "~/MAPTor_NET/Misc/Sadanandam_gene_expr_subtypes.csv", header=TRUE, row.names=1, sep="\t",stringsAsFactors = F)
balanced.centroid_importance = sort(rowSums(abs(balanced.centroid)), decreasing = T)
balanced.centroid = balanced.centroid[ match(names(balanced.centroid_importance),rownames(balanced.centroid)),]

#pure_data = read.table("~/Koop_Klinghammer/Data/35S.14_03_2018.normalized.tsv", sep ="\t", header = T, row.names = 1)
pure_data = expr_raw

### Preparation

rownames(pure_data) = str_to_upper(rownames(pure_data))
rownames(pure_data) = str_replace_all(rownames(pure_data), pattern = "\\_", "")
rownames(pure_data) = str_replace_all(rownames(pure_data), pattern = "-", "")
rownames(balanced.centroid) = str_to_upper(rownames(balanced.centroid))
rownames(balanced.centroid) = str_replace_all(rownames(balanced.centroid), pattern = "\\_", "")
rownames(balanced.centroid) = str_replace_all(rownames(balanced.centroid), pattern = "-", "")
colnames(pure_data) = str_replace_all(colnames(pure_data), pattern = "^X", "")

col_var = apply(as.matrix(pure_data),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 2)
row_var = apply(as.matrix(pure_data),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 1)

pure_data = pure_data[row_var > 0, col_var > 0]
dim(pure_data)

meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_data = meta_info[match(colnames(pure_data), meta_info$Name),]

###

expr = pure_data
source("~/Koop_Klinghammer/Scripts/Classification_scripts.R")
table( rownames(expr) %in% rownames(balanced.centroid) )

centroid_genes = rownames(balanced.centroid)
genes_matching_centroid = rownames(expr)[which(  (rownames(expr) %in% rownames(balanced.centroid) ) ) ]
genes_matching_centroid = c( genes_matching_centroid, rep("",length(centroid_genes) - length(genes_matching_centroid)) )
genes_not_matching_centroid = rownames(expr)[which(!(  (rownames(expr) %in% rownames(balanced.centroid) ) )) ]
genes_not_matching_centroid = c(genes_not_matching_centroid, rep("",length(centroid_genes) - length(genes_not_matching_centroid)) )

genes_present = data.frame(
  "Matching_to_centroid" = genes_matching_centroid,
  "Not_matching_to_centroid" = genes_not_matching_centroid,
  "Centroid_genes" = centroid_genes
)
#write.table(genes_present, "~/Koop_Klinghammer/Results/First_results_11S_1C/Centroid_genes_And_not_centroid_genes.tsv",sep ="\t",quote = F, row.names = F)

### centroid classification

pub_cor <<- matrix( as.double(), ncol = length( colnames( balanced.centroid )  ) )
expr2bc = centroid2expr( balanced.centroid[,], expr )
colnames(expr2bc$correlation) = c("Sample","Subtype","Correlation","P_value")
class_data = as.data.frame(expr2bc$correlation)

meta_match = match( class_data$Sample, meta_info$Name, nomatch = 0 )
meta_info$Subtype_Sadanandam[meta_match] = as.character( class_data$Subtype )
meta_info$Significance_Sadanandam[meta_match] = as.double( as.character( class_data$P_value ) )
#meta_info$Sig =meta_info$P_value < 0.05 

write.table(meta_info,"~/MAPTor_NET/Misc/Meta_information.tsv",sep ="\t",quote =F,row.names =F)
