library("stringr")
library("grid")     ## Need to attach (and not just load) grid package
draw_colnames_45 <- function (coln, gaps, ...) {coord = pheatmap:::find_coordinates(length(coln), gaps);x = coord$coord - 0.5 * coord$size;  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...));  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

mir_data = read.table("~/MAPTor_NET/Misc/Sadanandam_miR_stubypes.csv", sep ="\t", header = T, quote = "\"")
mir_data = as.data.frame(mir_data)

exp_data = read.table("~/MAPTor_NET/Misc/Sadanandam_gene_expr_subtypes.csv", sep ="\t", header = T, quote = "\"")
mir_data = as.data.frame(exp_data)

match_genes = intersect(mir_data$HGNC_gene, exp_data$HGNC_gene)
mir_match = match(match_genes, mir_data$HGNC_gene, nomatch = 0)
exp_match = match(match_genes, exp_data$HGNC_gene, nomatch = 0)
merge_data = data.frame(
    "HGNC_genes" = match_genes
)
mir_data = cbind(merge_data, mir_data[mir_match,-1], exp_data[exp_match,-1] )

normals = which(str_detect(pattern = "Normal", colnames(mir_data)))
MLP = which(str_detect(pattern = "MLP", colnames(mir_data)))
#MLP2 = which(str_detect(pattern = "MLP.2", colnames(mir_data)))
insulinoma = which(str_detect(pattern = "Insulinoma", colnames(mir_data)))
intermediate = which(str_detect(pattern = "Intermediate", colnames(mir_data)))

centroids = data.frame(
  "Genes" = mir_data$HGNC_gene,
  "Norm"  = rowMeans(mir_data[,normals]),
  "Insulinoma"  = rowMeans(mir_data[,insulinoma]),
  "Intermediate"  = rowMeans(mir_data[,intermediate]),
  #"MLP1"  = rowMeans(mir_data[,MLP1]),
  "MLP"  = rowMeans(mir_data[,MLP])
)

#selection = which( rowSums( abs(centroids[,c(-1)])) >  4 )
#pheatmap::pheatmap(centroids[ selection ,c(-1)], labels_row = centroids$Genes[selection])
#pheatmap::pheatmap(cor(centroids[,-1]))

z_expr = scale(expr_raw_raw)
z_expr[1:5,1:5]

matcher = match( as.character( centroids$Genes ), as.character( rownames(z_expr) ), nomatch = 0)

res = cor( z_expr[matcher,], centroids[matcher!= 0,-1]   )
#pheatmap(res)

result <<- data.frame(
  "Sample" = character(),
  "Subtype" = character(),
  "P_value" = double()
)

for ( sample in colnames(z_expr)){
  
    sample_col = which( colnames(z_expr)==sample)
    correlations = cor( z_expr[matcher,sample_col], centroids[matcher!= 0,-1]   )
    max_cor = which.max( correlations  )+1
    
    p_value = cor.test( z_expr[matcher,sample_col], centroids[ matcher != 0,max_cor])
    result <<- data.frame(
      "Sample" = c( as.character( result$Sample ), as.character( sample) ),
      "Subtype" = c( as.character( result$Subtype ), as.character( colnames( centroids)[ max_cor] ) ),
      "P_value" = c( result$P_value, p_value$p.value)
    )
}

meta_match = match(result$Sample, meta_info$ Name)

meta_info$Subtype_Sad = rep("", nrow(meta_info))
meta_info$P_value_Sad = rep("", nrow(meta_info))
meta_info$Subtype_Sad[meta_match] = as.character(result$Subtype)
meta_info$P_value_Sad[meta_match] = as.double( result$P_value )

meta_info$P_value_Sad[ as.double( meta_info$P_value_Sad ) > 0.05] = "Not_sig"
meta_info$P_value_Sad[ meta_info$P_value_Sad != "Not_sig" ] = "Sig"

write.table(meta_info, "~/MAPTor_NET/Misc/Meta_information.tsv", sep ="\t", quote =F, row.names = F)

### centroid
 
library("pamr")

classifier = rep( "", ncol(mir_data) )
classifier[normals] = "normal"
classifier[MLP] = "MLP"
#classifier[MLP2] = "MLP2"
classifier[insulinoma] = "insulinoma"
classifier[intermediate] = "intermediate"
classifier = classifier[-1]

train_mat = data.frame(mir_data[,-1])
rownames(train_mat) = mir_data$HGNC_gene
colnames(train_mat) = 1:ncol(train_mat)

my_data = list(
  x = as.matrix(train_mat),
  y = classifier,
  genenames = as.character ( mir_data$HGNC_gene ),
  geneid = as.character( mir_data$HGNC_gene)
)

fit = pamr.train(
  my_data,
  n.threshold = 100#,
  #ngroup.survival = 5
)
min_threshold = tail(which( fit$errors == min(fit$errors)),1)

centroids = as.data.frame( fit$centroids,ncol=5)
centroids$means = rowMeans( centroids )
centroids = centroids[order( centroids$means, decreasing = T),]
centroids = centroids[,-ncol(centroids)]

centroids = data.frame(
  "Genes" = rownames( centroids ),
  "Norm"  = centroids[ ,colnames(centroids ) == "normal"],
  "Insulinoma"  = centroids[ ,colnames(centroids ) == "insulinoma"],
  "Intermediate"  = centroids[ ,colnames(centroids ) == "intermediate"],
  #"MLP1"  = centroids[ ,colnames(centroids ) == "MLP1"],
  "MLP"  = centroids[ ,colnames(centroids ) == "MLP"]
)

write.table(centroids,"~/MAPTor_NET/Misc/Centroid_table_merge_c_e.tsv", sep ="\t", row.names = F, quote = F)

#### plots

expr = mir_data[,-1]
groups = rep("",ncol(expr))

groups[normals-1] = "Norm"
groups[intermediate-1] = "Intermediate"
groups[insulinoma-1] = "Insulinoma"
groups[MLP-1] = "MLP"

cor_mat = cor(expr);pcr = prcomp(t(cor_mat))
p = ggbiplot::ggbiplot(
  pcr,
  alpha = 1,
  obs.scale = .5,
  groups = as.character(groups),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F#, labels = rownames(df)
)
p = p + guides(color=guide_legend(title="Study")) + guides(shape=guide_legend(title="Grading"))+ guides(size=guide_legend(title="MKI67"))
p = p + geom_point( aes( shape = as.factor(df$Grading),color = as.factor(df$Study), size = mki67_exp_scaled ) )
p

mir_data = data.frame(
    "Name" = colnames(expr),
    "Group" = groups
)
rownames(mir_data) = mir_data$Name

dds_raw = DESeq2::DESeqDataSetFromMatrix(
  countData = round(expr,0),
  colData = mir_data,
  design = ~ Groups,
  tidy = F
)
dds_raw = DESeq2::estimateSizeFactors(dds_raw)

dds_n = DESeq2::varianceStabilizingTransformation(dds_raw)
