source("~/MAPTor_NET//Scripts/Update_meta_info_table.R")
#sig_thr = as.double(meta_data$P_value) > 5*10^(-2)
#meta_data$Classification[sig_thr] = "Not_sig"

meta_data$Insulin = as.double(meta_data$Insulin)
meta_data$GCG = as.double(meta_data$GCG)
meta_data$SST = as.double(meta_data$SST)
meta_data$PPY = as.double(meta_data$PPY)
#pd <- new("AnnotatedDataFrame", data = meta_data);rownames(pd) = meta_data$Name;sce <- newSCESet( tpmData = as.matrix(count_data), phenoData = pd)
pd <- new("AnnotatedDataFrame", data = meta_data);rownames(pd) = meta_data$Name;sce <- newSCESet( tpmData = as.matrix(c2), phenoData = pd)
sig_match = match(meta_data$Name[!(meta_data$Classification %in% "Not_sig")],colnames(count_data),nomatch = 0)
only_sig_count = count_data[,sig_match]
#cor_mat_sig = cor(only_sig_count); saveRDS(cor_mat, "~/MAPTor_Net_RNA_data/Results/8685/cor_mat_tpm_Muraro_noly_sig.RData")
cor_mat_sig = readRDS("~/MAPTor_Net_RNA_data/Results/8685/cor_mat_tpm_Muraro_noly_sig.RData")
#pcr_only_sig = prcomp(t(cor_mat_sig)); saveRDS(cor_mat_sig, "~/MAPTor_Net_RNA_data/Results/8685/pcr_mat_muraro_only_sig.RData")
pcr_only_sig = readRDS("~/MAPTor_Net_RNA_data/Results/8685/pcr_mat_muraro_only_sig.RData")

#pdf("~/MAPTor_Net_RNA_data/Results/8685/Iterative_classification_results/ggBiplot_Pca_Muraro_only_sig.pdf")
ggbiplot::ggbiplot(
  pcr_only_sig,
  obs.scale = 1,
  var.scale = 1, 
  labels.size = 4,
  alpha = 1,
  groups = meta_data$Classification[meta_data$Classification != "Not_sig"],
  ellipse = TRUE, 
  circle = TRUE,
  var.axes = F
)
#dev.off()

### additional vis
#c2 = count_data[rownames(count_data)%in% final_centroid_list$name,]
cor_mat_2 = cor(c2)

cor_mat = cor(count_data)
pd <- new("AnnotatedDataFrame", data = meta_data)
rownames(pd) = meta_data$SampleName
df = data.frame( meta_data )
df = df[, c("Study","Subtype")]

rownames(df) = colnames(cor_mat)

r = sample(1:ncol(cor_mat),300,replace = F)
pheatmap::pheatmap(
  cor_mat_2[r,r],
  annotation_col = df,
  show_rownames = F,
  show_colnames = F,
  #annotation_colors = aka3,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
)

###
prep_mat = scale((count_data))
#prep_mat = prep_mat[rownames(prep_mat) %in% c("ENSG00000115263","ENSG00000188175","ENSG00000254647","ENSG00000108849","ENSG00000157005"),]
prep_mat = prep_mat[rownames(prep_mat) %in% rownames(balanced.centroid)[1:20],]

prep_mat = data.frame(prep_mat)
colnames(prep_mat) = str_replace_all(colnames(prep_mat),pattern = "\\.","_")

hgnc_list = as.character( mapIds(
  org.Hs.eg.db,
  keys=row.names(prep_mat),
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals="first"
) )
rownames(prep_mat) = hgnc_list

labels = matrix(rep(colnames(prep_mat),nrow(prep_mat)), ncol = nrow(prep_mat))

vis_mat = data.frame(
  "Sample" = as.vector(apply(labels, FUN = c, MARGIN = 1)),
  "Gene" = as.factor(rep(rownames(prep_mat), ncol(prep_mat))),
  "TPM" = as.vector(apply(prep_mat, FUN = c, MARGIN = 2))
)

library("ggplot2")
vis_mat = vis_mat[vis_mat$Sample != "SRR3617197",]
samples = unique(vis_mat$Sample)

pl <- lapply( 1:length(samples),
  function(.x){
    sample = samples[.x]
    v_mat = vis_mat[ vis_mat$Sample == sample ,]
    expr_plot = qplot( data=v_mat,x=Gene, y = TPM, fill = factor(Gene))+geom_bar(stat="identity")
    expr_plot = expr_plot + theme( axis.text.x = element_text(angle = 330, vjust = 0.5) )
    expr_plot = expr_plot + ylab("") + xlab(sample) + theme(legend.position="none")
    expr_plot
  }
)
gridExtra::grid.arrange(grobs=pl, ncol=5,nrow = 5)

cor_mat = cor(m_t)
df = data.frame( meta_data );
df = df[, c("SST","PPY","GCG","Insulin") ]
rownames(df) = meta_data$Name
df$Insulin = as.double(df$Insulin);df$GCG = as.double(df$GCG); df$PPY = as.double(df$PPY); df$SST = as.double(df$SST)

prep_mat = cor((count_data));
prep_mat = prep_mat[rownames(prep_mat) %in% rownames(balanced.centroid)[1:30],]
hgnc_list = as.character( mapIds(
  org.Hs.eg.db,
  keys=row.names(prep_mat),
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals="first"
) )
rownames(prep_mat) = hgnc_list

pheatmap::pheatmap(
  prep_mat,
  annotation_col = df,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
)
