### load_data

final_centroid_list = read.table(
  "~/MAPTor_NET//Results/Classification/Centroids_IT6.tsv",
  sep = "\t",
  header = T,
  stringsAsFactors = F
)

alpha_marker = final_centroid_list$name[final_centroid_list$Alpha > .2]; beta_marker = final_centroid_list$name[final_centroid_list$Beta > .2]; gamma_marker = final_centroid_list$name[final_centroid_list$Gamma > .2]; delta_marker = final_centroid_list$name[final_centroid_list$Delta > .2]
pancreasMarkers = list("Alpha" = alpha_marker,"Beta" = beta_marker,"Gamma" = gamma_marker,"Delta" = delta_marker)

# count data baron

#eislet_ori = readRDS("~/MAPTor_NET//RNA_seq_raw_data/islet-eset.rds")
#count_data = exprs(eislet_ori)
#count_data = count_data[,pData(eislet_ori)$cellType %in% c("alpha","beta","gamma","delta")]
#hgnc_list = as.character( mapIds(org.Hs.eg.db,keys=rownames(count_data),column="ENSEMBL",keytype="SYMBOL",multiVals="first") )
#hgnc_list_uni = unique(hgnc_list)
#source("~/MAPTor_NET//Scripts/Variance_selection.R")
#count_data = variance_filtered_mat
#count_data[1:5,1:5]

#count_data = count_data[! is.na(rownames(count_data)),]
#count_data = count_data[match(  unique(rownames(count_data)), rownames(count_data) ),]
#

#eislet = new("ExpressionSet", exprs=(count_data))
#fData(eislet) = data.frame(rownames(count_data))
#source("~/MAPTor_NET/Scripts/Update_meta_info_table.R")
#pData(eislet) = meta_data[ match(colnames(exprs(eislet)),meta_data$Name),]

#

#p_data = pData(eislet_ori)
#p_data = p_data[rownames(p_data) %in% colnames(count_data),]
#p_data$Subtype = p_data$cellType
#p_data$Name = colnames(count_data)
#p_data = p_data[,c("Name","Subtype")]
#p_data$Subtype[p_data$Subtype == "pericyte"] = "Pericyte";p_data$Subtype[p_data$Subtype == "schwann"] = "Schwann";p_data$Subtype[p_data$Subtype == "other"] = "Other";p_data$Subtype[p_data$Subtype == "mast"] = "Mast";p_data$Subtype[p_data$Subtype == "macrophage"] = "Macrophage";p_data$Subtype[p_data$Subtype == "epsilon"] = "Epsilon";p_data$Subtype[p_data$Subtype == "endothelial"] = "Endothelial";p_data$Subtype[p_data$Subtype == "alpha"] = "Alpha";p_data$Subtype[p_data$Subtype == "beta"] = "Beta";p_data$Subtype[p_data$Subtype == "gamma"] = "Gamma";p_data$Subtype[p_data$Subtype == "delta"] = "Delta";p_data$Subtype[p_data$Subtype == "acinar"] = "Acinar";p_data$Subtype[p_data$Subtype == "stellate"] = "Stellate";p_data$Subtype[p_data$Subtype == "ductal"] = "Ductal"
#pData(eislet) = p_data

# marker init baron

#pancreasMarkers = readRDS("~/MAPTor_NET/Misc/Marker_genes_baron.list")
#data(pancreasMarkers)
#alpha_marker = hgnc_list = as.character( mapIds(org.Hs.eg.db,keys=as.character(unlist(pancreasMarkers["alpha"])),column="ENSEMBL",keytype="SYMBOL",multiVals="first") )
#pancreasMarkers$Alpha = pancreasMarkers$alpha;pancreasMarkers$Beta = pancreasMarkers$beta;pancreasMarkers$Gamma = pancreasMarkers$gamma;pancreasMarkers$Delta = pancreasMarkers$delta;
#pancreasMarkers = pancreasMarkers[c("Alpha","Beta","Gamma","Delta")]

# count_data based training

source("~/MAPTor_NET/Scripts/Update_meta_info_table.R")
count_data = count_data[,meta_data$Subtype %in% c("Alpha","Beta","Gamma","Delta")]
source("~/MAPTor_NET/Scripts/Update_meta_info_table.R")
eislet = new("ExpressionSet", exprs=as.matrix(count_data))
fData(eislet) = data.frame(rownames(count_data))
pData(eislet) = meta_data[ match( colnames(exprs(eislet)), meta_data$Name ) , ]

count_data[1:5,1:5]
table(meta_data$Subtype)

eset = new("ExpressionSet", exprs=(as.matrix(bam_counts)));
fData(eset) = data.frame(rownames(bam_counts))

# run
pancreasMarkers$Alpha = pancreasMarkers$Alpha[(pancreasMarkers$Alpha %in% intersect( rownames(bam_counts),rownames(count_data)))]
pancreasMarkers$Beta = pancreasMarkers$Beta[(pancreasMarkers$Beta %in% intersect( rownames(bam_counts),rownames(count_data)))]
pancreasMarkers$Gamma = pancreasMarkers$Gamm[(pancreasMarkers$Gamma %in% intersect( rownames(bam_counts),rownames(count_data)))]
pancreasMarkers$Delta = pancreasMarkers$Delta[(pancreasMarkers$Delta %in% intersect( rownames(bam_counts),rownames(count_data)))]

B <- bseqsc_basis(eislet, pancreasMarkers, clusters = 'Subtype', samples = 'Name', ct.scale = FALSE)
fit <- bseqsc_proportions(eset, B, verbose = TRUE, log = F)

###

res_coeff = t(fit$coefficients);res_coeff[ is.na(res_coeff) ] = 0.0
res_cor   = fit$stats;res_cor[ is.na(res_cor) ] = 0.0
res_coeff[ res_cor[,"Correlation"] <= .3, ] = .0
res_cor[ res_cor[,"Correlation"] <= .3, "Correlation" ] = 0.0

pheatmap::pheatmap(
  t(res_coeff),
  #annotation_col = df[c("INS","GCG","PPY","SST","MKI67","Correlation","MEN1","ATRX","RICTOR","MYC","DAXX","CDKN2A","Sad_KNN5","Study")],
  #annotation_colors = aka3,
  annotation_legend = T,
  show_colnames = T,
  show_rownames = T,
  cluster_rows = T,
  cluster_cols = T
)

###

meta_match = match( colnames(bam_counts), meta_info$Name, nomatch = 0 ) 
res_coeff_match = match( rownames(res_coeff), meta_info$Name, nomatch = 0 )
res_cor_match = match( rownames(res_cor), meta_info$Name, nomatch = 0 ) 

meta_info$Alpha = rep("",nrow(meta_info))
meta_info$Beta = rep("",nrow(meta_info))
meta_info$Gamma = rep("",nrow(meta_info))
meta_info$Delta = rep("",nrow(meta_info))
meta_info$Res_cor = rep("",nrow(meta_info))
meta_info$Alpha[res_coeff_match] = as.data.frame(res_coeff)$Alpha
meta_info$Beta[res_coeff_match] = as.data.frame(res_coeff)$Beta
meta_info$Gamma[res_coeff_match] = as.data.frame(res_coeff)$Gamma
meta_info$Delta[res_coeff_match] = as.data.frame(res_coeff)$Delta
meta_info$Res_cor[res_cor_match] = as.data.frame(res_cor)$Correlation

#write.table(meta_info, "~/MAPTor_NET/Misc/Meta_information.tsv",sep="\t", row.names = F, quote = F)