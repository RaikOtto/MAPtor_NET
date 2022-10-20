library("DESeq2")
library("stringr")
library("limma")

#expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Discovery_Cohort.S64.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Discovery_Cohort.S81.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X\\.)", "")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X)", "")
expr_raw[1:5,1:5]
dim(expr_raw)

meta_data = meta_info[colnames(expr_raw),]
design <- model.matrix(~0 + as.factor(meta_data$NET_NEC_UMAP))

#DGE = edgeR::DGEList(expr_raw)
#v = limma::voom(DGE,design = NULL)
#expr_raw = v$E

### DESeq2 normalization
meta_data$Study = factor(meta_data$Study)
meta_data$NET_NEC_UMAP = factor(meta_data$NET_NEC_UMAP)

dds_raw = DESeq2::DESeqDataSetFromMatrix(
  countData = round(expr_raw,0),
  colData = meta_data,
  design = ~ 0 + NET_NEC_UMAP,
  tidy = FALSE
)

dds_raw = estimateSizeFactors(dds_raw)
dds <- DESeq(dds_raw)
vsd = vst(dds_raw, blind=FALSE)
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, vsd$Study)
assay(vsd) <- mat
expr_raw = counts_batch_corrected <- assay(vsd)

expr_raw = assay(DESeq2::varianceStabilizingTransformation(dds_raw))
expr_raw[1:5,1:5]

### new dds_raw

dds <- DESeq(
  dds_raw,
  test = c("Wald", "LRT"),
  fitType = c("parametric", "local", "mean", "glmGamPoi"),
  betaPrior = FALSE
)

res = results(
  dds,
  contrast = c("NET_NEC_UMAP","NEC","NET")
)

res = res[!is.na(res$pvalue),]
summary(res$log2FoldChange)
res = res[order(res$pvalue,decreasing = F),]
sum(res$pvalue<= 0.05)

write.table(res, "~/Downloads/MAPTor-NET_plots_23_02_2022/Dif_exp_discovery_NET_NEC_PCA.tsv",sep ="\t", quote = F)
#write.table(expr_raw, "~/MAPTor_NET/BAMs_new/Publication_datasets/Discovery_Cohort.S64.HGNC.DESeq2.tsv",sep ="\t", quote = F)
res_mll1 = res[rownames(res) %in% sad_genes,]
res_hox =  res[grep(rownames(res),pattern = "^HOX|^MEN1"),]
res = rbind(res_mll1,res_hox)
#write.table(res, "~/Downloads/Dif_exp.discovery.NET_NEC.MLL1_HOX.tsv",sep ="\t", quote = F)
#res_mll1 = res[rownames(res) %in% sad_genes,]
#res_hox =  res[grep(rownames(res),pattern = "^HOX|^MEN1"),]
#res = rbind(res_mll1,res_hox)
#write.table(res, "~/MAPTor_NET/Misc/Differentially_expressed.validation.tsv",sep ="\t", quote = F)

resLFC <- lfcShrink(dds,res = res)
resLFC


###

row_var = apply(round(expr_raw,0), FUN = var, MARGIN = 1)
summary(row_var)
hist(row_var)
expr_raw = expr_raw[row_var > 20,]
dim(expr_raw)

###

meta_info_maptor = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info_maptor) = meta_info_maptor$Sample
colnames(meta_info_maptor) = str_replace(colnames(meta_info_maptor),pattern = "\\.","_")

###

#expr_raw = read.table("~/Deko_Projekt/Data/Cancer_Pancreas_Bulk_Array/Alvarez/GSE98894.raw.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = read.table("~/Deko_Projekt/Data/Cancer_Pancreas_Bulk_Array/GSE73338/GSE73338.all.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = read.table("~/Deko_Projekt/Data/Cancer_Pancreas_Bulk_Array/GSE73339/GSE73339.all.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = read.table("~/Deko_Projekt/Data/Cancer_Pancreas_Bulk_Array/Sato.S35.Ninon.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = read.table("~/MAPTor_NET/BAMs_new/Groetzinger/Groetzinger_new.S40.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = read.table("~/MAPTor_NET/BAMs_new/Scarpa/Scarpa_new.S29.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = read.table("~/MAPTor_NET/BAMs_new/Master/Master_new.S34.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = read.table("~/Deko_Projekt/Data/Cancer_Pancreas_Bulk_Array/Diedisheim.S66.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
dim(expr_raw)

no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match] = str_replace(colnames(expr_raw)[no_match], pattern = "^X","")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match] = paste("X",colnames(expr_raw)[no_match],sep ="")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[which(no_match)]
meta_data = meta_info[colnames(expr_raw),]

dim(expr_raw)
expr_raw[,meta_data[(meta_data$Histology_Primary == "Pancreatic"),"Sample"]]
expr_raw = expr_raw[,meta_data[(meta_data$Histology_Primary == "Pancreatic"),"Sample"]]
dim(expr_raw)

DGE = edgeR::DGEList(expr_raw)
v = limma::voom(DGE,design = NULL)
expr_raw = v$E

#write.table(expr_raw_export, "~/MAPTor_NET/BAMs_new/Publication_datasets/Discovery_Cohort.S64.DESeq2.tsv",sep ="\t", quote = F) #1
#write.table(expr_raw, "~/MAPTor_NET/BAMs_new/Visualizations/Missiaglia.tsv",sep ="\t", quote = F) #2
#write.table(expr_raw, "~/MAPTor_NET/BAMs_new/Visualizations/Sadanandam.S29.tsv",sep ="\t", quote = F) #3 
#write.table(expr_raw, "~/MAPTor_NET/BAMs_new/Visualizations/Sato_S35.DESeq2.tsv",sep ="\t", quote = F) #4
#write.table(expr_raw, "~/MAPTor_NET/BAMs_new/Visualizations/Charite_S29.DESeq2.tsv",sep ="\t", quote = F) #5
#write.table(expr_raw, "~/MAPTor_NET/BAMs_new/Visualizations/Scarpa.S29.HGNC.DESeq2.VOOM.tsv",sep ="\t", quote = F) #6
#write.table(expr_raw, "~/MAPTor_NET/BAMs_new/Visualizations/Master_S20.tsv",sep ="\t", quote = F) #7
#write.table(counts_batch_corrected, "~/MAPTor_NET/BAMs_new/Publication_datasets/Discovery_Cohort.S64.DESeq2.vsd.tsv",sep ="\t", quote = F) #8


