source("~/MAPTor_NET/Scripts/Library_init.R")

bulk_files_i = list.files(
  "~/MAPTor_NET/BAMs/Groetzinger/",
  pattern = ".bam", full.names = T
)
bulk_files_i = c(bulk_files_i, list.files(
  "~/MAPTor_NET/BAMs/CCLs/",
  pattern = ".bam", full.names = T
))
bulk_files_i = c(bulk_files_i, list.files(
  "~/MAPTor_NET/BAMs/KNIH/",
  pattern = "_6.hg19.bam", full.names = T
))

bam_files_bulk_data = Rsubread::featureCounts(
  files = bulk_files_i,
  annot.ext = normalizePath("~/MAPTor_NET/Misc/Homo_sapiens.GRCh37.85.gtf"), 
  isGTFAnnotationFile = T,
  isPairedEnd = T,
  nthreads = 1
)
#saveRDS(bam_files_bulk_data, "~/MAPTor_NET/BAMs/Groetzinger/Groetzinger_56.RData")

if(F){
bam_file_data = readRDS("~/MAPTor_NET/BAMs/Scarpa/Scarpa_29.RData")

g_len = read.table("~/MAPTor_NET//Misc/Homo_sapiens.GRCh37.85.gtf", sep ="\t", stringsAsFactors = F)
g_len = g_len[ g_len$V3 == "gene",]

gene_ids= as.character(unlist(sapply(g_len$V9, FUN=function(vec){return(as.character(unlist(str_split(vec,pattern=";")))[1])})))
gene_ids = str_replace(gene_ids, pattern = "gene_id ", "")
len = g_len$V5 - g_len$V4 + 1

g_len = data.frame( "ensembl_id" = gene_ids, "gene_length" = len)

row_var = apply(bam_counts, FUN = var, MARGIN = 1)
col_var = apply(bam_counts, FUN = var, MARGIN = 2)
bam_counts = bam_counts[row_var > 0,col_var > 0]

bam_counts = bam_counts[ rownames(bam_counts) %in% g_len$ensembl_id,]
g_match = match(rownames(bam_counts), g_len$ensembl_id,nomatch = 0)
g_len = as.integer(g_len$gene_length[g_match])

#count_data[count_data == 0] = 1
sce_n <- scater::newSCESet( countData = bam_counts)
sce_n <- scater::normaliseExprs( sce_n, method = "RLE")

#### TPM
bam_counts = scater::calculateTPM(sce_n, effective_length = g_len, calc_from = "counts")
bam_counts[1:5,1:5]
summary(bam_counts[,1])

#write.table(bam_counts,"~/MAPTor_NET/BAMs/RNA_Seq_85_samples_read_counts.deseq2_normalized.tpm.tsv",sep="\t",row.names=T,quote=F)
pheatmap(cor(bam_counts), annotation_col = df["Study"])

library("vsn")
meanSdPlot(log(as.matrix(bam_counts)), ranks = FALSE)
###

hgnc_list = as.character( mapIds(
  org.Hs.eg.db,
  keys=row.names(count_data),
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals="first"
) )

centroid_t = read.table("~/MAPTor_Net_RNA_data/Results/Classification/Centroids_IT6.tsv",sep="\t",header = T,stringsAsFactors = F) 
#genes_of_interest = c("FGFR1","FGFR2","FGFR3","BEX3","EGFR","NGFR","IGF1R")
genes_of_interest = centroid_t$id[1:20]

hgnc_list_map = which( hgnc_list %in% genes_of_interest)
exp_data = count_data[hgnc_list_map,]

res = data.frame(
  #"HGNC_symbol" = hgnc_list,
  "HGNC_symbol" = hgnc_list[hgnc_list_map],
  "ENSEMBL_id" = rownames(count_data)[hgnc_list_map],
  #"Counts" = count_data
  "TPM" = exp_data
)
res[ res == 0] = 1
res[3:ncol(res)] = log(res[3:ncol(res)])
res[3:ncol(res)] = round(res[3:ncol(res)],0)
rownames(res) = res$HGNC_symbol

res2 = data.frame(
  "HGNC_symbol" = hgnc_list,
  "ENSEMBL_id" = rownames(count_data),
  "TPM" = count_data
)
c_m = count_data
c_m [ c_m == 0] = 1
res2 = as.matrix(apply(log(c_m),FUN=scale,MARGIN = 2), ncol = ncol(count_data))
#res2 = log(c_m)
res2 = res2[hgnc_list_map,]
rownames(res2) = hgnc_list[hgnc_list_map]
pheatmap::pheatmap(res2)

res = res[ !is.na(res$HGNC_symbol),]
res

write.table(res,"~/MAPTor_Net_RNA_data/Misc/Raw_counts_Bulk_rna_cut_10.tsv",sep="\t",row.names=F,quote=F)
#write.table(res,"~/MAPTor_Net_RNA_data/Misc/TPM_Bulk_rna_cut_10.tsv",sep="\t",row.names=F,quote=F)

res_quant = apply(res[,3:ncol(res)], MARGIN = 2, FUN = function(vec){
  quantile = ecdf(vec)
  return(round(quantile(vec) * 100))
})
res_quant = cbind(as.character(res$HGNC_symbol),as.character(res$ENSEMBL_id),res_quant)
res_quant[1:5,1:5]

#write.table(res_quant,"~/MAPTor_Net_RNA_data/Misc/Quantile_TPM_Bulk_rna_cut_10.tsv",sep="\t",row.names=F,quote=F)
}

