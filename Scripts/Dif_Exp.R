library("DESeq2")
library("stringr")
library("limma")

expr_raw = read.table("~/MAPTor_NET/BAMs/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
dim(expr_raw)
c("HOXD1","HOXD12") %in% rownames(expr_raw)

meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info$Name = str_replace(meta_info$Name, pattern = "^X", "")
rownames(meta_info) = meta_info$Name
expr_raw = expr_raw[,colnames(expr_raw) %in%  meta_info$Name]
meta_data = meta_info[ colnames(expr_raw), ]
#meta_data = meta_data %>% filter( Included %in% c("Yes"))
rownames(meta_data) = meta_data$Name
expr_raw  = expr_raw[,meta_data$Name]
row_var = apply(expr_raw, FUN = var, MARGIN = 1)
#expr_raw = expr_raw[row_var > 1,]
meta_data = meta_info[colnames(expr_raw),]
expr_raw = expr_raw[ , meta_data$Grading != ""]
expr_raw = expr_raw[ , meta_data$Included == "Yes"]
meta_data = meta_info[colnames(expr_raw),]
expr_raw[1:5,1:5]
dim(expr_raw)

#expr_raw = expr_raw + -1*min(expr_raw)
c("HOXD1","HOXD12") %in% rownames(expr_raw)

dds_raw = DESeq2::DESeqDataSetFromMatrix(
  countData = round(expr_raw,0),
  colData = meta_data,
  design = ~ as.factor(NEC_NET_PCA),# + as.factor(Study),
  tidy = F
)
dds_raw = estimateSizeFactors(dds_raw)
res = DESeq(dds_raw)
expr_raw = assay(res)

res_t = results(res)
summary(res_t$pvalue < 0.05)

res_t = data.frame(res_t)
res_t["MEN1",]

#write.table(res_t,"~/MAPTor_NET/Results/Publication/Dif_exp_Deseq_2_HGNC_raw.tsv",sep="\t",quote =F , row.names = T)
design = model.matrix( ~ 0 + meta_data$NEC_NET_PCA ) # contrast
colnames(design) = c("NEC","NET")
cont.matrix = makeContrasts( contrast = NEC - NET,  levels = design )

fit = lmFit( expr_raw[ ,  ], design )
fit = contrasts.fit( fit, cont.matrix )
fit = eBayes( fit )
#volc_all = topTable(fit, adjust="fdr", sort.by="B", number = 50000)
volc_all = topTable( fit, coef = "contrast", number  = nrow(expr_raw), adjust  ="none", p.value = 1, lfc = 0)
volc_all = volc_all[order(volc_all$P.Value,decreasing = F),]

#write.table(volc_all,"~/MAPTor_NET/Results/Dif_exp_Master_original.tsv",sep ="\t",quote = F)
