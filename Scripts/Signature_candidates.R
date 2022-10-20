
misc_folder = "~/MAPTor_Net_RNA_data/Misc/"
maptor_folder = "~/MAPTor_Net_RNA_data/"

source( paste0( maptor_folder, "/Scripts/Library_init.R"))
source( paste0( maptor_folder, "/Scripts/Pre_processing.R"))

source(paste0(maptor_folder,"Scripts/Data_polishing.R"))

### alt ###
count_data = read.table("~/MAPTor_Net_RNA_data/Results/count_data.tsv",sep ="\t", header = T)
colnames(count_data) = str_replace(colnames(count_data),pattern = "^X","")
grep(colnames(count_data),pattern = "^X", value =T)
source(paste0(maptor_folder, "/Scripts/Update_meta_info_table.R"))
dim(count_data)
dim(meta_data)
#


pd <- new("AnnotatedDataFrame", data = meta_data)
rownames(pd) = meta_data$SampleName
count_data[count_data == 0] = 1
sce_n <- newSCESet( countData = count_data, phenoData = pd)

sce_n <- normaliseExprs( sce_n, method = "TMM")
count_data = scater::counts(sce_n)

var.fit <- trendVar(sce_n, trend="loess", use.spikes=FALSE, span=0.2)
var.out <- decomposeVar(sce_n, var.fit)

hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] 
nrow(hvg.out)

set.seed(100)
var.cor <- correlatePairs(sce_n, subset.row=rownames(hvg.out))
sig.cor <- var.cor$FDR <= 0.05
chosen <- unique(c(var.cor$gene1[sig.cor], var.cor$gene2[sig.cor]))
norm.exprs <- exprs(sce_n)[chosen,,drop=FALSE]
heat.vals <- norm.exprs - rowMeans(norm.exprs)


library("ggplot2")
sig_cands = read.table("~/MAPTor_Net_RNA_data/Results/Dif_exp/Signature_candidates.tsv",sep ="\t",header =T, stringsAsFactors = T)

for  ( gene in as.character(hgnc_list_chosen)){
  
  print(gene)
  
  map_id = match(gene , hgnc_list, nomatch = 0)
  exp_data = count_data[map_id,]
  data_mat = as.data.frame(cbind(exp_data, as.factor(meta_data$Subtype)))
  colnames(data_mat) = c("Exp","Subtype")
  
  data_mat$Subtype <- factor(data_mat$Subtype, levels = c("1", "2", "5","3","4","6"))
  
  vio_all = ggplot( data_mat, aes( x = factor(Subtype), y = log2(Exp)  ) ) 
  vio_all = vio_all + geom_violin(alpha=0.5, aes(col = factor(Subtype) ) )
  vio_all = vio_all + geom_jitter(alpha=0.25, aes( ), position = position_jitter(width = 0.05))
  #vio_all = vio_all + theme(legend.position=c(.5, .1),legend.direction="horizontal")
  vio_all = vio_all + scale_x_discrete(
    gene,
    labels = c(
      "Alpha","Beta","Gamma","Delta","Epsilon","Unknown")) 
  vio_all = vio_all + labs(y = "Log2 scRNA Expression", x = "")
  vio_all = vio_all + theme(legend.position="none")
  
  jpg_name = paste0( collapse = "", c("~/MAPTor_Net_RNA_data/Results/Signature_candidates/",as.character(gene),".jpeg"  ) )
  ggsave(jpg_name, vio_all)
}
