library("grid")
library("stringr")
library("dplyr")
library("viridis")
library("ggplot2")
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

expr_raw = read.table(
  file="~/MAPTor_NET/BAMs_new/Master/Master_new.S34.HGNC.DESeq2.VOOM.tsv",
  sep="\t", stringsAsFactors =  F, header = T)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
dim(expr_raw)

grep(rownames(expr_raw), pattern = "HOXB", value = T)

meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
#meta_info$Sample = str_replace(meta_info$Sample, pattern = "^X", "")
rownames(meta_info) = meta_info$Sample
expr_raw = expr_raw[,colnames(expr_raw) %in%  meta_info$Sample]
meta_data = meta_info[colnames(expr_raw),]
dim(expr_raw)
y_lim_val = 16

### HOX plots ###

### NKX

genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/BAMs/Hox_Data/NKX.csv",sep ="\t", stringsAsFactors = F, header = F)[,1]
genes_of_interest_hgnc_t[ ! (genes_of_interest_hgnc_t %in% rownames(expr_raw))]
contained_genes = length(genes_of_interest_hgnc_t[  (genes_of_interest_hgnc_t %in% rownames(expr_raw))])

expr = expr_raw[ genes_of_interest_hgnc_t[genes_of_interest_hgnc_t %in% rownames(expr_raw)],]
vis_mat = reshape2::melt(as.matrix(expr));
colnames(vis_mat) = c("Gene","Sample","Expression")
vis_mat$Gene = rep(rownames(expr),ncol(expr))
vis_mat$Gene = factor(
  as.character(vis_mat$Gene),
  levels= genes_of_interest_hgnc_t[  (genes_of_interest_hgnc_t %in% rownames(expr_raw))])

vis_mat$NEC_NET_Ori = meta_data[vis_mat$Sample,"NEC_NET_Ori"]
genes = names(table(vis_mat$Gene)); length(genes)

vis_mat_filt = vis_mat %>% filter( Gene %in% genes[1:6])
vis_mat_filt = vis_mat %>% filter( Gene %in% genes[7:length(genes)])

if (min(vis_mat_filt$Expression) < 0)
  vis_mat_filt$Expression = vis_mat_filt$Expression + -1 * min(vis_mat_filt$Expression)
p = ggplot( data = vis_mat_filt,aes( x = Gene, y = log2(Expression+1), fill = NEC_NET_Ori ))
p = p +  geom_boxplot( )
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")+ xlab("") + ylab("")
p = p + scale_fill_manual(values=c("red","blue"))
p = p + expand_limits(x = c(0,6)) 

p = p + ylim(0,y_lim_val)
p = p + geom_hline(yintercept = 3.3, linetype= "dashed")

genes_of_interest_hgnc_t[ ! (genes_of_interest_hgnc_t %in% rownames(expr_raw))]
svg("~/MAPTor_NET/Results/Publication/Figures/Figure_7_a_NKX_A.svg", width = 5, height = 5)
svg("~/MAPTor_NET/Results/Publication/Figures/Figure_7_a_NKX_B.svg", width = 5, height = 5)
p 
dev.off()

### HOXA

genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/BAMs/Hox_Data//HOXA.csv",sep ="\t", stringsAsFactors = F, header = F)[,1]
genes_of_interest_hgnc_t[ ! (genes_of_interest_hgnc_t %in% rownames(expr_raw))]
contained_genes = length(genes_of_interest_hgnc_t[  (genes_of_interest_hgnc_t %in% rownames(expr_raw))])
expr = expr_raw[ genes_of_interest_hgnc_t[genes_of_interest_hgnc_t %in% rownames(expr_raw)],]

vis_mat = reshape2::melt(as.matrix(expr))
colnames(vis_mat) = c("Gene","Sample","Expression")
vis_mat$Gene = rep(rownames(expr),ncol(expr))
vis_mat$NEC_NET_Ori = meta_data[vis_mat$Sample,"NEC_NET_Ori"]
vis_mat$Gene = factor(
  as.character(vis_mat$Gene),
  levels= c("HOXA1","HOXA2","HOXA3","HOXA4","HOXA5","HOXA6","HOXA7","HOXA9","HOXA10","HOXA11","HOXA13"))

genes = names(table(vis_mat$Gene)); length(genes)

vis_mat_filt = vis_mat %>% filter( Gene %in% genes[1:6])
vis_mat_filt = vis_mat %>% filter( Gene %in% genes[7:length(genes)])
if (min(vis_mat_filt$Expression) < 0)
  vis_mat_filt$Expression = vis_mat_filt$Expression + -1 * min(vis_mat_filt$Expression)
p = ggplot( data = vis_mat_filt,aes( x = Gene, y = log2(Expression+1), fill = NEC_NET_Ori ))
p = p +  geom_boxplot( )
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")+ xlab("") + ylab("")
p = p + scale_fill_manual(values=c("red","blue"))
p = p + expand_limits(x = c(0,6)) 

p = p + geom_hline(yintercept = 3.3, linetype= "dashed")
p = p + ylim(0,y_lim_val)
#p = p + scale_y_continuous( limits=c(0,y_lim_val))

svg("~/MAPTor_NET/Results/Publication/Figures/Figure_7_a_HOXA_A.svg", width = 5, height = 5)
svg("~/MAPTor_NET/Results/Publication/Figures/Figure_7_a_HOXA_B.svg", width = 5, height = 5)
p 
dev.off()

### HOXB

genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/BAMs/Hox_Data//HOXB.csv",sep ="\t", stringsAsFactors = F, header = F)[,1]
genes_of_interest_hgnc_t[ ! (genes_of_interest_hgnc_t %in% rownames(expr_raw))]
contained_genes = length(genes_of_interest_hgnc_t[  (genes_of_interest_hgnc_t %in% rownames(expr_raw))])
expr = expr_raw[ genes_of_interest_hgnc_t[genes_of_interest_hgnc_t %in% rownames(expr_raw)],]

vis_mat = reshape2::melt(as.matrix(expr))
colnames(vis_mat) = c("Gene","Sample","Expression")
vis_mat$Gene = rep(rownames(expr),ncol(expr))
vis_mat$NEC_NET_Ori = meta_data[vis_mat$Sample,"NEC_NET_Ori"]
vis_mat$Gene = factor(
  as.character(vis_mat$Gene),
  levels= c("HOXB1","HOXB2","HOXB3","HOXB4","HOXB5","HOXB6","HOXB7","HOXB8","HOXB9","HOXB13"))
genes = names(table(vis_mat$Gene)); length(genes)

vis_mat_filt = vis_mat %>% filter( Gene %in% genes[1:6])
vis_mat_filt = vis_mat %>% filter( Gene %in% genes[7:length(genes)])
if (min(vis_mat_filt$Expression) < 0)
  vis_mat_filt$Expression = vis_mat_filt$Expression + -1 * min(vis_mat_filt$Expression)

p = ggplot( data = vis_mat_filt,aes( x = Gene, y = log2(Expression+1), fill = NEC_NET_Ori ))
p = p +  geom_boxplot()
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")+ xlab("") + ylab("")
p = p + scale_fill_manual(values=c("red","blue"))
p = p + expand_limits(x = c(0,6)) 

p = p + ylim(0,y_lim_val)
p = p + geom_hline(yintercept = 3.3, linetype= "dashed")

genes_of_interest_hgnc_t[ ! (genes_of_interest_hgnc_t %in% rownames(expr_raw))]
svg("~/MAPTor_NET/Results/Publication/Figures/Figure_7_a_HOXB_A.svg", width = 5, height = 5)
svg("~/MAPTor_NET/Results/Publication/Figures/Figure_7_a_HOXB_B.svg", width = 5, height = 5)
p 
dev.off()

### HOXC

genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/BAMs/Hox_Data//HOXC.csv",sep ="\t", stringsAsFactors = F, header = F)[,1]
genes_of_interest_hgnc_t[ ! (genes_of_interest_hgnc_t %in% rownames(expr_raw))]
contained_genes = length(genes_of_interest_hgnc_t[  (genes_of_interest_hgnc_t %in% rownames(expr_raw))])
expr = expr_raw[ genes_of_interest_hgnc_t[genes_of_interest_hgnc_t %in% rownames(expr_raw)],]

vis_mat = reshape2::melt(as.matrix(expr))
colnames(vis_mat) = c("Gene","Sample","Expression")
vis_mat$Gene = rep(rownames(expr),ncol(expr))
vis_mat$NEC_NET_Ori = meta_data[vis_mat$Sample,"NEC_NET_Ori"]
vis_mat$Gene = factor(
  as.character(vis_mat$Gene),
  levels= c("HOXC4","HOXC5","HOXC6","HOXC8","HOXC9","HOXC10","HOXC11","HOXC12","HOXC13"))
genes = names(table(vis_mat$Gene)); length(genes)

vis_mat_filt = vis_mat %>% filter( Gene %in% genes[1:6])
vis_mat_filt = vis_mat %>% filter( Gene %in% genes[7:length(genes)])
if (min(vis_mat_filt$Expression) < 0)
  vis_mat_filt$Expression = vis_mat_filt$Expression + -1 * min(vis_mat_filt$Expression)
p = ggplot( data = vis_mat_filt,aes( x = Gene, y = log2(Expression+1), fill = NEC_NET_Ori ))
p = p +  geom_boxplot()
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")+ xlab("") + ylab("")
p = p + scale_fill_manual(values=c("red","blue"))
p = p + expand_limits(x = c(0,6)) 

p = p + geom_hline(yintercept = 3.3, linetype= "dashed")
p = p + ylim(0,y_lim_val)
genes_of_interest_hgnc_t[ ! (genes_of_interest_hgnc_t %in% rownames(expr_raw))]
svg("~/MAPTor_NET/Results/Publication/Figures/Figure_7_a_HOXC_A.svg", width = 5, height = 5)
svg("~/MAPTor_NET/Results/Publication/Figures/Figure_7_a_HOXC_B.svg", width = 5, height = 5)
p 
dev.off()

### HOXD

genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/BAMs/Hox_Data//HOXD.csv",sep ="\t", stringsAsFactors = F, header = F)[,1]
genes_of_interest_hgnc_t[ ! (genes_of_interest_hgnc_t %in% rownames(expr_raw))]
contained_genes = length(genes_of_interest_hgnc_t[  (genes_of_interest_hgnc_t %in% rownames(expr_raw))])
expr = expr_raw[ genes_of_interest_hgnc_t[genes_of_interest_hgnc_t %in% rownames(expr_raw)],]

vis_mat = reshape2::melt(as.matrix(expr))
colnames(vis_mat) = c("Gene","Sample","Expression")
vis_mat$NEC_NET_Ori = meta_data[vis_mat$Sample,"NEC_NET_Ori"]
vis_mat$Gene = factor(
  as.character(vis_mat$Gene),
  levels= genes_of_interest_hgnc_t[  (genes_of_interest_hgnc_t %in% rownames(expr_raw))])
genes = names(table(vis_mat$Gene)); length(genes)

vis_mat_filt = vis_mat %>% filter( Gene %in% genes[1:6])
vis_mat_filt = vis_mat %>% filter( Gene %in% genes[7:length(genes)])
if (min(vis_mat_filt$Expression) < 0)
  vis_mat_filt$Expression = vis_mat_filt$Expression + -1 * min(vis_mat_filt$Expression)
p = ggplot( data = vis_mat_filt,aes( x = Gene, y = log2(Expression+1), fill = NEC_NET_Ori ))
p = p +  geom_boxplot()
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")+ xlab("") + ylab("")
p = p + scale_fill_manual(values=c("red","blue"))
p = p + expand_limits(x = c(0,6)) 

p = p + ylim(0,y_lim_val)
p = p + geom_hline(yintercept = 3.3, linetype= "dashed")
genes_of_interest_hgnc_t[ ! (genes_of_interest_hgnc_t %in% rownames(expr_raw))]
svg("~/MAPTor_NET/Results/Publication/Figures//Figure_7_a_HOXD_A.svg", width = 5, height = 5)
svg("~/MAPTor_NET/Results/Publication/Figures//Figure_7_a_HOXD_B.svg", width = 5, height = 5)
p
dev.off()

### HMT_complex

genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/BAMs/Hox_Data//HMT_complex.csv",sep ="\t", stringsAsFactors = F, header = F)[,1]
genes_of_interest_hgnc_t[ ! (genes_of_interest_hgnc_t %in% rownames(expr_raw))]
contained_genes = length(genes_of_interest_hgnc_t[  (genes_of_interest_hgnc_t %in% rownames(expr_raw))])
matcher = match(genes_of_interest_hgnc_t, rownames(expr_raw))
expr = expr_raw[matcher ,]

vis_mat = reshape2::melt(as.matrix(expr))
colnames(vis_mat) = c("Gene","Sample","Expression")
vis_mat$NEC_NET_Ori = meta_data[vis_mat$Sample,"NEC_NET_Ori"]
vis_mat$Gene = factor(
  as.character(vis_mat$Gene),
  levels= genes_of_interest_hgnc_t[  (genes_of_interest_hgnc_t %in% rownames(expr_raw))])
genes = names(table(vis_mat$Gene)); length(genes)

vis_mat_filt = vis_mat %>% filter( Gene %in% genes[1:6])
vis_mat_filt = vis_mat %>% filter( Gene %in% genes[7:length(genes)])
if (min(vis_mat_filt$Expression) < 0)
  vis_mat_filt$Expression = vis_mat_filt$Expression + -1 * min(vis_mat_filt$Expression)

p = ggplot( data = vis_mat_filt,aes( x = Gene, y = log2(Expression+1), fill = NEC_NET_Ori ))
p = p +  geom_boxplot( )
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")+ xlab("") + ylab("")
p = p + scale_fill_manual(values=c("red","blue"))
p = p + expand_limits(x = c(0,6)) 

p = p + geom_hline(yintercept = 3.3, linetype= "dashed")
p = p + ylim(0,y_lim_val)

genes_of_interest_hgnc_t[ ! (genes_of_interest_hgnc_t %in% rownames(expr_raw))]
svg("~/MAPTor_NET/Results/Publication/Figures/Figure_7_a_HMT_A.svg", width = 5, height = 5)
svg("~/MAPTor_NET/Results/Publication/Figures//Figure_7_a_HMT_B.svg", width = 5, height = 5)
p  
dev.off()

### MEN1

genes_of_interest_hgnc_t = "MEN1"
genes_of_interest_hgnc_t[ ! (genes_of_interest_hgnc_t %in% rownames(expr_raw))]
contained_genes = length(genes_of_interest_hgnc_t[  (genes_of_interest_hgnc_t %in% rownames(expr_raw))])
genes = genes_of_interest_hgnc_t[genes_of_interest_hgnc_t %in% rownames(expr_raw)]
expr = expr_raw[ genes,]

vis_mat = reshape2::melt(as.matrix(expr))
colnames(vis_mat) = c("Gene","Sample","Expression")
vis_mat$NEC_NET_Ori = meta_data[vis_mat$Sample,"NEC_NET_Ori"]
vis_mat$Gene = factor(
  as.character(vis_mat$Gene),
  levels= genes_of_interest_hgnc_t[  (genes_of_interest_hgnc_t %in% rownames(expr_raw))])
genes = names(table(vis_mat$Gene)); length(genes)

if (min(vis_mat$Expression) < 0)
  vis_mat$Expression = vis_mat$Expression + -1 * min(vis_mat$Expression)

p = ggplot( data = vis_mat,aes( x = Gene, y = log2(Expression+1), fill = NEC_NET_Ori ))
p = p +  geom_boxplot( )
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")+ xlab("") + ylab("")
p = p + scale_fill_manual(values=c("red","blue"))
p = p + expand_limits(x = c(0,6)) 

p = p + ylim(0,y_lim_val)
p = p + geom_hline(yintercept = 3.3, linetype= "dashed")

genes_of_interest_hgnc_t[ ! (genes_of_interest_hgnc_t %in% rownames(expr_raw))]
svg("~/MAPTor_NET/Results/Publication/Figures//Figure_7_a_MEN1.svg", width = 5, height = 5)
p 
dev.off()
