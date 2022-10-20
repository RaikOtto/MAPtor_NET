library("ggplot2")
library("stringr")
library("umap")
library("dplyr")
library("grid")
library("ggrepel")

#dif_exp = read.table("~/MAPTor_NET/Misc/Discovery_no_batch_correction_differential_expression.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)

dif_exp = read.table("~/MAPTor_NET/Misc/Discovery_batch_correction_differential_expression.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
dif_exp$log2FoldChange = dif_exp$log2FoldChange*-1

source("~/MAPTor_NET/Misc//Visualization_colors.R")

i = 76
genes_of_interest_hgnc_t = read.table("~/MAPTor_NET//Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
genes_of_interest_hgnc_t$V1[i]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
sad_genes = sad_genes[sad_genes != ""]

expr = dif_exp[sad_genes,]
expr = expr[!is.na(expr$baseMean),]
diffexpressed_ori = diffexpressed = expr$pvalue
diffexpressed[(diffexpressed_ori < 0.05) & (abs(expr$log2FoldChange) >= .25)] = "Higher expressed NEC-like"
diffexpressed[(diffexpressed_ori < 0.05) & ( expr$log2FoldChange) < -.25 ] = "Higher expressed NET-like"
diffexpressed[(diffexpressed_ori >= 0.05) | (abs(expr$log2FoldChange) < .25)] = "Not Differentially expressed"

p = ggplot(
  data = expr,
  aes(
    x = log2FoldChange,
    y = -log10(pvalue),
    col = diffexpressed,
    label = rownames(expr)
  )
)
p = p + geom_point() + theme_minimal()
p = p + geom_vline(xintercept=c(-.25, .25), col="red") + geom_hline(yintercept=-log10(0.05), col="red")
p = p + scale_color_manual(values=c("red","blue", "black"))
p = p  + geom_text_repel() + theme(legend.position="top")

#svg("~/Downloads/Volcano_plot_discovery_batch_correct_colors_adjusted.svg", width = 10, height = 10)
print(p)
dev.off()

# 132502
