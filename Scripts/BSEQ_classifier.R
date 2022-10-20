source( "~/MAPTor_NET/Scripts/Library_init.R")
library("bseqsc")
library("limma")
library("grid")     ## Need to attach (and not just load) grid package
library("pheatmap")

meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
draw_colnames_45 <- function (coln, gaps, ...) {coord = pheatmap:::find_coordinates(length(coln), gaps);x = coord$coord - 0.5 * coord$size;  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...));  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

###

bam_counts = read.table("~/MAPTor_NET/BAMs/TPMs.81_Samples_13_11_2017.Groetzinger_Scarpa.tsv",sep ="\t", header = T)
colnames(bam_counts) = str_replace(colnames(bam_counts),pattern = "^X","")
colnames(bam_counts) = str_replace(colnames(bam_counts),pattern = "\\.","_")

#condition = (meta_info$Included == "Yes") & (meta_info$NEC_NET %in% c("NEC","NET")) #& (meta_info$Grading %in% c("G3"))
#include_list  = meta_info$Name[ condition ]
#exclude_list  = meta_info$Name[ ! condition ]
#bam_counts = bam_counts[, ( colnames(bam_counts) %in% include_list )]

dim(bam_counts)
bam_counts[1:5,1:5]
colnames(bam_counts)

count_data = bam_counts
source("~/MAPTor_NET/Scripts/Update_meta_info_table.R")
rownames(meta_data) = meta_data$Name
### normalization

obj = bam_counts_raw = bam_counts
row_var = apply(bam_counts, FUN = var, MARGIN = 1)
col_var = apply(bam_counts, FUN = var, MARGIN = 2)
table(row_var == 0)
table(col_var == 0)
bam_counts = bam_counts[row_var != 0,col_var != 0]
bam_counts = bam_counts[rowSums(bam_counts) >= 1,]
dim(bam_counts)

dds_raw = DESeq2::DESeqDataSetFromMatrix(
  countData = round(bam_counts,0),
  colData = meta_data,
  design = ~ NEC_NET,
  tidy = F
)
dds_raw = DESeq2::estimateSizeFactors(dds_raw)

dds_n = DESeq2::varianceStabilizingTransformation(dds_raw)
library("DESeq2")
bam_counts = assay( dds_n )

vsn::meanSdPlot( bam_counts)
