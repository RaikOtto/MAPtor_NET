library("stringr")
library("bedr")
library("grid")     ## Need to attach (and not just load) grid package

meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
draw_colnames_45 <- function (coln, gaps, ...) {coord = pheatmap:::find_coordinates(length(coln), gaps);x = coord$coord - 0.5 * coord$size;  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...));  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

## merge data together

raw_files = list.files("/home/janniklas/RNA_Seq/workflow_groetz/results.hg19/", full.names = T, pattern = ".vcf")
raw_files_names = list.files("/home/janniklas/RNA_Seq/workflow_groetz/results.hg19/", full.names = F, pattern = ".vcf")

sample_names = c()
mutation_count_y_chrom = c()

for ( i in 1:length(raw_files)){
    
    vcf_h = read.vcf( raw_files[i], verbose = F )
    print( nrow(vcf_h$vcf[ vcf_h$vcf$CHROM == "chrY", ]) )
    sample_names = c( sample_names, str_replace( raw_files_names[i], pattern = ".filtered.hg19.vcf", "") )
    mutation_count_y_chrom = c(mutation_count_y_chrom, nrow(vcf_h$vcf[ vcf_h$vcf$CHROM == "chrY", ]))
}

match_names =  meta_info$Name[ match( sample_names, meta_info$Name_2, nomatch = 0) ]
Mean_count = mean(mutation_count_y_chrom)
above_mean = rep("TRUE", length(mutation_count_y_chrom))
above_mean[ mutation_count_y_chrom < Mean_count ] = "FALSE"

y_cov_frame = data.frame(
  "Sample_name" = match_names,
  "Count_Muts_Ychrom" = mutation_count_y_chrom,
  "Above_mutation_mean" = above_mean
)
y_cov_frame = y_cov_frame[ order(mutation_count_y_chrom, decreasing = T),]
WriteXLS::WriteXLS(y_cov_frame, "~/MAPTor_NET/Results/Mutation_data/Y_chrom_mutations.xlsx")
###

mut_t = read.table("~/MAPTor_NET/Results/Mutation_data/mutation_matrix.txt", sep = "\t", header = T, stringsAsFactors = F)
mut_t[1:5,1:5]

raw_genes = str_split( rownames(mut_t), pattern = "_" )
genes = str_replace( sapply( raw_genes, FUN = function(vec){ return(head(vec,1)) } ), pattern = "p.", "" )
mutations = sapply( raw_genes, FUN = function(vec){ return(tail( as.character(vec),1) ) } )

row_means = sort(rowMeans(mut_t), decreasing = T)
summary(row_means)
mut_t_f = mut_t[ row_means >= 0.8,]
genes_f = genes[ row_means >= 0.8]
mutations_f = mutations[ row_means >= 0.8]

paste(genes_f, mutations_f, sep ="__")
sort(table(genes_f), decreasing = T)

dim(mut_t_f)
pheatmap::pheatmap( cor(mut_t_f) , show_rownames = F)
