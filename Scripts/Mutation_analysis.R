library(reshape2)
library(ggplot2)
library("stringr")
library("grid")     ## Need to attach (and not just load) grid package
draw_colnames_45 <- function (coln, gaps, ...) {coord = pheatmap:::find_coordinates(length(coln), gaps);x = coord$coord - 0.5 * coord$size;  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...));  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv", sep ="\t", header = T, stringsAsFactors = F)
hgnc_genes = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Stem_signatures.gmt", sep ="\t", header = F, stringsAsFactors = F)
hgnc_genes = str_trim( str_to_upper( hgnc_genes[nrow(hgnc_genes),3:ncol(hgnc_genes)] ))
hgnc_genes = hgnc_genes[hgnc_genes != ""]

## merge data together

raw_files = list.files("~/MAPTor_NET/VCFs/Annotated/Groetzinger/", full.names = T, pattern = ".tsv")
raw_files_names = list.files("~/MAPTor_NET/VCFs/Annotated/Groetzinger/", full.names = F, pattern = ".tsv")
raw_files_names = str_replace(raw_files_names, pattern = ".GATK.RNA_seq.hg38.vcf.annot.filt.tsv", "")

raw_files =  c( raw_files, list.files("~/MAPTor_NET/VCFs/Annotated/Scarpa/", full.names = T, pattern = ".tsv") )
raw_files_names = c( raw_files_names, list.files("~/MAPTor_NET/VCFs/Annotated/Scarpa/", full.names = F, pattern = ".tsv"))
raw_files_names = str_replace(raw_files_names, pattern = ".GATK.RNA_seq.hg38.vcf.annot.filt.tsv", "")

raw_files_names = str_replace_all(raw_files_names, pattern = "-", "_")
meta_match = match(raw_files_names, meta_info$Name_2, nomatch = 0)

raw_files_names = as.character(meta_info$Name[meta_match])
raw_files_names = str_replace(raw_files_names, pattern = ".GATK.RNA_seq.hg38.vcf.annot.filt.tsv", "")

mut_mat <<- cbind("",read.table(raw_files[1], sep ="\t", skip=0, comment.char = "?", header = T,nrows = 1, fill = T, quote = ""))
colnames_mut_mat = colnames(mut_mat)
colnames_mut_mat[1] = "Sample"
mut_mat <<- mut_mat[0,]

for ( i in 1:length(raw_files)){
  print(c(i,raw_files_names[i],raw_files[i]))
  i_mat   = read.table(raw_files[i],sep="\t",stringsAsFactors = F, fill = F, header = T, quote = "")
  i_mat   = i_mat[ (as.double(i_mat$ExAC_ALL) < 0.01 ) | (i_mat$ExAC_ALL == ".")  ,]
  i_mat   = cbind( rep(raw_files_names[i], nrow(i_mat)), i_mat )
  i_mat   = i_mat[ i_mat$avsnp147 != "rs2943772",]
  #i_mat   = i_mat[ i_mat$Gene.refGene == "KMT2A",]
  print(dim(i_mat))
  print(dim(mut_mat))
  mut_mat <<- rbind(mut_mat,i_mat)
}
colnames(mut_mat) = colnames_mut_mat
dim(mut_mat)
mut_mat[1:5,]

#write.table(mut_mat,"~/MAPTor_NET/Results/VCFs/Mutation_matrix.Groetzinger_Scarpa.raw.tsv",sep="\t",row.names= F, quote = F)

###

#mut_mat = read.table("~/MAPTor_NET/Results/VCFs/Mutation_matrix.Groetzinger_Scarpa.raw.tsv", sep ="\t", stringsAsFactors = F, header =T, fill = F)
mut_mat = mut_mat[!is.na(mut_mat$Sample),]
uni_samples = unique(as.character(mut_mat$Sample))
uni_genes = unique(mut_mat$Gene.refGene)

gene_mat <<- matrix(as.integer(),ncol = length(uni_samples))
colnames(gene_mat) = uni_samples

for ( i in 1:length(uni_genes)){
  
  #print(i)
  gene_match  = which(mut_mat$Gene.refGene %in% uni_genes[i])
  
  #if (mut_mat$cosmic82[i] == "."){
  #    mut_vec = rep(0,length(uni_samples))
  #} else {
  mut_samples = as.character(mut_mat$Sample[gene_match])
  mut_vec = rep(0,length(uni_samples))
  mut_vec[colnames(gene_mat) %in% mut_samples] = 1
  #}
  gene_mat <<- rbind( gene_mat, mut_vec)
  
  rownames(gene_mat)[nrow(gene_mat)] = uni_genes[i]
}
colnames(gene_mat)

row_sums = sort(rowSums(gene_mat), decreasing = T)
#summary(row_sums)
#sort((row_sums), decreasing = T)

#write.table(gene_mat,"~/MAPTor_NET/Results/VCFs/Gene_Mutation_matrix.Groetzinger_Scarpa.tsv",sep="\t",row.names= T, quote = F)

# Visualization

gene_mat = read.table("~/MAPTor_NET/Results/VCFs/Gene_Mutation_matrix.Groetzinger_Scarpa.tsv",sep="\t", header  = T,stringsAsFactors = F)
colnames(gene_mat) = str_replace(colnames(gene_mat) , pattern = "^X","")
include_list = meta_info$Name[meta_info$Included == "Yes"]
gene_mat = gene_mat[,which( colnames(gene_mat) %in% include_list )]
dim(gene_mat)
gene_mat[1:5,1:5]

#pdf("~/MAPTor_NET/Results/Mutation_data/pNET_panel_protein_CDS_mutations.pdf",onefile = F, paper="a4r",width = 28, height = 18)
new_rowname = sapply( rownames(gene_mat), FUN = function(name){return(head(as.character(unlist(str_split(name,pattern =";"))),1))} )
rownames(gene_mat) = new_rowname

rownames(meta_info) = meta_info$Name
for (gene in c("MEN1","MYC","RICTOR","ATRX","DAXX","RHEB","CDKN2A","TP53","PIK3CA","RB1","KRAS")){
  meta_info[,gene] = rep("-1",nrow(meta_info))  
  meta_info[colnames(gene_mat),gene] = as.character( gene_mat[gene,])
}

#write.table(meta_info,"~/MAPTor_NET/Misc/Mutation_info.tsv", sep ="\t", row.names = F, quote = F)

## tmp

gene_mat = read.table("~/MAPTor_NET/Misc/Mutation_info.tsv",sep="\t", header = T, row.names = 1)
colnames(gene_mat) = str_replace(pattern = "^X", colnames(gene_mat), "")
gene_mat = gene_mat[,( colnames(gene_mat) %in% include_list )]
dim(gene_mat)
gene_mat[1:5,1:5]

dim(gene_mat)
row_var = apply(gene_mat, FUN = var, MARGIN = 1)
col_var = apply(gene_mat, FUN = var, MARGIN = 2)
table(row_var == 0)
table(col_var == 0)
gene_mat = gene_mat[row_var != 0,col_var != 0]
dim(gene_mat)

# vis gene_mat_flt

create_vis_mat = function ( i_mat, type_analysis ){
  
  cormat <- round(cor(t( i_mat )),2)
  head(cormat)
  dim(cormat)
  
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  upper_tri <- get_upper_tri(cormat)
  upper_tri
  
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  
  # reorder
  
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  cormat <- reorder_cormat(cormat)
  upper_tri <- get_upper_tri(cormat)
  dim(upper_tri)
  
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  out_put_mat = melted_cormat[abs( melted_cormat$value) != 1.0,]
  
  switch(type_analysis,"Panel" = { 
    text_size = 10;
    #write.table(out_put_mat, "~/MAPTor_NET/Results/VCFs/Panel.gene_mutation_correlations.tsv",sep ="\t", quote = F)
    #write.table(i_mat, "~/MAPTor_NET/Results/VCFs/Panel.mutated_genes.tsv",sep ="\t", quote = F)
    #pdf("~/MAPTor_NET/Results/VCFs/Panel.Samples_to_gene_mutations.Panel_genes.pdf",onefile = F, paper="a4r",width = 28, height = 18)
    pheatmap::pheatmap(
      i_mat,
      annotation_col = meta_data[c("NEC_NET","Histology","Grading","Study")],
      annotation_colors = aka3,
      show_rownames = T,
      cluster_rows = F,
      cluster_cols = F
    )
    #dev.off()
  },"Sad" = {
    text_size = 6.5
    #write.table(out_put_mat, "~/MAPTor_NET/Results/VCFs/Sad_221.gene_mutation_correlations.tsv",sep ="\t", quote = F)
    #write.table(i_mat, "~/MAPTor_NET/Results/VCFs/Sad_221.mutated_genes.tsv",sep ="\t", quote = F)
    #pdf("~/MAPTor_NET/Results/VCFs/Sad_221.Samples_to_gene_mutations.pdf",onefile = F, paper="a4r",width = 28, height = 18)
    pheatmap::pheatmap(
      i_mat,
      annotation_col = meta_data[c("NEC_NET","Grading","Study")],
      annotation_colors = aka3,
      fontsize = 5,  
      show_rownames = T
    )
    #dev.off()
  },"All" = {
    text_size = 4
    #write.table(out_put_mat, "~/MAPTor_NET/Results/VCFs/All_200_highest.gene_mutation_correlations.tsv",sep ="\t", quote = F)
    #write.table(i_mat, "~/MAPTor_NET/Results/VCFs/All_200_highest.mutated_genes.tsv",sep ="\t", quote = F)
    #pdf("~/MAPTor_NET/Results/VCFs/All_200_highest.Samples_to_gene_mutations.pdf",onefile = F, paper="a4r",width = 28, height = 18)
    pheatmap::pheatmap(
      i_mat,
      annotation_col = df[c("NEC_NET","Grading","Study")],
      annotation_colors = aka3,
      show_rownames = T
    )
    #dev.off()
  })
  text_size = 10
  ggheatmap = ggplot(melted_cormat, aes(Var2, Var1, fill = value))+ geom_tile(color = "white") + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = text_size),axis.text.y= element_text(size=text_size)) + coord_fixed() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal") + guides(fill = guide_colorbar(barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5))
  ggheatmap
  print(ggheatmap)
}

### panel genes

genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Stem_signatures.gmt",sep="\t",header = F,stringsAsFactors = F)
#genes_of_interest_hgnc = as.character( unique( c(genes_of_interest_hgnc_t[18,3:ncol(genes_of_interest_hgnc_t)], "RB1", "NRAS") ) )
genes_of_interest_hgnc = as.character( genes_of_interest_hgnc_t[ which(genes_of_interest_hgnc_t$V1 == "pNET_Genes_of_interest_2"),3:ncol(genes_of_interest_hgnc_t)] )
genes_of_interest_hgnc = genes_of_interest_hgnc[ genes_of_interest_hgnc != ""]
#genes_of_interest_hgnc = unique(c(genes_of_interest_hgnc, "RB1", "NRAS"))

i_mat <<- gene_mat[ rownames(gene_mat) %in% hgnc_genes, ]
colnames(i_mat) = str_replace(colnames(i_mat), pattern = "^X", "")
i_mat <<- i_mat[ order( rowSums(i_mat), decreasing = T ), ]
for (i in 1:nrow(i_mat)){
  i_mat <<- i_mat[ ,order( i_mat[nrow(i_mat) - i + 1,], decreasing = T ) ]  
}

i_mat[1:5,1:5]
dim(i_mat)

#pdf("~/MAPTor_NET/Results/VCFs/Upper_Triangle_Correlation.Panel_genes.pdf",onefile = F, paper="a4r",width = 28, height = 18)
    create_vis_mat(gene_mat_flt, "Panel")
#dev.off()

meta_match = match(colnames(gene_mat_flt),meta_info$Name,nomatch = 0)
df = meta_info[meta_match,]
rownames(df) = df$Name
df$Status=df$Subtyp_2
df$Tissue=df$Subtype;df$Tissue[df$Tissue == "Beta"] = "CCL";df$Tissue[df$Tissue == "Undefined"] = "CCL";df$Tissue[df$Tissue == "Delta"] = "CCL"
df$Status[df$Study=="Scarpa"] = "Cancer"
df$Tissue[df$Tissue == ""] = "Pancreas"

row_var = apply(gene_mat_flt, FUN = var, MARGIN = 1)
col_var = apply(gene_mat_flt, FUN = var, MARGIN = 2)
table(row_var == 0)
table(col_var == 0)
gene_mat_flt = gene_mat_flt[row_var != 0,col_var != 0]
gene_mat_flt = gene_mat_flt[rowSums(gene_mat_flt) >= 1,]


#pdf("~/MAPTor_NET/Results/VCFs/Correlation.Samples.Panel_genes.pdf",onefile = F, paper="a4r",width = 28, height = 18)
pheatmap::pheatmap(
  cor((gene_mat_flt)),
  annotation_col = df[c("NEC_NET","Grading","Location","Histology","Study")],
  annotation_colors = aka3,
  show_rownames = F
)
#dev.off()

### 221 Sadanandam genes

gene_mat_flt = gene_mat[rownames(gene_mat) %in% genes_of_interest_sad,]
gene_mat_flt[1:5,1:5]
dim(gene_mat_flt)

#pdf("~/MAPTor_NET/Results/VCFs/Upper_Triangle_Correlation.221_Sad_genes.pdf",onefile = F, paper="a4r",width = 28, height = 18)
    create_vis_mat(gene_mat_flt,"Sad")
#dev.off()

meta_match = match(colnames(gene_mat_flt),meta_info$Name,nomatch = 0)
df = meta_info[meta_match,]
rownames(df) = df$Name
df$Status=df$Subtyp_2
df$Tissue=df$Subtype;df$Tissue[df$Tissue == "Beta"] = "CCL";df$Tissue[df$Tissue == "Undefined"] = "CCL";df$Tissue[df$Tissue == "Delta"] = "CCL"
df$Status[df$Study=="Scarpa"] = "Cancer"
df$Tissue[df$Tissue == ""] = "Pancreas"

row_var = apply(gene_mat_flt, FUN = var, MARGIN = 1)
col_var = apply(gene_mat_flt, FUN = var, MARGIN = 2)
table(row_var == 0)
table(col_var == 0)
gene_mat_flt = gene_mat_flt[row_var != 0,col_var != 0]
gene_mat_flt = gene_mat_flt[rowSums(gene_mat_flt) >= 1,]


#pdf("~/MAPTor_NET/Results/VCFs/Correlation.Samples.221_Sad_genes.pdf",onefile = F, paper="a4r",width = 28, height = 18)
pheatmap::pheatmap(
  cor((gene_mat_flt)),
  annotation_col = df[c("NEC_NET","Grading","Study")],
  annotation_colors = aka3,
  show_rownames = T
)
#dev.off()


### all genes

cor_mat = cor(t(gene_mat))genes_of_interest_hgnc
hist(rowMeans(cor_mat))
row_max = apply(cor_mat, MARGIN = 1, FUN = function(vec){return(max(abs(vec)[abs(vec) < .9]))})
sort_vec = sort( row_max, decreasing = T)
selection_vec = which( names( sort_vec[1:200] ) %in% rownames(gene_mat) )

gene_mat_flt = gene_mat[ selection_vec,]
dim(gene_mat_flt)

#pdf("~/MAPTor_NET/Results/VCFs/Upper_Triangle_Correlation.200_highest_cor_unsupervised_genes.pdf",onefile = F, paper="a4r",width = 28, height = 18)
  create_vis_mat(gene_mat_flt,"All")
#dev.off()

meta_match = match(colnames(gene_mat_flt),meta_info$Name,nomatch = 0)
df = meta_info[meta_match,]
rownames(df) = df$Name
df$Status=df$Subtyp_2
df$Tissue=df$Subtype;df$Tissue[df$Tissue == "Beta"] = "CCL";df$Tissue[df$Tissue == "Undefined"] = "CCL";df$Tissue[df$Tissue == "Delta"] = "CCL"
df$Status[df$Study=="Scarpa"] = "Cancer"
df$Tissue[df$Tissue == ""] = "Pancreas"

row_var = apply(gene_mat_flt, FUN = var, MARGIN = 1)
col_var = apply(gene_mat_flt, FUN = var, MARGIN = 2)
table(row_var == 0)
table(col_var == 0)
gene_mat_flt = gene_mat_flt[row_var != 0,col_var != 0]
gene_mat_flt = gene_mat_flt[rowSums(gene_mat_flt) >= 1,]

###

genes_of_interest_hgnc = c( genes_of_interest_hgnc, "GRB2" )
meta_data = cbind(meta_data, t(gene_mat_flt[,match(colnames(gene_mat_flt), meta_data$Name)]))

res_t <<- matrix(as.integer(), ncol = 2)
for( i in 1:nrow(gene_mat)){
  nec = aggregate( as.integer( gene_mat[i,]), by = list(meta_data$NEC_NET), FUN = sum)$x[1]
  net = aggregate( as.integer( gene_mat[i,]), by = list(meta_data$NEC_NET), FUN = sum)$x[2]
  res_t <<- rbind( res_t, 
    matrix( c(nec,net),ncol = 2  )
  )
  rownames(res_t)[i] = rownames(gene_mat)[i]
}
colnames(res_t) = c("NEC", "NET")
res_t = as.data.frame(res_t)

dif_vec = res_t$NEC - res_t$NET
res_t_f = res_t[ order(as.integer( dif_vec), decreasing = T), ]
genes_of_interest_hgnc = rownames(res_t_f)[1:50]

#write.table(res_t_f,"~/MAPTor_NET/Results/Mutation_data/Genes_differentially_mutated_NEC_NET.tsv",quote = F, sep ="\t")

dd = data.frame( cbind( meta_data$MEN1, meta_data$NEC_NET) )
table(dd)
