library("stringr")
library("ggplot2")

s_match = match( as.character( colnames(expr_raw)), as.character( meta_info$Name), nomatch = 0)
meta_data = meta_info[s_match,]
rownames(meta_data) = meta_data$Name

b_match = match(rownames(expr), rownames(balanced.centroid), nomatch = 0) 
cent_clust = t( balanced.centroid[ b_match,] )
clust_data = expr[b_match != 0,]

### GOI

GOI = c("SPRR3","KRT13","KRT1","CCL19","CXCL9","CCR7","VNN2")

dd = as.data.frame( t( expr[rownames(expr) %in% GOI,] ) )
dd$Subtype = meta_data$Subtype
dd2 = reshape2::melt( dd, value.name = "Subtype" )
colnames( dd2 )  = c( "Subtype", "Gene", "Exp" )

ggplot( dd2, aes( Gene, Exp, fill = Subtype ) ) + geom_boxplot( position = "dodge" ) + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "top")

###

data = table(meta_data[,c("Subtype","Best_response")])
fisher.test(data)

###

data = aggregate(meta_data$Chemozyklen, by = list(meta_data$Subtype), FUN = c)
t.test( 
  as.double( as.character(unlist(data[ data$Group.1 == "BA",2] ) ) ),
  as.double( as.character(unlist(data[ data$Group.1 == "MS",2] ) ) )
)

library(ggplot2)

vis_mat = meta_data
vis_mat
colnames(vis_mat) = c("Name","MEN1_exp","MEN1_mt_AF")
vis_mat = reshape2::melt(vis_mat, id = c("Name","MEN1_exp","MEN1_mt_AF"))
vis_mat$Name = factor(vis_mat$Name, levels = vis_mat$Name[order(vis_mat$MEN1_exp)] )

vis_mat$MEN1_mt_AF[vis_mat$MEN1_mt_AF < .5] = 0
vis_mat$MEN1_mt_AF[ (.5 <= vis_mat$MEN1_mt_AF) & ( vis_mat$MEN1_mt_AF < .8)] = .5
vis_mat$MEN1_mt_AF[vis_mat$MEN1_mt_AF >= .8] = 1.0

# Basic barplot
label_vec = meta_data$NEC_NET
label_vec = label_vec[order(meta_data$MEN1)]
label_vec[label_vec == "NET"] = "T"
label_vec[label_vec != "T"] = "C"
col_vec = label_vec
col_vec[col_vec == "T"] = "darkgreen"
col_vec[col_vec != "darkgreen"] = "brown"

p = ggplot( data = vis_mat)
p = p + geom_bar(aes( x = Name, y = MEN1_exp, fill = MEN1_mt_AF ),stat="identity", colour="black")
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + scale_fill_gradientn(colours = c("white","yellow","red"), breaks = c(0.0,.5,1.0))
p = p + annotate("text", x=1:42,y = 5.7,parse=TRUE, label = label_vec, color = col_vec, size = 4.5 )
p = p + xlab("") + ylab("MEN1 expression in log TPM") + theme(legend.position = "top")
p

vis_mat = reshape2::melt(meta_data[,c("Subtype","Best_response")])
ggplot( vis_mat, aes( Subtype, Best_response ) ) + geom_bar(aes( x = Subtype, y = Best_response, fill = Best_response),stat="identity" ) + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "top")

sub_vec = meta_data$Subtype[s_match]
r_mat[1:5,1:5]

exp_mat = apply(
  r_mat[ , !( colnames(r_mat) %in% c("Sample", "OS") )], 
  FUN = function(vec){
     return(
        aggregate(vec, by = list(sub_vec) , c)
     )
  },
  MARGIN = 1
)


#


library(ggplot2)

meta_match = match( colnames(pure_data), meta_data$Name )
meta_data$CXCL10 = pure_data[ rownames(pure_data) == "CXCL10" ]
meta_data$CXCL9 = pure_data[ rownames(pure_data) == "CXCL9" ]
meta_data$STAT1 = pure_data[ rownames(pure_data) == "STAT1" ]
meta_data$IFNG = pure_data[ rownames(pure_data) == "IFNG" ]
meta_data$E6 = pure_data[ rownames(pure_data) == "E6 (HPV16)" ]
meta_data$CDKN2A = pure_data[ rownames(pure_data) == "CDKN2A" ]
meta_data$AREG = pure_data[ rownames(pure_data) == "AREG" ]

#meta_data$Subtype == "MS"
vis_mat = meta_data[ ,c("Name","Subtype","CDKN2A")]
vis_mat
#colnames(vis_mat) = c("Name","MEN1_exp","MEN1_mt_AF")
vis_mat = reshape2::melt(vis_mat, id = c("Name","Subtype","CDKN2A"))
vis_mat$Name = factor(vis_mat$Name, levels = vis_mat$Name[order(vis_mat$CDKN2A)] )

# Basic barplot
col_vec = meta_data$OS
col_vec[col_vec > mean(col_vec)] = "darkgreen"
col_vec[col_vec != "darkgreen"] = "brown"
label_vec = meta_data$OS
label_vec[label_vec > mean(label_vec)] = "A"
label_vec[label_vec != "A"] = "B"

#aggregate( vis_mat, by = list(vis_mat$Subtype), FUN =  )
p = ggplot( data = vis_mat)
p = p + geom_boxplot(aes( x = Subtype, fill = Subtype ),stat="identity", colour="black")
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#p = p + scale_fill_gradientn(colours = c("white","yellow","red"), breaks = c(0.0,.5,1.0))
#p = p + annotate("text", x=1:ncol(pure_data),y = 15,parse=TRUE, label = label_vec, size = 4.5, angle = 0, color = col_vec )
p = p + xlab("") + ylab("CDKN2A") + theme(legend.position = "top")
p

###


library("tidyr")
library("ggplot2")
library("stringr")

# Daten werden geladen, kann dauern

res_count = pure_data
# DEFINE GENES OF INTEREST HERE

des_genes = c("AREG","EGFR")

#pan_genes = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Stem_signatures.gmt",sep ="\t", header = F,stringsAsFactors = F)
#pan_genes[4,3:ncol(pan_genes)] = str_to_upper(pan_genes[4,3:ncol(pan_genes)])
#write.table(pan_genes,"~/MAPTor_NET/BAMs/Kallisto_three_groups/Stem_signatures.gmt",sep="\t",quote=F,row.names = F, col.names = F)

data_mat = pure_data[ rownames(pure_data) %in% des_genes,]

# DEFINE SAMPLE COHORTS HERE

cohort_one = meta_data$Name[meta_data$Subtype == "BA"]
cohort_two = meta_data$Name[meta_data$Subtype == "CL"]
cohort_three = meta_data$Name[meta_data$Subtype == "MS"]

###

cohort_one_matrix = data_mat[,colnames(data_mat) %in% cohort_one]
gather_matrix_one <- reshape2:::melt.matrix(as.matrix(cohort_one_matrix))
gather_matrix_one$Cohort = "BA"
cohort_two_matrix = data_mat[,colnames(data_mat) %in% cohort_two]
gather_matrix_two <- reshape2:::melt.matrix(as.matrix(cohort_two_matrix))
gather_matrix_two$Cohort = "CL"
cohort_three_matrix = data_mat[,colnames(data_mat) %in% cohort_three]
gather_matrix_three <- reshape2:::melt.matrix(as.matrix(cohort_three_matrix))
gather_matrix_three$Cohort = "MS"

gather_matrix = rbind(gather_matrix_one,gather_matrix_two,gather_matrix_three)

colnames(gather_matrix) = c("Gene","Sample","Expression","Subtype")

### Plot

tidy_dist = ggplot( data = gather_matrix, aes( Gene, Expression ) )
tidy_dist = tidy_dist  + geom_boxplot(aes(fill = Subtype))
tidy_dist = tidy_dist + ylab("Expression") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
tidy_dist = tidy_dist + guides(fill=guide_legend(title="Subtypes"))
tidy_dist

### Pheatmap

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

vis_mat = data_mat[rownames(data_mat) %in% ensembl_list,]
rownames(vis_mat) = des_genes
pheatmap::pheatmap(
  vis_mat,
  annotation_col = meta_data[,c("Histology","Study")]
)

qplot(meta_data$OS, geom="histogram", bins = 8)
q = qplot(meta_data$OS, geom="histogram", bins = 8)
