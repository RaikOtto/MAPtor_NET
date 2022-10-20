library("ggpubr")
library("stringr")
library("reshape2")
library("ggplot2")
library("dplyr")
library("grid")

#443 407 444 579 450 409 452 PNET08

#draw_colnames_45 <- function (coln, gaps, ...) {
#  coord = pheatmap:::find_coordinates(length(coln), gaps)
#  x = coord$coord - 0.5 * coord$size
#  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
#  return(res)}
#assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

meta_info = read.table("~/MAPTor_NET/Misc/Meta_information_FLT3L.tsv",sep = "\t",header = T,stringsAsFactors = F)
#meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

expr_raw = read.table("~/MAPTor_NET/BAMs_new/RepSet.S54.Flt3L.HGNC.DESeq2.tsv",sep="\t", stringsAsFactors =  F, header = T)
#expr_raw = read.table("~/MAPTor_NET/BAMs_new/RepSet.S54.Flt3L.HGNC.DESeq2.tsv",sep="\t", stringsAsFactors =  F, header = T)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X\\.)", "")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^(X)", "")
expr_raw[1:5,1:5]
#expr_raw = expr_raw[,str_detect( colnames(expr_raw), pattern = "_", negate = T)]
dim(expr_raw)

meta_info$Sample[which(!(meta_info$Sample %in% colnames(expr_raw)))]
expr_raw = expr_raw[,colnames(expr_raw) %in% meta_info$Sample]

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
matcher = match(colnames(expr_raw),meta_info$Sample,nomatch = 0)
colnames(expr_raw)[matcher == 0]
meta_data = meta_info[colnames(expr_raw),]
#"132502" %in% colnames(expr_raw)

source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
#genes_of_interest_hgnc_t = read.table("~/Deko_Projekt/Misc//Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t = read.table("~/MAPTor_NET//Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1

liver_genes = genes_of_interest_hgnc_t[70,3:ncol(genes_of_interest_hgnc_t)]
i = 58 #THBD 38 # active dendrites 58
genes_of_interest_hgnc_t[i,1]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
sad_genes[!(sad_genes %in% rownames(expr_raw))]
#sad_genes = sad_genes[!(sad_genes %in% liver_genes)]
length(sad_genes)

expr_mat = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw));colnames(expr_mat) = colnames(expr_raw);rownames(expr_mat) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
#expr_mat = expr_mat[,meta_data[meta_data$NEC_NET %in% "NEC","Sample"]]
expr = expr_mat
expr[1:5,1:5]
dim(expr)

###
correlation_matrix = cor(expr)
pcr = prcomp(t(correlation_matrix))

meta_data$FLT3LG = log(as.double(expr_raw["FLT3LG",rownames(meta_data)]))
meta_data[meta_data$Study == "Riemer","Study"] = "Charite"
meta_data$average_expression = apply(expr_mat,MARGIN = 2, FUN =mean)

meta_data$THBD_expression = log(meta_data$average_expression+1)
meta_data$DC_signature = log(meta_data$average_expression+1)

#svg(filename = "~/Downloads/Heatmap_DC.svg", width = 10, height = 10)
p  =pheatmap::pheatmap(
  correlation_matrix,
  #expr_mat,
  annotation_col = meta_data[,c("DC_signature","FLT3LG","Grading","Study")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  #treeheight_col = 0,
  treeheight_row = 0,
  legend = T,
  fontsize_col = 7,
  clustering_method = "average"
)
dev.off()

cohort_vec_dendrit = meta_data$average_expression
cohort_vec_dendrit[meta_data$average_expression > mean(meta_data$average_expression) + .1] = "high"
cohort_vec_dendrit[ cohort_vec_dendrit != "high" ] = "low"

cohort_vec_flt3lg = as.double(expr_raw["FLT3LG",])
cohort_vec_flt3lg[cohort_vec_flt3lg > mean(as.double(expr_raw["FLT3LG",])) ] = "high"
cohort_vec_flt3lg[ cohort_vec_flt3lg != "high" ] = "low"

# change colors
col_vec_flt3lg = cohort_vec_flt3lg
col_vec_flt3lg[col_vec_flt3lg == "high"] = "#FF4C33"
col_vec_flt3lg[col_vec_flt3lg != "#FF4C33"] = "#33ACFF"

col_vec_dendrite = cohort_vec_dendrit
col_vec_dendrite[col_vec_dendrite == "high"] = "orange"
col_vec_dendrite[col_vec_dendrite != "orange"] = "green"

size_vec_real = as.double(expr_raw["FLT3LG",]) #meta_data$average_expression^2*.5
size_vec_real = size_vec_real *1.2
size_vec_background = (size_vec_real)  * 1.4

p = ggbiplot::ggbiplot(
  pcr,
  obs.scale =.75,
  #labels = meta_data$Sample,
  #labels.size = 3,
  groups = cohort_vec_dendrit,
  ellipse = FALSE,
  circle = TRUE,
  var.axes = F
)
p = p + geom_point( size = size_vec_background, aes( color = "black"))
p = p + geom_point(size = size_vec_real,  aes( color = as.character(cohort_vec_flt3lg)))
p = p + stat_ellipse( linetype = 1,aes(color=as.character(col_vec_flt3lg)),level=.66, type ="t",size=1.25)
p = p + stat_ellipse( linetype = 2,aes(color=col_vec_dendrite),level=.66, type ="t",size=1.25)
p = p + scale_color_manual( values = c("#33ACFF","#FF4C33","black","darkgreen","#FF4C33","#33ACFF","orange"), name = "Gene expression act. dendritic cells" ) + theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=13))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p

### correlation flt3lg and immuno-score

estimate_result = read.table("~/MAPTor_NET/BAMs_new/RepSet_S57.HGNC.ESTIMATE.RESULTS.tsv",sep ="\t", header =T,as.is = T,skip = 2)
colnames(estimate_result) = str_replace_all(colnames(estimate_result), pattern = "^X","")
estimate_result = as.data.frame(t(estimate_result))
colnames(estimate_result) = estimate_result[1,]
estimate_result = estimate_result[c(-1,-2),]
e_mat = as.data.frame(matrix(as.double(as.character(unlist(estimate_result))),ncol = ncol(estimate_result),nrow = nrow(estimate_result)))
colnames(e_mat) = colnames(estimate_result)
rownames(e_mat) = rownames(estimate_result)
e_mat

meta_data = meta_info[rownames(e_mat),]
e_mat = e_mat[which(rownames(e_mat) %in% colnames(expr_raw)),]
expr_FLT3LG = expr_raw["FLT3LG",rownames(e_mat)]
cor(as.double(expr_FLT3LG), as.double(meta_data$ImmuneScore))
plot(as.double(expr_FLT3LG), as.double(meta_data$ImmuneScore))

vis_mat = as.data.frame(cbind(as.double(expr_FLT3LG),as.double(meta_data$ImmuneScore)))
colnames(vis_mat) = c("FLT3LG","ImmuneScore")

smooth_plot = ggplot(vis_mat,aes(x=ImmuneScore, y=FLT3LG)) + 
  geom_point()+
  geom_smooth(method=lm, color="black") +theme_classic()

library("gridExtra")
xdensity <- ggplot(vis_mat, aes(ImmuneScore)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
xdensity

# Marginal density plot of y (right panel)
ydensity <- ggplot(vis_mat, aes(FLT3LG)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )
grid.arrange(xdensity, blankPlot, smooth_plot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
