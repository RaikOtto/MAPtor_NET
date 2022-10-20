###install

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)

### contamination SM Figures

library("estimate")
library(stringr)

meta_info = read.table("~/Deko_Projekt//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

#expr_raw = read.table("~/Deko_Projekt/Data/Bench_data/Alvarez.S104.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = read.table("~/MAPTor_NET/BAMs_new/Publication_datasets/Discovery_Cohort.S64.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
#ncol = ncol(expr_raw)
#nrow = nrow(expr_raw)
#rownames = rownames(expr_raw)
#colnames = colnames(expr_raw)
#expr_raw = as.double(as.character(as.matrix(expr_raw)))
#expr_raw = round(expr_raw,1)
#expr_raw[expr_raw < 1 ] = 1
#expr_raw = matrix(expr_raw, ncol = ncol, nrow = nrow)
gct_mat = cbind( rownames(expr_raw),rownames(expr_raw), expr_raw )
new_line = rep("",ncol(gct_mat))
new_line[1] = "#1.2"
new_line_2 = rep("",ncol(gct_mat))
new_line_2[1] = nrow(expr_raw)
new_line_2[2] = ncol(gct_mat) - 2
new_line_3 = rep("",ncol(expr_raw))
new_line_3 = c("NAME","Description",new_line_3)
new_line_3[3:ncol(gct_mat)] = colnames(expr_raw)
gct_mat = matrix(as.character(gct_mat),ncol = ncol(expr_raw) + 2)
gct_mat = rbind(new_line,new_line_2,new_line_3,gct_mat)


gct_mat[1:5,1:5]
dim(gct_mat)

write.table(gct_mat,"~/MAPTor_NET/Results/Estimate_Discovery_S64.gct",sep="\t",quote = F, row.names = FALSE)

#setwd("~/R//estimate/extdata/")
estimateScore( 
  "~/MAPTor_NET/Results/Estimate_Discovery_S64.gct",
  "~/MAPTor_NET/Results/Estimate_Discovery_S64.estimate.tsv")

#estimate_mat = as.data.frame(read.table("~/MAPTor_NET/BAMs_new/Alvarez.S212.HGNC.estimate.tsv",row.names = 1, header = T, sep = "\t"))
estimate_mat = as.data.frame(read.table("~/MAPTor_NET/BAMs_new/RepSet_S84.HGNC.estimate.tsv",row.names = 1, header = T, sep = "\t"))
rownames(estimate_mat) = str_replace(rownames(estimate_mat), pattern = "^X", "")
estimate_mat[1:5,]
no_match = rownames(estimate_mat) %in% meta_info$Sample == F
rownames(estimate_mat)[no_match] = paste("X",rownames(estimate_mat)[no_match],sep ="")
no_match = rownames(estimate_mat) %in% meta_info$Sample == F
estimate_mat[no_match,]

meta_info[match(rownames(estimate_mat),rownames(meta_info)),"StromalScore"] = as.double(estimate_mat$StromalScore)
meta_info[match(rownames(estimate_mat),rownames(meta_info)),"ImmuneScore"] = as.double(estimate_mat$ImmuneScore)
meta_info[match(rownames(estimate_mat),rownames(meta_info)),"ESTIMATEScore"] = as.double(estimate_mat$ESTIMATEScore)
meta_info[match(rownames(estimate_mat),rownames(meta_info)),"TumorPurity"] = as.double(estimate_mat$TumorPurity)

#write.table(meta_info,"~/Deko_Projekt//Misc/Meta_information.tsv",sep="\t",quote = F, row.names = F)
