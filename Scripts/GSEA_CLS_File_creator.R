library("stringr")

meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
sort(table(meta_info$Sample),decreasing = TRUE)

path_transcriptome_file = "~/MAPTor_NET/BAMs_new/Publication_datasets/Fröhling.S34.HGNC.DESeq2.VOOM.tsv"
expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, as.is = F,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
nomatch = match(colnames(expr_raw),meta_info$Sample,nomatch =  0)
colnames(expr_raw)[ nomatch == 0] = str_c("X",colnames(expr_raw)[ nomatch == 0])

table(match(colnames(expr_raw), meta_info$Sample, nomatch = 0) != 0)

meta_data = meta_info[colnames(expr_raw),]
candidates = meta_data$Sample[ 
  #meta_data$Site_of_primary %in% c("Pancreatic")
]
expr = expr_raw[,candidates]
meta_data = meta_info[colnames(expr),]
dim(meta_data)

### cls file

cohort_vec = meta_data$NET_NEC_UMAP
#cohort_vec = meta_data$Cohort

cohort_vec = str_replace(cohort_vec," ","_")
table(cohort_vec)

first_entity = unique(cohort_vec)[1]
second_entity = unique(cohort_vec)[2]

first_line = as.character(c(ncol(expr),length(unique(cohort_vec)),1))
first_line[4:(length(cohort_vec))] = ""
second_line = as.character(c("#",first_entity,second_entity))
second_line[4:(length(cohort_vec))] = ""

cls_file = rbind(first_line,second_line,cohort_vec)

write.table(cls_file, "~/MAPTor_NET/GSEA/Validation_Cohort.S34.DESeq2.tsv",row.names = F, quote = F, sep ="\t",col.names = FALSE)
write.table(cls_file, "~/MAPTor_NET/GSEA/Validation_Cohort.S34.DESeq2.tsv.cls",row.names = F, quote = F, sep ="\t",col.names = FALSE)

### gct file

first_line = "#1.0"
first_line[2:(ncol(expr)+2)] = ""
second_line = as.character(c(nrow(expr),ncol(expr)))
second_line[3:(ncol(expr)+2)] = ""
third_line = c("Name","Description",colnames(expr))
gct_matrix = cbind(rownames(expr),rownames(expr),expr)

output_gct = rbind(first_line,second_line, third_line, gct_matrix)
output_gct[1:5,1:5]

write.table(output_gct, "~/MAPTor_NET/GSEA/Validation_Cohort.S34.DESeq2.gct",row.names = F, quote = F, sep ="\t",col.names = FALSE)
write.table(output_gct, "~/MAPTor_NET/GSEA/Validation_Cohort.S34.DESeq2.gct.tsv",row.names = F, quote = F, sep ="\t",col.names = FALSE)

####

percentages_t = as.data.frame(read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.with_NENs.tsv", header = TRUE, sep ="\t"))
percentages_t = percentages_t %>% dplyr::filter( Model != "Endocrine_only")
matcher = match(colnames(expr),percentages_t$Sample)
cohort_vec = cohort_vec_numeric = as.double(rowSums(percentages_t[ matcher, c("Ductal","Acinar")]))
cohort_vec[cohort_vec_numeric > mean(cohort_vec_numeric)] = "Exocrine-like-high"
cohort_vec[cohort_vec_numeric <= mean(cohort_vec_numeric)] = "Exocrine-like-low"

#write.table(output_gct, "~/Deko_Projekt/GSEA/Charite_Diedisheim_Riemer_Scarpa_S134.Exocrine.HGNC.cls.tsv",row.names = F, quote = F, sep ="\t",col.names = FALSE)
#write.table(output_gct, "~/Deko_Projekt/GSEA/RepSet_S51.Exocrine.HGNC.cls",row.names = F, quote = F, sep ="\t",col.names = FALSE)
