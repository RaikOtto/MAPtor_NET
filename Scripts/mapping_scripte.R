library("stringr")

mapping_t=read.table("~/MAPTor_NET/Misc/mapping_table_ega_real_name.tsv",sep ="\t",header = T)
mapping_t = mapping_t[mapping_t$Mask != "",]
data_t = read.table("~/Deko_Projekt/Results/Submission/Supplement/SM_Table_3_Deconvolution_P_values.tsv",sep ="\t",header = T)

rownames(mapping_t) = mapping_t$Sample
matcher = match(as.character(data_t$Sample),as.character(mapping_t$Sample),nomatch = 0)
data_t$mask = mapping_t[data_t$Sample,"Mask"]

write.table(data_t,"~/Deko_Projekt/Results/Submission/Supplement/SM_Table_3_Deconvolution_P_values.mask.tsv",sep = "\t",quote = F)
