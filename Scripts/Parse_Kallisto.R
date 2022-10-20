### prep
library("sleuth")
library("biomaRt")
library("grid")
library("stringr")

#meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
"
tx2gene <- function(){
mart <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')

  t2g <- biomaRt::getBM(
    attributes = c(
        #'ensembl_transcript_id_version', 
        #'ensembl_gene_id',
        #'external_transcript_name',
        'illumina_humanht_12_v4',
        'external_gene_name'),
    mart = mart)
  t2g <- dplyr::rename(t2g,
     target_id = ensembl_transcript_id_version,
     ens_gene = ensembl_gene_id, 
     ext_gene = external_gene_name
      )
  return(t2g)
}
t2g <- tx2gene()
"

mart = biomaRt::useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'ensembl.org'
)
t2g = biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name"),
    mart = mart, useCache = FALSE)

t2g = dplyr::rename(
  t2g,
  target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id,
  ext_gene = external_gene_name)

#
Study_var = "Sato"

#meta_info_master = meta_info[meta_info$Study == Study_var,]

i_file_names_old = list.files( "~/MAPTor_NET/BAMs_new/Master/", recursive = T, pattern ="abundance.h5",full = F)
#i_file_names_old = list.files( "~/Deko_Projekt/Data/JGA/", recursive = T, pattern ="abundance.h5",full = F)
i_file_names_old = str_replace_all(i_file_names_old,pattern = "/abundance.h5","")
i_file_names_old = str_replace_all(i_file_names_old,"-","_")

i_files = list.files( "~/MAPTor_NET/BAMs_new/Master/", recursive = T, pattern ="abundance.h5",full = T)
i_file_names = list.files( "~/MAPTor_NET/BAMs_new/Master/", recursive = T, pattern ="abundance.h5",full = F)
#i_files = list.files( "~/Deko_Projekt/Data/JGA/", recursive = T, pattern ="abundance.h5",full = T)
#i_file_names = list.files( "~/Deko_Projekt/Data/JGA/", recursive = T, pattern ="abundance.h5",full = F)
i_file_names = str_replace_all(i_file_names,pattern = "/abundance.h5","")
i_file_names = str_replace_all(i_file_names,"-","_")

i_file_names_old[!(i_file_names_old %in% i_file_names)]

#i_file_names = str_replace_all(i_file_names,pattern = paste(study,"",sep = "/"),"")
matcher = match( i_file_names , meta_info$Raw_name )
i_file_names = meta_info$Sample[ matcher]

study_vec = rep(Study_var,length(i_file_names))#meta_info$Study[ match(i_file_names, meta_info$Raw_data_name ) ]
s2c = data.frame(path=i_files, sample=i_file_names,  Study = study_vec,stringsAsFactors = F)

so = sleuth_prep(s2c[,], target_mapping = t2g, aggregation_column = 'ens_gene', num_cores = 1, gene_mode = T)

# TPMs
so_1_21=so
#so_34 = so

dd = sleuth_to_matrix(so, "obs_norm","tpm")
dd[1:5,1:5]
dim(dd)

write.table(dd,"~/MAPTor_NET/BAMs_new/Sato_S21.ENSG.tsv",sep="\t",row.names = T,quote =F)

###

expr_raw[1:5,1:5]
rownames(expr_raw) = expr_raw[,1]
#expr_raw = expr_raw[,-1]
matcher = match(rownames(expr_raw),t2g$target_id,nomatch = 0)
expr_raw = expr_raw[matcher != 0,]
rownames(t2g) = t2g$target_id
mapping_t = t2g[rownames(expr_raw),]
gene_ids = mapping_t$ens_gene
hgnc_list = mapping_t$ext_gene

###

col_names = colnames(expr_raw)
row_names = expr_raw[-1,1]
expr_raw = expr_raw[,-1]
colnames(expr_raw) = col_names[1:3]
expr_raw[1:5,]
expr_raw = expr_raw[1:length(row_names),]
rownames(expr_raw) = row_names

####
library("tximport")

meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample

library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
k = keys(txdb, keytype = "TXNAME")
tx2gene = select(txdb, k, "GENEID", "TXNAME","HGNC_SYMBOL")
tx2gene = tx2gene[ !is.na(tx2gene[,2]), ]

library("RMariaDB") 
library("GenomicFeatures") 

txdb <-makeTxDbFromUCSC(genome="hg38", tablename="refGene")
k <- keys(txdb, keytype="TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

input_folder = "~/Downloads/quant/"

file_paths = list.files(input_folder,recursive = TRUE,pattern = "quant.sf",full.names = TRUE)

file_names = list.files(input_folder,recursive = TRUE,pattern = "quant.sf",full.names = FALSE)

file_names = str_replace(file_names,pattern = "/quant.sf","")
file_names = str_replace(file_names,pattern = "-","_")
file_names = str_replace(file_names,pattern = "transcripts_quant_","")
file_names = str_replace(file_names,pattern = "_R$","")

matcher = match(file_names, meta_info$Raw_name, nomatch = 0)
#matcher = match(file_names, meta_info$Raw_name, nomatch = 0)
meta_data = meta_info[matcher,]
file_names = meta_data$Sample

names(file_paths) = file_names
txi.salmon <- tximport(
  file_paths,
  type = "salmon",
  tx2gene = t2g,
  countsFromAbundance = "scaledTPM",
  ignoreTxVersion = TRUE
  )

head(txi.salmon$counts)
expr_raw = txi.salmon$counts

library("biomaRt")
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
res <- getBM(attributes=c("hgnc_symbol", "entrezgene_id"),mart = ensembl)
res = as.data.frame(res)
res = res[!(is.na(res[,2])),]
matcher = match(rownames(expr_raw), res$entrezgene_id, nomatch = 0)
hgnc_list = res[matcher,1]

### variance selection