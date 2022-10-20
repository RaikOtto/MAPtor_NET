library("stringr")

path_stem = "~/MAPTor_NET/RAW/Groetzinger/"

meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep ="\t",header = T)
meta_info = meta_info[meta_info$Mask != "",]
i_samples = as.character(meta_info$Raw_data_name)
o_samples = as.character(meta_info$Mask)

files_full = list.files("~/MAPTor_NET/RAW/Groetzinger/",full.names = T,recursive = T,pattern = ".gz")
files_full = str_replace_all(files_full, pattern ="-","_")
files_short = as.character(sapply(
  files_full,
  FUN = function(vec){
    return(tail(as.character(unlist(str_split(vec,pattern = "/"))),1))})
)

for ( i  in 1:length(files_full)) {
    
    match_index = which( str_detect( files_full[i], i_samples ) == TRUE )
    replacement = o_samples[match_index]
    
    if ( length(replacement) == 0)
        next
    
    real_name = str_replace_all(files_full[i], pattern = "\\_LR\\_", "-LR-")
    real_name = str_replace_all(real_name, pattern = "/AS\\_", "/AS-")
    remainder = str_replace(files_short[i],pattern = i_samples[match_index], replacement = "")
    new_file = paste(c(path_stem, replacement,remainder), collapse = "")
    
    print(replacement)
    print (i)
    
    # file.rename(real_name,new_file)
}

