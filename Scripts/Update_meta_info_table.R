meta_info = read.table(
  "~/MAPTor_NET//Misc/Meta_information.tsv",
  sep = "\t",
  header = T,
  stringsAsFactors = F
)

row_names = rownames(count_data)
col_map = match( colnames(count_data), meta_info$Name , nomatch = 0 )
table( col_map == 0 )
meta_data = meta_info[ col_map,]

count_data = count_data[, colnames(count_data) %in% meta_data$Name]
meta_data  = meta_data[ meta_data$Name %in% colnames(count_data),]

raw_names = as.character( sapply(FUN=str_split,as.character( unlist(colnames(count_data)) ), "/")) 
clean_name = as.character( sapply(raw_names,FUN = function(name){ 
  return(
    str_replace(  
      tail(
        unlist(
          as.character( unlist(name) )
        ),
        1
      )
      , ".hg19.bam",
      ""
    )
  )  
}))

colnames(count_data) = as.character(clean_name)
meta_data$SampleName = as.character(clean_name)
