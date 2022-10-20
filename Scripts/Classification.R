### init

centroid2expr = function( centroid, vd ){
  
  gene.sig = intersect( rownames( centroid ), rownames( vd ) )
  vd = t( scale( t( vd[ gene.sig, ] ) ) )
  centroid = centroid[ gene.sig, ]
  vclass = c()
  vcor = c()
  
  for( i in 1:ncol( vd ) ){

    d = vd[,i]
    c.cor = c()
    pv = c()
    
    for( j in colnames( centroid ) ){

      centroidj = centroid[ , j ]
      corj = cor.test( centroidj, d, use = "complete", method = "pearson" )
      c.cor[ j ] = corj$estimate
      pv[j]=corj$p.value
    }
    
    maxk = which.max(c.cor)
    pub_cor <<- rbind(pub_cor,c.cor)
    group = names( maxk )
    vcor = rbind( vcor, c( colnames( vd )[ i ], group, c.cor[ maxk ], pv[ maxk ] ) )
    
    if( pv[ maxk ] < .05 ){
      vclass[colnames(vd)[i]]=group
    }
  }
  
  
  return( list( overlap.gene = gene.sig, cluster = vclass, correlation = vcor ) )
}

final_centroid_list = read.table(
  "~/MAPTor_NET//Results/Classification/Centroids_IT3.tsv",
  sep = "\t",
  header = T,
  stringsAsFactors = F
)
balanced.centroid = as.data.frame(final_centroid_list)[,c(-1,-2, -1 * ncol(final_centroid_list))]
rownames(balanced.centroid) = as.character(final_centroid_list[,2])
balanced.centroid

variance_selection = function(  elem, input_matrix ){
  
  mapping = rownames(input_matrix) == as.character( elem )
  exprs_vals = matrix( log2( as.double( unlist( input_matrix[ mapping, ] ) ) ), ncol = dim(input_matrix)[2] )
  exprs_var  = apply( exprs_vals, MARGIN = 1, FUN = var  )
  selector   = which( exprs_var == max( exprs_var )  )[1]
  
  max_var_gene = matrix(
    as.double( exprs_vals[ selector, ] ) , 
    ncol = dim( input_matrix )[2]
  )
  colnames(max_var_gene) = colnames(input_matrix)
  
  return( max_var_gene )
}

m_t = count_data
colnames(m_t) = str_replace(colnames(m_t),pattern = "^X","")
colnames(m_t) = str_replace(colnames(m_t),pattern = "\\.","_")
m_t = m_t[rowMeans(m_t) >= 1,]
m_t = log(m_t)# LOG
m_t[1:5,1:5]
genes = rownames(m_t)
genes = as.character(sapply(genes, FUN = function(gene){return(  as.character(unlist(str_split(gene,"__")))[1] )}))

row_var = apply( m_t, FUN = var, MARGIN = 1 )
col_var = apply( m_t, FUN = var, MARGIN = 2 )
m_t = m_t[row_var != 0, col_var != 0]
dim(m_t)

ensembl_list = as.character( mapIds(
  org.Hs.eg.db,
  keys = genes,
  column="ENSEMBL",
  keytype="SYMBOL",
  multiVals="first"
) )

m_t = m_t[!is.na(ensembl_list),]
ensembl_list = ensembl_list[!is.na(ensembl_list)]
rownames(m_t) = ensembl_list

count_data_centroid = m_t
count_data_centroid[1:5,1:5]

en_r = rownames(balanced.centroid)
type_c = colnames(balanced.centroid)
balanced.centroid = t(apply( balanced.centroid, FUN = function(vec){return((as.double(unlist(vec))))} , MARGIN = 1 ))
balanced.centroid = as.data.frame(balanced.centroid)
colnames(balanced.centroid) = type_c
rownames(balanced.centroid) = en_r

# filter for unwanted samples

pub_cor <<- matrix( as.double(), ncol = ncol( balanced.centroid ) )
expr2bc = centroid2expr( balanced.centroid[1:200,], count_data_centroid )
groups = as.character( expr2bc$correlation[,2] )
groups = str_replace(groups, pattern = "-score", "")
colnames(expr2bc$correlation) = c("ID","Group","Cor","PValue")

p_value = as.double(as.character(unlist(expr2bc$correlation[,4])))
q_value = p.adjust(p_value,"BH")
correlation = expr2bc$correlation[,3]
expr2bc$correlation = cbind( expr2bc$correlation , q_value)
expr2bc$correlation[,4] = p_value

groups[p_value < 0.05] = "not_sig"

classification_t    = data.frame(
  "Sample_name" = colnames(count_data_centroid),
  "Literature" =  as.character(meta_info$Subtype)[match( as.character(colnames(count_data_centroid)),as.character( meta_info$Name) )],
  "Classification" =  as.character(groups),
  "Cor" = as.character(correlation),
  "P_value" = as.character(p_value),
  "Q_value" = as.character(q_value)
)
classification_t[1:5,]

write.table(
    classification_t,
    file =  "~/MAPTor_NET//Results/Classification_Lawlor.tsv",
    sep = "\t",
    row.names = F,
    quote = F
)

info_map = match( classification_t$Sample_name, meta_info$Name , nomatch = 0 )
meta_info$Classification[info_map] =as.character( classification_t$Classification[ info_map != 0])
meta_info$P_value[info_map] =as.character( classification_t$P_value[ info_map != 0])
write.table(meta_info, "~/MAPTor_NET//Misc/Meta_information.tsv",sep ="\t", row.names =T, quote =F)

source("~/MAPTor_NET///Scripts/Update_meta_info_table.R")

### extended part for roc curve

centroid2expr_all = function( centroid, vd ){
  
  result_matrix <<- matrix(as.character(), ncol = ncol(centroid))
  
  gene.sig = intersect( rownames( centroid ), rownames( vd ) )
  vd = t( scale( t( vd[ gene.sig, ] ) ) )
  centroid = centroid[ gene.sig, ]
  vclass = c()
  vcor = c()
  
  for( i in 1:ncol( vd ) ){
    
    d = vd[,i]
    c.cor = c()
    pv = c()
    
    for( j in colnames( centroid ) ){
      
      centroidj = centroid[ , j ]
      corj = cor.test( centroidj, d, use = "complete", method = "pearson" )
      c.cor[ j ] = corj$estimate
      pv[j]=corj$p.value
    }
    
    result_matrix = rbind( result_matrix, pv )
  }
  rownames(result_matrix) = colnames(count_data_centroid)
  
  return( result_matrix )
}

pub_cor <<- matrix( as.double(), ncol = ncol( balanced.centroid ) )
expr2bc_all = centroid2expr_all( balanced.centroid[1:200,], count_data_centroid )
expr2bc_all = cbind(classification_t[,1:3], expr2bc_all)
expr2bc_all[1:5,]

write.table(
  classification_t,
  file =  "~/MAPTor_Net_RNA_data/Results/Classification/Classification_muraro_all_p_values.tsv",
  sep = "\t",
  row.names = F,
  quote = F
)
