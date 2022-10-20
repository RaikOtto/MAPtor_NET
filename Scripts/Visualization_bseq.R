meta_info$Location = str_replace_all(meta_info$Location, pattern = "-", "_")

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
colnames(bam_counts) = str_replace(colnames(bam_counts), pattern = "^X", "")
genes_of_interest_hgnc_t = read.table("~/MAPTor_NET/Misc/CPDB_PR.tsv", sep = "\t", header = T, stringsAsFactors = F)[-c(11,10,4,1,2),]

meta_match = match(colnames(bam_counts),meta_info$Name, nomatch = 0)
col_map = match( colnames(bam_counts), meta_info$Name, nomatch = 0)
df = meta_info[col_map,];df[is.na(df)] = -1;rownames(df) = meta_info$Name[col_map];df = df[,colnames(df) != "Name"]
df$Status=df$Subtyp_2;
#df$Tissue=df$Subtype;df$Tissue[df$Tissue == "Beta"] = "CCL";df$Tissue[df$Tissue == "Undefined"] = "CCL";df$Tissue[df$Tissue == "Delta"] = "CCL"
df_map = match(rownames(df),colnames(bam_counts));expr_raw_scaled = scale(expr_raw)
df$MEN1_exp = expr_raw_scaled[ rownames(expr_raw_scaled) == "MEN1", df_map]
df$MEN1_exp_log = as.double(bam_counts[ rownames(bam_counts) == "ENSG00000133895", df_map])

df$ABL1_exp = expr_raw_scaled[ rownames(expr_raw_scaled) == "ABL1", df_map]
df$ABL1_exp_log = as.double(bam_counts[ rownames(bam_counts) == "ENSG00000097007", df_map])
df$ABL2_exp = expr_raw_scaled[ rownames(expr_raw_scaled) == "ABL2", df_map]
df$ABL2_exp_log = as.double(bam_counts[ rownames(bam_counts) == "ENSG00000143322", df_map])

df_map = match(rownames(df), colnames(bam_counts))
df$Marker_genes = log(colSums( ( 
  bam_counts[ which( rownames(bam_counts) %in% 
                       c("ENSG00000254647","ENSG00000115263","ENSG00000108849","ENSG00000157005"))
              ,df_map] ) ) )


df$Status = meta_info$Subtyp_2[col_map]
#df$Tissue=df$Subtype;df$Tissue[df$Tissue == "Beta"] = "CCL";df$Tissue[df$Tissue == "Undefined"] = "CCL";df$Tissue[df$Tissue == "Delta"] = "CCL";df$Tissue[df$Tissue == ""] = "Pancreas";df$Tissue[df$Tissue == ""] = "Pancreas"
df$Chemotherapy[df$Chemotherapy == ""] = "Missing"
df$Chemotherapy[df$Chemotherapy == "Unknown"] = "Missing"

aka3 = list(
    Chemotherapy = c("Missing" = "gray", "No" =  "green", "Yes" = "red" ),
    Histology   = c(
        Pancreatic_NEN = "Red",
        Colorectal_NEN = "Orange",
        Small_intestinal_NEN = "Yellow",
        Gastric_NEN = "purple",
        Liver = "Darkgreen",
        CUP = "pink"),
    Location = c(
        Primary = "red",
        Liver_Met = "Darkgreen",
        Normal = "White",
        Lymph_node_Met = "Yellow",
        Peritoneum_Met = "Black",
        Spleen_Met = "Cyan",
        Lung_Met = "Blue"),
    NEC_NET = c(NEC= "Brown", NET = "darkgreen"),
    Study = c(Groetzinger = "Green", Scarpa = "purple"),
    MKI67 = c(high = "White", medium = "gray", low = "black"),
    Purity = c(high = "White", medium = "gray", low = "Blue"),
    Correlation = c(high = "Red", medium = "Yellow", low = "Green"),
    FGFR3 = c(high = "White", medium = "gray", low = "Black"),
    FGFR4 = c(high = "White", medium = "gray", low = "Black"),
    MEN1 = c("0" = "White", "1" = "White", "-1" = "Black"),
    RICTOR = c("0" = "Gray", "1" = "White", "-1" = "Black"),
    MYC = c("0" = "Gray", "1" = "White", "-1" = "Black"),
    ATRX = c("0" = "Gray", "1" = "White", "-1" = "Black"),
    DAXX = c("0" = "Gray", "1" = "White", "-1" = "Black"),
    CDKN2A = c("0" = "Gray", "1" = "White", "-1" = "Black"),
    Cellularity = c(high = "White", medium = "gray", low = "Black"),
    Marker_Genes = c(high = "White", medium = "Yellow", low = "Red"),
    Functionality = c( Unknown = "White",Functional = "green", Non_Functional="red"),
    Grading = c( G1 = "Green",G2 = "Yellow", G3 = "Red"),
    Included = c(Yes = "green", No = "red"),
    Chemotherapy = c( "1" = "red", "0" = "green", Missing = "gray"),
    Significance_Sadanandam = c(Sig = "green", Not_sig = "red"),
    Subtype_Sadanandam = c(Norm = "darkgreen", Insulinoma = "Blue", MLP = "Orange", Intermediate = "brown")
)

grading_col_vec = df$Grading
grading_col_vec[grading_col_vec == "No_tumor"] = rep("white", sum(grading_col_vec=="No_tumor"))
grading_col_vec[grading_col_vec == "G1"] = rep("green", sum(grading_col_vec=="G1"))
grading_col_vec[grading_col_vec == "G2"] = rep("yellow", sum(grading_col_vec=="G2"))
grading_col_vec[grading_col_vec == "G3"] = rep("red", sum(grading_col_vec=="G3"))
