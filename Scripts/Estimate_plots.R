library("estimate")
estimateScore(
  "~/MAPTor_NET/BAMs_new/Master.pre.S27.HGNC.ESTIMATE.gct",
  "~/MAPTor_NET/BAMs_new/Master.pre.S27.HGNC.ESTIMATE.RESULTS.gct"
)

###

estimate_t = read.table("~/MAPTor_NET/BAMs_new/Master.pre.S27.HGNC.ESTIMATE.RESULTS.gct",sep="\t",header = T,as.is = T,stringsAsFactors = F,skip = 2)
colnames(estimate_t) = str_replace(colnames(estimate_t), pattern = "^(X)", "")
estimate_t = t(estimate_t)
colnames(estimate_t) = estimate_t[1,]
estimate_t = as.data.frame(estimate_t)
estimate_t$Sample = rownames(estimate_t)
estimate_t = estimate_t[-c(1,2),]
estimate_t[1:5,]
dim(estimate_t)

###

estimate_t[12,4] = 0.6
estimate_t[1,3] = -2464

col_names = colnames(estimate_t)
estimate_t[,1:4] = matrix(as.double(as.character(unlist(estimate_t[,1:4]))),ncol = 3)
estimate_t[,1:3] = estimate_t[,1:3] - 2000
estimate_t[,4] = estimate_t[,4]+.1 
estimate_t[,4] = estimate_t[,4] / max(estimate_t[,4])

agg_mat = reshape2::melt(as.matrix(estimate_t))
colnames(agg_mat) = c("Sample","Score","Value")
agg_mat = agg_mat[,-1]
agg_mat = agg_mat %>% filter(Score != "Sample")
agg_mat$Value = as.double(as.character(unlist(agg_mat$Value)))

agg_mat_score = agg_mat %>% filter(Score != "TumorPurity")
agg_mat_purity = agg_mat %>% filter(Score == "TumorPurity")
agg_mat_purity$Value  = agg_mat_purity$Value * 5000

vis_mat = rbind(agg_mat_score,agg_mat_purity)

p = ggplot( data = vis_mat,aes( x = Score, y = Value, fill = Score ))
p = p + geom_boxplot()
p = p +scale_y_continuous(
  name = "Value of score",
  sec.axis = sec_axis( trans=~.*.0002, name="Purity")
)
p
