x = eislet
clusters = "Subtype"
samples = "SampleName"
limit = TRUE
ct.scale = FALSE
markers = pancreasMarkers
function (x, markers, clusters, samples, ct.scale = TRUE, limit = TRUE) 
{
  if (ct.scale) {
    x <- cpm_cell_type(x, clusters = clusters, samples = samples)
  }
  ids <- intersect(unlist(markers), rownames(x))
  x <- x[ids, , drop = FALSE]
  clusters <- as.character(xbioc::pVar(x, clusters))
  samples <- as.character(xbioc::pVar(x, samples))
  res <- sapply(unique(clusters), function(ct) {
    rowMeans(sapply(unique(samples), function(sid) {
      y <- exprs(x)[, clusters %in% ct & samples %in% sid, 
                    drop = FALSE]
      rowMeans(y)
    }), na.rm = TRUE)
  })
  if (limit) 
    res <- res[, names(markers)]
  res
}