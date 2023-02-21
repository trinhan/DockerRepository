DMR_analysis <- function(
  sample1_beta,
  sample1_name,
  sample2_beta,
  sample2_name,
  row_ranges,
  min_sample_per_probe,
  plot_dir
) {

  print(paste0("Detecting differentially methylated probes between ",sample1_name, 
    " and ", sample2_name, "..." ))

	# fetch beta dfs:
  unfilt_mtx <- list(sample1_beta = as.matrix(sample1_beta),
    sample2_beta = as.matrix(sample2_beta) )
  
  # remove rows with < min_sample_per_probe samples not NA:
  filt_mtx <- lapply(unfilt_mtx, function(x) {
    return(x[apply(x, 1, function(y) {
          length(which(!is.na(y))) >= min_sample_per_probe
        }), ])
  })
  
  # keep only rows in both, and merge:
  filt_mtx$sample1_beta <- filt_mtx$sample1_beta[
    rownames(filt_mtx$sample1_beta) %in% rownames(filt_mtx$sample2_beta), ]
  filt_mtx$sample2_beta <- filt_mtx$sample2_beta[
    rownames(filt_mtx$sample2_beta) %in% rownames(filt_mtx$sample1_beta), ]
  beta_mtx <- do.call("cbind", filt_mtx)
  
  # adjust other objects for new se:
  row_ranges <- row_ranges[match(rownames(beta_mtx), names(row_ranges))]
  
  dm_se <- SummarizedExperiment(assays = list(counts = beta_mtx),
    rowRanges = row_ranges, colData = DataFrame(
      tissue = gsub("^.*_", "", colnames(beta_mtx)) ))
  
  res <- list(
    DMR = TCGAanalyze_DMC(
      data = dm_se,
      groupCol = "tissue",
      group1 = sample1_name,
      group2 = sample2_name,
      plot.filename = paste0(
      	plot_dir, sample1_name, "_vs_", sample2_name,"_DMR_volcano.png" ),
      p.cut = 1 ),
    filtered_probes = nrow(beta_mtx)
  )

}


do_PCA <- function(pca_in, col_vec, shape_vec) {
   
  pca_in$tissue <- rownames(pca_in)
  pca_in$tissue <- gsub("_.*_", "_", pca_in$tissue)

  # create source/type columns:
  pca_in$source <- gsub("_.*$", "", pca_in$tissue)
  pca_in$type <- gsub("^.*_", "", pca_in$tissue)

  # order by type:
  pca_in$type <- factor(pca_in$type, 
    levels = c(unique(grep("MPNST", pca_in$type, value = T)), 
    "pNF", "blood", "nonmalig", "malig" ))
  pca_in <- pca_in[order(pca_in$type),]

  pca_in$tissue <- factor(pca_in$tissue, levels = unique(pca_in$tissue))

  # adjust sizes of vectors:
  shape_vec <- shape_vec[seq_along(levels(pca_in$tissue))]
  size_vec <- rep(1.5, length(levels(pca_in$tissue)))
  col_vec <- col_vec[levels(pca_in$tissue)]

  # reverse order for plotting:
  pca_in$type <- factor(pca_in$type, 
    levels = c("blood", "nonmalig", "malig", 
      unique(grep("MPNST", pca_in$type, value = T)), "pNF" ))
  pca_in <- pca_in[order(pca_in$type),]

  # plot PCA, adjusting transparency if needed:
  p <- ggplot(pca_in, aes(x=PC1, y=PC2, shape=source, color=tissue, size=tissue ))
  p <- p + geom_point()
  p <- p + scale_color_manual(values = col_vec, labels = levels(pca_in$tissue))
  p <- p + scale_shape_manual(values=shape_vec)
  p <- p + scale_size_manual(values=size_vec)
  p <- p + theme_cowplot(12)

  return(p)
  
}


DMR_analysis <- function(
  sample1_beta,
  sample1_name,
  sample2_beta,
  sample2_name,
  row_ranges,
  min_sample_per_probe,
  plot_dir
) {

  print(
    paste0(
      "Detecting differentially methylated probes between ",
      sample1_name, " and ", sample2_name, "..."
    )
  )

  # fetch beta dfs:
  unfilt_mtx <- list(
    sample1_beta = as.matrix(sample1_beta),
    sample2_beta = as.matrix(sample2_beta) )
  
  # remove rows with < min_sample_per_probe samples not NA:
  filt_mtx <- lapply(unfilt_mtx, function(x) {
    return(
      x[apply(x, 1, function(y) {
          length(which(!is.na(y))) >= min_sample_per_probe
        }), ] )
  })
  
  beta_mtx <- do.call("cbind", filt_mtx)
  
  # adjust other objects for new se:
  row_ranges <- row_ranges[names(row_ranges) %in% rownames(beta_mtx)]
  
  DM_se <- SummarizedExperiment(
    assays = list(counts = beta_mtx),
    rowRanges = row_ranges,
    colData = DataFrame(
      tissue = c(rep(sample1_name, ncol(filt_mtx$sample1_beta)),
        rep(sample2_name, ncol(filt_mtx$sample2_beta)) )))
  
  DMR <- TCGAanalyze_DMC(
    data = DM_se,
    groupCol = "tissue",
    group1 = sample1_name,
    group2 = sample2_name,
    plot.filename = paste0(
      plot_dir, "/", sample1_name, "_vs_", sample2_name,"_DMR_volcano.png" ),
    p.cut = 0.05 )

  return(DMR)

}