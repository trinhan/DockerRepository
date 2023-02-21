record_probes <- function(dmr, current_record, label) {
  dmr <- split(dmr, dmr$status)
  new_record <- data.frame(filt = c(paste0(label, "_hyper"), paste0(label, "_hypo")),
    no_probes = c(nrow(dmr$hypermethylated), nrow(dmr$hypomethylated)) )
  if (current_record == "none") {
    print(new_record)
    return(new_record)
  } else {
    new_record <- rbind(current_record, new_record)
    print(new_record)
  }
}


calc_stats <- function(
  beta, 
  dmr,  
  data_source) {

  source("/home/jamtor/base.R")

  rownames(dmr) <- dmr$probe_id

  # keep dm probes:
  dm_beta <- beta[rownames(beta) %in% dmr$probe_id, ]
  
  # calculate stats for DM probes:
  stats_df <- as.data.frame(t(apply(dm_beta, 1, function(y) {
    quantiles <- quantile(y, c(0.1, 0.25, 0.75, 0.9), na.rm = TRUE)
    return(c(
      median = median(y, na.rm = TRUE),
      mean = mean(y, na.rm = TRUE),
      sd = sd(y, na.rm = TRUE),
      quantile = quantiles[1],
      quantile = quantiles[2],
      quantile = quantiles[3],
      quantile = quantiles[4] ))
  })))
  # round to 3 decimals:
  stats_df <- round2(stats_df, 3)
  # remove probes with NA:
  stats_df <- stats_df[!is.na(stats_df$median),]
  # add data_source type to stats_df colnames:
  colnames(stats_df) <- paste0(data_source, "_", colnames(stats_df))

  return(stats_df)

}


plot_DMR_record <- function(record_df, labels, cols, log_convert = F) {

  plot_df <- data.frame(
    stage = rep(labels, each = 2),
    type = paste0(gsub("^.*_", "", record_df$filt), "methylated"),
    number = record_df$no_probes )

  plot_df$stage <- factor(plot_df$stage, levels = unique(plot_df$stage))
  no_sets <- nrow(plot_df)/2

  # plot:
  p <- ggplot(plot_df, aes(x=stage, y=number, fill=type))
  p <- p + geom_bar(stat="identity", position = "dodge")
  p <- p + geom_text(
    aes(label=number), 
    position=position_dodge(width=0.9), 
    vjust=-0.25,
    fontface = "bold")
  p <- p + scale_fill_manual(values=rep(cols, no_sets))
  p <- p + ylab("No. of probes")
  p <- p + theme_cowplot(12)
  p <- p + theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.title.x=element_blank()
  )
  if (log_convert) {
    # convert to log:
    p <- p + scale_y_continuous(trans='log10')
  }

  return(p)
}


do_boxplot <- function(beta, fsize, xlab_angle = 45, combine_nonmalig = TRUE, 
  probe_title_info = FALSE, dmr_table = NULL, rowname_to_title = FALSE ) {

  library(reshape2)
  library(ggplot2)

  beta$probe_id <- rownames(beta)
  beta_df <- melt(beta, 
    id = "probe_id", 
    variable.name = "tissue", 
    value.name = "beta" )

  # format source/tissue:
  beta_df$tissue <- gsub("_.*_", "_", beta_df$tissue)

  # fetch number of samples of each tissue:
  nvals <- sapply(split(beta_df$tissue, beta_df$tissue), length)/
    length(unique(beta_df$probe_id))
  # order nvals:
  nvals <- nvals[match(unique(beta_df$tissue), names(nvals))]

  if (!combine_nonmalig) {

    # remove tissue types with < 5 samples:
    beta_df <- beta_df[!(beta_df$tissue %in% names(nvals)[nvals < 5]),]
    # order columns:
    beta_df$probe_id <- factor(beta_df$probe_id, levels = beta$probe_id)
    beta_df$tissue <- factor(beta_df$tissue, levels = unique(beta_df$tissue))
    levels(beta_df$tissue) = gsub("_", " ", 
      paste0(names(nvals), " (", nvals, ")") )

  } else {
    # order columns:
    beta_df$probe_id <- factor(beta_df$probe_id, levels = beta$probe_id)
    beta_df$tissue <- factor(beta_df$tissue, levels = unique(beta_df$tissue)) 
    levels(beta_df$tissue) = gsub("_", " ", 
      paste0(names(nvals), " (", nvals, ")") )
  }
  
  # plot median boxplot:
  p <- ggplot(beta_df, aes(x=probe_id, y=beta, fill=tissue))
  p <- p + geom_boxplot()
  p <- p + theme_classic()
  p <- p + scale_fill_manual(
    values = plot_cols$V1[seq_along(levels(beta_df$tissue))] )
  p <- p + ylab("Median beta value (1st/3rd quartiles, 95% CI)")
  if (probe_title_info) {
    dmr_row <- dmr_table[dmr_table$probe_id == unique(beta_df$probe_id),]
    ptitle <- paste(dmr_row$probe_id, 
      paste0(dmr_row$chr, ":", dmr_row$genomic_coord), dmr_row$gene, sep = ", ")
    p <- p + ggtitle(gsub(", none_annotated", "", ptitle))
    p <- p + theme(
      plot.title = element_text(hjust = 0.5, size = fsize), 
      text = element_text(size = fsize),
      axis.text.x = element_blank(), axis.title.x = element_blank(),
      axis.ticks.x=element_blank(), axis.title.y = element_text(size = 22) )
  } else if (rowname_to_title) {
    p <- p + ggtitle(as.character(beta_df$probe_id)[1])
    p <- p + theme(
      plot.title = element_text(hjust = 0.5, size = fsize),
      text = element_text(size=fsize),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(), axis.ticks.x=element_blank(), 
      axis.title.y = element_text(size = 22) ) 
  } else {
    p <- p + theme(
      text = element_text(size=fsize),
      axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1),
      axis.title.x = element_blank(), axis.ticks.x=element_blank() ) 
    dmr_row <- dmr_table[match(unique(beta_df$probe_id), dmr_table$probe_id),]
    dmr_row$gene <- paste0("(", gsub(";.*$", "", dmr_row$gene), ")")
    dmr_row$gene <- gsub("\\(none_annotated\\)", "", dmr_row$gene)
    p <- p + scale_x_discrete(labels=paste0(unique(beta_df$probe_id), " ", dmr_row$gene)) 
  }
  return(p)
  
}

test_meta_association <- function(meta_df, hyper_mean, bp_cols, fsize=8) {

  # add mean hypermethylated beta values for each sample:
  print("Dimensions of metadata before merging with means: ")
  print(dim(meta_df))
  meta_df <- merge(meta_df, hyper_mean, by="ID")
  print("Dimensions of metadata after merging with means: ")
  print(dim(meta_df))
  
  for (i in 2:(ncol(meta_df)-1)) {
    
    category <- colnames(meta_df)[i]
    print(category)
    
    # prepare df:
    anova_df <- meta_df[,c(1, i)]
    
    # remove any categories with only 1 sample:
    remove <- sapply(split(anova_df, anova_df[,2]), function(x) nrow(x) == 1)
    remove <- names(remove)[which(remove)]
    if (length(remove) > 0) {
      anova_df[,2][anova_df[,2] == remove & !is.na(anova_df[,2])] <- NA
    }
    
    if (length(unique(anova_df[,2][!is.na(anova_df[,2])])) > 2) {

      # bind with mean beta values:
      anova_df <- cbind(anova_df, as.data.frame(meta_df$Mean.beta))
      colnames(anova_df)[3] <- "Mean.beta"
      anova_df[,2] <- factor(anova_df[,2])
      levels(anova_df[,2])
      
      # Compute the analysis of variance
      res.aov <- aov(anova_df$Mean.beta ~ anova_df[,2])
      # Summary of the analysis
      res.summary <- summary(res.aov)
      print(res.summary)
      
      if (res.summary[[1]]$`Pr(>F)`[1] < 0.05) {

        # record variable if significant:
        if (!(exists("sig_vars"))) {
          sig_vars <- category
        } else {
          sig_vars[i] <- category
        }

        # Tukey pairwise comparisons:
        res.tukey <- TukeyHSD(res.aov)
        tukey_df <- as.data.frame(res.tukey$`anova_df[, 2]`)
        
        # compare pairwise using t-test:
        comp <- as.data.frame(t(combn(levels(anova_df[,2]), 2)))
        comp <- lapply(split(comp, rownames(comp)), function(x) t(x)[,1])
        
        bp <- ggboxplot(anova_df, x=category, y="Mean.beta", 
                        color=category, palette=bp_cols,
                        ylab="Mean.beta", xlab=category )
        bp <- bp + stat_compare_means(comparisons=comp, method="t.test",
          size=fsize-21 )
        bp <- bp + rotate_x_text(45)
        bp <- bp + theme(text = element_text(size=fsize))
        bp <- bp + scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1.4))
        
      } else {
        # visualise on boxplot:
        bp <- ggboxplot(anova_df, x=category, y="Mean.beta", 
          color=category, palette=bp_cols,
          ylab="Mean.beta", xlab=category )
        bp <- bp + annotate("text", 
          label=paste0("ANOVA p.val = ", round(res.summary[[1]]$`Pr(>F)`[1], 2)), 
          size=fsize-21, x=1.5, y=0.9 )
        bp <- bp + rotate_x_text(45)
        bp <- bp + theme(text = element_text(size=fsize))
        bp <- bp + scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1.01))
        
      }
      
    } else if (length(unique(anova_df[,2][!is.na(anova_df[,2])])) == 2) {
      
      # bind with mean beta values:
      t_df <- cbind(anova_df, as.data.frame(meta_df$Mean.beta))
      colnames(t_df)[3] <- "Mean.beta"
      t_df[,2] <- factor(t_df[,2])
      levels(t_df[,2])
      
      # conduct t-test:
      t_res <- t.test(split(t_df, t_df[,2])[[1]]$Mean.beta, 
        split(t_df, t_df[,2])[[2]]$Mean.beta )

      # record variable if significant:
      if (t_res$p.value < 0.05) {
        if (!(exists("sig_vars"))) {
          sig_vars <- category
        } else {
          sig_vars[i] <- category
        }
      }

      # conduct t-test and plot:
      comp <- as.data.frame(t(combn(levels(t_df[,2]), 2)))
      comp <- lapply(split(comp, rownames(comp)), function(x) t(x)[,1])
      
      bp <- ggboxplot(t_df, x=category, y="Mean.beta", 
                      color=category, palette=bp_cols,
                      ylab="Mean.beta", xlab=category )
      bp <- bp + stat_compare_means(comparisons=comp, method="t.test",
        size=fsize-21)
      bp <- bp + rotate_x_text(45)
      bp <- bp + theme(text = element_text(size=fsize))
      bp <- bp + scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1.01))
      
    } else {
      bp <- NA
    }
    
    if (i==2) {
      bps <- list(bp)
      names(bps) <- category
    } else {
      bps[[i-1]] <- bp
      names(bps)[i-1] <- category
    }

  }
  bps <- bps[!is.na(bps)]

  if (exists("sig_vars")) {
    return(list(sig_vars=sig_vars[!is.na(sig_vars)], bps=bps))
  } else {
    return(list(sig_vars="none", bps=bps))
  }

}
