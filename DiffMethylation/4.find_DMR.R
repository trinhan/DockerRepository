
norm_group <- "pNF"
min_tissue_samples <- 5
min_sample_per_probe <- 5
DM_source <- "koelsche_MPNST"

# define and create directories:
home_dir <- "/share/ScratchGeneral/jamtor"
#home_dir <- "/Users/torpor/clusterHome/"
project_dir <- file.path(home_dir, "projects/MPNST_ctMet")

script_dir <- file.path(project_dir, "scripts")
func_dir <- file.path(script_dir, "functions")
col_dir <- file.path(home_dir, "R/colour_palettes")
ref_dir <- file.path(project_dir, "refs")

raw_dir <- file.path(project_dir, "raw_files")
in_dir <- file.path(raw_dir, "Rdata")
nonmalig_in <- file.path(raw_dir, "TCGA_normal")

result_dir <- file.path(project_dir, "results")
Robject_dir <- file.path(result_dir, "beta_values/Rdata")
dir.create(Robject_dir, recursive = T)
plot_dir <- file.path(result_dir, "beta_values/plots")
dir.create(plot_dir)

out_path <- file.path(result_dir, "DMR")
out_Robject_dir <- file.path(out_path, "Rdata")
dir.create(out_Robject_dir, recursive = T)
table_dir <- file.path(out_path, "tables")
dir.create(table_dir, recursive = T)

# define tissues:
tissue_key <- data.frame(
  code = c(
    "BLCA", "blood", "BRCA", 
    "CESC", "CHOL", "COAD", 
    "ESCA", "HNSC", "KIRP", 
    "LIHC", "LUAD", "LUSC", 
    "mixed", "MPNST", "pNF", 
    "PAAD", "PRAD", "READ", 
    "SARC", "STAD", "THCA", 
    "THYM" ),
  name = c(
    "bladder", "blood", "breast",
    "cervical", "bile duct", "colon",
    "esophagus", "head/neck", "kidney",
    "liver", "lung", "lung",
    "mixed", "MPNST", "pNF", 
    "pancreas", "prostate", "rectum", 
    "sarcoma", "stomach", "thyroid", 
    "thymus" ))

# define assay types:
assay_types <- data.frame(
  Source = c(
    "Koelsche et. al., Nat Comms 2021",  "Koelsche et. al., Nat Comms 2021 validation", 
    "TCGA_MPNST", "Lyskjaer et. al., J. Path. 2020", "GEO", "TCGA" ),
  Assay = c("EPIC/450k", "EPIC", "450k", "EPIC", rep("450k", 2)) )


####################################################################################
### 0. Load packages, functions and colours ###
####################################################################################

library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
library(cowplot)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tibble)

# source functions:
source(file.path(func_dir, "4.find_DMR_functions.R"))

cols <- c(TCGA_MPNST = "black", koelsche.reference_MPNST = "blue", 
  koelsche.validation_MPNST = "#1B9E77", lyskjaer_MPNST = "#E7298A", 
  koelsche.reference_pNF = "#A6761D",   GEO_blood = "#F8766D", 
  TCGA_nonmalig = "#58B9DB", TCGA_malig = "#E6AB02")

shape_vec <- c(16, 17, 13, 15, 23, 8, 25, 7, 22, 24)


####################################################################################
### 1. Load data ###
####################################################################################

print(paste0("Loading data..."))

MPNST_beta <- readRDS(file.path(in_dir, "MPNST_beta.rds"))
pNF_beta <- readRDS(file.path(in_dir, "pNF_beta.rds"))
blood_beta <- readRDS(file.path(in_dir, "blood_beta.rds"))
nonmalig_beta <- readRDS(file.path(in_dir, "TCGA_nonmalig_beta.rds"))
TCGA_malig_beta <- readRDS(file.path(in_dir, "TCGA_malig_beta.rds"))

# separate koelsche validation samples:
MPNST_beta$koelsche_val_MPNST <- MPNST_beta$koelsche_MPNST[
  , grep("koelsche_v", colnames(MPNST_beta$koelsche_MPNST))]
MPNST_beta$koelsche_MPNST <- MPNST_beta$koelsche_MPNST[
  , grep("koelsche_v", colnames(MPNST_beta$koelsche_MPNST), invert=T)]

# merge nonmalig data:
nonmalig_merged <- do.call("cbind", nonmalig_beta)
colnames(nonmalig_merged) <- gsub("^.*\\.", "", colnames(nonmalig_merged))
colnames(nonmalig_merged) <- sub("_[^_]+$", "_nonmalig", colnames(nonmalig_merged))


if (!file.exists(file.path(ref_dir, "annot_450k.rds"))) {
  
  library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  annotdf_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  annotdf_450k <- annotdf_450k[grep("cg", annotdf_450k$Name),]

  annot_450k <- GRanges(
    seqnames=annotdf_450k$chr,
    ranges = IRanges(start=annotdf_450k$pos, width=1),
    strand = annotdf_450k$strand )

  mcols(annot_450k) <- subset(annotdf_450k, select = c(Name, CpG_rs, Islands_Name,
    Relation_to_Island, UCSC_RefGene_Name, UCSC_RefGene_Accession, Enhancer, 
    Regulatory_Feature_Name, Regulatory_Feature_Group ))
  colnames(mcols(annot_450k))[1] <- c("Probe_ID")
  names(annot_450k) <- annot_450k$Probe_ID

  saveRDS(annot_450k, file.path(ref_dir, "annot_450k.rds"))
  
} else {
  annot_450k <- readRDS(file.path(ref_dir, "annot_450k.rds"))
}


####################################################################################
### 2. All sample PCA plots ###
####################################################################################

# combine data:
all_beta <- c(MPNST_beta, list(pNF = pNF_beta, blood = blood_beta, 
    TCGA_malig = TCGA_malig_beta, nonmalig = nonmalig_merged) )
common_probes <- Reduce(intersect, lapply(all_beta, rownames))
all_beta <- lapply(all_beta, function(x) x[common_probes,])
all_beta <- do.call("cbind", all_beta)
colnames(all_beta) <- gsub("^.*?\\.", "", colnames(all_beta))

# perform PCA on all samples:
if (!file.exists(file.path(Robject_dir, "methylation_pca.rds"))) {
  # create PCA with no source subgroups:
  pca_mtx <- t(na.omit(all_beta))
  # perfom pca and save for report:
  pca_res <- prcomp(pca_mtx, scale. = TRUE)
  # plot variances:
  png(file.path(plot_dir, "methylation_pca_weights.png"))
    plot(pca_res)
  dev.off()
  saveRDS(pca_res, file.path(Robject_dir, "methylation_pca.rds"))
} else {
  pca_res <- readRDS(file.path(Robject_dir, "methylation_pca.rds"))
}

# relabel Koelsche validation samples:
rownames(pca_res$x) <- gsub("koelscheval", "koelsche.validation", rownames(pca_res$x))
rownames(pca_res$x) <- gsub("koelsche_", "koelsche.reference_", rownames(pca_res$x))

pca_plot <- do_PCA(pca_in = as.data.frame(pca_res$x), col_vec = cols, 
  shape_vec = shape_vec )

png(file.path(plot_dir, "methylation_pca.png"), height = 7, width = 10, res = 300, units = "in")
  plot(pca_plot)
dev.off()
pdf(file.path(plot_dir, "methylation_pca.pdf"), height = 7, width = 10)
  plot(pca_plot)
dev.off()

saveRDS(pca_plot, file.path(Robject_dir, "methylation_pca_plot.rds"))


####################################################################################
### 3. Prepare and summarise DM analysis data ###
####################################################################################

print("Preparing DM analysis...")

if (!file.exists(
  paste0(out_Robject_dir, "/unfiltered_MPNST_vs_", norm_group, "_DMR.rds") )) {
  
  # assign normal group and fetch probe coordinates:
  if (norm_group == "nonmalig") {
    normal_df <- do.call("cbind", nonmalig_beta)
  } else {
    normal_df <- pNF_beta
  }
  norm_probe_coords <- annot_450k[annot_450k$Probe_ID %in% rownames(normal_df)]

  # subset MPNST data if needed:
  MPNST_df <- MPNST_beta[DM_source][[1]]

  # create sample summary:
  MPNST_no <- sapply(unique(gsub("_.*$", "", colnames(MPNST_df))), function(x) {
    return(ncol(MPNST_df[,grep(x, colnames(MPNST_df))]))
  })
  names(MPNST_no) <- gsub("koelsche", "Koelsche et. al., Nat Comms 2021", 
    names(MPNST_no) )
  MPNST_no <- MPNST_no[order(names(MPNST_no))]

  # record MPNST numbers not to be used in DMR:
  alt_MPNST <- MPNST_beta[names(MPNST_beta) != DM_source]
  alt_MPNST_no <- sapply(alt_MPNST, ncol)
  alt_MPNST_no <- alt_MPNST_no[order(names(alt_MPNST_no))]

  pNF_no <- ncol(pNF_beta)
  nonmalig_no <- sapply(nonmalig_beta, function(x) {
    return(ncol(x))
  })

  # merge lung sample numbers:
  nonmalig_no <- c(
    nonmalig_no, sum(nonmalig_no[names(nonmalig_no) == "lung"]) )
  nonmalig_no <- nonmalig_no[names(nonmalig_no) != "lung"]
  names(nonmalig_no)[length(nonmalig_no)] <- "lung"

  # order non malig tissues:
  nonmalig_order <- c( "liver", "kidney", 
    sort(names(nonmalig_no)[!(names(nonmalig_no) %in% c("liver", "kidney"))]) )

  sample_summary <- data.frame(
    Tissue = c(rep("MPNST", length(c(MPNST_no, alt_MPNST_no))), "pNF", "blood", nonmalig_order),
    Source = c(names(MPNST_no), paste0(names(MPNST_no), " validation"), 
      "Lyskjaer et. al., J. Path. 2020", "TCGA_MPNST", "Koelsche et. al., Nat Comms 2021", 
      "GEO", rep("TCGA", length(nonmalig_order)) ),
    Sample_number = c(MPNST_no, alt_MPNST_no, pNF_no, ncol(blood_beta), nonmalig_no[nonmalig_order]) )
  sample_summary <- merge(sample_summary, assay_types, by="Source")
  sample_summary$Assay[sample_summary$Tissue == "pNF"] <- "450k"

  # adjust order to put all MPNST at top:
  sample_summary$Tissue <- factor(
    sample_summary$Tissue, levels = c("MPNST", "pNF", "blood", nonmalig_order) )
  sample_summary <- sample_summary[order(sample_summary$Tissue),]

  saveRDS(sample_summary, file.path(Robject_dir, "sample_summary.rds"))


  ####################################################################################
  ### 4. Identify MPNST vs clt DMRs ###
  ####################################################################################

  print("Finding DMR...")

  # keep common probes:
  MPNST_df <- MPNST_df[rownames(MPNST_df) %in% rownames(normal_df),]
  normal_df <- normal_df[rownames(normal_df) %in% rownames(MPNST_df),]
  
  # find DMRs between normal and MPNST groups: 
  DMR <- DMR_analysis(
    sample1_beta = MPNST_df, 
    sample1_name = "MPNST", 
    sample2_beta = normal_df, 
    sample2_name = norm_group,
    row_ranges = norm_probe_coords,
    min_sample_per_probe,
    plot_dir )

  saveRDS(DMR, 
    paste0(out_Robject_dir, "/unfiltered_MPNST_vs_", norm_group, "_DMR.rds") )

} else {
  DMR <- readRDS(
    paste0(out_Robject_dir, "/unfiltered_MPNST_vs_", norm_group, "_DMR.rds") )
}

print("Formatting and saving...")
  
# create summary dataframe:
DMR_df <- subset(DMR, select=-status)
colnames(DMR_df) <- c("MPNST_mean_beta", paste0(norm_group, "_mean_beta"),
  paste0("mean_diff_MPNST_vs_ctl"), "pval", "pval_adj" )

# round to 3 decimals:
DMR_df[,1:3] <- round(DMR_df[,1:3], 3)

# add dm_status:
DMR_df$status <- "hypermethylated"
DMR_df$status[DMR_df[,3] < 0] <- "hypomethylated"

# add assoc genes and coords:
dm_coords <- subset(
  as.data.frame(annot_450k[names(annot_450k) %in% rownames(DMR_df)]),
  select = c(seqnames, start, UCSC_RefGene_Name) )
colnames(dm_coords) <- c("chr", "genomic_coord", "gene")
dm_coords$gene[dm_coords$gene == ""] <- "none_annotated"
annot_DMR <- merge(DMR_df, dm_coords, by=0)
colnames(annot_DMR)[colnames(annot_DMR) == "Row.names"] <- "probe_id"
annot_DMR <- subset(annot_DMR, select = c(probe_id, chr, genomic_coord, gene, 
  status,mean_diff_MPNST_vs_ctl, pval, pval_adj ))

# sort:
annot_DMR <- arrange(annot_DMR, status, pval_adj, 
  desc(mean_diff_MPNST_vs_ctl) )

# add control group name and write:
colnames(annot_DMR) <- gsub("ctl", norm_group, colnames(annot_DMR))
write.table(annot_DMR, 
  paste0(table_dir, "/unfiltered_sig_MPNST_vs_", norm_group, "_DMR.tsv"),
  sep = "\t", quote = F, row.names = F, col.names = T )



