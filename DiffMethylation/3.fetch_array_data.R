library(GEOquery)
library(dplyr)
library(GenomicRanges)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(minfi)
library(tibble)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

home_dir <- "/share/ScratchGeneral/jamtor"
project_dir <- file.path(home_dir, "projects/MPNST_ctMet")
raw_dir <- file.path(project_dir, "raw_files")
ref_dir <- file.path(project_dir, "refs")
genome_dir <- file.path(project_dir, "genome")
Robject_dir <- file.path(raw_dir, "Rdata")

dir.create(Robject_dir)


####################################################################################
### 0. Define functions
####################################################################################

matched_meth_exp <- function(project, n = NULL){
  # get primary solid tumor samples: DNA methylation
  message("Download DNA methylation information")
  met450k <- GDCquery(project = project,
    data.category = "DNA methylation",
    platform = "Illumina Human Methylation 450",
    legacy = TRUE, 
    sample.type = c("Primary Tumor"))
  met450k.tp <-  met450k$results[[1]]$cases
  
  # Get patients with samples in both platforms
  patients <- unique(substr(met450k.tp,1,12))
  if(!is.null(n)) {
    patients <- patients[1:n] } # get only n samples
  names(patients) <- rep(project, length(patients))
  return(patients)
}

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
    "THYM"
  ),
  name = c(
    "bladder", "blood", "breast",
    "cervical", "bile duct", "colon",
    "esophagus", "head/neck", "kidney",
    "liver", "lung", "lung",
    "mixed", "MPNST", "pNF", 
    "pancreas", "prostate", "rectum", 
    "sarcoma", "stomach", "thyroid", 
    "thymus"
  )
)


####################################################################################
### 1. Download TCGA MPNST data
####################################################################################

# fetch samples:
TCGA_meta <- read.table(file.path(ref_dir, "formatted_TCGA_MPNST_metadata.tsv"), 
  sep = "\t", header=T )
TCGA_id <- TCGA_meta$ID

if (!file.exists(file.path(raw_dir, "TCGA_MPNST/TCGA_SARC_450K_hg19.rda"))) {

  # download all sarcoma data from TCGA:
  TCGAbiolinks:::getGDCprojects()$project_id
  
  cohort_name <- "TCGA-SARC"
  
  query.met <- GDCquery(
    project = cohort_name,
    data.category = "DNA methylation", 
    platform = "Illumina Human Methylation 450", 
    sample.type = "Primary Tumor",
    legacy = TRUE
  )
  
  GDCdownload(query.met)
    
  data <- GDCprepare(
    query = query.met,
    save = TRUE,
    save.filename = file.path(raw_dir, "TCGA_MPNST/TCGA_SARC_450K_hg19.rda"),
    summarizedExperiment = TRUE
  )

} else {
  load(file.path(raw_dir, "TCGA_MPNST/TCGA_SARC_450K_hg19.rda"))
}

# isolate beta values:
sarc_beta <- data@assays@data@listData[[1]]

# isolate MPNST samples:
for (i in seq_along(TCGA_id)) {

  print(paste0("Fetching ", TCGA_id[i], "..."))

  if (i==1) {
    MPNST_beta <- list(TCGA_MPNST = data.frame(sarc_beta[, grep(TCGA_id[i], colnames(sarc_beta))]))
    colnames(MPNST_beta$TCGA_MPNST)[i] <- TCGA_id[i]
  } else {
    MPNST_beta$TCGA_MPNST <- cbind(MPNST_beta,data.frame(sarc_beta[
      , grep(TCGA_id[i], colnames(sarc_beta)) ]))
    colnames(MPNST_beta$TCGA_MPNST)[i] <- TCGA_id[i]
  }
}

# delete of replace dashes with underscores for colnames and add 'MPNST':
colnames(MPNST_beta$TCGA_MPNST) <- gsub("^.*\\.", "", colnames(MPNST_beta$TCGA_MPNST))
colnames(MPNST_beta$TCGA_MPNST) <- gsub("-", "", colnames(MPNST_beta$TCGA_MPNST))
colnames(MPNST_beta$TCGA_MPNST) <- paste0(
  gsub("TCGA", "TCGA_", colnames(MPNST_beta$TCGA_MPNST)), "_MPNST")


####################################################################################
### 2. Download  Koelsche et. al., Nat. Comms. 2021 data
####################################################################################

if (!file.exists(file.path(Robject_dir, "koelsche_beta.rds"))) {

  koelsche <- list(
    beta_450k = read.table(
      file.path(raw_dir, "koelsche_2021/GSE140686_GPL13534_matrix_processed.txt"),
      sep = "\t",
      header = T ),
    beta_epic = read.table(
      file.path(raw_dir, "koelsche_2021/GSE140686_GPL21145_matrix_processed.txt"),
      sep = "\t",
      header = T ) )

  array_type <- lapply(koelsche, function(x) colnames(x)[grep("Detection", colnames(x), invert=T)])
  saveRDS(array_type, file.path(Robject_dir, "sample_array_types.rds"))

  # fetch ids of samples to keep:
  koelsche_meta <- read.table(file.path(ref_dir, "formatted_koelsche_MPNST_metadata.tsv"),
    sep="\t", header=T )
  koelsche_meta <- split(koelsche_meta, koelsche_meta$Revised.diagnosis)
  names(koelsche_meta) <- c("MPNST", "pNF")
  pNF_id <- gsub("_pNF", "", gsub(" ", ".", koelsche_meta$pNF$ID))
  MPNST_id <- gsub("_MPNST", "", gsub(" ", ".", koelsche_meta$MPNST$ID))
  
  for (i in 1:length(koelsche)) {

    print(i)
  
    # add and format colnames to match ids:
    kbeta <- koelsche[[i]] %>%
      column_to_rownames("ID_REF")
    colnames(kbeta) <- gsub("REFERENCE_SAMPLE\\.", "koelsche_", colnames(kbeta))
    colnames(kbeta) <- gsub("VALIDATION_SAMPLE\\.", "koelsche_v", colnames(kbeta))
    
    # fetch indices of pNF samples:
    pNF_keep <- which(colnames(kbeta) %in% pNF_id)
    # fetch indices of p-values:
    pNF_pkeep <- pNF_keep+1
    # fetch beta values:
    pNF_temp <- kbeta[,pNF_keep]
  
    if (ncol(pNF_temp) > 0) {
      # add source to colnames:
      colnames(pNF_temp) <- paste(gsub("^.*_SAMPLE\\.", "", colnames(pNF_temp)), 
        "pNF", sep = "_" )
      # check p values are all < 0.05:
      range(kbeta[,pNF_pkeep])
    }
    
    # fetch indices of MPNST samples:
    MPNST_keep <- which(colnames(kbeta) %in% MPNST_id)
    # fetch indices of p-values:
    MPNST_pkeep <- MPNST_keep+1
    # fetch beta values:
    MPNST_temp <- kbeta[,MPNST_keep]
    # add source to colnames:
    colnames(MPNST_temp) <- paste0(colnames(MPNST_temp), "_MPNST")

    # check p values are all < 0.05:
    range(kbeta[,MPNST_pkeep])
    
    if (i==1) {
      koelsche_beta <- list(MPNST = MPNST_temp, pNF = pNF_temp)
    } else {
  
      # keep probes already in MPNST_beta only:
      MPNST_filt <- MPNST_temp[
        rownames(MPNST_temp) %in% rownames(koelsche_beta$MPNST),]
      koelsche_beta$MPNST <- koelsche_beta$MPNST[
        rownames(koelsche_beta$MPNST) %in% rownames(MPNST_filt),]
      koelsche_beta$MPNST <- cbind(koelsche_beta$MPNST, MPNST_filt )
  
      if (ncol(pNF_temp) > 0) {
        # keep probes already in pNF_beta only:
        pNF_filt <- pNF_temp[
          rownames(pNF_temp) %in% rownames(koelsche_beta$pNF),]
        koelsche_beta$pNF <- koelsche_beta$pNF[
          rownames(koelsche_beta$pNF) %in% rownames(pNF_filt),]
        koelsche_beta$pNF <- cbind(
          koelsche_beta$pNF, pNF_filt )
      }
  
    }
  
  }
  saveRDS(koelsche_beta, file.path(Robject_dir, "koelsche_beta.rds"))

} else {
  koelsche_beta <- readRDS(file.path(Robject_dir, "koelsche_beta.rds"))
}

# add koelsche to MPNST:
MPNST_beta$koelsche_MPNST <- koelsche_beta$MPNST

# save RDS:
saveRDS(koelsche_beta$pNF, file.path(Robject_dir, "pNF_beta.rds"))


####################################################################################
### 3. Load and bind Lyskjaer J. Path 2020 data
####################################################################################

MPNST_beta$lyskjaer_MPNST <- readRDS(
  file.path(raw_dir, "lyskjaer_2020/beta/beta_values.rds"))

# keep common probes only:
common_probes <- Reduce(intersect, lapply(MPNST_beta, rownames))
MPNST_beta <- lapply(MPNST_beta, function(x) x[common_probes,])

# order MPNST_beta alphabetically:
MPNST_beta <- MPNST_beta[order(names(MPNST_beta))]

# save RDS:
saveRDS(MPNST_beta, file.path(Robject_dir, "MPNST_beta.rds"))


####################################################################################
### 4. Load and format TCGA non-malignant data ###
####################################################################################

print("Loading non-malignant beta values...")
  
if (!file.exists(file.path(Robject_dir, "nonmalig_se.rds"))) {
  
  # load nonmalig se objects:
  nonmalig_se_files <- list.files(file.path(raw_dir, "TCGA_normal/"), pattern = "rda")
  
  for (i in 1:length(nonmalig_se_files)) {
    
    print(paste0("Loading ", nonmalig_se_files[i], "..."))
    
    load(file.path(raw_dir, "TCGA_normal/", nonmalig_se_files[i]))
    
    if (i==1) {
      nonmalig_se <- list(data)
    } else {
      nonmalig_se[[i]] <- data
    }
    
    names(nonmalig_se)[i] <- gsub(".rda", "", nonmalig_se_files[i])
    
  }
  
  saveRDS(nonmalig_se, file.path(Robject_dir, "nonmalig_se.rds"))
  
} else {
  nonmalig_se <- readRDS(file.path(Robject_dir, "nonmalig_se.rds"))
}

# isolate beta matrices from se objects:
nonmalig_beta <- lapply(nonmalig_se, function(x) assays(x)[[1]])

# make all df:
nonmalig_beta <- lapply(nonmalig_beta, as.data.frame)

# rename elements:
names(nonmalig_beta) <- paste0(
  gsub("_.*$|^.*-", "", names(nonmalig_beta)), "_nonmalig" )

# change tissue codes to names:
names(nonmalig_beta) <- tissue_key$name[
    match(gsub("_nonmalig", "", names(nonmalig_beta)), tissue_key$code) ]

# change dashes to underscores in colnames:
for (i in seq_along(nonmalig_beta)) {
  colnames(nonmalig_beta[[i]]) <- gsub("-", "", colnames(nonmalig_beta[[i]]))
  colnames(nonmalig_beta[[i]]) <- gsub("TCGA", "TCGA_", colnames(nonmalig_beta[[i]]))
  colnames(nonmalig_beta[[i]]) <- paste0(colnames(nonmalig_beta[[i]]), "_", 
    names(nonmalig_beta)[i])
}

# save RDS:
saveRDS(nonmalig_beta, file.path(Robject_dir, "TCGA_nonmalig_beta.rds"))


####################################################################################
### 6. Fetch blood methylation data
####################################################################################

# load GEO blood data and fetch beta matrix:
blood_se <- readRDS(file.path(raw_dir, "GEO_blood/SE_GEO_Meth_Data.rds"))
blood_df <- assays(blood_se)[[1]][-1,]

# add source and tissue type to colnames:
colnames(blood_df) <- paste0("GEO_", colnames(blood_df), "_blood")

saveRDS(blood_df, file.path(Robject_dir, "blood_beta.rds"))
 
# save metadata:
write.table(as.data.frame(blood_se@colData), 
  file.path(ref_dir, "blood_metadata.tsv"), sep = "\t", row.names = F, col.names = T)


####################################################################################
### 5. Download malignant TCGA data for PCA comparison
####################################################################################

malig_projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-PRAD")

if (!file.exists(file.path(raw_dir, "TCGA_MPNST/TCGA_malig_450K_hg19.rda"))) {
  
  TCGA_malig <- lapply(malig_projects, matched_meth_exp, 30)
  
  samples <- do.call("c", TCGA_malig)
  
  query.met <- GDCquery(project = malig_projects,
    data.category = "DNA methylation",
    platform = "Illumina Human Methylation 450",
    legacy = TRUE, 
    barcode = samples )
  
  GDCdownload(query.met)
  
  data <- GDCprepare(query = query.met,
      save = TRUE,
      save.filename = file.path(raw_dir, "TCGA_MPNST/TCGA_malig_450K_hg19.rda"),
      summarizedExperiment = TRUE )

  saveRDS(samples, file.path(Robject_dir, "TCGA_malig_samples.rds"))

} else {
  load(file.path(raw_dir, "TCGA_MPNST/TCGA_malig_450K_hg19.rda"))
  samples <- readRDS(file.path(Robject_dir, "TCGA_malig_samples.rds"))
}

# isolate beta values:
TCGA_malig_beta <- data@assays@data@listData[[1]]

# add project to colnames:
colnames(TCGA_malig_beta) <- gsub("-", "", colnames(TCGA_malig_beta))
colnames(TCGA_malig_beta) <- paste0(gsub("-", "_", names(samples)), 
  gsub("TCGA", "", colnames(TCGA_malig_beta)), "_malig" )

# save RDS:
saveRDS(TCGA_malig_beta, file.path(Robject_dir, "TCGA_malig_beta.rds"))



