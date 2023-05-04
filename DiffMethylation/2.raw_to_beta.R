library(minfi)
library(naturalsort)

home_dir <- "/share/ScratchGeneral/jamtor"
project_dir <- file.path(home_dir, "projects/MPNST_ctMet")
ref_dir <- file.path(project_dir, "refs")
raw_dir <- file.path(project_dir, "raw_files")
in_dir <- file.path(raw_dir, "lyskjaer_2020/raw")
table_dir <- file.path(raw_dir, "lyskjaer_2020/tables")
plot_dir <- file.path(raw_dir, "lyskjaer_2020/plots")
out_dir <- file.path(raw_dir, "lyskjaer_2020/beta")

dir.create(table_dir)
dir.create(plot_dir)
dir.create(out_dir)


####################################################################################
### Load metadata and process raw values to beta ###
####################################################################################

lysk_meta <- read.table(file.path(ref_dir, "formatted_lyskjaer_MPNST_metadata.tsv"), 
  sep = "\t", header = T )

setwd(in_dir)

# split those idat files with different numbers of probes and preprocess in parallel:
meta_sub1 <- lysk_meta[lysk_meta$ID %in% c("lyskjaer_2_MPNST", 
  "lyskjaer_3_MPNST", "lyskjaer_4_MPNST", "lyskjaer_5_MPNST", "lyskjaer_6_MPNST",
  "lyskjaer_8_MPNST", "lyskjaer_9_MPNST", "lyskjaer_11_MPNST", "lyskjaer_13_MPNST",
  "lyskjaer_17_MPNST","lyskjaer_20_MPNST", "lyskjaer_21_MPNST","lyskjaer_23_MPNST",  
  "lyskjaer_25_MPNST", "lyskjaer_26_MPNST", "lyskjaer_30_MPNST", "lyskjaer_33_MPNST", 
  "lyskjaer_40_MPNST","lyskjaer_42_MPNST", "lyskjaer_43_MPNST", "lyskjaer_45_MPNST", 
  "lyskjaer_46_MPNST", "lyskjaer_49_MPNST", "lyskjaer_50_MPNST", "lyskjaer_52_MPNST", 
  "lyskjaer_53_MPNST", "lyskjaer_54_MPNST", "lyskjaer_56_MPNST", "lyskjaer_57_MPNST", 
  "lyskjaer_81_MPNST", "lyskjaer_84_MPNST", "lyskjaer_97_MPNST", "lyskjaer_176_MPNST" ), ]

meta_sub2 <- lysk_meta[
  lysk_meta$ID %in% c("lyskjaer_1_MPNST", "lyskjaer_14_MPNST", "lyskjaer_15_MPNST", 
    "lyskjaer_18_MPNST", "lyskjaer_19_MPNST", "lyskjaer_27_MPNST", "lyskjaer_29_MPNST", 
    "lyskjaer_32_MPNST", "lyskjaer_34_MPNST", "lyskjaer_35_MPNST", "lyskjaer_37_MPNST", 
    "lyskjaer_39_MPNST", "lyskjaer_58_MPNST", "lyskjaer_59_MPNST", "lyskjaer_61_MPNST", 
    "lyskjaer_85_MPNST", "lyskjaer_96_MPNST", "lyskjaer_101_MPNST", "lyskjaer_103_MPNST",
    "lyskjaer_105_MPNST", "lyskjaer_107_MPNST", "lyskjaer_131_MPNST", "lyskjaer_132_MPNST", 
    "lyskjaer_136_MPNST", "lyskjaer_138_MPNST", "lyskjaer_139_MPNST", "lyskjaer_141_MPNST", 
    "lyskjaer_143_MPNST", "lyskjaer_146_MPNST", "lyskjaer_153_MPNST" ), ]

RGsets <- list(
  set1 = read.metharray.exp(targets = meta_sub1),
  set2 = read.metharray.exp(targets = meta_sub2) )

# check probe numbers:
sapply(RGsets, function(x) {
  print(dim(x@assays@data@listData$Green))
  print(dim(x@assays@data@listData$Red))
})

# combine sets:
RGset <- combineArrays(RGsets$set1, RGsets$set2, 
  outType = "IlluminaHumanMethylationEPIC" )
# check probe number:
print(dim(RGset@assays@data@listData$Green))

# detect probes which have failed for any sample:
detP <- detectionP(RGset)
to_remove <- apply(detP, 1, function (y) any(y > 0.01))

# Preprocess and normalise:
mpnst <- preprocessNoob(RGset)

# check probe number:
print(dim(mpnst@assays@data@listData$Meth))

# remove failed probes and convert to beta values:
meth <- mpnst[!rownames(mpnst) %in% names(which(to_remove)),]
beta_vals <- getBeta(meth)

# check probe number:
dim(beta_vals)

# replace colnames with IDs and histo groups:
colnames(beta_vals) <- lysk_meta$ID[match(colnames(beta_vals), lysk_meta$Basename)]

# save beta values:
saveRDS(beta_vals, file.path(out_dir, "beta_values.rds"))


