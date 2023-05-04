
library(naturalsort)

#home_dir <- "/Users/torpor/clusterHome"
home_dir <- "/share/ScratchGeneral/jamtor"
project_dir <- file.path(home_dir, "projects/MPNST_ctMet")
ref_dir <- file.path(project_dir, "refs")


####################################################################################
### 1. Load and format Koelsche metadata ###
####################################################################################

koelsche_ref_meta <- read.table(
  file.path(ref_dir, "koelsche_reference_metadata.tsv"), sep = "\t", 
  header = T )
koelsche_ref_meta$Set <- "reference"

koelsche_val_meta <- read.table(
  file.path(ref_dir, "koelsche_validation_metadata.tsv"), sep = "\t", 
  header = T )
koelsche_val_meta$ID <- koelsche_val_meta$X
koelsche_val_meta$Set <- "validation"

koelsche_samples <- read.table(
  file.path(ref_dir, "koelsche_relevant_samples.tsv"), sep = "\t", header = T )

koelsche_meta <- merge(koelsche_ref_meta, koelsche_val_meta, all=T)
koelsche_meta <- koelsche_meta[koelsche_meta$ID %in% koelsche_samples$ID,]
koelsche_meta$Diagnosis[is.na(koelsche_meta$Diagnosis)] <- 
  "Malignant peripheral nerve sheath tumour"

koelsche_meta$ID <- gsub("REFERENCE_SAMPLE ", "koelsche_", koelsche_meta$ID)
koelsche_meta$ID <- gsub("VALIDATION_SAMPLE ", "koelsche_v", koelsche_meta$ID)
koelsche_meta$ID[koelsche_meta$Diagnosis != "Neurofibroma (plexiform)"] <- 
  paste0(koelsche_meta$ID[koelsche_meta$Diagnosis != "Neurofibroma (plexiform)"], "_MPNST")
koelsche_meta$ID[koelsche_meta$Diagnosis == "Neurofibroma (plexiform)"] <- 
  paste0(koelsche_meta$ID[koelsche_meta$Diagnosis == "Neurofibroma (plexiform)"], "_pNF")

koelsche_meta$Diagnosis[!is.na(koelsche_meta$V12.2_MaxCalDiag_MCF)] <- 
  koelsche_meta$V12.2_MaxCalDiag_MCF[!is.na(koelsche_meta$V12.2_MaxCalDiag_MCF)]
koelsche_meta$Diagnosis[grep("MPNST", koelsche_meta$Diagnosis)] <- "MPNST"
koelsche_meta$Diagnosis <- gsub("Malignant peripheral nerve sheath tumour", "MPNST",
  koelsche_meta$Diagnosis )
koelsche_meta$Diagnosis <- gsub("Neurofibroma \\(plexiform\\)", "pNF",
  koelsche_meta$Diagnosis )

koelsche_meta$Manifestation[koelsche_meta$Manifestation == ""] <- NA

koelsche_meta$Batch.year <- gsub("^.*-", "20", koelsche_meta$Batch)
koelsche_meta$Batch.year[!is.na(koelsche_meta$Batch.date)] <- 
  gsub("^.*-", "20", koelsche_meta$Batch.date[!is.na(koelsche_meta$Batch.date)])
koelsche_meta$Batch.year[koelsche_meta$Batch.year == ""] <- NA

koelsche_meta$Site[!is.na(koelsche_meta$Location)] <- koelsche_meta$Location[
  !is.na(koelsche_meta$Location)]
koelsche_meta$Site <- gsub(" ", "_", gsub(", .*$", "", koelsche_meta$Site))
koelsche_meta$Site[koelsche_meta$Site == ""] <- NA

koelsche_meta$Methylation.Class.Name[
  grep("MPNST-like", koelsche_meta$Methylation.Class.Name)] <- "SARC"
koelsche_meta$Methylation.Class.Name[
  grep("sheath", koelsche_meta$Methylation.Class.Name)] <- "MPNST"
koelsche_meta$Methylation.Class.Name[
  !is.na(koelsche_meta$V12.2_MaxCalDiag_MCF)] <- koelsche_meta$V12.2_MaxCalDiag_MCF[
    !is.na(koelsche_meta$V12.2_MaxCalDiag_MCF)]
koelsche_meta$Methylation.Class.Name[koelsche_meta$Methylation.Class.Name == 
  "methylation class malignant peripheral nerve sheath tumour"] <- "MPNST"
koelsche_meta$Methylation.Class.Name[koelsche_meta$Methylation.Class.Name == 
  "methylation class sarcoma (MPNST-like)"] <- "SARC"
koelsche_meta$Methylation.Class.Name[koelsche_meta$Methylation.Class.Name == 
  "SARC (MPNST-like)"] <- "SARC"
koelsche_meta$Methylation.Class.Name[koelsche_meta$Methylation.Class.Name == 
  "methylation class neurofibroma (plexiform)"] <- "pNF"

koelsche_meta$DNA[!is.na(koelsche_meta$DNA.preparation)] <- koelsche_meta$DNA.preparation[
  !is.na(koelsche_meta$DNA.preparation)]

koelsche_meta$Supplier <- gsub(" ", "\\.", koelsche_meta$Supplier)
koelsche_meta$Supplier[
  grep("bingen", koelsche_meta$Supplier)] <- "Other.Neuropathology"
koelsche_meta$Supplier[!is.na(koelsche_meta$Supplier.study)] <- 
  gsub(" ", ".", koelsche_meta$Supplier.study[!is.na(koelsche_meta$Supplier.study)])

koelsche_meta$Cellularity <- factor(
  koelsche_meta$Tumour.cell.content..absolute. < 0.5)
levels(koelsche_meta$Cellularity) <- c("> 0.5", "< 0.5")

koelsche_meta$Institutional.diagnosis[
  koelsche_meta$Institutional.diagnosis == "Malignant peripheral nerve sheath tumour"] <- 
  "MPNST"
koelsche_meta$Orig.diagnosis <- gsub(
  ", ", "_", koelsche_meta$Institutional.diagnosis )

koelsche_meta$Orig.diagnosis <- gsub("Sarcoma_NOS", "sarcoma_NOS", 
    koelsche_meta$Orig.diagnosis)
koelsche_meta$Orig.diagnosis[grep("pNF", koelsche_meta$ID)] <- "pNF"
koelsche_meta$Diagnosis[grep("pNF", koelsche_meta$ID)] <- "pNF"

# load array types:
koelsche_array <- read.table(file.path(ref_dir, "koelsche_array_types.tsv"),
  sep = "\t", header=T)

koelsche_meta <- merge(koelsche_meta, koelsche_array, by="ID", all.x=T)

koelsche_meta <- subset(koelsche_meta, select = c(ID, Diagnosis, Set, 
  Orig.diagnosis, V12.2_MaxCalScore_MCF, V12.2.Result, Validation, 
  Manifestation, Site, Cellularity, Methylation.Class.Name, 
  Supplier, DNA, Array, Batch.year ))
colnames(koelsche_meta) <- c("ID", "Diagnosis", "Set", "Orig.diagnosis",
  "Classification.score", "Classification.outcome", "Validation", "Stage",
  "Tissue.type", "Cellularity", "Methylation.class", "Supplier",
  "Preservation.method", "Array", "Batch.year" )

write.table(koelsche_meta, 
  file.path(ref_dir, "formatted_koelsche_MPNST_metadata.tsv"),
  sep="\t", row.names=F, col.names=T )


####################################################################################
### 2. Load and format Lyskjaer metadata ###
####################################################################################

# read and format sample sheet:
meta <- read.table(file.path(ref_dir, "lyskjaer_sample_manifest.txt"), header = T, 
  sep = "\t")
meta <- subset(meta, select = c("Source.Name", "Characteristics.age.", 
  "Characteristics.sex.", "Characteristics.organism.part.", 
  "Characteristics.disease.", "Array.Data.File"))
meta$Array <- sapply(strsplit(meta$Array.Data.File, "_"), function(x) x[2])
meta$Slide <- gsub("_.*$", "", meta$Array.Data.File)
colnames(meta) <- c("ID", "Age", "Sex", "Site", "Disease", "Basename", 
  "Array", "Slide" )
meta$Basename <- sub("([0-9]+_[A-Za-z0-9]+).*", "\\1", meta$Basename)
meta <- meta[!duplicated(meta$Basename),]
meta$ID <- gsub("Case_", "Case", meta$ID)

# create age groups:
meta$Age.group <- factor(meta$Age > 25)
levels(meta$Age.group) <- c("<25 y.o.", ">25 y.o.")

# add additional MPNST metadata:
paper_meta <- read.table(
  file.path(ref_dir, "lyskjaer_metadata.tsv"), sep = "\t", header = T)

# format columns:
paper_meta$ID <- paste0("Case", paper_meta$Case)

paper_meta$Diagnosis <- paper_meta$Tumour.subtype
paper_meta$Revised.histo.group <- sapply(
  strsplit(paper_meta$Histopathological.groups, " "), function(x) x[3] )
paper_meta$Diagnosis <- gsub("Sarcoma, NOS", "sarcoma_NOS", 
    paper_meta$Diagnosis)
  paper_meta$Diagnosis <- gsub("MPNST RT-associated", "RT_assoc_MPNST", 
    paper_meta$Diagnosis)
  paper_meta$Diagnosis <- gsub("MPNST ", "MPNST", paper_meta$Diagnosis)

paper_meta$Nerve.origin <- gsub("Y\\. ", "", paper_meta$Arising.from.a.Nerve)
paper_meta$Nerve.origin <- gsub("y", "Y", paper_meta$Nerve.origin)
paper_meta$Grade <- paste0("grade.", paper_meta$Grade)

paper_meta$H3K27me3 <- paste0("H3K27me3.", paper_meta$H3K27me3)

paper_meta$SUZ12.Genomics <- gsub("Normal", "normal", 
  paper_meta$SUZ12.Genomics )
paper_meta$SUZ12.Genomics <- gsub("LOH and.*$", "LOH_and_mutation",
  paper_meta$SUZ12.Genomics )
paper_meta$SUZ12.Genomics[
  !(grepl("normal|LOH", paper_meta$SUZ12.Genomics)) & 
  !is.na(paper_meta$SUZ12.Genomics) ] <- "mutation"

paper_meta$EED.Genomics <- gsub("Normal", "normal", 
  paper_meta$EED.Genomics )
paper_meta$EED.Genomics <- gsub("LOH and.*$", "LOH_and_mutation",
  paper_meta$EED.Genomics )
paper_meta$EED.Genomics <- gsub("LOH .*$", "LOH", paper_meta$EED.Genomics)
paper_meta$EED.Genomics[
  !(grepl("normal|LOH", paper_meta$EED.Genomics)) & 
  !is.na(paper_meta$EED.Genomics) ] <- "mutation"

paper_meta$EZH2.Genomics <- gsub("Normal", "normal", 
  paper_meta$EZH2.Genomics )

paper_meta$TP53.Genomics <- gsub("Normal", "normal", 
  paper_meta$TP53.Genomics )
paper_meta$TP53.Genomics <- gsub("LOH .*$", "LOH_and_mutation",
  paper_meta$TP53.Genomics )
paper_meta$TP53.Genomics[
  !(grepl("normal|LOH", paper_meta$TP53.Genomics)) & 
  !is.na(paper_meta$TP53.Genomics) ] <- "mutation"

paper_meta$EZH2.Genomics <- gsub("Normal", "normal", 
  paper_meta$EZH2.Genomics )
paper_meta$EZH2.Genomics <- gsub("LOH and.*$", "LOH_and_mutation",
  paper_meta$EZH2.Genomics )
paper_meta$EZH2.Genomics <- gsub("LOH .*$", "LOH", paper_meta$EZH2.Genomics)
paper_meta$EZH2.Genomics[
  !(grepl("normal|LOH", paper_meta$EZH2.Genomics)) & 
  !is.na(paper_meta$EZH2.Genomics) ] <- "mutation"

paper_meta$KDM6B.Genomics <- gsub("Normal", "normal", 
  paper_meta$KDM6B.Genomics )
paper_meta$KDM6B.Genomics <- gsub("LOH .*$", "LOH_and_mutation", 
  paper_meta$KDM6B.Genomics )

paper_meta$SSTR2.Amplification <- gsub("YES", "Y", 
  paper_meta$SSTR2.Amplification )
paper_meta$SSTR2.Amplification <- gsub("NO", "N", 
  paper_meta$SSTR2.Amplification )

paper_meta$INI1 <- gsub("Retained", "retained", paper_meta$INI.1.status)
paper_meta$INI1 <- gsub("Retained", "retained", paper_meta$INI.1.status)

# indicate NF1 if NF1 phenotype or mutation recorded, for filtering:
paper_meta$NF1.assoc <- "N"
paper_meta$NF1.assoc[paper_meta$NF1..Phenotype. == "Y" | paper_meta$NF1.germline.mutation == "Y" |
  !(paper_meta$NF1.Somatic.mutation == "Normal" | 
  is.na(paper_meta$NF1.Somatic.mutation))] <- "Y"

# annotate samples with NF1 mutation status:
paper_meta$NF1.somatic.mutation <- "none.detected"
paper_meta$NF1.somatic.mutation[!(paper_meta$NF1.Somatic.mutation == "Normal" | 
  is.na(paper_meta$NF1.Somatic.mutation))] <- "mutation.detected"

colnames(paper_meta) <- gsub("NF1.germline", "NF1.Germline", colnames(paper_meta))
paper_meta$NF1.germline.mutation <- "none.detected"
paper_meta$NF1.germline.mutation[!(paper_meta$NF1.Germline.mutation == "N" | 
  paper_meta$NF1.Germline.mutation == "N " |
  is.na(paper_meta$NF1.Germline.mutation))] <- "mutation.detected"

paper_meta$Background.NF <- gsub("Neurofibroma", "present", 
  paper_meta$Background.component )
paper_meta$Background.NF[
  grep("[a,A]typical", paper_meta$Background.component)] <- "atypical"
paper_meta$Background.NF[
  grep("[p,P]lexiform", paper_meta$Background.component)] <- "plexiform"

# subset and merge metadata:
paper_meta <- subset(paper_meta, select = c(ID, Diagnosis, Grade, Nerve.origin, 
  Revised.histo.group, NF1.assoc, Background.NF, NF1.somatic.mutation, 
  NF1.germline.mutation, NF1..Phenotype., H3K27me3, SUZ12.Genomics, 
  EED.Genomics, EZH2.Genomics, TP53.Genomics, KDM6B.Genomics,
  SSTR2.Amplification, SSTR2.IHC.intensity, INI1 ))
meta <- merge(meta, paper_meta, by="ID", all.x=T)

# add methylation group data, keeping only those in first metadata table (additional
# samples in second are not MPNST):
meth_meta <- read.table(file.path(ref_dir, "lyskjaer_methylation_groups.tsv"), 
  sep = "\t", header = T )
meth_meta$ID <- paste0("Case", meth_meta$Case)
meth_meta <- subset(meth_meta, select = c(ID, Histopathological.groups, MeGroup))
colnames(meth_meta) <- c("ID", "Orig.histo.group", "MeGroup")
all_meta <- merge(meta, meth_meta, by="ID", all.x=T)

# add other metadata:
gd_meta <- read.table(file.path(ref_dir, "lyskjaer_gd_groups.tsv"), 
  sep = "\t", header = T )
gd_meta$ID <- paste0("Case", gd_meta$Case)
gd_meta$Genome.doubling <- "N"
gd_meta$Genome.doubling[gd_meta$X1xGD] <- "1x"
gd_meta$Genome.doubling[gd_meta$X2xGD] <- "2x"
gd_meta <- subset(gd_meta, select = c(ID, Genome.doubling, Dead, Survival.Days))

# stratify survival days into groups:
gd_meta$Survival.Days <- factor(gd_meta$Survival.Days > 1500)
levels(gd_meta$Survival.Days) <- c("< 1500", "> 1500")

all_meta <- merge(all_meta, gd_meta, by="ID", all.x=T)

# add RT induced column:
all_meta$RT.induced <- "N"
all_meta$RT.induced[grepl("RT", all_meta$Diagnosis)] <- "Y"

# keep non-epithelioid MPNST and benign nerve sheath tumour histo groups only:
dim(all_meta)
filt_meta <- all_meta[all_meta$Revised.histo.group %in% c("0", "1A", "1B", "2"),]
dim(filt_meta)

# change histo group names:
filt_meta$Revised.histo.group <- factor(filt_meta$Revised.histo.group)
levels(filt_meta$Revised.histo.group) <- c("benign_or_grade1_MPNST", "classical_MPNST", 
  "classicalhetero_MPNST", "nonclassical_MPNST")

filt_meta$Orig.histo.group <- gsub(" \\(.*$", "", filt_meta$Orig.histo.group)
filt_meta$Orig.histo.group <- gsub(" $", "", filt_meta$Orig.histo.group)
filt_meta$Orig.histo.group <- factor(filt_meta$Orig.histo.group)
levels(filt_meta$Orig.histo.group) <- c("benign_or_grade1_MPNST", "classical_MPNST", 
  "classicalhetero_MPNST", "nonclassical_MPNST", "sarcoma_NOS", "UPS")

# order by case no:
filt_meta <- filt_meta[naturalorder(filt_meta$ID),]
filt_meta$ID <- paste("lyskjaer", gsub("Case", "", filt_meta$ID), "MPNST", sep="_")
colnames(filt_meta) <- gsub("\\.\\.Phenotype\\.", "\\.phenotype", 
  colnames(filt_meta))
colnames(filt_meta) <- gsub("Dead", "Deceased", colnames(filt_meta))
colnames(filt_meta) <- gsub("Survival.Days", "Survival.rate.days", 
  colnames(filt_meta))
colnames(filt_meta) <- gsub("\\.Genomics", "", colnames(filt_meta))
colnames(filt_meta) <- gsub("Amplification", "amp", colnames(filt_meta))

# subset:
filt_meta <- subset(filt_meta, select = -c(Disease, Age))

write.table(filt_meta, 
  file.path(ref_dir, "formatted_lyskjaer_MPNST_metadata.tsv"), sep = "\t", 
  quote = F,
  row.names = F,
  col.names = T )


####################################################################################
### 3. Load and format TCGA metadata ###
####################################################################################

# load TCGA metadata:
TCGA_meta <- read.table(
  file.path(ref_dir, "sarc_tcga_clinical_data.tsv"), sep = "\t", header = T )

# filter for columns with more than 2 non-NA entries:
TCGA_filt <- TCGA_meta[,
  which(apply(TCGA_meta, 2, function(x) length(unique(x[!is.na(x)]))) > 1)]

TCGA_filt$ID2 <- gsub("-", "", 
  paste0(gsub("TCGA-", "TCGA_", TCGA_filt$Sample.ID), "_MPNST") )

TCGA_filt$Sex <- gsub("Female", "female", gsub("Male", "male", TCGA_filt$Sex))

TCGA_filt$Age.group <- factor(TCGA_filt$Diagnosis.Age < 25)
levels(TCGA_filt$Age.group) <- c(">25 y.o.", "<25 y.o.")

TCGA_filt$Age.group2 <- factor(TCGA_filt$Diagnosis.Age < 50)
levels(TCGA_filt$Age.group2) <- c(">50 y.o.", "<50 y.o.")

TCGA_filt$Days.till.collection <- factor(
  TCGA_filt$Days.to.Sample.Collection. < 600)
levels(TCGA_filt$Days.till.collection) <- c("> 600", "< 600")

TCGA_filt$Disease.status <- gsub("[0,1]:|\\/", "", TCGA_filt$Disease.Free.Status)

TCGA_filt$Frac.genome.alt <- factor(TCGA_filt$Fraction.Genome.Altered < 0.5)
levels(TCGA_filt$Frac.genome.alt) <- c("> 0.5", "< 0.5")

TCGA_filt$Mutation.count <- factor(TCGA_filt$Mutation.Count < 30)
levels(TCGA_filt$Mutation.count) <- c("> 30", "< 30")

TCGA_filt$Lesion.length <- TCGA_filt$Nte.lesion.radiologic.length
TCGA_filt$Lesion.length[is.na(TCGA_filt$Lesion.length)] <- 
  TCGA_filt$Nte.lesion.pathologic.length[is.na(TCGA_filt$Lesion.length)]
TCGA_filt$Lesion.length <- factor(as.numeric(TCGA_filt$Lesion.length) < 10)
levels(TCGA_filt$Lesion.length) <- c("> 10", "< 10")

TCGA_filt$Lesion.depth <- TCGA_filt$Nte.lesion.radiologic.depth
TCGA_filt$Lesion.depth[is.na(TCGA_filt$Lesion.depth)] <- 
  TCGA_filt$Nte.lesion.pathologic.depth[is.na(TCGA_filt$Lesion.depth)]
TCGA_filt$Lesion.depth <- sapply(
  strsplit(TCGA_filt$Lesion.depth, "\\|"), function(x) x[1])
TCGA_filt$Lesion.depth <- factor(as.numeric(TCGA_filt$Lesion.depth) < 5)
levels(TCGA_filt$Lesion.depth) <- c("> 5", "< 5")

TCGA_filt$Months.survival <- factor(
  as.numeric(TCGA_filt$Overall.Survival..Months.) < 25)
levels(TCGA_filt$Months.survival) <- c("> 25", "< 25")

TCGA_filt$Survival <- gsub("[0,1]:", "", TCGA_filt$Overall.Survival.Status)

TCGA_filt$Necrosis <- "extensive"
TCGA_filt$Necrosis[grep("no necrosis", TCGA_filt$Tumor.total.necrosis)] <- "none"
TCGA_filt$Necrosis[grep("focal", TCGA_filt$Tumor.total.necrosis)] <- "focal"

TCGA_filt$SUZ12 <- c(rep("normal", 8), "mutation")

TCGA_meta <- subset(TCGA_filt, select = c(ID2, Sex, Age.group, Age.group2,
  Frac.genome.alt, Mutation.count, Lesion.length, Lesion.depth, Tissue.Source.Site,
  Mpnst.nf.familial.or.sporadic, SUZ12, Necrosis, Person.Neoplasm.Status, Disease.status, 
  Disease.Surgical.Margin.Status, Survival, Months.survival, 
  Adjuvant.Postoperative.Pharmaceutical.Therapy.Administered.Indicator, 
  Did.patient.start.adjuvant.postoperative.radiotherapy., 
  New.Neoplasm.Event.Post.Initial.Therapy.Indicator, Days.till.collection ))

colnames(TCGA_meta) <- c("ID", "Sex", "Age.group", "Age.group2", "Frac.genome.alt", 
  "Mutation.count", "Lesion.length", "Lesion.depth", "Tissue.source.site",
  "Familial.or.sporadic", "SUZ12", "Necrosis", "Tumour.status", "Disease.status", 
  "Surgical.margin", "Survival.status", "Survival.rate.months",
  "ACT", "ART", "Post.therapy.neoplasm", "Days.till.collection" )

write.table(TCGA_meta, 
  file.path(ref_dir, "formatted_TCGA_MPNST_metadata.tsv"), sep = "\t", 
  quote = F,
  row.names = F,
  col.names = T )

