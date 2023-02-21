install.packages("R.utils")
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE, adhoc = TRUE)
bam <- args$bam
out <- args$out
int <- args$int
tag <- args$tag
roi <- args$roi

BiocManager::install("GenomicAlignments", update = FALSE)
library(GenomicAlignments)

param <- ScanBamParam(what = "flag")
if (!is.null(tag)) bamTag(param) <- tag
if (!is.null(roi)) bamWhich(param) <- GRanges(roi)

ga <- readGAlignments(bam, use.names = TRUE, param = param)
gr <- granges(ga, use.mcols = TRUE)

if (!is.null(tag) && all(is.na(mcols(gr)[[tag]]))) stop(paste(tag, "tag not found"))

## add whether alignment is for first or second read in template
mcols(gr)$R1 <- as.logical(bitwAnd(gr$flag, 0x40))
mcols(gr)$R2 <- as.logical(bitwAnd(gr$flag, 0x80))

## only keep reads with both mates
keep <- intersect(names(gr[gr$R1]), names(gr[gr$R2]))
gr <- gr[names(gr) %in% keep]

## only keep reads with nonsplit alignments
keep <- names(which(lengths(split(gr, names(gr))) == 2L))
gr <- gr[names(gr) %in% keep]

R1 <- gr[gr$R1]
R2 <- gr[gr$R2]
R2 <- R2[names(R1)]
if (!is.null(tag)) stopifnot(identical(mcols(R1)[[tag]], mcols(R2)[[tag]]))

## only keep reads with both mates on the same chromosome
## NB currently not checking F/R orientation
keep <- which(seqnames(R1) == seqnames(R2))
R1 <- R1[keep]
R2 <- R2[keep]

## ignore R1/R2 read orientation consistent with Picard MarkDuplicates
ignore_orientation <- TRUE

## label unique molecules based on mate start positions and tag (if applicable)
if (ignore_orientation) {

    grl <- pc(as(R1, "GRangesList"), as(R2, "GRangesList"))
    grl <- range(grl, ignore.strand = TRUE)
    stopifnot(all(lengths(grl) == 1L))
    molecules <- as.character(unlist(grl))

} else {

    molecules <- paste(flank(R1, -1), flank(R2, -1))

}
if (!is.null(tag)) molecules <- paste(molecules, mcols(R1)[[tag]])

set.seed(0)

N <- length(molecules)
r <- seq_len(floor(N/int)) * int
r <- c(int/10, r, N)
m <- sapply(r, function (x) { length(unique(sample(molecules, x))) })

df <- data.frame(
    reads = r,
    molecules = m,
    reads_per_molecule = signif(r/m, 2),
    stringsAsFactors = FALSE)

write.table(df, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE, file = out)
