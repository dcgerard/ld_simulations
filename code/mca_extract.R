#####################
## Extact SNPs from McAllister Data
#####################

library(stringr)
library(VariantAnnotation)
anodf <- read.csv("./data/gerardii/McAllister_Miller_Locality_Ploidy_Info.csv")
fl <-"./data/gerardii/McAllister.Miller.all.mergedRefGuidedSNPs.vcf.gz"

## choose arbitrary region
start_loc <- 1000000
end_loc <- 2000000
chr1_gr <- GRanges("1", IRanges(start = start_loc, end = end_loc))
param <- ScanVcfParam(which = chr1_gr)
compressVcf <- bgzip(fl, tempfile())
idx <- indexTabix(compressVcf, "vcf")
tab <- TabixFile(compressVcf, idx)
mca <- readVcf(tab, "1", param)

## Keep only biallelic snps
which_ba <- sapply(alt(mca), length) == 1
mca <- mca[which_ba, ]

## Remove SNPs with low MAF
which_maf <- info(mca)$AF > 0.05 & info(mca)$AF < 0.95
stopifnot(length(table(sapply(which_maf, length))) == 1)
which_maf <- unlist(which_maf)
mca <- mca[which_maf, ]

## Extract read-count matrices
DP <- geno(mca)$DP
AD <- geno(mca)$AD
stopifnot(length(table(sapply(AD, length))) == 2)

get_elem <- function(x, num) {
  if (length(x) < num) {
    return(NA)
  } else {
    return(x[[num]])
  }
}

refmat <- sapply(AD, get_elem, num = 1)
dim(refmat) <- dim(AD)
dimnames(refmat) <- dimnames(AD)
altmat <- sapply(AD, get_elem, num = 2)
dim(altmat) <- dim(AD)
dimnames(altmat) <- dimnames(AD)

## Remove snps with high missingness
goodsnp <- rowMeans(is.na(DP)) < 0.5

DP <- DP[goodsnp, ]
refmat <- refmat[goodsnp, ]
altmat <- altmat[goodsnp, ]

## remove individuals with high missingness
goodind <- colMeans(is.na(DP)) < 0.5

DP <- DP[, goodind]
refmat <- refmat[, goodind]
altmat <- altmat[, goodind]
stopifnot(all(DP == refmat + altmat, na.rm = TRUE))

## split individuals based on ploidy
sixind <- anodf$Individual[anodf$Ploidy.Level == 6]
nonind <- anodf$Individual[anodf$Ploidy.Level == 9]
candidate <- str_split_fixed(colnames(DP), pattern = ":", n = 4)[, 1]
stopifnot(candidate %in% anodf$Individual)

which_six <- candidate %in% sixind
which_non <- candidate %in% nonind

sizemat_six <- DP[, which_six]
refmat_six <- refmat[, which_six]

sizemat_non <- DP[, which_non]
refmat_non <- refmat[, which_non]

saveRDS(object = sizemat_six, file = "./output/mca/sizemat_hex.RDS")
saveRDS(object = refmat_six, file = "./output/mca/refmat_hex.RDS")
saveRDS(object = sizemat_non, file = "./output/mca/sizemat_non.RDS")
saveRDS(object = refmat_non, file = "./output/mca/refmat_non.RDS")
