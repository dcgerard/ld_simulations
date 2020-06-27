#####################
## Extact SNPs from McAllister Data
#####################

library(stringr)
library(VariantAnnotation)
anodf <- read.csv("./data/gerardii/McAllister_Miller_Locality_Ploidy_Info.csv")
fl <-"./data/gerardii/McAllister.Miller.all.mergedRefGuidedSNPs.vcf.gz"

## choose arbitrary region
chlist <- list(chr1_gr = GRanges("1", IRanges(start = 7000000, end = 7100000)),
               chr2_gr = GRanges("10", IRanges(start = 7000000, end = 7100000)))

compressVcf <- bgzip(fl, tempfile())
idx <- indexTabix(compressVcf, "vcf")
tab <- TabixFile(compressVcf, idx)
for (i in seq_along(chlist)) {
  param <- ScanVcfParam(which = chlist[[i]])
  mca <- readVcf(tab, as.character(i), param)

  ## Keep only biallelic snps
  which_ba <- sapply(alt(mca), length) == 1
  mca <- mca[which_ba, ]

  ## Remove SNPs with low MAF
  which_maf <- info(mca)$AF > 0.1 & info(mca)$AF < 0.9
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

  if (i == 1) {
    sizemat_f <- DP
    refmat_f <- refmat
    locdf_f <- data.frame(snp = rownames(DP), loc = i)
  } else {
    sizemat_f <- rbind(sizemat_f, DP)
    refmat_f <- rbind(refmat_f, refmat)

    locdf <- data.frame(snp = rownames(DP), loc = i)
    locdf_f <- rbind(locdf_f, locdf)
  }
}

## Remove snps with high missingness
goodsnp <- rowMeans(sizemat_f, na.rm = TRUE) >= 3
sizemat_f <- sizemat_f[goodsnp, ]
refmat_f <- refmat_f[goodsnp, ]
locdf_f <- locdf_f[goodsnp, ]

## remove individuals with high missingness
goodind <- str_split_fixed(colnames(sizemat_f), pattern = ":", n = 4)[, 1] %in% anodf$Individual
sizemat_f <- sizemat_f[, goodind]
refmat_f <- refmat_f[, goodind]

## split individuals based on ploidy
sixind <- anodf$Individual[anodf$Ploidy.Level == 6]
nonind <- anodf$Individual[anodf$Ploidy.Level == 9]
candidate <- str_split_fixed(colnames(sizemat_f), pattern = ":", n = 4)[, 1]
stopifnot(candidate %in% anodf$Individual)

which_six <- candidate %in% sixind
which_non <- candidate %in% nonind

sizemat_six <- sizemat_f[, which_six]
refmat_six <- refmat_f[, which_six]

sizemat_non <- sizemat_f[, which_non]
refmat_non <- refmat_f[, which_non]

## Remove duplicated rows
which_bad_six <- duplicated(sizemat_six) & duplicated(refmat_six)
sizemat_six <- sizemat_six[!which_bad_six, ]
refmat_six  <- refmat_six[!which_bad_six, ]

which_bad_non <- duplicated(sizemat_non) & duplicated(refmat_non)
sizemat_non <- sizemat_non[!which_bad_non, ]
refmat_non  <- refmat_non[!which_bad_non, ]

locdf_f <- locdf_f[!which_bad_non, ]

saveRDS(object = sizemat_six, file = "./output/mca/sizemat_hex.RDS")
saveRDS(object = refmat_six, file = "./output/mca/refmat_hex.RDS")
saveRDS(object = sizemat_non, file = "./output/mca/sizemat_non.RDS")
saveRDS(object = refmat_non, file = "./output/mca/refmat_non.RDS")
write.csv(x = locdf_f, file = "./output/mca/locdf.csv", row.names = FALSE)
