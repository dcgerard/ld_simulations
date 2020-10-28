############################
##  Extract Uitdewilligen SNPS
############################


library(vcfR)
library(tidyverse)
uitvcf <- read.vcfR(file = "./data/NewPlusOldCalls.headed.vcf")
uitdf <- read_csv2("./data/CSV-file S1 - Sequence variants filtered DP15.csv",
                   na = "#N/A")

## Filter for a single contig
set.seed(1)
nsnp <- 40
region1_start <- sample(seq_len(nrow(uitdf)), 1)
region1_end <- region1_start + nsnp - 1

region2_start <- sample(seq_len(nrow(uitdf)), 1)
region2_end <- region2_start + nsnp - 1

# subset
subuitdf <- uitdf[c(region1_start:region1_end, region2_start:region2_end), ]
subuitdf$region <- c(rep(1, nsnp), rep(2, nsnp))

which_keep <- getFIX(uitvcf)[, "ID"] %in% subuitdf$`Variant name`
stopifnot(sum(which_keep) == nrow(subuitdf))


refmat  <- extract.gt(uitvcf, element = "RA")
class(refmat) <- "numeric"
sizemat <- extract.gt(uitvcf, element = "DP")
class(sizemat) <- "numeric"

subref <- refmat[which_keep, ]
subsize <- sizemat[which_keep, ]

stopifnot(rownames(subref) == subuitdf$`Variant name`)
stopifnot(rownames(subsize) == subuitdf$`Variant name`)

write_csv(x = subuitdf, file = "./output/uit/uit_suc.csv")
saveRDS(object = subref, file = "./output/uit/refmat_suc.RDS")
saveRDS(object = subsize, file = "./output/uit/sizemat_suc.RDS")
