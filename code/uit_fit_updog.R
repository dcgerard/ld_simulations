library(updog)
library(tidyverse)

# Number of threads to use for multithreaded computing. This must be
# specified in the command-line shell; e.g., to use 8 threads, run
# command
#
#  R CMD BATCH '--args nc=8' mouthwash_sims.R
#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nc <- 1
} else {
  eval(parse(text = args[[1]]))
}
cat(nc, "\n")

refmat <- readRDS("./output/uit/refmat_suc.RDS")
sizemat <- readRDS("./output/uit/sizemat_suc.RDS")
ploidy <- 4

mout <- multidog(refmat  = refmat,
                 sizemat = sizemat,
                 ploidy  = ploidy,
                 model   = "norm",
                 nc      = nc)

# Filter out monoallelic SNPs
moutsub <- filter_snp(mout, pmax(Pr_0, Pr_1, Pr_2, Pr_3, Pr_4) < 0.90)

uitdf <- read_csv("./output/uit/uit_suc.csv")
subuitdf <- filter(uitdf, uitdf$`Variant name` %in% moutsub$snpdf$snp)
stopifnot(nrow(subuitdf) == nrow(moutsub$snpdf))

write_csv(x = subuitdf, file = "./output/uit/subuit_suc.csv")
saveRDS(object = moutsub, file = "./output/uit/uit_updog_fit.RDS")
