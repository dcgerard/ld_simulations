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

mout$inddf %>%
  group_by(snp) %>%
  summarize(freq = mean(geno) / 4) %>%
  filter(freq < 0.9, freq > 0.1) ->
  keepdf

moutsub <- filter_snp(mout, snp %in% keepdf$snp)

saveRDS(object = moutsub, file = "./output/uit/uit_updog_fit.RDS")
