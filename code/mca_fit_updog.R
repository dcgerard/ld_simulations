#################################
## Fit updog on McAllister Snps
#################################

library(updog)

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

refhex <- readRDS("./output/mca/refmat_hex.RDS")
sizehex <- readRDS("./output/mca/sizemat_hex.RDS")
mouthex <- multidog(refmat  = refhex,
                    sizemat = sizehex,
                    ploidy  = 6,
                    model   = "norm",
                    nc      = nc,
                    bias_init = 1,
                    seq = 0.001,
                    od = 0.01,
                    update_bias = FALSE,
                    update_od = FALSE,
                    update_seq = FALSE)

# mouthex <- filter_snp(x = mouthex, expr = (alpha > 0.1) & (alpha < 0.9))

saveRDS(object = mouthex, file = "./output/mca/updog_fits_hex.RDS")

refnon <- readRDS("./output/mca/refmat_non.RDS")
sizenon <- readRDS("./output/mca/sizemat_non.RDS")
moutnon <- multidog(refmat  = refnon,
                    sizemat = sizenon,
                    ploidy  = 9,
                    model   = "norm",
                    nc      = nc,
                    bias_init = 1,
                    seq = 0.001,
                    od = 0.01,
                    update_bias = FALSE,
                    update_od = FALSE,
                    update_seq = FALSE)

# moutnon <- filter_snp(x = moutnon, expr = (alpha > 0.1) & (alpha < 0.9))

saveRDS(object = moutnon, file = "./output/mca/updog_fits_non.RDS")
