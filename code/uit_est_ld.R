###############
## Estimate LD
###############

library(updog)
library(ldsep)

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

mout <- readRDS("./output/uit/uit_updog_fit.RDS")

genolike_array <- format_multidog(x = mout, varname = paste0("logL_", 0:4))
postmean_mat <- format_multidog(x = mout, varname = "postmean")
genomat <- format_multidog(x = mout, varname = "geno")
K <- 4

ld_hap_genolike <- mldest(geno = genolike_array,
                          K    = K,
                          nc   = nc,
                          type = "hap",
                          se   = TRUE)

saveRDS(object = ld_hap_genolike, file = "./output/uit/ldest_hap_genolike.RDS")

ld_hap_geno <- mldest(geno = genomat,
                      K    = K,
                      nc   = nc,
                      type = "hap",
                      se   = TRUE)

saveRDS(object = ld_hap_geno, file = "./output/uit/ldest_hap_geno.RDS")

ld_comp_genolike <- mldest(geno  = genolike_array,
                           K     = K,
                           nc    = nc,
                           type  = "comp",
                           model = "norm",
                           se    = TRUE)

saveRDS(object = ld_comp_genolike, file = "./output/uit/ldest_comp_genolike.RDS")

ld_comp_genolike_flex <- mldest(geno  = genolike_array,
                                K     = K,
                                nc    = nc,
                                type  = "comp",
                                model = "flex",
                                se    = FALSE)

saveRDS(object = ld_comp_genolike_flex, file = "./output/uit/ldest_comp_genolike_flex.RDS")

ld_comp_geno <- mldest(geno = postmean_mat,
                       K    = K,
                       nc   = nc,
                       type = "comp",
                       se   = TRUE)

saveRDS(object = ld_comp_geno, file = "./output/uit/ldest_comp_geno.RDS")
