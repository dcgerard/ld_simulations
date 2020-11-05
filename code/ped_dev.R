##################
## Plots demonstrating Deviation from HWE
## Also generate true correlation matrix
##################

library(corrplot)
library(xtable)
source("./code/ped_funs.R")
ploidy <- 4
seed <- 1
nind <- 10000
chrome_size <- 100
snp_to_keep <- c(50, 51, 60, 70, 80, 90, 99, 100)
ngen <- 10
alpha1 <- 0.9

set.seed(seed)
## generate parameter file
gen_par(ploidy = ploidy,
        seed = seed,
        file = "./PedigreeSim/sim.par")

## generate pedigree
gen_ped(nind = nind,
        ngen = ngen,
        file = "./PedigreeSim/sim.ped")

## generate map file
gen_map(snp_loc = snp_to_keep,
        file = "./PedigreeSim/sim.map")

## generate founder genotype
gen_sim(nind = nind,
        snp_loc = snp_to_keep,
        alpha1 = alpha1,
        file = "./PedigreeSim/sim.gen")

pref_pair_vec <- c(1/3, 2/3, 1)
quadprop_vec <- c(0, 1/3, 2/3)
paramdf <- expand.grid(pref_pair = pref_pair_vec, quadprop = quadprop_vec)

## holds the true correlation matrix
r2_true_list <- list()

## holds the genotype frequencies
hwe_df <- data.frame(pp = rep(NA_real_, 2 * length(snp_to_keep) + 2),
                     qq = rep(NA_real_, 2 * length(snp_to_keep) + 2),
                     loc = rep(NA_real_, 2 * length(snp_to_keep) + 2),
                     `0` = rep(NA_real_, 2 * length(snp_to_keep) + 2),
                     `1` = rep(NA_real_, 2 * length(snp_to_keep) + 2),
                     `2` = rep(NA_real_, 2 * length(snp_to_keep) + 2),
                     `3` = rep(NA_real_, 2 * length(snp_to_keep) + 2),
                     `4` = rep(NA_real_, 2 * length(snp_to_keep) + 2))
hwe_df[1, ] <- c(0, 0, NA, dbinom(0:ploidy, ploidy, 0.5))
ygen1 <- dbinom(x = 0:(ploidy/2), size = ploidy/2, prob = 0.1)
ygen2 <- dbinom(x = 0:(ploidy/2), size = ploidy/2, prob = 0.9)
hwe_df[2, ] <- c(1, 0, NA, convolve(x = ygen1, y = rev(ygen2), type = "open"))

for (index in seq_len(nrow(paramdf))) {
  pref_pair <- paramdf$pref_pair[[index]]
  quadprop <- paramdf$quadprop[[index]]

  ## generate chromosome file
  gen_chrom(length = chrome_size,
            centromere = round(chrome_size / 2),
            prefPairing = pref_pair,
            quadrivalents = quadprop,
            file = "./PedigreeSim/sim.chrom")

  ## Run PedigreeSim
  system("cd PedigreeSim; java -jar PedigreeSim.jar sim.par")

  dose <- read.table("./PedigreeSim/sim_out_alleledose.dat", header = TRUE)
  genomat <- as.matrix(dose[paste0("F", ngen-1, "_", seq_len(nind))])
  rownames(genomat) <- dose$marker

  r2_true_list[[index]] <- cor(t(genomat))^2

  pdf(file = paste0("./output/ped/true_r/heatmap_pp",
                    round(100 * pref_pair),
                    "_qq",
                    round(100 * quadprop),
                    ".pdf"),
      width = 3,
      height = 3,
      family = "Times")

  corrplot(corr = r2_true_list[[index]],
           method = "color",
           type = "upper",
           diag = FALSE)
  dev.off()

  ## Get genotype distributions
  geno50 <- c(genomat[1, ], 0:ploidy) ## subtract off one from each group in table later
  geno100 <- c(genomat[length(snp_to_keep), ], 0:ploidy) ## subract off one from each group in table later

  hwe_df[index * 2 + 1, 1] <- paramdf$pref_pair[[index]]
  hwe_df[index * 2 + 2, 1] <- paramdf$pref_pair[[index]]

  hwe_df[index * 2 + 1, 2] <- paramdf$quadprop[[index]]
  hwe_df[index * 2 + 2, 2] <- paramdf$quadprop[[index]]

  hwe_df[index * 2 + 1, 3] <- 50
  hwe_df[index * 2 + 2, 3] <- 100

  hwe_df[index * 2 + 1, 4:8] <- prop.table(table(geno50) - 1) ## since added quasi counts
  hwe_df[index * 2 + 2, 4:8] <- prop.table(table(geno100) - 1) ## since added quasi counts
}

write.csv(x = hwe_df, file = "./output/ped/true_r/genodist.csv", row.names = FALSE)
saveRDS(object = r2_true_list, file = "./output/ped/true_r/true_r_list.RDS")
