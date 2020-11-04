####################
## Simulations using Pedigree Sim
####################

library(updog)
library(ldsep)
source("./code/ped_funs.R")

## Set simulation parameter settings
ploidy        <- 4
nind          <- 1000                        ## population size
ngen          <- 10                          ## number of generations
chrome_size   <- 100                         ## size of chromosome in cM
prob_hom_pair <- c(1/3, 2/3, 1)              ## Probability that homologues pair
pref_pair     <- (3 * prob_hom_pair - 1) / 2 ## PedigreeSim preferential pairing parameter
quadprop      <- c(0, 1/3, 2/3)              ## Probability of quadrivalent pairing. 2/3 is max
alpha1        <- 0.9                         ## allele frequency of founder population on subgenome 1
alpha2        <- 1 - alpha1                  ## allele frequency of founder population on subgenome 2
niter         <- 100                         ## number of iterations per scenario
seq           <- 0.01                        ## sequencing error rate
od            <- 0.01                        ## overdispersion parameters
bias          <- 1                           ## allele bias
readdepth     <- 10                          ## read depth
nsamp         <- 100                         ## sample size

# which pairs of SNPs to compare.
snp_to_keep <- c(50, 51, 60, 70, 80, 90, 99, 100)

simdf <- expand.grid(seed        = seq_len(niter),
                     pref_pair   = pref_pair,
                     quadprop    = quadprop,
                     ploidy      = ploidy,
                     nind        = nind,
                     ngen        = ngen,
                     chrome_size = chrome_size,
                     alpha1      = alpha1,
                     seq         = seq,
                     bias        = bias,
                     od          = od,
                     readdepth   = readdepth,
                     nsamp       = nsamp)

## r^2 with A50
simdf[paste0("true_A50_A", snp_to_keep[-1])] <- NA_real_
simdf[paste0("hapgeno_A50_A", snp_to_keep[-1])] <- NA_real_
simdf[paste0("compmean_A50_A", snp_to_keep[-1])] <- NA_real_
simdf[paste0("haplike_A50_A", snp_to_keep[-1])] <- NA_real_
simdf[paste0("complikeflex_A50_A", snp_to_keep[-1])] <- NA_real_
simdf[paste0("complikenorm_A50_A", snp_to_keep[-1])] <- NA_real_

## r^2 with A100
simdf[paste0("true_A100_A", snp_to_keep[-length(snp_to_keep)])] <- NA_real_
simdf[paste0("hapgeno_A100_A", snp_to_keep[-length(snp_to_keep)])] <- NA_real_
simdf[paste0("compmean_A100_A", snp_to_keep[-length(snp_to_keep)])] <- NA_real_
simdf[paste0("haplike_A100_A", snp_to_keep[-length(snp_to_keep)])] <- NA_real_
simdf[paste0("complikeflex_A100_A", snp_to_keep[-length(snp_to_keep)])] <- NA_real_
simdf[paste0("complikenorm_A100_A", snp_to_keep[-length(snp_to_keep)])] <- NA_real_

for (index in seq_len(nrow(simdf))) {
  set.seed(simdf$seed)
  cat("Iteration:", index, "\n")

  ## generate parameter file
  gen_par(ploidy = simdf$ploidy[[index]],
          seed = simdf$seed[[index]],
          file = "./PedigreeSim/sim.par")

  ## generate pedigree
  gen_ped(nind = simdf$nind[[index]],
          ngen = simdf$ngen[[index]],
          file = "./PedigreeSim/sim.ped")

  ## generate chromosome file
  gen_chrom(length = simdf$chrome_size[[index]],
            centromere = round(simdf$chrome_size[[index]] / 2),
            prefPairing = simdf$pref_pair[[index]],
            quadrivalents = simdf$quadprop[[index]],
            file = "./PedigreeSim/sim.chrom")

  ## generate map file
  gen_map(snp_loc = snp_to_keep,
          file = "./PedigreeSim/sim.map")

  ## generate founder genotype
  gen_sim(nind = simdf$nind[[index]],
          snp_loc = snp_to_keep,
          alpha1 = simdf$alpha1[[index]],
          file = "./PedigreeSim/sim.gen")

  ## Run PedigreeSim
  system("cd PedigreeSim; java -jar PedigreeSim.jar sim.par")

  ## Extract doses
  dose <- read.table("./PedigreeSim/sim_out_alleledose.dat", header = TRUE)
  genomat <- as.matrix(dose[paste0("F", ngen-1, "_", seq_len(nind))])
  rownames(genomat) <- dose$marker

  ## calculate "true" correlation matrix
  r2_true <- cor(t(genomat))^2

  ## subsample SNPs/individuals
  ind_to_keep <- sort(sample(seq_len(nind), size = simdf$nsamp[[index]]))
  sub_geno <- genomat[, ind_to_keep]

  ## generate read counts
  sizemat <- matrix(simdf$readdepth[[index]],
                    nrow = length(snp_to_keep),
                    ncol = simdf$nsamp[[index]])
  refmat <- matrix(NA_real_,
                   nrow = length(snp_to_keep),
                   ncol = simdf$nsamp[[index]])
  for (j in seq_along(snp_to_keep)) {
    refmat[j, ] <- rflexdog(sizevec = sizemat[j, ],
                            geno = sub_geno[j, ],
                            ploidy = simdf$ploidy[[index]],
                            seq = simdf$seq[[index]],
                            bias = simdf$bias[[index]],
                            od = simdf$od[[index]])
  }
  dimnames(sizemat) <- dimnames(sub_geno)
  dimnames(refmat) <- dimnames(sub_geno)

  ## fit updog
  trash <- capture.output({
    mout <- multidog(refmat = refmat,
                     sizemat = sizemat,
                     ploidy = simdf$ploidy[[index]],
                     seq = 0.01,
                     bias_init = 1,
                     od = 0.01,
                     update_seq = FALSE,
                     update_bias = FALSE,
                     update_od = FALSE,
                     model = "norm")
  })

  geno_mode <- format_multidog(x = mout, varname = "geno")
  geno_mean <- format_multidog(x = mout, varname = "postmean")
  geno_like <- format_multidog(x = mout, varname = paste0("logL_", 0:simdf$ploidy[[index]]))

  ## fit ldsep
  hap_geno      <- mldest(geno = geno_mode, K = simdf$ploidy[[index]], type = "hap", se = FALSE)
  comp_mean     <- mldest(geno = geno_mean, K = simdf$ploidy[[index]], type = "comp", se = FALSE)
  hap_like      <- mldest(geno = geno_like, K = simdf$ploidy[[index]], type = "hap", se = FALSE)
  comp_likeflex <- mldest(geno = geno_like, K = simdf$ploidy[[index]], type = "comp", model = "flex", se = FALSE)
  comp_likenorm <- mldest(geno = geno_like, K = simdf$ploidy[[index]], type = "comp", model = "norm", se = FALSE)

  ## compare to "true" r2
  r2_hap_geno <- format_lddf(obj = hap_geno, element = "r2")
  r2_comp_mean <- format_lddf(obj = comp_mean, element = "r2")
  r2_hap_like <- format_lddf(obj = hap_like, element = "r2")
  r2_comp_likeflex <- format_lddf(obj = comp_likeflex, element = "r2")
  r2_comp_likenorm <- format_lddf(obj = comp_likenorm, element = "r2")

  ## Add to data frame
  simdf[index, paste0("true_A50_A", snp_to_keep[-1])]         <- r2_true[1, -1]
  simdf[index, paste0("hapgeno_A50_A", snp_to_keep[-1])]      <- r2_hap_geno[1, -1]
  simdf[index, paste0("compmean_A50_A", snp_to_keep[-1])]     <- r2_comp_mean[1, -1]
  simdf[index, paste0("haplike_A50_A", snp_to_keep[-1])]      <- r2_hap_like[1, -1]
  simdf[index, paste0("complikeflex_A50_A", snp_to_keep[-1])] <- r2_comp_likeflex[1, -1]
  simdf[index, paste0("complikenorm_A50_A", snp_to_keep[-1])] <- r2_comp_likenorm[1, -1]

  simdf[index, paste0("true_A100_A", snp_to_keep[-length(snp_to_keep)])]         <- r2_true[-length(snp_to_keep), length(snp_to_keep)]
  simdf[index, paste0("hapgeno_A100_A", snp_to_keep[-length(snp_to_keep)])]      <- r2_hap_geno[-length(snp_to_keep), length(snp_to_keep)]
  simdf[index, paste0("compmean_A100_A", snp_to_keep[-length(snp_to_keep)])]     <- r2_comp_mean[-length(snp_to_keep), length(snp_to_keep)]
  simdf[index, paste0("haplike_A100_A", snp_to_keep[-length(snp_to_keep)])]      <- r2_hap_like[-length(snp_to_keep), length(snp_to_keep)]
  simdf[index, paste0("complikeflex_A100_A", snp_to_keep[-length(snp_to_keep)])] <- r2_comp_likeflex[-length(snp_to_keep), length(snp_to_keep)]
  simdf[index, paste0("complikenorm_A100_A", snp_to_keep[-length(snp_to_keep)])] <- r2_comp_likenorm[-length(snp_to_keep), length(snp_to_keep)]
}

write.csv(x = simdf, file = "./output/ped/ped_sim_out.csv", row.names = FALSE)
