################
## Run simulations comparing MLE to moment-based LD estimators
################

set.seed(1)
library(updog)
library(ldsep)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)

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

## Parameters of simulations setting ------------------------------------------
nind   <- c(100) # number of individuals
size   <- c(1, 5, 10, 50, 100) # depth per individual
ploidy <- c(2, 4, 6, 8)

# major allele frequencies and D
mean_prop_seq <- list(c(mu1 = 0.5, mu2 = 0.5),
                      c(mu1 = 0.5, mu2 = 0.75),
                      c(mu1 = 0.9, mu2 = 0.9))
se_prop_seq <- c(0.5, 1)
cor_seq <- c(0, 0.5, 0.9)

bias <- 1 # allele bias
od   <- 0.01 # overdispersion
seq  <- 0.01 # sequencing error

itermax <- 200

get_ldvec <- function(mu1, mu2, sigma11, sigma12, sigma22, K) {
  sigma <- matrix(c(sigma11, sigma12, sigma12, sigma22), nrow = 2, ncol = 2)
  L <- t(chol(sigma))
  par <- c(mu1, mu2, L[lower.tri(L, diag = TRUE)])
  ldvec <- ldsep:::ld_from_pbnorm(par = par, K = K)
  return(c(ldvec[c("D", "r2", "r", "z", "Dprime", "Dprimeg")],
         l11 = L[1, 1],
         l21 = L[2, 1],
         l22 = L[2, 2]))
}

paramdf <- expand.grid(seed   = seq_len(itermax),
                       nind   = nind,
                       size   = size,
                       mean_prop = mean_prop_seq,
                       se_prop = se_prop_seq,
                       cor    = cor_seq,
                       bias   = bias,
                       od     = od,
                       seq    = seq,
                       ploidy = ploidy)
paramdf <- cbind(paramdf, do.call(rbind, paramdf$mean_prop))
paramdf %>%
  select(-mean_prop) %>%
  mutate(mu1 = mu1 * ploidy,
         mu2 = mu2 * ploidy,
         sigma11 = (se_prop * ploidy) ^ 2,
         sigma22 = sigma11,
         sigma12 = cor * sqrt(sigma11) * sqrt(sigma22)) %>%
  select(-cor, -se_prop) %>%
  rowwise() %>%
  mutate(ldvec = list(get_ldvec(mu1 = mu1,
                                mu2 = mu2,
                                sigma11 = sigma11,
                                sigma12 = sigma12,
                                sigma22 = sigma22,
                                K = ploidy))) ->
  paramdf

paramdf <- cbind(paramdf, do.call(rbind, paramdf$ldvec))
paramdf <- select(paramdf, -ldvec)

# populate columns to be filled ----------------------------------------------
paramdf %>%
  mutate(mle_D_est           = NA_real_,
         mle_D_se            = NA_real_,
         mle_Dprime_est      = NA_real_,
         mle_Dprime_se       = NA_real_,
         mle_r2_est          = NA_real_,
         mle_r2_se           = NA_real_,
         mle_r_est           = NA_real_,
         mle_r_se            = NA_real_,
         mle_z_est           = NA_real_,
         mle_z_se            = NA_real_,
         mle_pab_est         = NA_real_,
         mle_pAb_est         = NA_real_,
         mle_paB_est         = NA_real_,
         mle_pAB_est         = NA_real_,
         mle_time            = NA_real_,
         gen_D_est           = NA_real_,
         gen_D_se            = NA_real_,
         gen_Dprime_est      = NA_real_,
         gen_Dprime_se       = NA_real_,
         gen_r2_est          = NA_real_,
         gen_r2_se           = NA_real_,
         gen_z_est           = NA_real_,
         gen_z_se            = NA_real_,
         gen_r_est           = NA_real_,
         gen_r_se            = NA_real_,
         gen_pab_est         = NA_real_,
         gen_pAb_est         = NA_real_,
         gen_paB_est         = NA_real_,
         gen_pAB_est         = NA_real_,
         gen_time            = NA_real_,
         mom_D_est           = NA_real_,
         mom_D_se            = NA_real_,
         mom_Dprime_est      = NA_real_,
         mom_Dprime_se       = NA_real_,
         mom_Dprimeg_est     = NA_real_,
         mom_Dprimeg_se      = NA_real_,
         mom_r2_est          = NA_real_,
         mom_r2_se           = NA_real_,
         mom_z_est           = NA_real_,
         mom_z_se            = NA_real_,
         mom_r_est           = NA_real_,
         mom_r_se            = NA_real_,
         mom_time            = NA_real_,
         com_D_est           = NA_real_,
         com_D_se            = NA_real_,
         com_Dprime_est      = NA_real_,
         com_Dprime_se       = NA_real_,
         com_Dprimeg_est     = NA_real_,
         com_Dprimeg_se      = NA_real_,
         com_r2_est          = NA_real_,
         com_r2_se           = NA_real_,
         com_z_est           = NA_real_,
         com_z_se            = NA_real_,
         com_r_est           = NA_real_,
         com_r_se            = NA_real_,
         com_time            = NA_real_,
         comnorm_D_est       = NA_real_,
         comnorm_D_se        = NA_real_,
         comnorm_Dprime_est  = NA_real_,
         comnorm_Dprime_se   = NA_real_,
         comnorm_Dprimeg_est = NA_real_,
         comnorm_Dprimeg_se  = NA_real_,
         comnorm_r2_est      = NA_real_,
         comnorm_r2_se       = NA_real_,
         comnorm_z_est       = NA_real_,
         comnorm_z_se        = NA_real_,
         comnorm_r_est       = NA_real_,
         comnorm_r_se        = NA_real_,
         comnorm_time        = NA_real_) ->
  paramdf

## shuffle order to equalize computation time across nodes
paramdf <- paramdf[sample(seq_len(nrow(paramdf))), ]

## Register workers ----------------------------------------------------------
if (nc == 1) {
  foreach::registerDoSEQ()
} else {
  cl = parallel::makeCluster(nc)
  doParallel::registerDoParallel(cl = cl)
  if (foreach::getDoParWorkers() == 1) {
    stop(paste0("nc > 1 but only one core registered from ",
                "foreach::getDoParWorkers()."))
  }
}

i <- 1
simdf <- foreach::foreach(i = seq_len(nrow(paramdf)),
                          .combine = rbind,
                          .export = c("ldest",
                                      "rflexdog",
                                      "flexdog",
                                      "%>%",
                                      "gather",
                                      "as_tibble",
                                      "mutate")) %dopar% {
                            set.seed(paramdf$seed[[i]])

                            # generate data -----------------------------------
                            mu <- c(paramdf$mu1[[i]], paramdf$mu2[[i]])
                            sigma <- matrix(c(paramdf$sigma11[[i]],
                                              paramdf$sigma12[[i]],
                                              paramdf$sigma12[[i]],
                                              paramdf$sigma22[[i]]),
                                            nrow = 2,
                                            ncol = 2)
                            distmat <- ldsep::pbnorm_dist(mu = mu,
                                                          sigma = sigma,
                                                          K = paramdf$ploidy[[i]],
                                                          log = FALSE)

                            genomat <- matrix(1:((paramdf$ploidy[[i]] + 1)^2),
                                              ncol = paramdf$ploidy[[i]] + 1,
                                              nrow = paramdf$ploidy[[i]] + 1)
                            colnames(genomat) <- 0:paramdf$ploidy[[i]]
                            rownames(genomat) <- 0:paramdf$ploidy[[i]]
                            genomat %>%
                              as_tibble(rownames = "g1") %>%
                              gather(-g1, key = "g2", value = "ind") %>%
                              mutate(g1 = as.numeric(g1),
                                     g2 = as.numeric(g2)) ->
                              genodf

                            indvec <- sample(x = genodf$ind,
                                             size = paramdf$nind[[i]],
                                             replace = TRUE,
                                             prob = c(distmat))

                            gA <- genodf$g1[indvec]
                            gB <- genodf$g2[indvec]
                            sizevec <- rep(paramdf$size[[i]], paramdf$nind[[i]])
                            refA <- updog::rflexdog(sizevec = sizevec,
                                                    geno    = gA,
                                                    ploidy  = paramdf$ploidy[[i]],
                                                    seq     = paramdf$seq[[i]],
                                                    bias    = paramdf$bias[[i]],
                                                    od      = paramdf$od[[i]])
                            refB <- updog::rflexdog(sizevec = sizevec,
                                                    geno    = gB,
                                                    ploidy  = paramdf$ploidy[[i]],
                                                    seq     = paramdf$seq[[i]],
                                                    bias    = paramdf$bias[[i]],
                                                    od      = paramdf$od[[i]])

                            ## Run updog --------------------------------------
                            foutA <- updog::flexdog(refvec      = refA,
                                                    sizevec     = sizevec,
                                                    ploidy      = paramdf$ploidy[[i]],
                                                    verbose     = FALSE,
                                                    bias_init   = paramdf$bias[[i]],
                                                    update_bias = FALSE,
                                                    seq         = paramdf$seq[[i]],
                                                    update_seq  = FALSE,
                                                    od          = paramdf$od[[i]],
                                                    update_od   = FALSE,
                                                    model       = "hw")
                            foutB <- updog::flexdog(refvec      = refB,
                                                    sizevec     = sizevec,
                                                    ploidy      = paramdf$ploidy[[i]],
                                                    verbose     = FALSE,
                                                    bias_init   = paramdf$bias[[i]],
                                                    update_bias = FALSE,
                                                    seq         = paramdf$seq[[i]],
                                                    update_seq  = FALSE,
                                                    od          = paramdf$od[[i]],
                                                    update_od   = FALSE,
                                                    model       = "hw")

                            ## Estimate ld ------------------------------------
                            tryCatch({
                              paramdf$mle_time[[i]] <- system.time(
                                ldmle <- ldsep::ldest(ga = foutA$genologlike,
                                                      gb = foutB$genologlike,
                                                      K = paramdf$ploidy[[i]],
                                                      type = "gam")
                              )[[3]]

                              paramdf$mle_D_est[[i]]      <- ldmle[["D"]]
                              paramdf$mle_D_se[[i]]       <- ldmle[["D_se"]]
                              paramdf$mle_Dprime_est[[i]] <- ldmle[["Dprime"]]
                              paramdf$mle_Dprime_se[[i]]  <- ldmle[["Dprime_se"]]
                              paramdf$mle_r2_est[[i]]     <- ldmle[["r2"]]
                              paramdf$mle_r2_se[[i]]      <- ldmle[["r2_se"]]
                              paramdf$mle_z_est[[i]]      <- ldmle[["z"]]
                              paramdf$mle_z_se[[i]]       <- ldmle[["z_se"]]
                              paramdf$mle_r_est[[i]]      <- ldmle[["r"]]
                              paramdf$mle_r_se[[i]]       <- ldmle[["r_se"]]
                              paramdf$mle_pab_est[[i]]    <- ldmle[["p_ab"]]
                              paramdf$mle_pAb_est[[i]]    <- ldmle[["p_Ab"]]
                              paramdf$mle_paB_est[[i]]    <- ldmle[["p_aB"]]
                              paramdf$mle_pAB_est[[i]]    <- ldmle[["p_AB"]]
                            }, error = function(e) NULL)

                            tryCatch({
                              paramdf$gen_time[[i]] <- system.time(
                                ldgen <- ldsep::ldest(ga = foutA$geno,
                                                      gb = foutB$geno,
                                                      K = paramdf$ploidy[[i]],
                                                      type = "gam")
                              )[[3]]

                              paramdf$gen_D_est[[i]]      <- ldgen[["D"]]
                              paramdf$gen_D_se[[i]]       <- ldgen[["D_se"]]
                              paramdf$gen_Dprime_est[[i]] <- ldgen[["Dprime"]]
                              paramdf$gen_Dprime_se[[i]]  <- ldgen[["Dprime_se"]]
                              paramdf$gen_r2_est[[i]]     <- ldgen[["r2"]]
                              paramdf$gen_r2_se[[i]]      <- ldgen[["r2_se"]]
                              paramdf$gen_z_est[[i]]      <- ldgen[["z"]]
                              paramdf$gen_z_se[[i]]       <- ldgen[["z_se"]]
                              paramdf$gen_r_est[[i]]      <- ldgen[["r"]]
                              paramdf$gen_r_se[[i]]       <- ldgen[["r_se"]]
                              paramdf$gen_pab_est[[i]]    <- ldgen[["p_ab"]]
                              paramdf$gen_pAb_est[[i]]    <- ldgen[["p_Ab"]]
                              paramdf$gen_paB_est[[i]]    <- ldgen[["p_aB"]]
                              paramdf$gen_pAB_est[[i]]    <- ldgen[["p_AB"]]
                            }, error = function(e) NULL)

                            tryCatch({
                              paramdf$mom_time[[i]] <- system.time(
                                ldmom <- ldsep::ldest(ga = foutA$postmean,
                                                      gb = foutB$postmean,
                                                      K = paramdf$ploidy[[i]],
                                                      type = "comp")
                              )[[3]]

                              paramdf$mom_D_est[[i]]       <- ldmom[["D"]]
                              paramdf$mom_D_se[[i]]        <- ldmom[["D_se"]]
                              paramdf$mom_Dprime_est[[i]]  <- ldmom[["Dprime"]]
                              paramdf$mom_Dprime_se[[i]]   <- ldmom[["Dprime_se"]]
                              paramdf$mom_Dprimeg_est[[i]] <- ldmom[["Dprimeg"]]
                              paramdf$mom_Dprimeg_se[[i]]  <- ldmom[["Dprimeg_se"]]
                              paramdf$mom_r2_est[[i]]      <- ldmom[["r2"]]
                              paramdf$mom_r2_se[[i]]       <- ldmom[["r2_se"]]
                              paramdf$mom_z_est[[i]]       <- ldmom[["z"]]
                              paramdf$mom_z_se[[i]]        <- ldmom[["z_se"]]
                              paramdf$mom_r_est[[i]]       <- ldmom[["r"]]
                              paramdf$mom_r_se[[i]]        <- ldmom[["r_se"]]
                            }, error = function(e) NULL)

                            tryCatch({
                              paramdf$com_time[[i]] <- system.time(
                                ldcom <- ldsep::ldest(ga = foutA$genologlike,
                                                      gb = foutB$genologlike,
                                                      K = paramdf$ploidy[[i]],
                                                      type = "comp",
                                                      model = "flex",
                                                      pen = 1,
                                                      se = FALSE)
                              )[[3]]

                              paramdf$com_D_est[[i]]       <- ldcom[["D"]]
                              paramdf$com_D_se[[i]]        <- ldcom[["D_se"]]
                              paramdf$com_Dprime_est[[i]]  <- ldcom[["Dprime"]]
                              paramdf$com_Dprime_se[[i]]   <- ldcom[["Dprime_se"]]
                              paramdf$com_Dprimeg_est[[i]] <- ldcom[["Dprimeg"]]
                              paramdf$com_Dprimeg_se[[i]]  <- ldcom[["Dprimeg_se"]]
                              paramdf$com_r2_est[[i]]      <- ldcom[["r2"]]
                              paramdf$com_r2_se[[i]]       <- ldcom[["r2_se"]]
                              paramdf$com_z_est[[i]]       <- ldcom[["z"]]
                              paramdf$com_z_se[[i]]        <- ldcom[["z_se"]]
                              paramdf$com_r_est[[i]]       <- ldcom[["r"]]
                              paramdf$com_r_se[[i]]        <- ldcom[["r_se"]]
                            }, error = function(e) NULL)

                            tryCatch({
                              paramdf$comnorm_time[[i]] <- system.time(
                                ldcomnorm <- ldsep::ldest(ga = foutA$genologlike,
                                                      gb = foutB$genologlike,
                                                      K = paramdf$ploidy[[i]],
                                                      type = "comp",
                                                      model = "norm")
                              )[[3]]

                              paramdf$comnorm_D_est[[i]]       <- ldcomnorm[["D"]]
                              paramdf$comnorm_D_se[[i]]        <- ldcomnorm[["D_se"]]
                              paramdf$comnorm_Dprime_est[[i]]  <- ldcomnorm[["Dprime"]]
                              paramdf$comnorm_Dprime_se[[i]]   <- ldcomnorm[["Dprime_se"]]
                              paramdf$comnorm_Dprimeg_est[[i]] <- ldcomnorm[["Dprimeg"]]
                              paramdf$comnorm_Dprimeg_se[[i]]  <- ldcomnorm[["Dprimeg_se"]]
                              paramdf$comnorm_r2_est[[i]]      <- ldcomnorm[["r2"]]
                              paramdf$comnorm_r2_se[[i]]       <- ldcomnorm[["r2_se"]]
                              paramdf$comnorm_z_est[[i]]       <- ldcomnorm[["z"]]
                              paramdf$comnorm_z_se[[i]]        <- ldcomnorm[["z_se"]]
                              paramdf$comnorm_r_est[[i]]       <- ldcomnorm[["r"]]
                              paramdf$comnorm_r_se[[i]]        <- ldcomnorm[["r_se"]]
                            }, error = function(e) NULL)

                            paramdf[i, , drop = FALSE]
                          }

if (nc > 1) {
  parallel::stopCluster(cl)
}

write.csv(simdf, "./output/comp/comp_sims_out.csv", row.names = FALSE)
