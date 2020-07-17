################
## Run simulations comparing MLE to moment-based LD estimators
################

set.seed(1)
library(updog)
library(ldsep)
library(dplyr)
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
rseq <- c(0, 0.5, 0.9) ## r has to be less than 0.57 in the 0.5/0.75 scenario
p0505  <- lapply(c(0, 0.5, 0.9) * sqrt(0.5^4), function(x) c("pA" = 0.5, "pB" = 0.5, "D" = x))
p05075 <- lapply(c(0, 0.5) * sqrt(0.5^2 * 0.25 * 0.75), function(x) c("pA" = 0.5, "pB" = 0.75, "D" = x))
p0909  <- lapply(c(0, 0.5, 0.9) * sqrt(0.9^2 * 0.1^2), function(x) c("pA" = 0.9, "pB" = 0.9, "D" = x))
papbd  <- c(p0505, p05075, p0909)

bias <- 1 # allele bias
od   <- 0.01 # overdispersion
seq  <- 0.01 # sequencing error

itermax <- 200

get_Dprime <- function(pab, pAb, paB, pAB) {
  pA <- pAb + pAB
  pB <- paB + pAB
  D  <- pAB - pA * pB
  Dprime <- ifelse(D < 0,
                   D / pmin(pA * pB, (1 - pA) * (1 - pB)),
                   D / pmin(pA * (1 - pB), (1 - pA) * pB))
  return(Dprime)
}

get_Dprimeg <- function(pab, pAb, paB, pAB, K) {
  qmat <- ldsep::get_prob_array(K = K, prob = c(pab, pAb, paB, pAB), log_p = FALSE)
  ldout <- ldsep::Dprime(qmat = qmat, type = "geno")
  return(ldout[["Dprime"]])
}

paramdf <- expand.grid(seed   = seq_len(itermax),
                       nind   = nind,
                       size   = size,
                       papd   = papbd,
                       bias   = bias,
                       od     = od,
                       seq    = seq,
                       ploidy = ploidy)
paramdf <- cbind(paramdf, do.call(rbind, paramdf$papd))
paramdf %>%
  select(-papd) %>%
  mutate(pAB = D + pA * pB,
         pAb = pA - pAB,
         paB = pB - pAB,
         pab = 1 - paB - pAb - pAB,
         r2  = round(D ^ 2 / (pA * (1 - pA) * pB * (1 - pB)), digits = 2),
         r   = sign(D) * sqrt(r2),
         z   = atanh(r),
         Dprime = get_Dprime(pab = pab, pAb = pAb, paB = paB, pAB = pAB)) %>%
  rowwise() %>%
  mutate(Dprimeg = get_Dprimeg(pab = pab, pAb = pAb, paB = paB, pAB = pAB, K = ploidy)) ->
  paramdf

## confirm parameters
stopifnot(abs(sort(rseq) - sort(unique(paramdf$r))) < 0.0001)
stopifnot(paramdf$pAB > 0, paramdf$pAB < 1,
          paramdf$paB > 0, paramdf$paB < 1,
          paramdf$pAb > 0, paramdf$pAb < 1,
          paramdf$pab > 0, paramdf$pab < 1)

# populate columns to be filled ----------------------------------------------
paramdf %>%
  mutate(mle_D_est            = NA_real_,
         mle_D_se             = NA_real_,
         mle_Dprime_est       = NA_real_,
         mle_Dprime_se        = NA_real_,
         mle_r2_est           = NA_real_,
         mle_r2_se            = NA_real_,
         mle_r_est            = NA_real_,
         mle_r_se             = NA_real_,
         mle_z_est            = NA_real_,
         mle_z_se             = NA_real_,
         mle_pab_est          = NA_real_,
         mle_pAb_est          = NA_real_,
         mle_paB_est          = NA_real_,
         mle_pAB_est          = NA_real_,
         mle_time             = NA_real_,
         gen_D_est            = NA_real_,
         gen_D_se             = NA_real_,
         gen_Dprime_est       = NA_real_,
         gen_Dprime_se        = NA_real_,
         gen_r2_est           = NA_real_,
         gen_r2_se            = NA_real_,
         gen_z_est            = NA_real_,
         gen_z_se             = NA_real_,
         gen_r_est            = NA_real_,
         gen_r_se             = NA_real_,
         gen_pab_est          = NA_real_,
         gen_pAb_est          = NA_real_,
         gen_paB_est          = NA_real_,
         gen_pAB_est          = NA_real_,
         gen_time             = NA_real_,
         mom_D_est            = NA_real_,
         mom_D_se             = NA_real_,
         mom_Dprime_est       = NA_real_,
         mom_Dprime_se        = NA_real_,
         mom_Dprimeg_est      = NA_real_,
         mom_Dprimeg_se       = NA_real_,
         mom_r2_est           = NA_real_,
         mom_r2_se            = NA_real_,
         mom_z_est            = NA_real_,
         mom_z_se             = NA_real_,
         mom_r_est            = NA_real_,
         mom_r_se             = NA_real_,
         mom_time             = NA_real_,
         com_D_est            = NA_real_,
         com_D_se             = NA_real_,
         com_Dprime_est       = NA_real_,
         com_Dprime_se        = NA_real_,
         com_Dprimeg_est      = NA_real_,
         com_Dprimeg_se       = NA_real_,
         com_r2_est           = NA_real_,
         com_r2_se            = NA_real_,
         com_z_est            = NA_real_,
         com_z_se             = NA_real_,
         com_r_est            = NA_real_,
         com_r_se             = NA_real_,
         com_time             = NA_real_,
         comnorm_D_est        = NA_real_,
         comnorm_D_se         = NA_real_,
         comnorm_Dprime_est   = NA_real_,
         comnorm_Dprime_se    = NA_real_,
         comnorm_Dprimeg_est  = NA_real_,
         comnorm_Dprimeg_se   = NA_real_,
         comnorm_r2_est       = NA_real_,
         comnorm_r2_se        = NA_real_,
         comnorm_z_est        = NA_real_,
         comnorm_z_se         = NA_real_,
         comnorm_r_est        = NA_real_,
         comnorm_r_se         = NA_real_,
         comnorm_time         = NA_real_) ->
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
                                      "flexdog")) %dopar% {
                            set.seed(paramdf$seed[[i]])

                            # generate data -----------------------------------
                            prob <- c(paramdf$pab[[i]],
                                      paramdf$pAb[[i]],
                                      paramdf$paB[[i]],
                                      paramdf$pAB[[i]])
                            hapmat <- stats::rmultinom(n = paramdf$nind[[i]],
                                                       size = paramdf$ploidy[[i]],
                                                       prob = prob)
                            gA <- hapmat[2, ] + hapmat[4, ]
                            gB <- hapmat[3, ] + hapmat[4, ]
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

                              paramdf$mom_D_est[[i]]        <- ldmom[["D"]]
                              paramdf$mom_D_se[[i]]         <- ldmom[["D_se"]]
                              paramdf$mom_Dprime_est[[i]]   <- ldmom[["Dprime"]]
                              paramdf$mom_Dprime_se[[i]]    <- ldmom[["Dprime_se"]]
                              paramdf$mom_Dprimeg_est[[i]]  <- ldmom[["Dprimeg"]]
                              paramdf$mom_Dprimeg_se[[i]]   <- ldmom[["Dprimeg_se"]]
                              paramdf$mom_r2_est[[i]]       <- ldmom[["r2"]]
                              paramdf$mom_r2_se[[i]]        <- ldmom[["r2_se"]]
                              paramdf$mom_z_est[[i]]        <- ldmom[["z"]]
                              paramdf$mom_z_se[[i]]         <- ldmom[["z_se"]]
                              paramdf$mom_r_est[[i]]        <- ldmom[["r"]]
                              paramdf$mom_r_se[[i]]         <- ldmom[["r_se"]]
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

                              paramdf$com_D_est[[i]]        <- ldcom[["D"]]
                              paramdf$com_D_se[[i]]         <- ldcom[["D_se"]]
                              paramdf$com_Dprime_est[[i]]   <- ldcom[["Dprime"]]
                              paramdf$com_Dprime_se[[i]]    <- ldcom[["Dprime_se"]]
                              paramdf$com_Dprimeg_est[[i]]  <- ldcom[["Dprimeg"]]
                              paramdf$com_Dprimeg_se[[i]]   <- ldcom[["Dprimeg_se"]]
                              paramdf$com_r2_est[[i]]       <- ldcom[["r2"]]
                              paramdf$com_r2_se[[i]]        <- ldcom[["r2_se"]]
                              paramdf$com_z_est[[i]]        <- ldcom[["z"]]
                              paramdf$com_z_se[[i]]         <- ldcom[["z_se"]]
                              paramdf$com_r_est[[i]]        <- ldcom[["r"]]
                              paramdf$com_r_se[[i]]         <- ldcom[["r_se"]]
                            }, error = function(e) NULL)

                            tryCatch({
                              paramdf$comnorm_time[[i]] <- system.time(
                                ldcomnorm <- ldsep::ldest(ga = foutA$genologlike,
                                                      gb = foutB$genologlike,
                                                      K = paramdf$ploidy[[i]],
                                                      type = "comp",
                                                      model = "norm")
                              )[[3]]

                              paramdf$comnorm_D_est[[i]]        <- ldcomnorm[["D"]]
                              paramdf$comnorm_D_se[[i]]         <- ldcomnorm[["D_se"]]
                              paramdf$comnorm_Dprime_est[[i]]   <- ldcomnorm[["Dprime"]]
                              paramdf$comnorm_Dprime_se[[i]]    <- ldcomnorm[["Dprime_se"]]
                              paramdf$comnorm_Dprimeg_est[[i]]  <- ldcomnorm[["Dprimeg"]]
                              paramdf$comnorm_Dprimeg_se[[i]]   <- ldcomnorm[["Dprimeg_se"]]
                              paramdf$comnorm_r2_est[[i]]       <- ldcomnorm[["r2"]]
                              paramdf$comnorm_r2_se[[i]]        <- ldcomnorm[["r2_se"]]
                              paramdf$comnorm_z_est[[i]]        <- ldcomnorm[["z"]]
                              paramdf$comnorm_z_se[[i]]         <- ldcomnorm[["z_se"]]
                              paramdf$comnorm_r_est[[i]]        <- ldcomnorm[["r"]]
                              paramdf$comnorm_r_se[[i]]         <- ldcomnorm[["r_se"]]
                            }, error = function(e) NULL)

                            paramdf[i, , drop = FALSE]
                          }

if (nc > 1) {
  parallel::stopCluster(cl)
}

write.csv(simdf, "./output/mle/mle_sims_out.csv", row.names = FALSE)
