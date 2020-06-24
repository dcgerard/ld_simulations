##########################
## Demonstrate Flexibility of proportional bivariate normal
##########################

library(ldsep)
library(tidyverse)
library(ggthemes)

#' KL-divergence between two matrices of probabilities.
#'
#' @param A A matrix of probabilities. A[i,j] is the probability of genotype
#'     i-1 in locus 1 and j-1 in locus 2. Should be the same dimension as B
#' @param B A matrix of probabilities. B[i,j] is the probability of genotype
#'     i-1 in locus 1 and j-1 in locus 2. Should be the same dimension as A
#'
#' @author David Gerard
kl_div <- function(A, B) {
  TOL <- sqrt(.Machine$double.eps)
  sum((A * log(A/B))[A > TOL])
}

#' Objective to minimize to find closest proportional biviariate normal to A
#'
#' @param A A matrix of probabilities. Col and row dimensions should be the same.
#' @param par A vector of length 5. par[1:2] is the mean vector. par[3:5]
#'     contains L[1, 1], L[2, 1], and L[2, 2] where sigma = L %*% t(L)
#'
#' @return The KL-divergence between the proportional normal distribution
#'     defined by par and A
#'
#' @example
#' K <- 6
#' A <- matrix(runif((K + 1)^2), nrow = K+1)
#' A <- A / sum(A)
#' par <- c(3, 4, 4, 1, 3)
#' obj(par, A)
#'
#' @author David Gerard
obj <- function(par, A) {
  stopifnot(nrow(A) == ncol(A))
  stopifnot(length(par) == 5)
  K <- nrow(A) - 1
  L <- matrix(data = c(par[[3]], par[[4]], 0, par[[5]]),
              byrow = FALSE,
              ncol = 2,
              nrow = 2)
  sigma <- L %*% t(L)
  mu <- par[1:2]
  distB <- pbnorm_dist(mu = mu, sigma = sigma, K = K, log = FALSE)
  kl_div(A = A, B = distB)
}

#' Generate distribution under HWE
#'
#' @param pA The reference allele frequency at locus 1.
#' @param pB The reference allele frequency at locus 2
#' @param D The LD coefficient
#' @param K The ploidy of the species.
#'
#' @author David Gerard
hwe_dist <- function(pA, pB, D, K) {
  stopifnot(D >= -pA * pB,
            D >= -(1 - pA) * (1 - pB),
            D <= pA * (1 - pB),
            D <= (1 - pA) * pB)
  pAB <- pA * pB + D
  pAb <- pA * (1 - pB) - D
  paB <- (1 - pA) * pB - D
  pab <- (1 - pA) * (1 - pB) + D

  distmat <- exp(ldsep:::get_prob_array(K = K, prob = c(pab, pAb, paB, pAB)))
  colnames(distmat) <- 0:K
  rownames(distmat) <- 0:K
  return(distmat)
}

# These are the same settings as in the simulations under HWE ----
K <- 6
simdf <- data.frame(pA = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9),
                    pB = c(0.5, 0.5, 0.5, 0.75, 0.75, 0.9, 0.9, 0.9),
                    r = c(0, 0.5, 0.9, 0, 0.5, 0, 0.5, 0.9))
simdf$D <- simdf$r * sqrt(simdf$pA * (1 - simdf$pA) * simdf$pB * (1 - simdf$pB))

# Find pbnorm dist closest to each hwe dist ----
for (i in seq_len(nrow(simdf))) {
  A <- hwe_dist(pA = simdf$pA[[i]], pB = simdf$pB[[i]], D = simdf$D[[i]], K = K)
  pg1 <- rowSums(A)
  pg2 <- colSums(A)
  mu_init <- c(sum(0:K * pg1),
               sum(0:K * pg2))
  sigma_init <- matrix(NA_real_, nrow = 2, ncol = 2)
  sigma_init[1, 1] <- sum((0:K - mu_init[[1]])^2 * pg1)
  sigma_init[2, 2] <- sum((0:K - mu_init[[2]])^2 * pg2)
  sigma_init[1, 2] <- sum(tcrossprod(0:K) * A) - prod(mu_init)
  sigma_init[2, 1] <- sigma_init[1, 2]
  L_init <- t(chol(sigma_init))
  par_init <- c(mu_init, L_init[lower.tri(L_init, diag = TRUE)])
  oout <- stats::optim(par = par_init,
                       fn = obj,
                       method = "L-BFGS-B",
                       lower = c(-Inf, -Inf, 0, -Inf, 0),
                       upper = rep(Inf, 5),
                       A = A)

  mu_final <- oout$par[1:2]
  L_final <- matrix(0, nrow = 2, ncol = 2)
  L_final[1, 1] <- oout$par[[3]]
  L_final[2, 1] <- oout$par[[4]]
  L_final[2, 2] <- oout$par[[5]]
  sigma_final <- L_final %*% t(L_final)

  B <- pbnorm_dist(mu = mu_final, sigma = sigma_final, K = K, log = FALSE)
  colnames(B) <- 0:K
  rownames(B) <- 0:K

  as_tibble(A, rownames = "g1") %>%
    gather(-g1, key = "g2", value = "probA") ->
    Adf

  as_tibble(B, rownames = "g1") %>%
    gather(-g1, key = "g2", value = "probB") ->
    Bdf

  probtemp <- left_join(Adf, Bdf, by = c("g1", "g2"))
  probtemp$scenario <- paste0("list(p[A]==",
                              simdf$pA[[i]],
                              ", p[B]==",
                              simdf$pB[[i]],
                              ", r==",
                              simdf$r[[i]],
                              ")")
  if (i == 1) {
    probdf <- probtemp
  } else {
    probdf <- bind_rows(probdf, probtemp)
  }
}

probdf %>%
  rename(HWE = "probA", Normal = "probB") %>%
  gather(-g1, -g2, -scenario, key = "Distribution", value = "Probability") ->
  probdflong
ggplot(probdflong, aes(x = g1, y = g2, size = Probability, pch = Distribution, color = Distribution)) +
  geom_point() +
  xlab("Dosage at Locus 1") +
  ylab("Dosage at Locus 2") +
  theme_bw() +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_colorblind() +
  facet_wrap(~scenario,
             labeller = label_parsed)  +
  theme(strip.background = element_rect(fill = "white"))->
  pl

ggsave(filename = "./output/compare_norm/normdist.pdf",
       plot = pl,
       height = 5,
       width = 6,
       family = "Times")

# Look at random pbnorm distributions ------
for (j in 1:9) {
  set.seed(j)
  mu <- runif(n = 2, min = 0, max = K)
  sigma <- rWishart(n = 1, df = 2, Sigma = diag(c(1, 1)))
  B <- pbnorm_dist(mu = mu, sigma = sigma[,,1], K = K, log = FALSE)
  colnames(B) <- 0:K
  rownames(B) <- 0:K

  as_tibble(B, rownames = "g1") %>%
    gather(-g1, key = "g2", value = "Probability") ->
    tempB
  tempB$seed <- j

  if (j == 1) {
    Bdf <- tempB
  } else {
    Bdf <- bind_rows(Bdf, tempB)
  }
}

ggplot(Bdf, aes(x = g1, y = g2, size = Probability)) +
  geom_point() +
  xlab("Dosage at Locus 1") +
  ylab("Dosage at Locus 2") +
  theme_bw() +
  facet_wrap(~seed)  +
  theme(strip.background = element_rect(fill = "white"))->
  pl

ggsave(filename = "./output/compare_norm/randnorm.pdf",
       plot = pl,
       height = 5,
       width = 6,
       family = "Times")
