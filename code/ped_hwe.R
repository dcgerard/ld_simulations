##################
## Explore effects of random mating
##################

library(updog)
library(xtable)

#' Update dosage probabilities given random mating
#'
#' Assuming random mating and autopolyploidy, and given current genotype
#' probabilities, this function will return the genotype probabilities in
#' the next generation.
#'
#' @param prob A vector of probabilities. \code{prob[i]} is the probability
#'     of dosage \code{i-1}.
#'
#' @return A vector of probabilities. Element \code{i} is the probability
#'     in the next generation of dosage \code{i-1}
update_prob <- function(prob) {
  stopifnot(abs(sum(prob) - 1) < sqrt(.Machine$double.eps))
  ploidy <- length(prob) - 1
  qarray <- updog::get_q_array(ploidy = ploidy)
  pdosage_mat <- outer(X = prob, Y = prob, FUN = `*`)
  apply(qarray, 3, function(x) sum(x * pdosage_mat))
}

ngen <- 10
pstart <- c(0, 0, 1, 0, 0)
ploidy <- length(pstart) - 1
alpha <- sum(0:ploidy * pstart) / ploidy
hwe_prob <- dbinom(x = 0:ploidy, size = ploidy, prob = alpha)

probmat <- matrix(NA_real_, nrow = ngen + 1, ncol = ploidy + 1)
probmat[1, ] <- pstart
probmat[ngen+1, ] <- hwe_prob

pnew <- pstart
for (index in 2:ngen) {
  pnew <- update_prob(pnew)
  probmat[index, ] <- pnew
}
colnames(probmat) <- 0:ploidy
rownames(probmat) <- c(paste0("Generation ", 1:ngen), "HWE")

print(xtable(probmat, digits = 4), file = "./output/ped/hwe_prob.txt")
