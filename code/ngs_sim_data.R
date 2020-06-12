##########################
## Simulations comparing ngsLD to ldsep
##########################

library(updog)
library(reshape2)
library(ldsep)
library(readr)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nc <- 1
} else {
  eval(parse(text = args[[1]]))
}

TOL <- sqrt(.Machine$double.eps)
rhap_given_loc <- function(locvec, freqvec) {
  stopifnot(length(locvec) == length(freqvec))
  hapvec <- rep(NA_real_, length = length(locvec))
  hapvec[[1]] <- sample(x = c(0, 1),
                        size = 1,
                        prob = c(1 - freqvec[[1]], freqvec[[1]]))
  for (i in 2:length(hapvec)) {
    corval <- exp(-abs(locvec[[i]] - locvec[[i-1]]))

    pA <- freqvec[[i-1]]
    pB <- freqvec[[i]]

    dval <- corval * sqrt(pA * (1 - pA) * pB * (1 - pB))

    dval <- min(c(dval, pA * (1 - pB), pB * (1 - pA)))

    pAB = dval + pA * pB
    pAb = pA - pAB
    paB = pB - pAB
    pab = 1 - paB - pAb - pAB

    stopifnot(pAB > -TOL, pAb > -TOL, paB > -TOL, pab > -TOL,
              pAB < 1 + TOL, pAb < 1 + TOL, paB < 1 + TOL, pab < 1 + TOL)

    pAB[pAB < 0] <- 0
    pAB[pAB > 1] <- 1

    paB[paB < 0] <- 0
    paB[paB > 1] <- 1

    pAb[pAb < 0] <- 0
    pAb[pAb > 1] <- 1

    pab[pab < 0] <- 0
    pab[pab > 1] <- 1

    if (hapvec[[i-1]] == 1) {
      pB_gA <- pAB / (pAb + pAB)
    } else {
      pB_gA <- paB / (pab + paB)
    }
    hapvec[[i]] <- sample(x = c(0, 1),
                          size = 1,
                          prob = c(1 - pB_gA, pB_gA))
  }
  return(hapvec)
}

gensim <- function(nind, locvec, freqvec) {
  nsnp <- length(locvec)
  genomat <- matrix(NA, nrow = nsnp, ncol = nind)

  for (i in seq_len(nind)) {
    hvec1 <- rhap_given_loc(locvec = locvec, freqvec = freqvec)
    hvec2 <- rhap_given_loc(locvec = locvec, freqvec = freqvec)
    genomat[, i] <- hvec1 + hvec2
  }

  return(genomat)
}

library(updog)
set.seed(1)
nind <- 100
nsnp <- 20
ploidy <- 2
bias <- 1
od <- 0.01
seq <- 0.01
size <- 10
freqvec <- stats::runif(n = nsnp, min = 0.1, max = 0.9)
locvec <- seq(0, 5, length = nsnp)
genomat <- gensim(nind = nind, locvec = locvec, freqvec = freqvec)

sizemat <- matrix(data = size, nrow = nsnp, ncol = nind)
refmat <- matrix(NA_real_, nrow = nsnp, ncol = nind)
for (i in seq_len(nsnp)) {
  refmat[i, ] <- updog::rflexdog(sizevec = sizemat[i, ],
                                 geno = genomat[i, ],
                                 ploidy = 2,
                                 seq = seq,
                                 bias = bias,
                                 od = od)
}

colnames(refmat) <- seq_len(nind)
colnames(sizemat) <- seq_len(nind)
rownames(refmat) <- seq_len(nsnp)
rownames(sizemat) <- seq_len(nsnp)

mout <- multidog(refmat = refmat,
                 sizemat = sizemat,
                 ploidy = 2,
                 model = "hw",
                 nc = nc,
                 bias_init = bias,
                 update_bias = FALSE,
                 od = od,
                 update_od = FALSE,
                 seq = seq,
                 update_seq = FALSE)

loglarray <- format_multidog(x = mout, varname = c("logL_0", "logL_1", "logL_2"))

saveRDS(object = loglarray, file = "./output/ngs_out/updog_format.RDS")

## format data for ngsLD
llmat <- dcast(melt(data = loglarray), formula = snp ~ ind + variable, value.var = "value")
llmat <- llmat[, -1]

write.table(x = llmat, file = "./output/ngs_out/llike.tsv", sep = '\t', row.names = FALSE)
gzip("./output/ngs_out/llike.tsv", destname = "./output/ngs_out/llike.tsv.gz")

posdf <- data.frame(chrome = 1, pos = seq_len(nsnp))
write.table(x = posdf, file = "./output/ngs_out/pos.tsv", sep = '\t', row.names = FALSE, col.names = FALSE)
gzip("./output/ngs_out/pos.tsv", destname = "./output/ngs_out/pos.tsv.gz")
