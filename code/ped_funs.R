###################
## Functions for generating data from PedigreeSim
## only for tetraploids
###################

## Generate sim.par

gen_par <- function(ploidy, seed, file) {

  partext <- paste0("PLOIDY = ", ploidy, "\n",
                    "MAPFUNCTION = HALDANE\n",
                    "MISSING = NA\n",
                    "PEDFILE = sim.ped\n",
                    "CHROMFILE = sim.chrom\n",
                    "MAPFILE = sim.map\n",
                    "FOUNDERFILE = sim.gen\n",
                    "OUTPUT = sim_out\n",
                    "NATURALPAIRING = 0\n",
                    "SEED = ", seed, "\n")
  cat(file = file, partext)
}

gen_par(ploidy = 4, seed = 1, file = "./PedigreeSim/sim.par")

# Function to generate sim.ped (pedigree file)
gen_ped <- function(nind, ngen, file) {
  peddf <- data.frame(Name = paste0("P", seq_len(nind)),
                      Parent1 = NA,
                      Parent2 = NA)
  for (index in seq_len(ngen - 1)) {
    if (index == 1) {
      tempdf <- data.frame(Name = paste0("F", index, "_", seq_len(nind)),
                           Parent1 = sample(paste0("P", seq_len(nind)), size = nind, replace = TRUE),
                           Parent2 = sample(paste0("P", seq_len(nind)), size = nind, replace = TRUE))
    } else {
      tempdf <- data.frame(Name = paste0("F", index, "_", seq_len(nind)),
                           Parent1 = sample(paste0("F", index - 1, "_", seq_len(nind)), size = nind, replace = TRUE),
                           Parent2 = sample(paste0("F", index - 1, "_", seq_len(nind)), size = nind, replace = TRUE))
    }
    peddf <- rbind(peddf, tempdf)
    # peddf$Name[(nind * index + 1):(nind * (index + 1))] <- tempdf$Name
    # peddf$Parent1[(nind * index + 1):(nind * (index + 1))] <- tempdf$Parent1
    # peddf$Parent2[(nind * index + 1):(nind * (index + 1))] <- tempdf$Parent2
  }
  write.table(x = peddf, file = file, row.names = FALSE, quote = FALSE)
}

# gen_ped(nind = nind, ngen = ngen, file = "./PedigreeSim/sim.ped")

## Generate sim.chrom
gen_chrom <- function(length, centromere, prefPairing, quadrivalents, file) {
  chromedf <- data.frame(chromosome = "A",
                         length = length,
                         centromere = centromere,
                         prefPairing = prefPairing,
                         quadrivalents = quadrivalents)
  write.table(x = chromedf, file = file, row.names = FALSE, quote = FALSE)
}

# gen_chrom(length = chrome_size,
#           centromere = round(chrome_size / 2),
#           prefPairing = 1,
#           quadrivalents = 0,
#           file = "./PedigreeSim/sim.chrom")

## Generate sim.map
gen_map <- function(snp_loc, file) {
  df <- data.frame(marker = paste0("A", snp_loc),
                   chromosome = "A",
                   position = snp_loc)
  write.table(df, file = file, row.names = FALSE, quote = FALSE)
}

# gen_map(snp_loc = c(50, 51, 60, 70, 80, 90, 99, 100), file = "./PedigreeSim/sim.map")


## Generate sim.gen
# 1 and 3 are homologues, 2 and 4 are homologues

#' generate parental dosages
#' Parents start in perfect LD.
#'
#' @param alpha1 The allele frequency of subgenome 1. Subgenome 2 has allele
#'     frequency 1 - alpha1. So the overall allele frequency is 0.5.
#'     This allows D to start at 1.
gen_sim <- function(nind, snp_loc, alpha1, file, ploidy = 4) {
  stopifnot(alpha1 >= 0, alpha1 <= 1)
  alpha1 <- max(alpha1, 1 - alpha1)

  nsnp <- length(snp_loc)

  n_flip <- 2 * nind * alpha1 - nind
  n_noflip <- nind - n_flip

  genmat_flip <- matrix(rep(rep(rep(c("A", "B"), ploidy/2), n_flip), nsnp), byrow = TRUE, nrow = nsnp, ncol = n_flip * ploidy)
  genmat_noflip <- matrix(rep(rep(rep(c("A", "B"), each = ploidy/2), n_noflip), nsnp), byrow = TRUE, nrow = nsnp, ncol = n_noflip * ploidy)

  genmat <- cbind(paste0("A", snp_loc), genmat_flip, genmat_noflip)
  colnames(genmat) <- c("marker", paste0(rep(paste0("P", seq_len(nind)), each = ploidy), "_", seq_len(ploidy)))

  write.table(x = genmat,
              file = file,
              row.names = FALSE,
              quote = FALSE)
}

# gen_sim(nind = nind, snp_loc = c(50, 51, 60, 70, 80, 90, 99, 100), alpha1 = alpha1, file = "./PedigreeSim/sim.gen")
