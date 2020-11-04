##################
## Plots demonstrating Deviation from HWE
## Also generate true correlation matrix
##################

## Generate under HWE ---------------------------------------------------------
library(corrplot)
library(xtable)
source("./code/ped_funs.R")
ploidy <- 4
seed <- 1
nind <- 10000
chrome_size <- 100
snp_to_keep <- c(50, 51, 60, 70, 80, 90, 99, 100)
ngen <- 10
pref_pair <- 1/3
quadprop <- 0
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

r2_true_hwe <- cor(t(genomat))^2

pdf(file = "./output/ped/hwe_heatmap.pdf", width = 3, height = 3, family = "Times")
corrplot(corr = r2_true_hwe[50:100, 50:100],
         method = "color",
         type = "upper",
         diag = FALSE,
         tl.pos = "n")
dev.off()

geno50 <- genomat[50, ]
geno100 <- genomat[100, ]


hwe_mat <- matrix(NA_real_, ncol = ploidy + 1, nrow = 3)
hwe_mat[1, ] <- dbinom(0:ploidy, ploidy, 0.5)
hwe_mat[2, ] <- prop.table(table(geno50))
hwe_mat[3, ] <- prop.table(table(geno100))
hwe_mat
rownames(hwe_mat) <- c("HWE", "50 cM", "100 cM")
colnames(hwe_mat) <- 0:ploidy

print(xtable(hwe_mat, digits = 4), file = "./output/ped/ped_hwe.txt")


## Generate under HWD pp -------------------------------------------------------
pref_pair <- 0
quadprop <- 2/3

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

r2_true_dr <- cor(t(genomat))^2

pdf(file = "./output/ped/hwe_violated_heatmap_dr.pdf", width = 3, height = 3, family = "Times")
corrplot(corr = r2_true_dr[50:100, 50:100],
         method = "color",
         type = "upper",
         diag = FALSE,
         tl.pos = "n")
dev.off()

geno50 <- genomat[50, ]
geno100 <- genomat[100, ]


hwe_mat <- matrix(NA_real_, ncol = ploidy + 1, nrow = 3)
hwe_mat[1, ] <- dbinom(0:ploidy, ploidy, 0.5)
hwe_mat[2, ] <- prop.table(table(geno50))
hwe_mat[3, ] <- prop.table(table(geno100))
hwe_mat
rownames(hwe_mat) <- c("HWE", "50 cM", "100 cM")
colnames(hwe_mat) <- 0:ploidy

print(xtable(hwe_mat, digits = 4), file = "./output/ped/ped_hwe_violated_dr.txt")


## Generate under HWD double reduction -----------------------------------------
pref_pair <- 1
quadprop <- 0

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

r2_true_pp <- cor(t(genomat))^2

pdf(file = "./output/ped/hwe_violated_heatmap_pp.pdf", width = 3, height = 3, family = "Times")
corrplot(corr = r2_true_pp[50:100, 50:100],
         method = "color",
         type = "upper",
         diag = FALSE,
         tl.pos = "n")
dev.off()

geno50 <- genomat[50, ]
geno100 <- genomat[100, ]

hwe_mat <- matrix(NA_real_, ncol = ploidy + 1, nrow = 4)
hwe_mat[1, ] <- dbinom(0:ploidy, ploidy, 0.5)

ygen1 <- dbinom(x = 0:(ploidy/2), size = ploidy/2, prob = 0.1)
ygen2 <- dbinom(x = 0:(ploidy/2), size = ploidy/2, prob = 0.9)
hwe_mat[2, ] <- convolve(x = ygen1, y = rev(ygen2), type = "open")

hwe_mat[3, ] <- prop.table(table(geno50))
hwe_mat[4, ] <- prop.table(table(geno100))
hwe_mat
rownames(hwe_mat) <- c("HWE Auto", "HWE Allo", "50 cM", "100 cM")
colnames(hwe_mat) <- 0:ploidy

print(xtable(hwe_mat, digits = 4), file = "./output/ped/ped_hwe_violated_pp.txt")
