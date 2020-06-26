########################
## Plot LD Estimates in McAllister Data
########################

library(ldsep)
library(tidyverse)
locdf <- read_csv("./output/mca/locdf.csv")

## hexaploid ----------
ldest_hap_geno <- readRDS(file = "./output/mca/ldest_hap_geno_hex.RDS")
pdf(file = "./output/mca/mca_plots/heat_mca_hap_geno.pdf",
    width = 4,
    height = 3,
    family = "Times")
plot(ldest_hap_geno, element = "r2", is.corr = TRUE, title = "")
dev.off()

ldest_hap_genolike <- readRDS(file = "./output/mca/ldest_hap_genolike_hex.RDS")
pdf(file = "./output/mca/mca_plots/heat_mca_hap_genolike.pdf",
    width = 4,
    height = 3,
    family = "Times")
plot(ldest_hap_genolike, element = "r2", is.corr = TRUE, title = "")
dev.off()

ldest_comp_geno <- readRDS(file = "./output/mca/ldest_comp_geno_hex.RDS")
pdf(file = "./output/mca/mca_plots/heat_mca_comp_geno.pdf",
    width = 4,
    height = 3,
    family = "Times")
plot(ldest_comp_geno, element = "r2", is.corr = TRUE, title = "")
dev.off()

ldest_comp_genolike <- readRDS(file = "./output/mca/ldest_comp_genolike_hex.RDS")
pdf(file = "./output/mca/mca_plots/heat_mca_comp_genolike.pdf",
    width = 4,
    height = 3,
    family = "Times")
plot(ldest_comp_genolike, element = "r2", title = "")
dev.off()

ldest_comp_genolike_flex <- readRDS(file = "./output/mca/ldest_comp_genolike_flex_hex.RDS")
ldest_comp_genolike_flex$r2[ldest_comp_genolike_flex$r2 > 1] <- NA
pdf(file = "./output/mca/mca_plots/heat_mca_comp_genolike_flex.pdf",
    width = 4,
    height = 3,
    family = "Times")
plot(ldest_comp_genolike_flex, element = "r2", title = "")
dev.off()

## shrinkage ------
rout0 <- ldshrink(obj = ldest_comp_genolike)
pdf(file = "./output/mca/mca_plots/shrink_mca_comp_genolike.pdf",
    width = 4,
    height = 3,
    family = "Times")
corrplot::corrplot(corr = rout0^2,
                   method = "color",
                   diag = FALSE,
                   type = "upper",
                   tl.pos = "n")
dev.off()

rout1 <- ldshrink(obj = ldest_hap_genolike)
pdf(file = "./output/mca/mca_plots/shrink_mca_hap_genolike.pdf",
    width = 4,
    height = 3,
    family = "Times")
corrplot::corrplot(corr = rout1^2,
                   method = "color",
                   diag = FALSE,
                   type = "upper",
                   tl.pos = "n")
dev.off()

rout2 <- ldshrink(obj = ldest_comp_geno, optmethod = "mixEM")
pdf(file = "./output/mca/mca_plots/shrink_mca_comp_geno.pdf",
    width = 4,
    height = 3,
    family = "Times")
corrplot::corrplot(corr = rout2^2,
                   method = "color",
                   diag = FALSE,
                   type = "upper",
                   tl.pos = "n")
dev.off()

rout3 <- ldshrink(obj = ldest_hap_geno, optmethod = "mixEM")
pdf(file = "./output/mca/mca_plots/shrink_mca_hap_geno.pdf",
    width = 4,
    height = 3,
    family = "Times")
corrplot::corrplot(corr = rout3^2,
                   method = "color",
                   diag = FALSE,
                   type = "upper",
                   tl.pos = "n")
dev.off()

## nonoploid -------------
ldest_hap_geno <- readRDS(file = "./output/mca/ldest_hap_geno_non.RDS")
plot(ldest_hap_geno, element = "r2")

ldest_hap_genolike <- readRDS(file = "./output/mca/ldest_hap_genolike_non.RDS")
plot(ldest_hap_genolike, element = "r2")

ldest_comp_geno <- readRDS(file = "./output/mca/ldest_comp_geno_non.RDS")
plot(ldest_comp_geno, element = "r2")

ldest_comp_genolike <- readRDS(file = "./output/mca/ldest_comp_genolike_non.RDS")
plot(ldest_comp_genolike, element = "r2")

ldest_comp_genolike_flex <- readRDS(file = "./output/mca/ldest_comp_genolike_flex_non.RDS")
ldest_comp_genolike_flex$r2[ldest_comp_genolike_flex$r2 > 1] <- NA
plot(ldest_comp_genolike_flex, element = "r2")


