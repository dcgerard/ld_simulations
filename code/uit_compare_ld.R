library(ldsep)
library(tidyverse)
library(GGally)

uit_info <- read_csv("./output/uit/uit_suc.csv")

ld_hap_geno <- readRDS("./output/uit/ldest_hap_geno.RDS")
ld_comp_geno <- readRDS("./output/uit/ldest_comp_geno.RDS")
ld_hap_genolike <- readRDS("./output/uit/ldest_hap_genolike.RDS")
ld_comp_genolike <- readRDS("./output/uit/ldest_comp_genolike.RDS")
ld_comp_genolike_flex <- readRDS("./output/uit/ldest_comp_genolike_flex.RDS")

pdf(file = "./output/uit/uit_fig/heat_hap_geno.pdf",
    height = 6,
    width = 6,
    family = "Times")
plot(ld_hap_geno,
     element = "r2",
     title = "")
dev.off()

pdf(file = "./output/uit/uit_fig/heat_hap_genolike.pdf",
    height = 6,
    width = 6,
    family = "Times")
plot(ld_hap_genolike,
     element = "r2",
     title = "")
dev.off()

pdf(file = "./output/uit/uit_fig/heat_comp_geno.pdf",
    height = 6,
    width = 6,
    family = "Times")
plot(ld_comp_geno,
     element = "r2",
     title = "")
dev.off()

pdf(file = "./output/uit/uit_fig/heat_comp_genolike_flex.pdf",
    height = 6,
    width = 6,
    family = "Times")
plot(ld_comp_genolike_flex,
     element = "r2",
     title = "")
dev.off()

pdf(file = "./output/uit/uit_fig/heat_comp_genolike.pdf",
    height = 6,
    width = 6,
    family = "Times")
plot(ld_comp_genolike,
     element = "r2",
     title = "")
dev.off()

ld_df <- tibble(`Haplotypic, Posterior\nMode Genotypes` = ld_hap_geno$r2,
                `Haplotypic, Genotype\nLikelihoods` = ld_hap_genolike$r2,
                `Composite, Posterior\nMean Genotypes` = ld_comp_geno$r2,
                `Composite,\nGeneral Categorical` = ld_comp_genolike_flex$r2,
                `Composite,\nProportional Normal` = ld_comp_genolike$r2)

my_hist <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_histogram(..., fill = "white", color = "black", bins = 15)
}

my_point <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(..., alpha = 1/4) +
    geom_abline(intercept = 0, slope = 1, lty = 2, col = 2)
}

ggpairs(ld_df,
        lower = list(continuous = my_point),
        diag = list(continuous = my_hist),
        upper = "blank") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) ->
  pl

ggsave(filename = "./output/uit/uit_fig/uit_pairs.pdf",
       plot = pl,
       height = 6.6,
       width = 6.6,
       family = "Times")

## Look at difference

ldhapgeno <- format_lddf(obj = ld_hap_geno, element = "r2")
stopifnot(uit_info$`Variant name`[uit_info$`Variant name` %in% rownames(ldhapgeno)] == rownames(ldhapgeno))

uitsub <- uit_info[match(rownames(ldhapgeno), uit_info$`Variant name`), ]
stopifnot(uitsub$`Variant name` == rownames(ldhapgeno))

which_contig <- unique(uitsub$Contig)
stopifnot(length(which_contig) == 2)
is1 <- uitsub$Contig == which_contig[[1]]
is2 <- uitsub$Contig == which_contig[[2]]

ldmat11_hap_geno <- ldhapgeno[is1, is1]
ldmat12_hap_geno <- ldhapgeno[is1, is2]
ldmat22_hap_geno <- ldhapgeno[is2, is2]

ldmat_hap_genolike <- format_lddf(obj = ld_hap_genolike, element = "r2")
ldmat11_hap_genolike <- ldmat_hap_genolike[is1, is1]
ldmat12_hap_genolike <- ldmat_hap_genolike[is1, is2]
ldmat22_hap_genolike <- ldmat_hap_genolike[is2, is2]

ldmat_comp_geno <- format_lddf(obj = ld_comp_geno, element = "r2")
ldmat11_comp_geno <- ldmat_comp_geno[is1, is1]
ldmat12_comp_geno <- ldmat_comp_geno[is1, is2]
ldmat22_comp_geno <- ldmat_comp_geno[is2, is2]

ldmat_comp_genolike_flex <- format_lddf(obj = ld_comp_genolike_flex, element = "r2")
ldmat11_comp_genolike_flex <- ldmat_comp_genolike_flex[is1, is1]
ldmat12_comp_genolike_flex <- ldmat_comp_genolike_flex[is1, is2]
ldmat22_comp_genolike_flex <- ldmat_comp_genolike_flex[is2, is2]

ldmat_comp_genolike <- format_lddf(obj = ld_comp_genolike, element = "r2")
ldmat11_comp_genolike  <- ldmat_comp_genolike[is1, is1]
ldmat12_comp_genolike  <- ldmat_comp_genolike[is1, is2]
ldmat22_comp_genolike  <- ldmat_comp_genolike[is2, is2]

df11 <- data.frame(hap_geno = ldmat11_hap_geno[upper.tri(ldmat11_hap_geno)],
                   hap_genolike = ldmat11_hap_genolike[upper.tri(ldmat11_hap_genolike)],
                   comp_geno = ldmat11_comp_geno[upper.tri(ldmat11_comp_geno)],
                   comp_genolike_flex = ldmat11_comp_genolike_flex[upper.tri(ldmat11_comp_genolike_flex)],
                   comp_genolike = ldmat11_comp_genolike[upper.tri(ldmat11_comp_genolike)])
df11 %>%
  transmute(`Haplotypic, Posterior\nMode Genotypes` = comp_genolike - hap_geno,
            `Haplotypic, Genotype\nLikelihoods` = comp_genolike - hap_genolike,
            `Composite, Posterior\nMean Genotypes` = comp_genolike - comp_geno,
            `Composite,\nGeneral Categorical` = comp_genolike - comp_genolike_flex) %>%
  gather(key = "Estimator", value = "Difference") %>%
  ggplot() +
  geom_boxplot(aes(x = Estimator, y = Difference)) +
  geom_hline(yintercept = 0, lty = 2, col = 2) +
  theme_bw() ->
  pl

ggsave(filename = "./output/uit/uit_fig/diff11.pdf",
       plot = pl,
       height = 3,
       width = 6,
       family = "Times")


df22 <- data.frame(hap_geno = ldmat22_hap_geno[upper.tri(ldmat22_hap_geno)],
                   hap_genolike = ldmat22_hap_genolike[upper.tri(ldmat22_hap_genolike)],
                   comp_geno = ldmat22_comp_geno[upper.tri(ldmat22_comp_geno)],
                   comp_genolike_flex = ldmat22_comp_genolike_flex[upper.tri(ldmat22_comp_genolike_flex)],
                   comp_genolike = ldmat22_comp_genolike[upper.tri(ldmat22_comp_genolike)])
df22 %>%
  transmute(`Haplotypic, Posterior\nMode Genotypes` = comp_genolike - hap_geno,
            `Haplotypic, Genotype\nLikelihoods` = comp_genolike - hap_genolike,
            `Composite, Posterior\nMean Genotypes` = comp_genolike - comp_geno,
            `Composite,\nGeneral Categorical` = comp_genolike - comp_genolike_flex) %>%
  gather(key = "Estimator", value = "Difference") %>%
  ggplot() +
  geom_boxplot(aes(x = Estimator, y = Difference)) +
  geom_hline(yintercept = 0, lty = 2, col = 2) +
  theme_bw() ->
  pl

ggsave(filename = "./output/uit/uit_fig/diff22.pdf",
       plot = pl,
       height = 3,
       width = 6,
       family = "Times")


df12 <- tibble(`Haplotypic,\nPosterior\nMode Genotypes` = c(ldmat12_hap_geno),
               `Haplotypic,\nGenotype\nLikelihoods` = c(ldmat12_hap_genolike),
               `Composite,\nPosterior\nMean Genotypes` = c(ldmat12_comp_geno),
               `Composite,\nGeneral\nCategorical` = c(ldmat12_comp_genolike_flex),
               `Composite,\nProportional\nNormal` = c(ldmat12_comp_genolike))
df12 %>%
  gather(key = "Estimator", value = "LD") %>%
  ggplot() +
  geom_boxplot(aes(x = Estimator, y = LD)) +
  geom_hline(yintercept = 0, lty = 2, col = 2) +
  theme_bw() ->
  pl

ggsave(filename = "./output/uit/uit_fig/box12.pdf",
       plot = pl,
       height = 3,
       width = 6,
       family = "Times")
