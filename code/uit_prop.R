## properties of SNPs

set.seed(1)
library(tidyverse)
library(updog)
mout <- readRDS("./output/uit/uit_updog_fit.RDS")

mout$inddf %>%
  group_by(snp) %>%
  summarize(maf = 1 - mean(postmean) / 4) ->
  mafdf

ggplot(mafdf, aes(x = maf)) +
  geom_histogram(fill = "white", color = "black", bins = 10) +
  theme_bw() +
  xlab("Alternative Allele Frequency") ->
  pl

ggsave(filename = "./output/uit/uit_fig/maf.pdf",
       plot = pl,
       height = 2,
       width = 4,
       family = "Times")

mout$inddf %>%
  group_by(snp) %>%
  summarize(meansize = mean(size)) %>%
  ggplot(aes(x = meansize)) +
  geom_histogram(fill = "white", color = "black", bins = 10) +
  theme_bw() +
  xlab("Mean Read Depth") -> pl

ggsave(filename = "./output/uit/uit_fig/readdepth.pdf",
       plot = pl,
       height = 2,
       width = 4,
       family = "Times")

## HWE

#' x A vector of counts. x[i] is the number of individuals with dosage i-1.
basic_hwe_test <- function(x) {
  ploidy <- length(x) - 1
  alphahat <- sum(x / sum(x) * 0:ploidy) / ploidy
  theoprop <- dbinom(x = 0:ploidy, size = ploidy, prob = alphahat)
  chisq.test(x = x, p = theoprop, simulate.p.value = TRUE)$p.value
}

mout$inddf %>%
  group_by(snp, geno) %>%
  count() %>%
  ungroup() %>%
  pivot_wider(id_cols = "snp", values_from = "n", names_from = "geno") %>%
  select(-`NA`) ->
  tabdf
tabdf[is.na(tabdf)] <- 0
tabdf %>%
  select(snp, `0`, `1`, `2`, `3`, `4`) ->
  tabdf

countmat <- as.matrix(tabdf[-1])
class(countmat) <- "numeric"
pval <- apply(countmat, 1, basic_hwe_test)


qplot(pval, bins = 10, fill = I("white"), color = I("black")) +
  theme_bw() +
  xlab("p-value") +
  ylab("Count") ->
  pl

ggsave(filename = "./output/uit/uit_fig/hwe_p_hist.pdf",
       plot = pl,
       height = 3,
       width = 4, family = "Times")
