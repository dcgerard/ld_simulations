## properties of SNPs

library(tidyverse)
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
