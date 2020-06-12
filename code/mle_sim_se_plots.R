library(tidyverse)
library(ggthemes)
simdf <- read_csv("../simulations/output/mle/mle_sims_out.csv")



# Haplotypic, Genotype Likelihood -----------------------------------------
simdf %>%
  group_by(nind, size, ploidy, pA, pB, z) %>%
  summarize(mean_z_se = mean(mle_z_se, na.rm = TRUE),
            se_z = sd(mle_z_est)) %>%
  mutate(`n and depth` = as.character(size == 1 & nind == 100),
         `n and depth` = recode(`n and depth`,
                                `TRUE` = "n = 100, depth = 1",
                                `FALSE` = "other")) %>%
  ggplot(aes(x = se_z, y = mean_z_se, color = `n and depth`, shape = `n and depth`)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  xlab("True Standard Error") +
  ylab("Mean of Estimated Standard Errors") +
  scale_color_colorblind() ->
  pl

ggsave(filename = "./output/fig/mle_se_est.pdf",
       plot = pl, height = 3, width = 6, family = "Times")

simdf %>%
  filter(nind == 100,
         pA == 0.5,
         pB == 0.5,
         r == 0) %>%
  ggplot(aes(sample = mle_z_est)) +
  facet_grid(size ~ ploidy) +
  geom_qq() +
  geom_qq_line(lty = 2, col = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) ->
  pl

ggsave(filename = "./output/fig/qq_nind100_pA50_pB50_r0.pdf",
       plot = pl,
       height = 8,
       width = 6,
       family = "Times")

simdf %>%
  filter(nind == 100,
         pA == 0.9,
         pB == 0.9,
         r == 0) %>%
  ggplot(aes(sample = mle_z_est)) +
  facet_grid(size ~ ploidy) +
  geom_qq() +
  geom_qq_line(lty = 2, col = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) ->
  pl

ggsave(filename = "./output/fig/qq_nind100_pA90_pB90_r0.pdf",
       plot = pl,
       height = 8,
       width = 6,
       family = "Times")


# Haplotypic, Posterior Modes ------------------------------------------
simdf %>%
  group_by(nind, size, ploidy, pA, pB, z) %>%
  summarize(mean_z_se = mean(gen_D_se, na.rm = TRUE),
            se_z = sd(gen_D_est)) %>%
  mutate(`n and depth` = as.character(size == 1),
         `n and depth` = recode(`n and depth`,
                                `TRUE` = "n = 100, depth = 1",
                                `FALSE` = "other")) %>%
  ggplot(aes(x = se_z, y = mean_z_se, color = `n and depth`, shape = `n and depth`)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  xlab("True Standard Error") +
  ylab("Mean of Estimated Standard Errors") +
  scale_color_colorblind() ->
  pl

pl

# Composite, Genotype Likelihoods ------------------------------------------
simdf %>%
  group_by(nind, size, ploidy, pA, pB, z) %>%
  summarize(mean_z_se = mean(com_D_se, na.rm = TRUE),
            se_z = sd(com_D_est)) %>%
  mutate(`n and depth` = as.character(size == 1),
         `n and depth` = recode(`n and depth`,
                                `TRUE` = "n = 100, depth = 1",
                                `FALSE` = "other")) %>%
  ggplot(aes(x = se_z, y = mean_z_se, color = as.factor(ploidy), shape = `n and depth`)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  xlab("True Standard Error") +
  ylab("Mean of Estimated Standard Errors") +
  scale_color_colorblind() ->
  pl

pl
