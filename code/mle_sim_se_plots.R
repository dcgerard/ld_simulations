library(tidyverse)
library(ggthemes)
library(latex2exp)
simdf <- read_csv("./output/mle/mle_sims_out.csv")




# Haplotypic, Genotype Likelihood -----------------------------------------
simdf %>%
  group_by(nind, size, ploidy, pA, pB, r) %>%
  summarize(meanse_d = mean(mle_D_se, na.rm = TRUE),
            se_d = sd(mle_D_est, na.rm = TRUE),
            meanse_r2 = mean(mle_r2_se, na.rm = TRUE),
            se_r2 = sd(mle_r2_est, na.rm = TRUE),
            meanse_dprime = mean(mle_Dprime_se, na.rm = TRUE),
            se_dprime = sd(mle_Dprime_est, na.rm = TRUE),
            meanse_z = mean(mle_z_se, na.rm = TRUE),
            se_z = sd(mle_z_est, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(`n and depth` = as.character(size == 1 & nind == 100),
         `n and depth` = recode(`n and depth`,
                                `TRUE` = "n = 100, depth = 1",
                                `FALSE` = "other")) %>%
  gather(starts_with("se"), starts_with("meanse"), key = "type_measure", value = "value") %>%
  separate(col = "type_measure", into = c("type", "measure")) %>%
  spread(key = "type", value = "value") %>%
  mutate(measure = recode(measure,
                          d = "D",
                          dprime = "D-prime",
                          r2 = "r-squared",
                          z = "z")) %>%
  ggplot(aes(x = se, y = meanse, color = `n and depth`, shape = `n and depth`)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  xlab("True Standard Error") +
  ylab("Mean of Estimated Standard Errors") +
  scale_color_colorblind() +
  facet_wrap(~measure,
             scales = "free") +
  theme(strip.background = element_rect(fill = "white")) ->
  pl

ggsave(filename = "./output/mle_se_plots/mle_se_est.pdf",
       plot = pl, height = 4, width = 6, family = "Times")

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

ggsave(filename = "./output/mle_se_plots/qq_nind100_pA50_pB50_r0.pdf",
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

ggsave(filename = "./output/mle_se_plots/qq_nind100_pA90_pB90_r0.pdf",
       plot = pl,
       height = 8,
       width = 6,
       family = "Times")

# Composite, Genotype Likelihoods ------------------------------------------
simdf %>%
  group_by(nind, size, ploidy, pA, pB, r) %>%
  summarize(meanse_d = mean(comnorm_D_se, na.rm = TRUE),
            se_d = sd(comnorm_D_est, na.rm = TRUE),
            meanse_r2 = mean(comnorm_r2_se, na.rm = TRUE),
            se_r2 = sd(comnorm_r2_est, na.rm = TRUE),
            meanse_dprime = mean(comnorm_Dprime_se, na.rm = TRUE),
            se_dprime = sd(comnorm_Dprime_est, na.rm = TRUE),
            meanse_z = mean(comnorm_z_se, na.rm = TRUE),
            se_z = sd(comnorm_z_est, na.rm = TRUE)) %>%
  ungroup() %>%
  gather(starts_with("se"), starts_with("meanse"), key = "type_measure", value = "value") %>%
  separate(col = "type_measure", into = c("type", "measure")) %>%
  spread(key = "type", value = "value") %>%
  mutate(measure = recode(measure,
                          d = "D",
                          dprime = "D-prime",
                          r2 = "r-squared",
                          z = "z")) %>%
  ggplot(aes(x = se, y = meanse)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  xlab("True Standard Error") +
  ylab("Mean of Estimated Standard Errors") +
  scale_color_colorblind() +
  facet_wrap(~measure,
             scales = "free") +
  theme(strip.background = element_rect(fill = "white")) ->
  pl

ggsave(filename = "./output/mle_se_plots/comnorm_se_est.pdf",
       plot = pl, height = 4, width = 4.5, family = "Times")
