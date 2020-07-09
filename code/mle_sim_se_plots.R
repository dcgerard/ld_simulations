library(tidyverse)
library(ggthemes)
library(latex2exp)
simdf <- read_csv("./output/mle/mle_sims_out.csv")

# Haplotypic, Genotype Likelihood -----------------------------------------
simdf %>%
  group_by(nind, size, ploidy, pA, pB, r) %>%
  summarize(meanse_d = median(mle_D_se, na.rm = TRUE),
            se_d = sd(mle_D_est, na.rm = TRUE),
            meanse_r2 = median(mle_r2_se, na.rm = TRUE),
            se_r2 = sd(mle_r2_est, na.rm = TRUE),
            meanse_dprime = median(mle_Dprime_se, na.rm = TRUE),
            se_dprime = sd(mle_Dprime_est, na.rm = TRUE),
            meanse_z = median(mle_z_se, na.rm = TRUE),
            se_z = sd(mle_z_est, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(`n and depth` = as.character(size == 1 & nind == 100),
         `n and depth` = recode(`n and depth`,
                                `TRUE` = "n = 100,\ndepth = 1",
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
  ylab("Median of Estimated\nStandard Errors") +
  scale_color_colorblind() +
  facet_wrap(~measure,
             scales = "free", nrow = 1) +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) ->
  pl

ggsave(filename = "./output/mle_se_plots/mle_se_est.pdf",
       plot = pl, height = 1.8, width = 6.5, family = "Times")

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
  summarize(meanse_d = median(comnorm_D_se, na.rm = TRUE),
            se_d = sd(comnorm_D_est, na.rm = TRUE),
            meanse_r2 = median(comnorm_r2_se, na.rm = TRUE),
            se_r2 = sd(comnorm_r2_est, na.rm = TRUE),
            meanse_dprime = median(comnorm_Dprime_se, na.rm = TRUE),
            se_dprime = sd(comnorm_Dprime_est, na.rm = TRUE),
            meanse_z = median(comnorm_z_se, na.rm = TRUE),
            se_z = sd(comnorm_z_est, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(depth = as.character(size == 1),
         depth = recode(depth,
                        `TRUE` = "1",
                        `FALSE` = "other")) %>%
  gather(starts_with("se"), starts_with("meanse"), key = "type_measure", value = "value") %>%
  separate(col = "type_measure", into = c("type", "measure")) %>%
  spread(key = "type", value = "value") %>%
  mutate(measure = recode(measure,
                          d = "D",
                          dprime = "D-prime",
                          r2 = "r-squared",
                          z = "z")) %>%
  ggplot(aes(x = se, y = meanse, color = depth, shape = depth)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  xlab("True Standard Error") +
  ylab("Median of Estimated\nStandard Errors") +
  scale_color_colorblind() +
  facet_wrap(~measure,
             scales = "free",
             nrow = 1) +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) ->
  pl

ggsave(filename = "./output/mle_se_plots/comnorm_se_est.pdf",
       plot = pl, height = 1.8, width = 6.5, family = "Times")
