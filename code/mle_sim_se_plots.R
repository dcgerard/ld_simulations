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
            se_dprime = sd(mle_Dprime_est, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Depth = as.character(size == 1),
         Depth = recode(Depth,
                                `TRUE` = "1",
                                `FALSE` = "Other")) %>%
  gather(starts_with("se"), starts_with("meanse"), key = "type_measure", value = "value") %>%
  separate(col = "type_measure", into = c("type", "measure")) %>%
  spread(key = "type", value = "value") %>%
  mutate(measure = parse_factor(measure, levels = c("d", "dprime", "r2")),
         measure = recode(measure,
                          d = "hat(D)[g][l]",
                          dprime = "paste(hat(D), minute)[g][l]",
                          r2 = "paste(hat(r)^2, {})[g][l]")) %>%
  ggplot(aes(x = se, y = meanse, color = Depth, shape = Depth)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  xlab("True Standard Error") +
  ylab("Median of Estimated\nStandard Errors") +
  scale_color_colorblind() +
  facet_wrap(~measure,
             scales = "free", nrow = 1,
             labeller = label_parsed) +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) ->
  pl

ggsave(filename = "./output/mle_se_plots/mle_se_est.pdf",
       plot = pl, height = 2.2, width = 6.5, family = "Times")

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
            se_dprime = sd(comnorm_Dprime_est, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Depth = as.character(size == 1),
         Depth = recode(Depth,
                        `TRUE` = "1",
                        `FALSE` = "other")) %>%
  gather(starts_with("se"), starts_with("meanse"), key = "type_measure", value = "value") %>%
  separate(col = "type_measure", into = c("type", "measure")) %>%
  spread(key = "type", value = "value") %>%
  mutate(measure = parse_factor(measure, levels = c("d", "dprime", "r2")),
         measure = recode(measure,
                          d = "hat(Delta)[p][n]",
                          dprime = "paste(hat(Delta), minute, {}[p][n])",
                          r2 = "paste(hat(rho)^2, {}[p][n])")) %>%
  ggplot(aes(x = se, y = meanse, color = Depth, shape = Depth)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  xlab("True Standard Error") +
  ylab("Median of Estimated\nStandard Errors") +
  scale_color_colorblind() +
  facet_wrap(~measure,
             scales = "free",
             nrow = 1,
             labeller = label_parsed) +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) ->
  pl

ggsave(filename = "./output/mle_se_plots/comnorm_se_est.pdf",
       plot = pl, height = 2.2, width = 6.5, family = "Times")
