library(tidyverse)
library(ggthemes)
library(latex2exp)
library(gridExtra)
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
  pl1

ggsave(filename = "./output/mle_se_plots/mle_se_est.pdf",
       plot = pl1, height = 2.2, width = 6.5, family = "Times")

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
  pl2

ggsave(filename = "./output/mle_se_plots/comnorm_se_est.pdf",
       plot = pl2, height = 2.2, width = 6.5, family = "Times")

## Combinded ----------------------------------------------------------------

pl1 <- pl1 + theme(axis.title.x = element_blank())
pl2 <- pl2 + theme(legend.position = "none")

lay <- rbind(c(1, 1),
             c(2, 3))

plblank <- ggplot() + ggplot2::geom_blank() + theme_minimal()


pdf(file = "./output/mle_se_plots/like_combined_se.pdf",
    width = 6.5,
    height = 4.5,
    family = "Times")
grid.arrange(pl1, pl2, plblank,
             layout_matrix = lay,
             widths = c(100, 16),
             heights = c(90, 100))
dev.off()

## All together -------------------------------------------------------------
simdf %>%
  select(seed,
         size,
         ploidy,
         pA,
         pB,
         r,
         starts_with("mle_"),
         starts_with("gen_"),
         starts_with("mom_"),
         starts_with("comnorm_")) %>%
  select(-ends_with("_time"),
         -contains("pab"),
         -contains("paB"),
         -contains("pAb"),
         -contains("pAB")) %>%
  gather(ends_with("_est"),
         ends_with("_se"),
         key = "method_estimator_sevest",
         value = "value") %>%
  filter(!is.na(value)) %>%
  separate(col = "method_estimator_sevest", into = c("Method", "Estimator", "sevest")) %>%
  spread(key = sevest, value = value) %>%
  group_by(size, ploidy, pA, pB, r, Method, Estimator) %>%
  summarize(sdest = sd(est, na.rm = TRUE),
            medse = median(se, na.rm = TRUE)) %>%
  mutate(Method = recode(Method,
                         comnorm = "pn",
                         gen = "g",
                         mle = "gl",
                         mom = "mom"),
         Estimator = recode(Estimator,
                            D = "D",
                            r = "r",
                            r2 = "r^2",
                            Dprime = "paste(D, minute)",
                            Dprimeg = "paste(Delta, minute)[g]")) ->
  sumdf


label_me <- function(labels) {
  label_parsed(labels, multi_line = FALSE)
}

sumdf %>%
  filter(size > 1) %>%
  group_by(Method, Estimator) %>%
  mutate(Estimator = parse_factor(Estimator),
         Method = parse_factor(Method)) %>%
  ggplot(aes(x = sdest, y = medse)) +
  facet_wrap( ~ Method + Estimator,
             labeller = label_me,
             scales = "free",
             drop = FALSE,
             ncol = 6) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = 2, lty = 2) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  xlab("True Standard Error") +
  ylab("Median of Estimated Standard Errors") ->
  pl

ggsave(filename = "./output/mle_se_plots/zscore_se.pdf",
       plot = pl,
       height = 4.5,
       width = 6.5,
       family = "Times")
