library(tidyverse)
library(ggthemes)
library(latex2exp)
simdf <- read_csv("./output/mle/mle_sims_out.csv")

simdf %>%
  select(pA, pB, nind) %>%
  distinct() ->
  combodf

get_Dprime <- function(pab, pAb, paB, pAB) {
  pA <- pAb + pAB
  pB <- paB + pAB
  D  <- pAB - pA * pB
  Dprime <- ifelse(D < 0,
                   D / pmin(pA * pB, (1 - pA) * (1 - pB)),
                   D / pmin(pA * (1 - pB), (1 - pA) * pB))
  return(Dprime)
}

for (index in seq_len(nrow(combodf))) {
  pAnow <- combodf$pA[[index]]
  pBnow <- combodf$pB[[index]]
  nindnow <- combodf$nind[[index]]

  simdf %>%
    mutate(Dprime = get_Dprime(pab = pab, pAb = pAb, paB = paB, pAB = pAB)) %>%
    filter(pA == pAnow, pB == pBnow, nind == nindnow) %>%
    select(size,
           ploidy,
           D,
           r2,
           z,
           Dprime,
           r,
           ends_with("D_est"),
           ends_with("r2_est"),
           ends_with("Dprime_est")) %>%
    pivot_longer(cols = c("mle_D_est",
                          "gen_D_est",
                          "mom_D_est",
                          "com_D_est",
                          "comnorm_D_est",
                          "mle_r2_est",
                          "gen_r2_est",
                          "mom_r2_est",
                          "com_r2_est",
                          "comnorm_r2_est",
                          "mle_Dprime_est",
                          "gen_Dprime_est",
                          "mom_Dprime_est",
                          "com_Dprime_est",
                          "comnorm_Dprime_est"),
                 names_to = "method_param",
                 values_to = "estimate")  %>%
    mutate(method_param = str_remove(method_param, "_est$")) %>%
    separate(method_param, into = c("method", "param")) %>%
    mutate(truth = case_when(param == "D" ~ D,
                             param == "r2" ~ r2,
                             param == "Dprime" ~ Dprime)) %>%
    mutate(truth = round(truth, digits = 2)) %>%
    mutate(method = recode(method,
                           "comnorm" = "Composite,\nProportional Normal",
                           "mom" = "Composite, Posterior\nMean Genotypes",
                           "com" = "Composite,\nGeneral Categorical",
                           "mle" = "Haplotypic, Genotype\nLikelihoods",
                           "gen" = "Haplotypic, Posterior\nMode Genotypes")) %>%
    group_by(size, ploidy, method, param, truth) %>%
    summarize(bias = mean(estimate - truth, na.rm = TRUE),
              biaslow = stats::quantile(estimate - truth, probs = 0.025, na.rm = TRUE),
              biashigh = stats::quantile(estimate - truth, probs = 0.975, na.rm = TRUE),
              se = sd(estimate, na.rm = TRUE),
              mse = mean((estimate - truth)^2), na.rm = TRUE) ->
    sumdf

  ## D plots ------------------------------------------------
  sumdf %>%
    filter(param == "D") %>%
    ggplot(aes(x = size, y = bias, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ truth) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_colorblind() +
    xlab("Read Depth") +
    ylab("D Bias") ->
    pl

  ggsave(filename = paste0("./output/mle_plots/D_bias_nind_",
                           nindnow,
                           "_pA_",
                           round(pAnow * 100),
                           "_pB_",
                           round(pBnow * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  sumdf %>%
    filter(param == "D") %>%
    ggplot(aes(x = size, y = se, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ truth) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    geom_line() +
    geom_point() +
    scale_color_colorblind() +
    xlab("Read Depth") +
    ylab("D Standard Error") +
    scale_y_log10() ->
    pl

  ggsave(filename = paste0("./output/mle_plots/D_se_nind_",
                           nindnow,
                           "_pA_",
                           round(pAnow * 100),
                           "_pB_",
                           round(pBnow * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  sumdf %>%
    filter(param == "D") %>%
    ggplot(aes(x = size, y = mse, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ truth) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    geom_line() +
    geom_point() +
    scale_color_colorblind() +
    scale_y_log10() +
    xlab("Read Depth") +
    ylab("D MSE") ->
    pl

  ggsave(filename = paste0("./output/mle_plots/D_mse_nind_",
                           nindnow,
                           "_pA_",
                           round(pAnow * 100),
                           "_pB_",
                           round(pBnow * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  ## R2 plots --------------------------------------------------
  sumdf %>%
    filter(param == "r2") %>%
    ggplot(aes(x = size, y = bias, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ truth) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_colorblind() +
    xlab("Read Depth") +
    ylab(TeX("$r^2$ Bias")) ->
    pl

  ggsave(filename = paste0("./output/mle_plots/r2_bias_nind_",
                           nindnow,
                           "_pA_",
                           round(pAnow * 100),
                           "_pB_",
                           round(pBnow * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  sumdf %>%
    filter(param == "r2") %>%
    ggplot(aes(x = size, y = se, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ truth) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    geom_line() +
    geom_point() +
    scale_color_colorblind() +
    xlab("Read Depth") +
    ylab(TeX("$r^2$ Standard Error")) +
    scale_y_log10() ->
    pl

  ggsave(filename = paste0("./output/mle_plots/r2_se_nind_",
                           nindnow,
                           "_pA_",
                           round(pAnow * 100),
                           "_pB_",
                           round(pBnow * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  sumdf %>%
    filter(param == "r2") %>%
    ggplot(aes(x = size, y = mse, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ truth) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    geom_line() +
    geom_point() +
    scale_color_colorblind() +
    xlab("Read Depth") +
    ylab(TeX("$r^2$ MSE")) +
    scale_y_log10() ->
    pl

  ggsave(filename = paste0("./output/mle_plots/r2_mse_nind_",
                           nindnow,
                           "_pA_",
                           round(pAnow * 100),
                           "_pB_",
                           round(pBnow * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")



  ## Dprime plots --------------------------------------------------
  sumdf %>%
    filter(param == "Dprime") %>%
    ggplot(aes(x = size, y = bias, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ truth) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_colorblind() +
    xlab("Read Depth") +
    ylab(TeX("$D'$ Bias")) ->
    pl

  ggsave(filename = paste0("./output/mle_plots/Dprime_bias_nind_",
                           nindnow,
                           "_pA_",
                           round(pAnow * 100),
                           "_pB_",
                           round(pBnow * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  sumdf %>%
    filter(param == "Dprime") %>%
    ggplot(aes(x = size, y = se, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ truth) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    geom_line() +
    geom_point() +
    scale_color_colorblind() +
    xlab("Read Depth") +
    ylab(TeX("$r^2$ Standard Error")) +
    scale_y_log10() ->
    pl

  ggsave(filename = paste0("./output/mle_plots/Dprime_se_nind_",
                           nindnow,
                           "_pA_",
                           round(pAnow * 100),
                           "_pB_",
                           round(pBnow * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  sumdf %>%
    filter(param == "Dprime") %>%
    ggplot(aes(x = size, y = mse, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ truth) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    geom_line() +
    geom_point() +
    scale_color_colorblind() +
    xlab("Read Depth") +
    ylab(TeX("$r^2$ MSE")) +
    scale_y_log10() ->
    pl

  ggsave(filename = paste0("./output/mle_plots/Dprime_mse_nind_",
                           nindnow,
                           "_pA_",
                           round(pAnow * 100),
                           "_pB_",
                           round(pBnow * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

}
