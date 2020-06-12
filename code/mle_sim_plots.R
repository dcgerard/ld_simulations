library(tidyverse)
library(ggthemes)
library(latex2exp)
simdf <- read_csv("./output/mle/mle_sims_out.csv")

simdf %>%
  select(pA, pB, nind) %>%
  distinct() ->
  combodf


for (index in seq_len(nrow(combodf))) {
  pAnow <- combodf$pA[[index]]
  pBnow <- combodf$pB[[index]]
  nindnow <- combodf$nind[[index]]

  simdf %>%
    filter(pA == pAnow, pB == pBnow, nind == nindnow) %>%
    select(size,
           ploidy,
           D,
           r2,
           z,
           r,
           ends_with("D_est"),
           ends_with("r2_est"),
           ends_with("z_est")) %>%
    pivot_longer(cols = c("mle_D_est",
                          "gen_D_est",
                          "mom_D_est",
                          "com_D_est",
                          "mle_r2_est",
                          "gen_r2_est",
                          "mom_r2_est",
                          "com_r2_est",
                          "mle_z_est",
                          "gen_z_est",
                          "mom_z_est",
                          "com_z_est"),
                 names_to = "method_param",
                 values_to = "estimate")  %>%
    mutate(method_param = str_remove(method_param, "_est$")) %>%
    separate(method_param, into = c("method", "param")) %>%
    mutate(truth = case_when(param == "D" ~ D,
                             param == "r2" ~ r2,
                             param == "z" ~ z)) %>%
    mutate(truth = round(truth, digits = 2)) %>%
    mutate(method = recode(method,
                           "mom" = "Composite, Posterior\nMean Genotypes",
                           "com" = "Composite, Genotype\nLikelihoods",
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

}
