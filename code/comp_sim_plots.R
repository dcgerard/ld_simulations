library(tidyverse)
library(ggthemes)
library(latex2exp)
simdf <- read_csv("./output/comp/comp_sims_out.csv")

simdf %>%
  mutate(mu1prop = round(mu1 / ploidy, digits = 2),
         mu2prop = round(mu2 / ploidy, digits = 2),
         sigma_prop = round(sigma11 / ploidy ^2, digits = 2),
         cor = round(sigma12 / sqrt(sigma11 * sigma22), digits = 2)) ->
  simdf

simdf %>%
  select(mu1prop, mu2prop, sigma_prop) %>%
  distinct() ->
  combodf

## legend stuff
colorvec <- ggthemes::colorblind_pal()(5)
names(colorvec) <- c("gen", "mle", "mom", "com", "comnorm")
labelvec_D <- c(gen = expression(hat(D)[g]),
                mle = expression(hat(D)[g][l]),
                mom = expression(hat(Delta)[m][o][m]),
                com = expression(hat(Delta)[g][c]),
                comnorm = expression(hat(Delta)[p][n]))

labelvec_r <- c(gen = expression(paste(hat(r)^2, {}[g])),
                mle = expression(paste(hat(r)^2, {}[g][l])),
                mom = expression(paste(hat(rho)^2, {}[m][o][m])),
                com = expression(paste(hat(rho)^2, {}[g][c])),
                comnorm = expression(paste(hat(rho)^2, {}[p][n])))

labelvec_Dprime <- c(gen = expression(paste(hat(D), minute, {}[g])),
                     mle = expression(paste(hat(D), minute, {}[g][l])),
                     mom = expression(paste(hat(Delta), minute, {}[m][o][m])),
                     com = expression(paste(hat(Delta), minute, {}[g][c])),
                     comnorm = expression(paste(hat(Delta), minute, {}[p][n])))

for (index in seq_len(nrow(combodf))) {
  mu1prop_now <- combodf$mu1prop[[index]]
  mu2prop_now <- combodf$mu2prop[[index]]
  sigmaprop_now <- combodf$sigma_prop[[index]]

  simdf %>%
    filter(mu1prop == mu1prop_now,
           mu2prop == mu2prop_now,
           sigma_prop == sigmaprop_now) %>%
    select(size,
           ploidy,
           D,
           r2,
           z,
           Dprime,
           r,
           cor,
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
    group_by(size, ploidy, method, param, truth, cor,) %>%
    summarize(bias = mean(estimate - truth, na.rm = TRUE),
              biaslow = stats::quantile(estimate - truth, probs = 0.025, na.rm = TRUE),
              biashigh = stats::quantile(estimate - truth, probs = 0.975, na.rm = TRUE),
              se = sd(estimate, na.rm = TRUE),
              mse = mean((estimate - truth)^2), na.rm = TRUE) %>%
    mutate(method = parse_factor(method, levels = c("gen", "mle", "mom", "com", "comnorm"))) ->
    sumdf

  ## D plots ------------------------------------------------
  sumdf %>%
    filter(param == "D") %>%
    ggplot(aes(x = size, y = bias, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ cor) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"), legend.text.align = 0) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_manual(values = colorvec, labels = labelvec_D) +
    scale_shape(labels = labelvec_D) +
    scale_linetype(labels = labelvec_D) +
    xlab("Read Depth") +
    ylab(TeX("$\\Delta$ Bias")) ->
    pl

  ggsave(filename = paste0("./output/comp/comp_plots/D_bias_mu1prop_",
                           round(mu1prop_now * 100),
                           "_mu2prop_",
                           round(mu2prop_now * 100),
                           "_sigmaprop_",
                           round(sigmaprop_now * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  sumdf %>%
    filter(param == "D") %>%
    ggplot(aes(x = size, y = se, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ cor) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"), legend.text.align = 0) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = colorvec, labels = labelvec_D) +
    scale_shape(labels = labelvec_D) +
    scale_linetype(labels = labelvec_D) +
    xlab("Read Depth") +
    ylab(TeX("$\\Delta$ Standard Error")) +
    scale_y_log10() ->
    pl

  ggsave(filename = paste0("./output/comp/comp_plots/D_se_mu1prop_",
                           round(mu1prop_now * 100),
                           "_mu2prop_",
                           round(mu2prop_now * 100),
                           "_sigmaprop_",
                           round(sigmaprop_now * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  sumdf %>%
    filter(param == "D") %>%
    ggplot(aes(x = size, y = mse, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ cor) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"), legend.text.align = 0) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = colorvec, labels = labelvec_D) +
    scale_shape(labels = labelvec_D) +
    scale_linetype(labels = labelvec_D) +
    scale_y_log10() +
    xlab("Read Depth") +
    ylab(TeX("$\\Delta$ MSE")) ->
    pl

  ggsave(filename = paste0("./output/comp/comp_plots/D_mse_mu1prop_",
                           round(mu1prop_now * 100),
                           "_mu2prop_",
                           round(mu2prop_now * 100),
                           "_sigmaprop_",
                           round(sigmaprop_now * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  ## R2 plots --------------------------------------------------
  sumdf %>%
    filter(param == "r2") %>%
    ggplot(aes(x = size, y = bias, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ cor) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"), legend.text.align = 0) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_manual(values = colorvec, labels = labelvec_r) +
    scale_shape(labels = labelvec_r) +
    scale_linetype(labels = labelvec_r) +
    xlab("Read Depth") +
    ylab(TeX("$\\rho^2$ Bias")) ->
    pl

  ggsave(filename = paste0("./output/comp/comp_plots/r2_bias_mu1prop_",
                           round(mu1prop_now * 100),
                           "_mu2prop_",
                           round(mu2prop_now * 100),
                           "_sigmaprop_",
                           round(sigmaprop_now * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  sumdf %>%
    filter(param == "r2") %>%
    ggplot(aes(x = size, y = se, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ cor) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"), legend.text.align = 0) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = colorvec, labels = labelvec_r) +
    scale_shape(labels = labelvec_r) +
    scale_linetype(labels = labelvec_r) +
    xlab("Read Depth") +
    ylab(TeX("$\\rho^2$ Standard Error")) +
    scale_y_log10() ->
    pl

  ggsave(filename = paste0("./output/comp/comp_plots/r2_se_mu1prop_",
                           round(mu1prop_now * 100),
                           "_mu2prop_",
                           round(mu2prop_now * 100),
                           "_sigmaprop_",
                           round(sigmaprop_now * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  sumdf %>%
    filter(param == "r2") %>%
    ggplot(aes(x = size, y = mse, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ cor) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"), legend.text.align = 0) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = colorvec, labels = labelvec_r) +
    scale_shape(labels = labelvec_r) +
    scale_linetype(labels = labelvec_r) +
    xlab("Read Depth") +
    ylab(TeX("$\\rho^2$ MSE")) +
    scale_y_log10() ->
    pl

  ggsave(filename = paste0("./output/comp/comp_plots/r2_mse_mu1prop_",
                           round(mu1prop_now * 100),
                           "_mu2prop_",
                           round(mu2prop_now * 100),
                           "_sigmaprop_",
                           round(sigmaprop_now * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")



  ## Dprime plots --------------------------------------------------
  sumdf %>%
    filter(param == "Dprime") %>%
    ggplot(aes(x = size, y = bias, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ cor) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"), legend.text.align = 0) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_manual(values = colorvec, labels = labelvec_Dprime) +
    scale_shape(labels = labelvec_Dprime) +
    scale_linetype(labels = labelvec_Dprime) +
    xlab("Read Depth") +
    ylab(TeX("$\\Delta'$ Bias")) ->
    pl

  ggsave(filename = paste0("./output/comp/comp_plots/Dprime_bias_mu1prop_",
                           round(mu1prop_now * 100),
                           "_mu2prop_",
                           round(mu2prop_now * 100),
                           "_sigmaprop_",
                           round(sigmaprop_now * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  sumdf %>%
    filter(param == "Dprime") %>%
    ggplot(aes(x = size, y = se, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ cor) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"), legend.text.align = 0) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = colorvec, labels = labelvec_Dprime) +
    scale_shape(labels = labelvec_Dprime) +
    scale_linetype(labels = labelvec_Dprime) +
    xlab("Read Depth") +
    ylab(TeX("$\\Delta'$ Standard Error")) +
    scale_y_log10() ->
    pl

  ggsave(filename = paste0("./output/comp/comp_plots/Dprime_se_mu1prop_",
                           round(mu1prop_now * 100),
                           "_mu2prop_",
                           round(mu2prop_now * 100),
                           "_sigmaprop_",
                           round(sigmaprop_now * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")

  sumdf %>%
    filter(param == "Dprime") %>%
    ggplot(aes(x = size, y = mse, color = method, lty = method, shape = method)) +
    facet_grid(ploidy ~ cor) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"), legend.text.align = 0) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = colorvec, labels = labelvec_Dprime) +
    scale_shape(labels = labelvec_Dprime) +
    scale_linetype(labels = labelvec_Dprime) +
    xlab("Read Depth") +
    ylab(TeX("$\\Delta'$ MSE")) +
    scale_y_log10() ->
    pl

  ggsave(filename = paste0("./output/comp/comp_plots/Dprime_mse_mu1prop_",
                           round(mu1prop_now * 100),
                           "_mu2prop_",
                           round(mu2prop_now * 100),
                           "_sigmaprop_",
                           round(sigmaprop_now * 100),
                           ".pdf"),
         plot = pl,
         height = 6,
         width = 6,
         family = "Times")
}
