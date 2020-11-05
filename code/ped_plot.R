######################
## Plot the results of ped_sim.R
######################

library(tidyverse)
simdf <- read_csv("./output/ped/ped_sim_out.csv")

simdf %>%
  pivot_longer(cols = starts_with(c("hapgeno_", "compmean_", "haplike_", "complikeflex_", "complikenorm_")),
               names_to = "method_snp",
               values_to = "ldest") %>%
  separate(col = method_snp, into = c("method", "snp1", "snp2"), sep = "_") %>%
  mutate(snp1 = parse_number(snp1),
         snp2 = parse_number(snp2),
         true = case_when(snp1 == 50 & snp2 == 51 ~ true_A50_A51,
                          snp1 == 50 & snp2 == 60 ~ true_A50_A60,
                          snp1 == 50 & snp2 == 70 ~ true_A50_A70,
                          snp1 == 50 & snp2 == 80 ~ true_A50_A80,
                          snp1 == 50 & snp2 == 90 ~ true_A50_A90,
                          snp1 == 50 & snp2 == 99 ~ true_A50_A99,
                          snp1 == 50 & snp2 == 100 ~ true_A50_A100,
                          snp1 == 100 & snp2 == 50 ~ true_A100_A50,
                          snp1 == 100 & snp2 == 51 ~ true_A100_A51,
                          snp1 == 100 & snp2 == 60 ~ true_A100_A60,
                          snp1 == 100 & snp2 == 70 ~ true_A100_A70,
                          snp1 == 100 & snp2 == 80 ~ true_A100_A80,
                          snp1 == 100 & snp2 == 90 ~ true_A100_A90,
                          snp1 == 100 & snp2 == 99 ~ true_A100_A99)) %>%
  select(-starts_with("true_")) %>%
  mutate(diff = ldest - true,
         quadprop = case_when(near(quadprop, 0) ~ "0",
                              near(quadprop, 1/3) ~ "1/3",
                              near(quadprop, 2/3) ~ "2/3"),
         pref_pair = case_when(near(pref_pair, 0) ~ "1/3",
                               near(pref_pair, 0.5) ~ "2/3",
                               near(pref_pair, 1) ~ "1"),
         quadprop = parse_factor(quadprop, levels = c("0", "1/3", "2/3")),
         pref_pair = parse_factor(pref_pair, levels = c("1/3", "2/3", "1")),
         method = parse_factor(method, levels = c("hapgeno", "haplike", "compmean", "complikeflex", "complikenorm"))) ->
  simlong

## make sure I didn't miss a pairwise comparison
stopifnot(!any(is.na(simlong$true)))


labelvec_r <- c(hapgeno = expression(paste(hat(r)^2, {}[g])),
                haplike = expression(paste(hat(r)^2, {}[g][l])),
                compmean = expression(paste(hat(rho)^2, {}[m][o][m])),
                complikeflex = expression(paste(hat(rho)^2, {}[g][c])),
                complikenorm = expression(paste(hat(rho)^2, {}[p][n])))

snpcomparedf <- data.frame(snp1 = c(rep(50, 3), rep(100, 3)),
                           snp2 = c(51, 60, 70, 80, 90, 99))

for (index in seq_len(nrow(snpcomparedf))) {
  simlong %>%
    filter(snp1 == snpcomparedf$snp1[[index]], snp2 == snpcomparedf$snp2[[index]]) %>%
    ggplot(aes(x = method, y = diff)) +
    geom_boxplot() +
    facet_grid(pref_pair ~ quadprop) +
    geom_hline(yintercept = 0, lty = 2, col = 2) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    xlab("Method") +
    ylab("Difference") +
    scale_x_discrete(labels = labelvec_r) ->
    pl

  ggsave(filename = paste0("./output/ped/simplots/diff_",
                           snpcomparedf$snp1[[index]],
                           "_",
                           snpcomparedf$snp2[[index]],
                           ".pdf"),
         plot = pl,
         height = 6.5,
         width = 7,
         family = "Times")
}


###################
## Plots of True R2
###################

simdf %>%
  select(pref_pair, quadprop, starts_with("true_")) %>%
  pivot_longer(cols = starts_with("true_"), names_to = c("snp1_snp2"), values_to = "ld") %>%
  mutate(snp1_snp2 = str_remove(snp1_snp2, "^true_")) %>%
  separate(snp1_snp2, into = c("snp1", "snp2")) %>%
  mutate(snp1 = parse_number(snp1),
         snp2 = parse_number(snp2),
         snps = paste0("(", snp1, ",", snp2, ")"),
         quadprop = case_when(near(quadprop, 0) ~ "0",
                              near(quadprop, 1/3) ~ "1/3",
                              near(quadprop, 2/3) ~ "2/3"),
         pref_pair = case_when(near(pref_pair, 0) ~ "1/3",
                               near(pref_pair, 0.5) ~ "2/3",
                               near(pref_pair, 1) ~ "1"),
         quadprop = parse_factor(quadprop, levels = c("0", "1/3", "2/3")),
         pref_pair = parse_factor(pref_pair, levels = c("1/3", "2/3", "1"))) %>%
  filter((snp1 == 50 & snp2 == 51) |
           (snp1 == 50 & snp2 == 60) |
           (snp1 == 50 & snp2 == 70) |
           (snp1 == 100 & snp2 == 80) |
           (snp1 == 100 & snp2 == 90) |
           (snp1 == 100 & snp2 == 99)) ->
  truelong

truelong %>%
  ggplot(aes(x = snps, y = ld)) +
  geom_boxplot() +
  facet_grid(pref_pair ~ quadprop) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("SNP Pairs") +
  ylab(expression(paste(textstyle(True), ~~rho^2))) +
  ylim(0, 1) ->
  pl

ggsave(filename = "./output/ped/simplots/true_r2_box.pdf",
       plot = pl,
       height = 6.5,
       width = 7,
       family = "Times")
