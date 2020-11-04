######################
## Plot the results of ped_sim.R
######################

library(tidyverse)
simdf <- read_csv("./output/ped/ped_sim_out.csv")

simdf %>%
  gather(hapgeno_A51:complikenorm_A90, key = "method_snp", value = "ld") %>%
  separate(col = method_snp, into = c("method", "SNP")) %>%
  mutate(true = case_when(str_detect(SNP, "A51") ~ true_A51,
                          str_detect(SNP, "A60") ~ true_A60,
                          str_detect(SNP, "A70") ~ true_A70,
                          str_detect(SNP, "A80") ~ true_A80,
                          str_detect(SNP, "A90") ~ true_A90)) %>%
  select(-contains("true_")) %>%
  filter(SNP == "A51", method == "hapgeno")
