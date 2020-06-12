#####################
## Plot ngsLD and ldsep estimates
#####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

ngsout <- read_tsv("./output/ngs_out/ngs_fit.tsv",
                   col_names = FALSE,
                   col_types = cols(
                     X1 = col_character(),
                     X2 = col_character(),
                     X3 = col_double(),
                     X4 = col_double(),
                     X5 = col_double(),
                     X6 = col_double(),
                     X7 = col_double()
                   ))
names(ngsout) <- c("pos1", "pos2", "dist", "r2mom", "D", "Dprime", "r2")
ldsepout <- read_tsv(file = "./output/ngs_out/lsep_out.tsv")

ngsout %>%
  separate(pos1, into = c("chrome1", "i")) %>%
  separate(pos2, into = c("chrome2", "j")) %>%
  select(-chrome1, -chrome2) %>%
  mutate(i = parse_number(i),
         j = parse_number(j))->
  ngsout

full <- left_join(ldsepout, ngsout, by = c("i", "j"))

full %>%
  rename(r2_ngsLD = r2.x,
         D_ngsLD = D.x,
         Dprime_ngsLD = Dprime,
         r2_ldsep = r2.y,
         r_ldsep = r,
         D_ldsep = D.y) ->
  full

full %>%
  ggplot(aes(x = D_ngsLD, y = D_ldsep)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  xlab("D ngsLD") +
  ylab("D ldsep") ->
  pl

ggsave(filename = "./output/fig/D_ngsld_ldsep.pdf",
       plot = pl,
       height = 3,
       width = 3,
       family = "Times")


# library(corrplot)
# nloci <- max(max(ldsepout$i), max(ldsepout$j))
# cormat_ldsep <- matrix(NA_real_, ncol = nloci, nrow = nloci)
# cormat_ldsep[as.matrix(ldsepout[, c("i", "j")])] <- ldsepout[["r2"]]
#
# cormat_ngsLD <- matrix(NA_real_, ncol = nloci, nrow = nloci)
# cormat_ngsLD[as.matrix(ngsout[, c("i", "j")])] <- ngsout[["r2"]]
#
# corrplot(corr = cormat_ngsLD, type = "upper")
# corrplot(corr = cormat_ldsep, type = "upper")

