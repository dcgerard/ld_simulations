####################
## Plot true r2
###################

library(tidyverse)
library(GGally)
library(corrplot)

pref_pair_vec <- c(1/3, 2/3, 1)
quadprop_vec <- c(0, 1/3, 2/3)
paramdf <- expand.grid(pref_pair = pref_pair_vec, quadprop = quadprop_vec)
r2_true_list <- readRDS("./output/ped/true_r/true_r_list.RDS")

pdf(file = "./output/ped/true_r/heatmap_pp_qq.pdf", height = 7, width = 7, family = "Times")
par(mfrow = c(3, 3))
for (index in seq_along(r2_true_list)) {
  corrplot(corr = r2_true_list[[index]],
           method = "color",
           type = "upper",
           diag = FALSE,
           title = paste0("quad = ", round(paramdf$quadprop[[index]], digits = 2), ", pp = ", round(paramdf$pref_pair[[index]], digits = 2)),
           mar=c(0,0,1,0))
}
dev.off()

trilist <- map(r2_true_list, ~.[upper.tri(.)])
names(trilist) <- paste0("quad = ",
                         round(paramdf$quadprop, digits = 2),
                         "\npp = ",
                         round(paramdf$pref_pair, digits = 2))

my_hist <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(..., color = "black") +
    xlim(0, 1)
}

my_point <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(..., alpha = 1/4) +
    geom_abline(intercept = 0, slope = 1, lty = 2, col = 2) +
    xlim(0, 1) +
    ylim(0, 1)
}

bind_rows(trilist) %>%
  ggpairs(lower = list(continuous = my_point),
          diag = list(continuous = my_hist),
          upper = "blank") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) ->
  pl

ggsave(filename = "./output/ped/true_r/truer_scatter.pdf",
       plot = pl,
       height = 7,
       width = 7,
       family = "Times")

