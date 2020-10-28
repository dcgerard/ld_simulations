library(ldsep)
library(tidyverse)
library(GGally)
library(gridExtra)
library(latex2exp)

uit_info <- read_csv("./output/uit/subuit_suc.csv")

ld_hap_geno <- readRDS("./output/uit/ldest_hap_geno.RDS")
ld_comp_geno <- readRDS("./output/uit/ldest_comp_geno.RDS")
ld_hap_genolike <- readRDS("./output/uit/ldest_hap_genolike.RDS")
ld_comp_genolike <- readRDS("./output/uit/ldest_comp_genolike.RDS")
ld_comp_genolike_flex <- readRDS("./output/uit/ldest_comp_genolike_flex.RDS")

pdf(file = "./output/uit/uit_fig/heat_hap_geno.pdf",
    height = 3,
    width = 3,
    family = "Times")
plot(ld_hap_geno,
     element = "r2",
     title = "")
dev.off()

pdf(file = "./output/uit/uit_fig/heat_hap_genolike.pdf",
    height = 3,
    width = 3,
    family = "Times")
plot(ld_hap_genolike,
     element = "r2",
     title = "")
dev.off()

pdf(file = "./output/uit/uit_fig/heat_comp_geno.pdf",
    height = 3,
    width = 3,
    family = "Times")
plot(ld_comp_geno,
     element = "r2",
     title = "")
dev.off()

pdf(file = "./output/uit/uit_fig/heat_comp_genolike_flex.pdf",
    height = 3,
    width = 3,
    family = "Times")
plot(ld_comp_genolike_flex,
     element = "r2",
     title = "")
dev.off()

pdf(file = "./output/uit/uit_fig/heat_comp_genolike.pdf",
    height = 3,
    width = 3,
    family = "Times")
plot(ld_comp_genolike,
     element = "r2",
     title = "")
dev.off()

ld_hap_geno %>%
  select(snpi, snpj) %>%
  left_join(uit_info[c("region", "Variant name")], by = c(snpi = "Variant name")) %>%
  left_join(uit_info[c("region", "Variant name")], by = c(snpj = "Variant name")) %>%
  mutate(within = if_else(region.x == region.y,
                          "within",
                          "between")) ->
  wbdf

ld_df <- tibble(hap_geno = ld_hap_geno$r2,
                hap_genolike = ld_hap_genolike$r2,
                comp_geno = ld_comp_geno$r2,
                comp_genolike = ld_comp_genolike_flex$r2,
                comp_genolike_norm = ld_comp_genolike$r2,
                `Within/Between` = wbdf$within)
names(ld_df) <- c(
  as.character(expression(paste(hat(r)^2, {}[g]))),
  as.character(expression(paste(hat(r)^2, {}[g][l]))),
  as.character(expression(paste(hat(rho)^2, {}[m][o][m]))),
  as.character(expression(paste(hat(rho)^2, {}[g][c]))),
  as.character(expression(paste(hat(rho)^2, {}[p][n]))),
  "Within/Between"
)

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

ld_df %>%
  filter(`Within/Between` == "within") %>%
  ggpairs(columns = seq_len(ncol(ld_df) - 1),
        lower = list(continuous = my_point),
        diag = list(continuous = my_hist),
        upper = "blank",
        labeller = label_parsed) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) ->
  pl

ggsave(filename = "./output/uit/uit_fig/uit_pairs.pdf",
       plot = pl,
       height = 6.6,
       width = 6.6,
       family = "Times")

ld_df %>%
  filter(`Within/Between` == "between") %>%
  ggpairs(columns = seq_len(ncol(ld_df) - 1),
          lower = list(continuous = my_point),
          diag = list(continuous = my_hist),
          upper = "blank",
          labeller = label_parsed) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) ->
  pl

ggsave(filename = "./output/uit/uit_fig/uit_pairs_between.pdf",
       plot = pl,
       height = 6.6,
       width = 6.6,
       family = "Times")

## Look at difference
# ldhapgeno <- format_lddf(obj = ld_hap_geno, element = "r2")
# stopifnot(uit_info$`Variant name`[uit_info$`Variant name` %in% rownames(ldhapgeno)] == rownames(ldhapgeno))
#
# uitsub <- uit_info[match(rownames(ldhapgeno), uit_info$`Variant name`), ]
# stopifnot(uitsub$`Variant name` == rownames(ldhapgeno))
#
# which_contig <- unique(uitsub$Contig)
# stopifnot(length(which_contig) == 2)
# is1 <- uitsub$Contig == which_contig[[1]]
# is2 <- uitsub$Contig == which_contig[[2]]
#
# ldmat11_hap_geno <- ldhapgeno[is1, is1]
# ldmat12_hap_geno <- ldhapgeno[is1, is2]
# ldmat22_hap_geno <- ldhapgeno[is2, is2]
#
# ldmat_hap_genolike <- format_lddf(obj = ld_hap_genolike, element = "r2")
# ldmat11_hap_genolike <- ldmat_hap_genolike[is1, is1]
# ldmat12_hap_genolike <- ldmat_hap_genolike[is1, is2]
# ldmat22_hap_genolike <- ldmat_hap_genolike[is2, is2]
#
# ldmat_comp_geno <- format_lddf(obj = ld_comp_geno, element = "r2")
# ldmat11_comp_geno <- ldmat_comp_geno[is1, is1]
# ldmat12_comp_geno <- ldmat_comp_geno[is1, is2]
# ldmat22_comp_geno <- ldmat_comp_geno[is2, is2]
#
# ldmat_comp_genolike_flex <- format_lddf(obj = ld_comp_genolike_flex, element = "r2")
# ldmat11_comp_genolike_flex <- ldmat_comp_genolike_flex[is1, is1]
# ldmat12_comp_genolike_flex <- ldmat_comp_genolike_flex[is1, is2]
# ldmat22_comp_genolike_flex <- ldmat_comp_genolike_flex[is2, is2]
#
# ldmat_comp_genolike <- format_lddf(obj = ld_comp_genolike, element = "r2")
# ldmat11_comp_genolike  <- ldmat_comp_genolike[is1, is1]
# ldmat12_comp_genolike  <- ldmat_comp_genolike[is1, is2]
# ldmat22_comp_genolike  <- ldmat_comp_genolike[is2, is2]
#
# df11 <- data.frame(hap_geno = ldmat11_hap_geno[upper.tri(ldmat11_hap_geno)],
#                    hap_genolike = ldmat11_hap_genolike[upper.tri(ldmat11_hap_genolike)],
#                    comp_geno = ldmat11_comp_geno[upper.tri(ldmat11_comp_geno)],
#                    comp_genolike_flex = ldmat11_comp_genolike_flex[upper.tri(ldmat11_comp_genolike_flex)],
#                    comp_genolike = ldmat11_comp_genolike[upper.tri(ldmat11_comp_genolike)])
#
# df11 %>%
#   transmute(hap_geno = comp_genolike - hap_geno,
#             hap_genolike = comp_genolike - hap_genolike,
#             comp_geno = comp_genolike - comp_geno,
#             comp_genolike_flex = comp_genolike - comp_genolike_flex) %>%
#   gather(key = "Estimator", value = "Difference") %>%
#   mutate(Estimator = parse_factor(Estimator,
#                                   levels = c("hap_geno",
#                                              "hap_genolike",
#                                              "comp_geno",
#                                              "comp_genolike_flex"))) %>%
#   ggplot() +
#   geom_boxplot(aes(x = Estimator, y = Difference), outlier.size = 0.3) +
#   scale_x_discrete(labels = c(
#     hap_geno = expression(paste(hat(r)^2, {}[g])),
#     hap_genolike = expression(paste(hat(r)^2, {}[g][l])),
#     comp_geno = expression(paste(hat(rho)^2, {}[m][o][m])),
#     comp_genolike_flex = expression(paste(hat(rho)^2, {}[g][c]))
#   )) +
#   geom_hline(yintercept = 0, lty = 2, col = 2) +
#   theme_bw() +
#   ylab(TeX("$\\hat{\\rho}^2_{pn}$ - Estimator")) +
#   ggtitle("(A)") ->
#   pl1
#
# ggsave(filename = "./output/uit/uit_fig/diff11.pdf",
#        plot = pl1,
#        height = 3,
#        width = 6,
#        family = "Times")
#
#
# df22 <- data.frame(hap_geno = ldmat22_hap_geno[upper.tri(ldmat22_hap_geno)],
#                    hap_genolike = ldmat22_hap_genolike[upper.tri(ldmat22_hap_genolike)],
#                    comp_geno = ldmat22_comp_geno[upper.tri(ldmat22_comp_geno)],
#                    comp_genolike_flex = ldmat22_comp_genolike_flex[upper.tri(ldmat22_comp_genolike_flex)],
#                    comp_genolike = ldmat22_comp_genolike[upper.tri(ldmat22_comp_genolike)])
# df22 %>%
#   transmute(hap_geno = comp_genolike - hap_geno,
#             hap_genolike = comp_genolike - hap_genolike,
#             comp_geno = comp_genolike - comp_geno,
#             comp_genolike_flex = comp_genolike - comp_genolike_flex) %>%
#   gather(key = "Estimator", value = "Difference") %>%
#   mutate(Estimator = parse_factor(Estimator,
#                                   levels = c("hap_geno",
#                                              "hap_genolike",
#                                              "comp_geno",
#                                              "comp_genolike_flex"))) %>%
#   ggplot() +
#   geom_boxplot(aes(x = Estimator, y = Difference), outlier.size = 0.3) +
#   scale_x_discrete(labels = c(
#     hap_geno = expression(paste(hat(r)^2, {}[g])),
#     hap_genolike = expression(paste(hat(r)^2, {}[g][l])),
#     comp_geno = expression(paste(hat(rho)^2, {}[m][o][m])),
#     comp_genolike_flex = expression(paste(hat(rho)^2, {}[g][c]))
#   )) +
#   geom_hline(yintercept = 0, lty = 2, col = 2) +
#   theme_bw() +
#   ylab(TeX("$\\hat{\\rho}^2_{pn}$ - Estimator")) +
#   ggtitle("(B)")  ->
#   pl2
#
# ggsave(filename = "./output/uit/uit_fig/diff22.pdf",
#        plot = pl2,
#        height = 3,
#        width = 6,
#        family = "Times")
#
#
# df12 <- tibble(hap_geno = c(ldmat12_hap_geno),
#                hap_genolike = c(ldmat12_hap_genolike),
#                comp_geno = c(ldmat12_comp_geno),
#                comp_genolike_flex = c(ldmat12_comp_genolike_flex),
#                comp_genolike_norm = c(ldmat12_comp_genolike))
# df12 %>%
#   gather(key = "Estimator", value = "LD") %>%
#   mutate(Estimator = parse_factor(Estimator,
#                                   levels = c("hap_geno",
#                                              "hap_genolike",
#                                              "comp_geno",
#                                              "comp_genolike_flex",
#                                              "comp_genolike_norm"))) %>%
#   ggplot() +
#   geom_boxplot(aes(x = Estimator, y = LD), outlier.size = 0.3) +
#   scale_x_discrete(labels = c(
#     hap_geno = expression(paste(hat(r)^2, {}[g])),
#     hap_genolike = expression(paste(hat(r)^2, {}[g][l])),
#     comp_geno = expression(paste(hat(rho)^2, {}[m][o][m])),
#     comp_genolike_flex = expression(paste(hat(rho)^2, {}[g][c])),
#     comp_genolike_norm = expression(paste(hat(rho)^2, {}[p][n]))
#   )) +
#   geom_hline(yintercept = 0, lty = 2, col = 2) +
#   theme_bw() +
#   ggtitle("(C)")  ->
#   pl3
#
# ggsave(filename = "./output/uit/uit_fig/box12.pdf",
#        plot = pl3,
#        height = 3,
#        width = 6,
#        family = "Times")
#
#
# pdf(file = "./output/uit/uit_fig/uit_box_combo.pdf",
#     width = 6.5,
#     height = 2.3,
#     family = "Times")
# gridExtra::grid.arrange(pl1, pl2, pl3, nrow = 1)
# dev.off()
