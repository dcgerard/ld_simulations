########
## Visualize the effects of preferential pairing
########

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggthemes)
ploidy <- 4
p1seq <- seq(0.1, 0.5, by = 0.1)
p2seq <- seq(0.9, 0.5, by = -0.1)

for (i in seq_along(p1seq)) {
  p1 <- p1seq[[i]]
  p2 <- p2seq[[i]]
  pbar <- (p1 + p2) / 2

  yauto <- dbinom(x = 0:ploidy, size = ploidy, prob = pbar)

  ygen1 <- dbinom(x = 0:(ploidy/2), size = ploidy/2, prob = p1)
  ygen2 <- dbinom(x = 0:(ploidy/2), size = ploidy/2, prob = p2)
  yallo <- convolve(x = ygen1, y = rev(ygen2), type = "open")

  if (i == 1) {
    probdf <- data.frame(x = 0:ploidy, yauto = yauto, yallo = yallo, p = paste("list(p[1]==", p1, ", p[2]==", p2, ")"))
  } else {
    tempdf <- data.frame(x = 0:ploidy, yauto = yauto, yallo = yallo, p = paste("list(p[1]==", p1, ", p[2]==", p2, ")"))
    probdf <- rbind(probdf, tempdf)
  }
}

problong <- pivot_longer(data = probdf,
                         cols = c("yauto", "yallo"),
                         names_to = "type",
                         values_to = "prob")
problong <- mutate(problong, Type = recode(type,
                                           "yauto" = "Autopolyploid",
                                           "yallo" = "Allopolyploid"),
                   x = x - 1/10 + (Type == "Autopolyploid") / 5)

ggplot(problong, aes(x = x, xend = x, y = 0, yend = prob, color = Type)) +
  geom_segment(lwd = 1.7) +
  facet_wrap(.~ p, labeller = label_parsed) +
  scale_color_colorblind() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  ylab("Pr(dosage)") +
  xlab("Dosage") ->
  pl

ggsave(filename = "./output/ped/allo_auto_dist.pdf",
       plot = pl,
       width = 6,
       height = 3,
       family = "Times")
