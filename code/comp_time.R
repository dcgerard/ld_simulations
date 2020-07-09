#######################
## Time comparisons for comp sims
#######################

library(tidyverse)
library(ggthemes)
simdf <- read_csv("./output/comp/comp_sims_out.csv")

simdf %>%
  select(seed, size, ploidy, mu1, mu2, sigma11, sigma22, sigma12, ends_with("time")) %>%
  mutate(p1 = round(mu1 / ploidy, digits = 2),
         p2 = round(mu2 / ploidy, digits = 2),
         sigmaprop = round(sigma11 / ploidy^2, digits = 2),
         covprop = round(sigma12 / sqrt(sigma11 * sigma22), digits = 2)) %>%
  select(size, ploidy, p1, p2, sigmaprop, covprop, ends_with("time")) %>%
  gather(ends_with("time"), key = "method", value = "time") %>%
  mutate(method = str_replace(method, "_time$", "")) %>%
  group_by(size, ploidy, p1, p2, sigmaprop, covprop, method) %>%
  summarize(mean = mean(time),
            lower = quantile(time, 0.025),
            upper = quantile(time, 0.975)) %>%
  group_by(method) %>%
  mutate(index = row_number()) %>%
  ungroup() %>%
  mutate(Ploidy = as.character(ploidy),
         method = parse_factor(method, levels = c("gen", "mle", "mom", "com", "comnorm")),
         method = recode(method,
                         gen = as.character(expression(list(hat(D)[g], paste(hat(D), minute, {}[g]), paste(hat(r)^2, {}[g])))),
                         mle = as.character(expression(list(hat(D)[g][l], paste(hat(D), minute, {}[g][l]), paste(hat(r)^2, {}[g][l])))),
                         mom = as.character(expression(list(hat(Delta)[m][o][m], paste(hat(Delta), minute, {}[m][o][m]), paste(hat(rho)^2, {}[m][o][m])))),
                         com = as.character(expression(list(hat(Delta)[g][c], paste(hat(Delta), minute, {}[g][c]), paste(hat(rho)^2, {}[g][c])))),
                         comnorm = as.character(expression(list(hat(Delta)[p][n], paste(hat(Delta), minute, {}[p][n]), paste(hat(rho)^2, {}[p][n]))))
         )) %>%
  ggplot() +
  geom_point(mapping = aes(x = index, y = mean, color = Ploidy, shape = Ploidy)) +
  facet_wrap(~method, labeller = label_parsed) +
  ylab("Time (seconds)") +
  xlab("Scenario") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_colorblind() ->
  pl

ggsave(filename = "./output/comp/comp_plots/comp_time.pdf",
       plot = pl,
       width = 6.5,
       height = 4.5,
       family = "Times")
