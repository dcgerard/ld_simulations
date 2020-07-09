############################
## Timing of MLE sims
############################

library(tidyverse)
library(ggthemes)
simdf <- read_csv("./output/mle/mle_sims_out.csv")

simdf %>%
  select(seed, size, ploidy, pA, pB, r,
         ends_with("time")) %>%
  gather(ends_with("time"), key = "method", value = "time") %>%
  mutate(method = str_replace(method, "_time$", "")) %>%
  group_by(size, ploidy, pA, pB, r, method) %>%
  summarise(mean = mean(time, na.rm = TRUE),
            lower = quantile(time, 0.025, na.rm = TRUE),
            upper = quantile(time, 0.975, na.rm = TRUE)) %>%
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

ggsave(filename = "./output/mle_plots/mle_time.pdf",
       plot = pl,
       width = 6.5,
       height = 4.5,
       family = "Times")
