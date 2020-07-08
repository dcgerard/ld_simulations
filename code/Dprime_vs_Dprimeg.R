###########################
## Compare Dprime and Dprimeg under HWE
###########################

library(ldsep)
library(tidyverse)
library(latex2exp)

nval <- 100
Kseq <- c(2, 4, 6, 8)
pAseq <- c(0.01, 0.25, 0.5, 0.75, 0.99)
pBseq <- pAseq
paramdf <- expand.grid(K = Kseq, pA = pAseq, pB = pBseq)
paramdf$data <- vector(mode = "list", length = nrow(paramdf))
paramdf <- as_tibble(paramdf)

for (index in seq_len(nrow(paramdf))) {
  pA <- paramdf$pA[[index]]
  pB <- paramdf$pB[[index]]
  K <- paramdf$K[[index]]

  Dmin <- max(-pA * pB, -(1 - pA) * (1 - pB))
  Dmax <- min(pA * (1 - pB), pB * (1 - pA))
  Dseq <- seq(Dmin, Dmax, length = nval)

  pAB_seq <- pA * pB + Dseq
  pAb_seq <- pA * (1 - pB) - Dseq
  paB_seq <- (1 - pA) * pB - Dseq
  pab_seq <- (1 - pA) * (1 - pB) + Dseq
  stopifnot(abs(pAB_seq + pAb_seq + paB_seq + pab_seq - 1) < 10^-6)

  Dprime_seq <- rep(NA_real_, times = nval)
  Dprime_seq[Dseq <= 0] <- Dseq[Dseq <= 0] / abs(Dmin)
  Dprime_seq[Dseq > 0] <- Dseq[Dseq > 0] / Dmax

  Dprimeg_seq <- rep(NA_real_, times = nval)
  for (i in seq_along(Dprimeg_seq)) {
    prob <- c(pab_seq[[i]], pAb_seq[[i]], paB_seq[[i]], pAB_seq[[i]])
    qmat <- get_prob_array(K = K, prob = prob, log_p = FALSE)
    Dprimeg_seq[[i]] <- Dprime(qmat = qmat, type = "geno")[["Dprime"]]
  }

  df <- data.frame(Dprime = Dprime_seq, Dprimeg = Dprimeg_seq)

  paramdf$data[[index]] <- df
}

paramdf <- unnest(paramdf, cols = data)

for (i in seq_along(Kseq)) {
  Know <- Kseq[[i]]
  paramdf %>%
    filter(K == Know) %>%
    mutate(pA = paste0("p[A]==", pA),
           pB = paste0("p[B]==", pB)) %>%
    ggplot(aes(x = Dprime, y = Dprimeg)) +
    geom_line() +
    facet_grid(pA ~ pB, labeller = label_parsed) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab(TeX("$D'$")) +
    ylab(TeX("$\\Delta'_g$")) +
    geom_abline(intercept = 0, slope = 1, col = 2, lty = 2) ->
    pl

  ggsave(filename = paste0("./output/dprime_v_dprimeg/dprimediff_", Know, ".pdf"),
         plot = pl,
         width = 6,
         height = 6,
         family = "Times")
}
