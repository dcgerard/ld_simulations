library(ldsep)
nind <- 100
nsnp <- 20
loglarray <- readRDS(file = "./output/ngs_out/updog_format.RDS")

ldsep_time <- system.time({
  ldsep_out <- mldest_genolike(genoarray = loglarray, nc = 1, pen = 1)
})

write.table(x = ldsep_out,
            file = "./output/ngs_out/lsep_out.tsv",
            sep = "\t",
            row.names = FALSE)

com <- paste0("./ngsLD/ngsLD",
              " --geno ./output/ngs_out/llike.tsv.gz",
              " --n_ind ", nind,
              " --n_sites ", nsnp,
              " --out ./output/ngs_out/ngs_fit.tsv",
              " --pos ./output/ngs_out/pos.tsv.gz",
              " --max_kb_dist 0",
              " --log_scale 1",
              " --n_threads 1")

ngsLD_time <- system.time({
  system(command = com)
})
