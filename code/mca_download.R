## Download McAllister data from dryad
library(rdryad)
mcdat <- dryad_download(dois = "10.5061/dryad.05qs7")

for (i in seq_along(mcdat[[1]])) {
  cmd <- paste0("cp ", mcdat[[1]][[i]], " ./data/gerardii")
  system(cmd)
}
