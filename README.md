
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduce the Results of Gerard (2020)

## Introduction

This repository contains the code and instructions needed to reproduce
all of the results from Gerard (2020).

If you find a bug, please create an
[issue](https://github.com/dcgerard/ld_simulations/issues).

## Instructions

You’ll need to have up-to-date versions of R and python3 to run these
simulatoins.

1.  Install the appropriate R packages
    
    ``` r
    install.packages("BiocManager")
    BiocManager::install(c("tidyverse",
                           "ggthemes",
                           "latex2exp",
                           "Rcpp",
                           "RcppArmadillo",
                           "foreach",
                           "doParallel",
                           "corrplot",
                           "devtools",
                           "vcfR",
                           "VariantAnnotation"))
    devtools::install_github("dcgerard/updog")
    devtools::install_github("dcgerard/ldsep")
    ```

2.  Install ngsLD (Fox et al. 2019) in the working directory of the
    simulation files: <https://github.com/fgvieira/ngsLD>
    
    ``` bash
    git clone https://github.com/fgvieira/ngsLD.git
    cd ngsLD
    make
    ```

3.  Download the following files from Uitdewilligen (2013), place them
    in the “data” folder, and extract them.
    
      - <https://doi.org/10.1371/journal.pone.0062355.s004>
      - <https://doi.org/10.1371/journal.pone.0062355.s007>
      - <https://doi.org/10.1371/journal.pone.0062355.s009>
      - <https://doi.org/10.1371/journal.pone.0062355.s010>
    
    You can do this with the following bash code from the top directory
    of the ld\_simulations
    repo.
    
    ``` bash
    wget --directory-prefix=data -nc https://doi.org/10.1371/journal.pone.0062355.s004
    mv ./data/journal.pone.0062355.s004 ./data/journal.pone.0062355.s004.gz
    wget --directory-prefix=data -nc https://doi.org/10.1371/journal.pone.0062355.s007
    mv ./data/journal.pone.0062355.s007 ./data/journal.pone.0062355.s007.gz 
    wget --directory-prefix=data -nc https://doi.org/10.1371/journal.pone.0062355.s009
    mv ./data/journal.pone.0062355.s009 ./data/journal.pone.0062355.s009.xls 
    wget --directory-prefix=data -nc https://doi.org/10.1371/journal.pone.0062355.s010
    mv ./data/journal.pone.0062355.s010 ./data/journal.pone.0062355.s010.xls 
    7z e ./data/journal.pone.0062355.s004.gz -o./data/
    7z e ./data/journal.pone.0062355.s007.gz -o./data/
    ```
    
    As long as you have 7z installed, `make` will attempt to download
    and extract the files for you.

4.  Download the data from McAllister and Miller (2016)
    <https://doi.org/10.5061/dryad.05qs7> and unzip it. I unzipped using
    7z:
    
    ``` bash
    7z e doi_10.5061_dryad.05qs7__v1.zip 
    ```

5.  Run `make`. This will run all of the simulations.
    
    ``` bash
    make
    ```

You may choose to run only part of the simulations `make mle`, `make
ngsLD`, `make uit`, `make mca`, `make norm`, or `make comp`.

6.  Get coffee/sweets. Running `make sims` should take a few hours. You
    should get some coffee\! Here is a list of some of my favorite
    places:
    
      - Washington, DC
          - [Colony
            Club](https://www.yelp.com/biz/colony-club-washington)
          - [Little Red
            Fox](https://www.yelp.com/biz/little-red-fox-washington)
          - [Three
            Fifty](https://www.yelp.com/biz/three-fifty-bakery-and-coffee-bar-washington)
            (for the baked goods, coffee is meh)
      - Chicago
          - [Sawada
            Coffee](https://www.yelp.com/biz/sawada-coffee-chicago)
          - [Plein Air
            Cafe](https://www.yelp.com/biz/plein-air-cafe-and-eatery-chicago-2)
      - Seattle
          - [Bauhaus
            Ballard](https://www.yelp.com/biz/bauhaus-ballard-seattle)
          - [Cafe
            Solstice](https://www.yelp.com/biz/cafe-solstice-seattle)
      - Columbus
          - [Yeah, Me
            Too](https://www.yelp.com/biz/yeah-me-too-columbus)
          - [Stauf’s Coffee
            Roasters](https://www.yelp.com/biz/staufs-coffee-roasters-columbus-2)

Note that I’ve also only tried this on Ubuntu.

# Session Information

    #> R version 4.0.2 (2020-06-22)
    #> Platform: x86_64-pc-linux-gnu (64-bit)
    #> Running under: Ubuntu 20.04 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
    #> 
    #> locale:
    #>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    #>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    #>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    #>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    #> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    #> 
    #> attached base packages:
    #> [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    #> [8] methods   base     
    #> 
    #> other attached packages:
    #>  [1] ldsep_0.0.0.9005            updog_2.0.1                
    #>  [3] VariantAnnotation_1.34.0    Rsamtools_2.4.0            
    #>  [5] Biostrings_2.56.0           XVector_0.28.0             
    #>  [7] SummarizedExperiment_1.18.1 DelayedArray_0.14.0        
    #>  [9] matrixStats_0.56.0          Biobase_2.48.0             
    #> [11] GenomicRanges_1.40.0        GenomeInfoDb_1.24.2        
    #> [13] IRanges_2.22.2              S4Vectors_0.26.1           
    #> [15] BiocGenerics_0.34.0         vcfR_1.11.0                
    #> [17] devtools_2.3.0              usethis_1.6.1              
    #> [19] corrplot_0.84               doParallel_1.0.15          
    #> [21] iterators_1.0.12            foreach_1.5.0              
    #> [23] RcppArmadillo_0.9.900.1.0   Rcpp_1.0.4.6               
    #> [25] latex2exp_0.4.0             ggthemes_4.2.0             
    #> [27] forcats_0.5.0               stringr_1.4.0              
    #> [29] dplyr_1.0.0                 purrr_0.3.4                
    #> [31] readr_1.3.1                 tidyr_1.1.0                
    #> [33] tibble_3.0.1                ggplot2_3.3.2              
    #> [35] tidyverse_1.3.0             BiocManager_1.30.10        
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] colorspace_1.4-1         ellipsis_0.3.1           rprojroot_1.3-2         
    #>  [4] fs_1.4.1                 rstudioapi_0.11          remotes_2.1.1           
    #>  [7] bit64_0.9-7              AnnotationDbi_1.50.0     fansi_0.4.1             
    #> [10] lubridate_1.7.9          xml2_1.3.2               codetools_0.2-16        
    #> [13] splines_4.0.2            knitr_1.29               pkgload_1.1.0           
    #> [16] jsonlite_1.7.0           broom_0.5.6              cluster_2.1.0           
    #> [19] dbplyr_1.4.4             compiler_4.0.2           httr_1.4.1              
    #> [22] backports_1.1.8          assertthat_0.2.1         Matrix_1.2-18           
    #> [25] cli_2.0.2                htmltools_0.5.0          prettyunits_1.1.1       
    #> [28] tools_4.0.2              gtable_0.3.0             glue_1.4.1              
    #> [31] GenomeInfoDbData_1.2.3   rappdirs_0.3.1           cellranger_1.1.0        
    #> [34] vctrs_0.3.1              ape_5.4                  nlme_3.1-147            
    #> [37] rtracklayer_1.48.0       pinfsc50_1.2.0           xfun_0.15               
    #> [40] ps_1.3.3                 testthat_2.3.2           rvest_0.3.5             
    #> [43] lifecycle_0.2.0          XML_3.99-0.3             MASS_7.3-51.6           
    #> [46] zlibbioc_1.34.0          scales_1.1.1             BSgenome_1.56.0         
    #> [49] hms_0.5.3                curl_4.3                 yaml_2.2.1              
    #> [52] memoise_1.1.0            biomaRt_2.44.1           RSQLite_2.2.0           
    #> [55] stringi_1.4.6            desc_1.2.0               permute_0.9-5           
    #> [58] GenomicFeatures_1.40.0   BiocParallel_1.22.0      pkgbuild_1.0.8          
    #> [61] rlang_0.4.6              pkgconfig_2.0.3          bitops_1.0-6            
    #> [64] evaluate_0.14            lattice_0.20-41          GenomicAlignments_1.24.0
    #> [67] bit_1.1-15.2             processx_3.4.2           tidyselect_1.1.0        
    #> [70] magrittr_1.5             R6_2.4.1                 generics_0.0.2          
    #> [73] DBI_1.1.0                pillar_1.4.4             haven_2.3.1             
    #> [76] withr_2.2.0              mgcv_1.8-31              RCurl_1.98-1.2          
    #> [79] modelr_0.1.8             crayon_1.3.4             BiocFileCache_1.12.0    
    #> [82] rmarkdown_2.3            progress_1.2.2           grid_4.0.2              
    #> [85] readxl_1.3.1             blob_1.2.1               callr_3.4.3             
    #> [88] vegan_2.5-6              reprex_0.3.0             digest_0.6.25           
    #> [91] openssl_1.4.1            munsell_0.5.0            viridisLite_0.3.0       
    #> [94] askpass_1.1              sessioninfo_1.1.1

## References

<div id="refs" class="references">

<div id="ref-fox2019ngsld">

Fox, Emma A, Alison E Wright, Matteo Fumagalli, and Filipe G Vieira.
2019. “ngsLD: evaluating linkage disequilibrium using genotype
likelihoods.” *Bioinformatics* 35 (19): 3855–6.
<https://doi.org/10.1093/bioinformatics/btz200>.

</div>

<div id="ref-gerard2020pairwise">

Gerard, David. 2020. “Pairwise Linkage Disequilibrium Estimation for
Polyploids.” *Unpublished Manuscript*.

</div>

<div id="ref-mcallister2016single">

McAllister, Christine A., and Allison J. Miller. 2016. “Single
Nucleotide Polymorphism Discovery via Genotyping by Sequencing to Assess
Population Genetic Structure and Recurrent Polyploidization in
*Andropogon Gerardii*.” *American Journal of Botany* 103 (7): 1314–25.
<https://doi.org/10.3732/ajb.1600146>.

</div>

<div id="ref-uitdewilligen2013next">

Uitdewilligen, Anne-Marie A. AND D’hoop, Jan G. A. M. L. AND Wolters.
2013. “A Next-Generation Sequencing Method for Genotyping-by-Sequencing
of Highly Heterozygous Autotetraploid Potato.” *PLOS ONE* 8 (5). Public
Library of Science: 1–14.
<https://doi.org/10.1371/journal.pone.0062355>.

</div>

</div>
