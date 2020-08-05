
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduce the Results of Gerard (2020)

[![DOI](https://img.shields.io/badge/doi-https://doi.org/10.6084/m9.figshare.12765803.v1-blue.svg?style=flat&labelColor=gainsboro&logoWidth=40&logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAAFAAAAAZCAYAAACmRqkJAAAKi0lEQVR4Ae3ZaVBUV97H8evuE0EfH32MmkcfoyAuGjXKgkvMaFRAFuiloemWvRuEXlgEBREXBYJiXAQUFeKocUniQiKogAJhQWwWENDEjLNYvjFLzUzNkplEZb5kTme6nCRjKlOpSZlb9SmL2%2Ffcuv3re87%2FnKP0TYfOcslqPMbt63xBKuh09MTxgi7HKT1Sj1TvKp%2BMkZB6%2FXT8c4AjUYPyVdfb7Qs6HTIJ8EHe7Ul%2B152CphDabRQ0uMr7%2FRQgh%2B8qU6%2FBiPDVGv0jq0uGE94b0ZZ3j%2B25MTetoMsh%2FWD91OBqT9%2Fsehd5EqGV17nKMzTqOHvaRMMLEp7qACfinq%2FW1BBx5ZxB13x5X3Jr1v%2Fz9pUcaHU63PiicjrhvXfNRbY1Th49Q6Y1vu6zyqSjzX3aVIgf4OkKToxhgxpd5OMzV0bYE4CRN1Chu34pnTfwnV03FiTlfzDRXBHo6dfgIq8sX6ByV6vjthGc0UdrrPPVGFQBxlSjzJQWENVUZkebceiLpyM8IZSx7O7Zl4JivUNMZX5h8Rt4%2B2L0llKfgu6JKa%2BXvpB5bZ48%2Ba3F6lil2pDkE2rODzCsU0VUnNFHNZQqdS3lx3Utl%2FMILQcfYt5TEeC1GSprgAq0XlgYGLQyxJTlr0uK0DVX7E5s2ZtOgHvLw5fLK9xVmcqguEj%2F2LXbwsvPBkZZKl4j5NcIKinaUsLbejFWZ7m8Do2cmwnb4cFqArRwx3TEYzi%2Bz7DTD0uhxnj8cAEWWUZK%2BTcdhh4pmTWUsW01Y1uCUmNY7Rtqzo5svJSS0poVXtg6yVj7sn9qunek3j8xPVXXeMFoaDkev6lDF7ene7Y5r2taNAXmEBXaP69zevaOjuUeeZ0zhzJuPsM5CdYvOhZVqBMhBqIVDt8zwGdQjR4of9AA%2BXJjUFpww7GodnHAQca4srDAWCXjW3pETal%2BbfumuOLKqSm17vIQtWr1Uu3JYy6JbXuXFbRN1R8pm5byxtG5CcdOz9EUVc7I5IeQEWQ7wWVwzwrsRn%2BbAFeiCxNsKv5Y9P03BFgjAlT90AGOQy2T47fObl00ocFZHl%2B2UGXw0RjzNUWHTPFthckHWh18al8KsGuaFigVVzlKuY%2BG9z37qvuoGlelpsJVldrgrFjbOE%2BeWe8uW18W84qCqc4s7tmCIgzI75hs%2FaJKNFu7rF%2BIIIhr%2BmIQ%2Btn8LQkDMQOeWAYnDHgsQI3NNU7W9j4h5t72o%2FEyvLEQ%2F%2Bu7ymzbOxbCAeOxAgtghz6YgOVYiufEOUlqu0M37ho%2BYn%2FnpJT8bsejVSt90uqdFdlGmV7hF7cuWXetNCShLX%2BI3nKhN%2ByvCs%2Bs6GQpWB33fzKNQR%2BqWr022yvc94q7spBCY%2Bbzkou6ZfJNPf89ZN%2FdidYHnIsKfIzjCMIc7MAwSJiMPFxGMcKQixGwx07R%2FiEe4CNsxFCbAJvwifj8LkIgYRHa8Lm47jNY8AokmMS5NryPh%2FijOB%2BOX4h7foEuyPHlisMtylJpzu1YspkQ36YbLqnx8F1X4abaqmYs9DGmLlrk4CE9XlHlKZskxfpt%2FUJLzyhV23dG%2BITF72fqo9njEaokwIu8lSbG1N4wx273CrP%2B%2BjniQVZhGrzQjlEioFIRcjDM6MIdjBVtHogvl4W9qIX8sTfwU5SgU%2FzdhdGYLcJ9BzvRID6vgx2SxN8PUI9KnIEWH4n7FuIo%2FoRfYV5vMMV4wHRFs%2BvG%2FKl05ZrDVdP11T7eulK3oNQcz%2FAXcj3DpMePjO44KetDL2lDh%2FmV1S3nNoeWnJb7RSXmMJl%2BI0GmH13rKs8lvEdQwfoWKmCxdmGbAEdgAW5jFiQhBb8WXSYTPSjGCBHaMPR5LMANkOCM%2B%2FgD3MS5Z8W1ElzwW3HNJCSI9tcw2ub%2BO8T5LPTBQBy1nusNcB7ztximI1sIsSSzXb04v3vyusJmx63nMufHXlV6LvpEShDd9x%2FHFYWXVPuSX7%2FD7zmpcjuWRupbyvaHnj8Z7BNsUFCArm70iTRcd5bFEN7oxwJs%2FpoA%2FwfBaLJ2Z2EFbmEsNKL7fYYPUI9DIqj%2Fsgkw0CasW%2BL6RbBDFI7gTZSKzz6Gk02AJ23G3QF4xybYU8INce6s5CJNlTyXhYwKv%2FRWMiEeimquzIhrPpGzuSNCsbvLec2%2Brpmh2e0yu%2FxOp96wv6p8X0xeIZW5Bo2%2F6ucdvb%2FdMWVDm8lX11pRpD16OJ6VyZsrQ8yK%2BVFJ9h4UhwEHDj5JgGE23UkSfoZujMMzSESNCPBT9KAFjqi2rcIYZRPgYmzDQ9xDLSz4%2FGsCPIE%2BNkWrTJy%2FhRrRthpVyJJExbnmG2I%2B6x%2BT%2FHxYyQkzQfJGlufpWy6bYlvPUEgu%2BHlHJA5boo7rE3blnBR7r6mv%2BvCBMYEag%2Faqsyr1%2BIk5a%2Fd2z9zGBDpZ31qulCWk9443Hfg5BuJJAgxAG0ZBEmS4DZ7RKIliMVi0d8UvRUCeuPoNAf4Z%2FmgV13pAwiwR3iffFKBQJM5noB%2F6Y5h45v7Wwf0cDtD1DlMIeiugWmZOy5Cv3RgjX7%2FF4GdMXasOjgurmqdafqpojltml9IjvOJ8NMu9lNL5gQmXdMu0BTefz8loMyoJvivs3VMZvhpjqaig%2FZ8gwJGYIsIKRh%2FY4wh%2Bg%2FGQoxYbREgZ%2BB3uww1V3xKgN%2BrwCNtF4Pvx8NveQCEYX%2BAukhCIYuHZLy%2FyDjHbJQfo7PTK1dEBWqPBX2vS%2B2hNW1XquDURypiwXStCjVWuyrSKQC%2FdoUaHtOT2HENoyal4b40x7rK7ylip9NIV3Jy0P6fD24fl3Ra6uoe3PNqOH2Pw3x%2FC8K8CHIU%2BIpQ7OI8yNOJ9TMJO%2FAU9Nn6PjRiGmm%2FpwgsRLQpKjwjuU%2Fz1CQK0R4G4T4%2FwCHWYKlmcA6xr4SA2EzobXeUa9vh21LgpdKxK8hqd5RsaXWS7S9YvlhU2O7ya3ekXrm%2B9lK3KzFH6a4y5V92Ve5hkM4d02EShMestZekE2IxZX7MWdkAgBtmsi9U2lXEwliAOK%2BGLTowThWIZkrEVSSKYgegPOUxwtFmdaBGLsRgg2qeKtosQDh2GYzbisUIEaPvcQ8T5VGzCKowBk2I3mTVALe4wd4tumKcoaZirSKte4RtVrvXwLrw%2BJXV%2F18Ts3BtLEmOaS0yRtRdMfpGJhTKNMbDJWR5V7eEbUNDtcIQAd1PJMwnuJl6E9KQHY7AAHkzQoBkj8B%2B%2FpTWQ4Maezne1P3x1esLBuqmB%2BbccNhJMGetbM%2BGZIi1V%2FoRyOXB77sKVWuPmrd4RBvYQm9ihVue%2F7xDPGljB50MoJmO%2By36gCGsQovCyCGwOarD9R7PLLXZOJjKZvse%2FDQQSvffG7F1rWrZPiLKUX2DPr1hbfHAKb0kDBSeTed5MQj94Pn1xBMvA%2B2IDYTAkcXzXANPRjHq04ACeFeH9aAIcBC3LOq%2FY5pPDeYtO4yRTmzUhbx9LozCEea8ybaHoxDNmVtPltxSVzxhCm3Asg4Tvs683Aa5wwkD8qP9XbgQqUbb6Tp09U5Os3rWiV4jZv2OuvxPdvht70RfST8fjATZd7P33OYzxZ%2FdF7FwcgqPU0yMR2vMYDulpDfBvw%2BGCdBePpq8AAAAASUVORK5CYII%3D)](https://doi.org/10.6084/m9.figshare.12765803.v1)

## Introduction

This repository contains the code and instructions needed to reproduce
all of the results from Gerard (2020).

Some of the plots from the simulations mentioned in the manuscript can
be found
[here](https://github.com/dcgerard/ld_simulations/tree/master/output/mle_plots)
and
[here](https://github.com/dcgerard/ld_simulations/tree/master/output/comp/comp_plots).

If you find a bug, please create an
[issue](https://github.com/dcgerard/ld_simulations/issues).

## Instructions

You will need to have R and GNU Make installed to run these simulations.
There are other dependencies for
[ngsLD](https://github.com/fgvieira/ngsLD) and some of the R packages I
used. Check their documentation if you are having trouble installing
them.

Note that I have only tried these out on R 4.0.2 using Ubuntu 20.04.

1.  Install the appropriate R packages
    
    ``` r
    install.packages("BiocManager")
    BiocManager::install(c("tidyverse",
                           "ggthemes",
                           "gridExtra",
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
    of the ld\_simulations repo.
    
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
    <https://doi.org/10.5061/dryad.05qs7> and unzip it to the
    ./data/gerardii directory. I unzipped using 7z:
    
    ``` bash
    7z e doi_10.5061_dryad.05qs7__v1.zip 
    ```

5.  Your data directory should now look like this:
    
        ./data/CSV-file S1 - Sequence variants filtered DP15.csv
        ./data/gerardii/doi_10.5061_dryad.05qs7__v1.zip
        ./data/gerardii/McAllister_et_al_2017_Data_from _Single_nucleotide_polymorphism.pdf
        ./data/gerardii/McAllister_Miller_Locality_Ploidy_Info.csv
        ./data/gerardii/McAllister.Miller.all.mergedRefGuidedSNPs.vcf.gz
        ./data/gerardii/McAllister.Miller.all.mergedUNEAKSNPs.vcf.gz
        ./data/journal.pone.0062355.s004.gz
        ./data/journal.pone.0062355.s007.gz
        ./data/journal.pone.0062355.s009.xls
        ./data/journal.pone.0062355.s010.xls
        ./data/NewPlusOldCalls.headed.vcf
        ./data/Workbook

6.  Run `make`. This will run all of the analyses.
    
    ``` bash
    make
    ```
    
    You may choose to run only part of the analyses via:
    
      - `make mle`: Simulations under HWE.
      - `make comp`: Simulations under violations from HWE.
      - `make norm`: Visualizing the proportional bivariate normal
        distribution.
      - `make ddiff`: Comparing, under HWE, D’ and the Δ’ that
        conditions on the marginal genotype distributions.
      - `make ngsLD`: Verifying that
        [`ldsep`](https://cran.r-project.org/package=ldsep) and
        [`ngsLD`](https://github.com/fgvieira/ngsLD) (Fox et al. 2019)
        provide the same results in diploids.
      - `make uit`: Real-data analysis using the data from Uitdewilligen
        (2013).
      - `make mca`: Real-data analysis using the data from McAllister
        and Miller (2016).

7.  Get coffee/sweets. Running `make sims` should take a few hours. You
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

# Session Information

    #> R version 4.0.2 (2020-06-22)
    #> Platform: x86_64-pc-linux-gnu (64-bit)
    #> Running under: Ubuntu 20.04.1 LTS
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
    #>  [1] ldsep_1.0.0                 updog_2.0.2                
    #>  [3] VariantAnnotation_1.34.0    Rsamtools_2.4.0            
    #>  [5] Biostrings_2.56.0           XVector_0.28.0             
    #>  [7] SummarizedExperiment_1.18.2 DelayedArray_0.14.1        
    #>  [9] matrixStats_0.56.0          Biobase_2.48.0             
    #> [11] GenomicRanges_1.40.0        GenomeInfoDb_1.24.2        
    #> [13] IRanges_2.22.2              S4Vectors_0.26.1           
    #> [15] BiocGenerics_0.34.0         vcfR_1.11.0                
    #> [17] devtools_2.3.1              usethis_1.6.1              
    #> [19] corrplot_0.84               doParallel_1.0.15          
    #> [21] iterators_1.0.12            foreach_1.5.0              
    #> [23] RcppArmadillo_0.9.900.2.0   Rcpp_1.0.5                 
    #> [25] latex2exp_0.4.0             gridExtra_2.3              
    #> [27] ggthemes_4.2.0              forcats_0.5.0              
    #> [29] stringr_1.4.0               dplyr_1.0.1                
    #> [31] purrr_0.3.4                 readr_1.3.1                
    #> [33] tidyr_1.1.1                 tibble_3.0.3               
    #> [35] ggplot2_3.3.2               tidyverse_1.3.0            
    #> [37] BiocManager_1.30.10        
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] colorspace_1.4-1         ellipsis_0.3.1           rprojroot_1.3-2         
    #>  [4] fs_1.5.0                 rstudioapi_0.11          remotes_2.2.0           
    #>  [7] bit64_4.0.2              AnnotationDbi_1.50.3     fansi_0.4.1             
    #> [10] lubridate_1.7.9          xml2_1.3.2               codetools_0.2-16        
    #> [13] splines_4.0.2            knitr_1.29               pkgload_1.1.0           
    #> [16] jsonlite_1.7.0           broom_0.7.0              cluster_2.1.0           
    #> [19] dbplyr_1.4.4             compiler_4.0.2           httr_1.4.2              
    #> [22] backports_1.1.8          assertthat_0.2.1         Matrix_1.2-18           
    #> [25] cli_2.0.2                htmltools_0.5.0          prettyunits_1.1.1       
    #> [28] tools_4.0.2              gtable_0.3.0             glue_1.4.1              
    #> [31] GenomeInfoDbData_1.2.3   rappdirs_0.3.1           cellranger_1.1.0        
    #> [34] vctrs_0.3.2              ape_5.4                  nlme_3.1-147            
    #> [37] rtracklayer_1.48.0       pinfsc50_1.2.0           xfun_0.16               
    #> [40] ps_1.3.3                 testthat_2.3.2           rvest_0.3.6             
    #> [43] lifecycle_0.2.0          XML_3.99-0.5             MASS_7.3-51.6           
    #> [46] zlibbioc_1.34.0          scales_1.1.1             BSgenome_1.56.0         
    #> [49] hms_0.5.3                curl_4.3                 yaml_2.2.1              
    #> [52] memoise_1.1.0            biomaRt_2.44.1           RSQLite_2.2.0           
    #> [55] stringi_1.4.6            desc_1.2.0               permute_0.9-5           
    #> [58] GenomicFeatures_1.40.1   BiocParallel_1.22.0      pkgbuild_1.1.0          
    #> [61] rlang_0.4.7              pkgconfig_2.0.3          bitops_1.0-6            
    #> [64] evaluate_0.14            lattice_0.20-41          GenomicAlignments_1.24.0
    #> [67] bit_4.0.3                processx_3.4.3           tidyselect_1.1.0        
    #> [70] magrittr_1.5             R6_2.4.1                 generics_0.0.2          
    #> [73] DBI_1.1.0                pillar_1.4.6             haven_2.3.1             
    #> [76] withr_2.2.0              mgcv_1.8-31              RCurl_1.98-1.2          
    #> [79] modelr_0.1.8             crayon_1.3.4             BiocFileCache_1.12.0    
    #> [82] rmarkdown_2.3            progress_1.2.2           grid_4.0.2              
    #> [85] readxl_1.3.1             blob_1.2.1               callr_3.4.3             
    #> [88] vegan_2.5-6              reprex_0.3.0             digest_0.6.25           
    #> [91] openssl_1.4.2            munsell_0.5.0            viridisLite_0.3.0       
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
Polyploids.” *bioRxiv*. <https://doi.org/10.1101/2020.08.03.234476>.

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
of Highly Heterozygous Autotetraploid Potato.” *PLOS ONE* 8 (5): 1–14.
<https://doi.org/10.1371/journal.pone.0062355>.

</div>

</div>
