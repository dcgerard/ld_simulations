
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduce the Results of Gerard (2021)

[![DOI](https://img.shields.io/badge/doi-10.6084/m9.figshare.12765803-blue.svg?style=flat&labelColor=whitesmoke&logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAAB8AAAAfCAYAAAAfrhY5AAAJsklEQVR42qWXd1DTaRrHf%2BiB2Hdt5zhrAUKz4IKEYu9IGiGFFJJQ0gkJCAKiWFDWBRdFhCQUF3UVdeVcRQEBxUI3yY9iEnQHb3bdW1fPubnyz%2F11M7lvEHfOQee2ZOYzPyDv%2B3yf9%2Fk95YX4fx%2BltfUt08GcFEuPR4U9hDDZ%2FVngIlhb%2FSiI6InkTgLzgDcgfvtnovhH4BzoVlrbwr55QnhCtBW4QHXnFrZbPBaQoBh4%2FSYH2EnpBEtqcDMVzB93wA%2F8AFwa23XFGcc8CkT3mxz%2BfXWtq9T9IQlLIXYEuHojudb%2BCM7Hgdq8ydi%2FAHiBXyY%2BLjwFlAEnS6Jnar%2FvnQVhvdzasad0eKvWZKe8hvDB2ofLZ%2FZEcWsh%2BhyIuyO5Bxs2iZIE4nRv7NWAb0EO8AC%2FWPxjYAWuOEX2MSXZVgPxzmRL3xKz3ScGpx6p6QnOx4mDIFqO0w6Q4fEhO5IzwxlSwyD2FYHzwAW%2BAZ4fEsf74gCumykwNHskLM7taQxLYjjIyy8MUtraGhTWdkfhkFJqtvuVl%2F9l2ZquDfEyrH8B0W06nnpH3JtIyRGpH1iJ6SfxDIHjRXHJmdQjLpfHeN54gnfFx4W9QRnovx%2FN20aXZeTD2J84hn3%2BqoF2Tqr14VqTPUCIcP%2B5%2Fly4qC%2BUL3sYxSvNj1NwsVYPsWdMUfomsdkYm3Tj0nbV0N1wRKwFe1MgKACDIBdMAhPE%2FwicwNWxll8Ag40w%2BFfhibJkGHmutjYeQ8gVlaN%2BjO51nDysa9TwNUFMqaGbKdRJZFfOJSp6mkRKsv0rRIpEVWjAvyFkxNOEpwvcAVPfEe%2Bl8ojeNTx3nXLBcWRrYGxSRjDEk0VlpxYrbe1ZmaQ5xuT0u3r%2B2qe5j0J5uytiZPGsRL2Jm32AldpxPUNJ3jmmsN4x62z1cXrbedXBQf2yvIFCeZrtyicZZG2U2nrrBJzYorI2EXLrvTfCSB43s41PKEvbZDEfQby6L4JTj%2FfIwam%2B4%2BwucBu%2BDgNK05Nle1rSt9HvR%2FKPC4U6LTfvUIaip1mjIa8fPzykii23h2eanT57zQ7fsyYH5QjywwlooAUcAdOh5QumgTHx6aAO7%2FL52eaQNEShrxfhL6albEDmfhGflrsT4tps8gTHNOJbeDeBlt0WJWDHSgxs6cW6lQqyg1FpD5ZVDfhn1HYFF1y4Eiaqa18pQf3zzYMBhcanlBjYfgWNayAf%2FASOgklu8bmgD7hADrk4cRlOL7NSOewEcbqSmaivT33QuFdHXj5sdvjlN5yMDrAECmdgDWG2L8P%2BAKLs9ZLZ7dJda%2BB4Xl84t7QvnKfvpXJv9obz2KgK8dXyqISyV0sXGZ0U47hOA%2FAiigbEMECJxC9aoKp86re5O5prxOlHkcksutSQJzxZRlPZmrOKhsQBF5zEZKybUC0vVjG8PqOnhOq46qyDTDnj5gZBriWCk4DvXrudQnXQmnXblebhAC2cCB6zIbM4PYgGl0elPSgIf3iFEA21aLdHYLHUQuVkpgi02SxFdrG862Y8ymYGMvXDzUmiX8DS5vKZyZlGmsSgQqfLub5RyLNS4zfDiZc9Edzh%2FtCE%2BX8j9k%2FqWB071rcZyMImne1SLkL4GRw4UPHMV3jjwEYpPG5uW5fAEot0aTSJnsGAwHJi2nvF1Y5OIqWziVCQd5NT7t6Q8guOSpgS%2Fa1dSRn8JGGaCD3BPXDyQRG4Bqhu8XrgAp0yy8DMSvvyVXDgJcJTcr1wQ2BvFKf65jqhvmxXUuDpGBlRvV36XvGjQzLi8KAKT2lYOnmxQPGorURSV0NhyTIuIyqOmKTMhQ%2BieEsgOgpc4KBbfDM4B3SIgFljvfHF6cef7qpyLBXAiQcXvg5l3Iunp%2FWv4dH6qFziO%2BL9PbrimQ9RY6MQphEfGUpOmma7KkGzuS8sPUFnCtIYcKCaI9EXo4HlQLgGrBjbiK5EqMj2AKWt9QWcIFMtnVvQVDQV9lXJJqdPVtUQpbh6gCI2Ov1nvZts7yYdsnvRgxiWFOtNJcOMVLn1vgptVi6qrNiFOfEjHCDB3J%2BHDLqUB77YgQGwX%2Fb1eYna3hGKdlqJKIyiE4nSbV8VFgxmxR4b5mVkkeUhMgs5YTi4ja2XZ009xJRHdkfwMi%2BfocaancuO7h%2FMlcLOa0V%2FSw6Dq47CumRQAKhgbOP8t%2BMTjuxjJGhXCY6XpmDDFqWlVYbQ1aDJ5Cptdw4oLbf3Ck%2BdWkVP0LpH7s9XLPXI%2FQX8ws%2Bj2In63IcRvOOo%2BTTjiN%2BlssfRsanW%2B3REVKoavBOAPTXABW4AL7e4NygHdpAKBscmlDh9Jysp4wxbnUNna3L3xBvyE1jyrGIkUHaqQMuxhHElV6oj1picvgL1QEuS5PyZTEaivqh5vUCKJqOuIgPFGESns8kyFk7%2FDxyima3cYxi%2FYOQCj%2F%2B9Ms2Ll%2Bhn4FmKnl7JkGXQGDKDAz9rUGL1TIlBpuJr9Be2JjK6qPzyDg495UxXYF7JY1qKimw9jWjF0iV6DRIqE%2B%2FeWG0J2ofmZTk0mLYVd4GLiFCOoKR0Cg727tWq981InYynvCuKW43aXgEjofVbxIqrm0VL76zlH3gQzWP3R3Bv9oXxclrlO7VVtgBRpSP4hMFWJ8BrUSBCJXC07l40X4jWuvtc42ofNCxtlX2JH6bdeojXgTh5TxOBKEyY5wvBE%2BACh8BtOPNPkApjoxi5h%2B%2FFMQQNpWvZaMH7MKFu5Ax8HoCQdmGkJrtnOiLHwD3uS5y8%2F2xTSDrE%2F4PT1yqtt6vGe8ldMBVMEPd6KwqiYECHDlfbvzphcWP%2BJiZuL5swoWQYlS%2Br7Yu5mNUiGD2retxBi9fl6RDGn4Ti9B1oyYy%2BMP5G87D%2FCpRlvdnuy0PY6RC8BzTA40NXqckQ9TaOUDywkYsudxJzPgyDoAWn%2BB6nEFbaVxxC6UXjJiuDkW9TWq7uRBOJocky9iMfUhGpv%2FdQuVVIuGjYqACbXf8aa%2BPeYNIHZsM7l4s5gAQuUAzRUoT51hnH3EWofXf2vkD5HJJ33vwE%2FaEWp36GHr6GpMaH4AAPuqM5eabH%2FhfG9zcCz4nN6cPinuAw6IHwtvyB%2FdO1toZciBaPh25U0ducR2PI3Zl7mokyLWKkSnEDOg1x5fCsJE9EKhH7HwFNhWMGMS7%2BqxyYsbHHRUDUH4I%2FAheQY7wujJNnFUH4KdCju83riuQeHU9WEqNzjsJFuF%2FdTDAZ%2FK7%2F1WaAU%2BAWymT59pVMT4g2AxcwNa0XEBDdBDpAPvgDIH73R25teeuAF5ime2Ul0OUIiG4GpSAEJeYW9wDTf43wfwHgHLKJoPznkwAAAABJRU5ErkJggg%3D%3D)](https://doi.org/10.6084/m9.figshare.12765803)

## Introduction

This repository contains the code and instructions needed to reproduce
all of the results from Gerard (2021).

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
[ngsLD](https://github.com/fgvieira/ngsLD),
[PedigreeSim](https://github.com/PBR/pedigreeSim), and some of the R
packages I used. Check their documentation if you are having trouble
installing them.

Note that I have only tried these out on R version 4.0.5 (2021-03-31)
running Ubuntu 20.04.2 LTS.

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

    -   <https://doi.org/10.1371/journal.pone.0062355.s004>
    -   <https://doi.org/10.1371/journal.pone.0062355.s007>
    -   <https://doi.org/10.1371/journal.pone.0062355.s009>
    -   <https://doi.org/10.1371/journal.pone.0062355.s010>

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

    -   `make mle`: Simulations under HWE.
    -   `make comp`: Simulations under violations from HWE.
    -   `make norm`: Visualizing the proportional bivariate normal
        distribution.
    -   `make ddiff`: Comparing, under HWE, D’ and the Δ’ that
        conditions on the marginal genotype distributions.
    -   `make ngsLD`: Verifying that
        [`ldsep`](https://cran.r-project.org/package=ldsep) and
        [`ngsLD`](https://github.com/fgvieira/ngsLD) (Fox et al. 2019)
        provide the same results in diploids.
    -   `make uit`: Real-data analysis using the data from
        Uitdewilligen (2013).
    -   `make mca`: Real-data analysis using the data from McAllister
        and Miller (2016).
    -   `make ped`: Simulations under *interpretable* violations from
        HWE, where data were simulated using
        [PedigreeSim](https://github.com/PBR/pedigreeSim) (Voorrips and
        Maliepaard 2012).

7.  Get coffee/sweets. Running `make sims` should take a few hours. You
    should get some coffee! Here is a list of some of my favorite
    places:

    -   Washington, DC
        -   [Doubles](https://www.yelp.com/biz/doubles-washington)
        -   [Little Red
            Fox](https://www.yelp.com/biz/little-red-fox-washington)
        -   [Three
            Fifty](https://www.yelp.com/biz/three-fifty-bakery-and-coffee-bar-washington)
            (for the baked goods, coffee is meh)
    -   Chicago
        -   [Sawada
            Coffee](https://www.yelp.com/biz/sawada-coffee-chicago)
        -   [Plein Air
            Cafe](https://www.yelp.com/biz/plein-air-cafe-and-eatery-chicago-2)
    -   Seattle
        -   [Bauhaus
            Ballard](https://www.yelp.com/biz/bauhaus-ballard-seattle)
        -   [Cafe
            Solstice](https://www.yelp.com/biz/cafe-solstice-seattle)
    -   Columbus
        -   [Yeah, Me
            Too](https://www.yelp.com/biz/yeah-me-too-columbus)
        -   [Stauf’s Coffee
            Roasters](https://www.yelp.com/biz/staufs-coffee-roasters-columbus-2)

# Session Information

    #> R version 4.0.5 (2021-03-31)
    #> Platform: x86_64-pc-linux-gnu (64-bit)
    #> Running under: Ubuntu 20.04.2 LTS
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
    #>  [1] ldsep_2.0.3                 updog_2.1.0                
    #>  [3] VariantAnnotation_1.36.0    Rsamtools_2.6.0            
    #>  [5] Biostrings_2.58.0           XVector_0.30.0             
    #>  [7] SummarizedExperiment_1.20.0 Biobase_2.50.0             
    #>  [9] GenomicRanges_1.42.0        GenomeInfoDb_1.26.4        
    #> [11] IRanges_2.24.1              S4Vectors_0.28.1           
    #> [13] MatrixGenerics_1.2.1        matrixStats_0.58.0         
    #> [15] BiocGenerics_0.36.0         vcfR_1.12.0                
    #> [17] devtools_2.3.2              usethis_2.0.1              
    #> [19] corrplot_0.84               doParallel_1.0.16          
    #> [21] iterators_1.0.13            foreach_1.5.1              
    #> [23] RcppArmadillo_0.10.2.2.0    Rcpp_1.0.6                 
    #> [25] latex2exp_0.5.0             gridExtra_2.3              
    #> [27] ggthemes_4.2.4              forcats_0.5.1              
    #> [29] stringr_1.4.0               dplyr_1.0.5                
    #> [31] purrr_0.3.4                 readr_1.4.0                
    #> [33] tidyr_1.1.3                 tibble_3.1.0               
    #> [35] ggplot2_3.3.3               tidyverse_1.3.0            
    #> [37] BiocManager_1.30.12        
    #> 
    #> loaded via a namespace (and not attached):
    #>   [1] readxl_1.3.1             backports_1.2.1          BiocFileCache_1.14.0    
    #>   [4] splines_4.0.5            listenv_0.8.0            BiocParallel_1.24.1     
    #>   [7] digest_0.6.27            htmltools_0.5.1.1        fansi_0.4.2             
    #>  [10] magrittr_2.0.1           memoise_2.0.0            BSgenome_1.58.0         
    #>  [13] cluster_2.1.1            remotes_2.3.0            globals_0.14.0          
    #>  [16] modelr_0.1.8             doFuture_0.12.0          askpass_1.1             
    #>  [19] prettyunits_1.1.1        colorspace_2.0-0         blob_1.2.1              
    #>  [22] rvest_1.0.0              rappdirs_0.3.3           haven_2.3.1             
    #>  [25] xfun_0.22                callr_3.6.0              crayon_1.4.1            
    #>  [28] RCurl_1.98-1.3           jsonlite_1.7.2           ape_5.4-1               
    #>  [31] glue_1.4.2               gtable_0.3.0             zlibbioc_1.36.0         
    #>  [34] DelayedArray_0.16.3      pkgbuild_1.2.0           scales_1.1.1            
    #>  [37] rngtools_1.5             DBI_1.1.1                viridisLite_0.3.0       
    #>  [40] progress_1.2.2           bit_4.0.4                httr_1.4.2              
    #>  [43] ellipsis_0.3.1           pkgconfig_2.0.3          XML_3.99-0.6            
    #>  [46] dbplyr_2.1.0             utf8_1.2.1               tidyselect_1.1.0        
    #>  [49] rlang_0.4.10             AnnotationDbi_1.52.0     munsell_0.5.0           
    #>  [52] cellranger_1.1.0         tools_4.0.5              cachem_1.0.4            
    #>  [55] cli_2.4.0                generics_0.1.0           RSQLite_2.2.5           
    #>  [58] broom_0.7.6              evaluate_0.14            fastmap_1.1.0           
    #>  [61] yaml_2.2.1               processx_3.5.1           knitr_1.31              
    #>  [64] bit64_4.0.5              fs_1.5.0                 doRNG_1.8.2             
    #>  [67] future_1.21.0            nlme_3.1-152             xml2_1.3.2              
    #>  [70] biomaRt_2.46.3           compiler_4.0.5           rstudioapi_0.13         
    #>  [73] curl_4.3                 testthat_3.0.2           reprex_2.0.0            
    #>  [76] stringi_1.5.3            ps_1.6.0                 GenomicFeatures_1.42.3  
    #>  [79] desc_1.3.0               lattice_0.20-41          Matrix_1.3-2            
    #>  [82] vegan_2.5-7              permute_0.9-5            vctrs_0.3.7             
    #>  [85] pillar_1.5.1             lifecycle_1.0.0          bitops_1.0-6            
    #>  [88] rtracklayer_1.50.0       R6_2.5.0                 parallelly_1.24.0       
    #>  [91] sessioninfo_1.1.1        codetools_0.2-18         MASS_7.3-53.1           
    #>  [94] assertthat_0.2.1         pkgload_1.2.0            openssl_1.4.3           
    #>  [97] rprojroot_2.0.2          withr_2.4.1              pinfsc50_1.2.0          
    #> [100] GenomicAlignments_1.26.0 GenomeInfoDbData_1.2.4   mgcv_1.8-33             
    #> [103] hms_1.0.0                grid_4.0.5               rmarkdown_2.7           
    #> [106] lubridate_1.7.10

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-fox2019ngsld" class="csl-entry">

Fox, Emma A, Alison E Wright, Matteo Fumagalli, and Filipe G Vieira.
2019. “<span class="nocase"><span class="nocase">ngsLD</span>:
evaluating linkage disequilibrium using genotype likelihoods</span>.”
*Bioinformatics* 35 (19): 3855–56.
<https://doi.org/10.1093/bioinformatics/btz200>.

</div>

<div id="ref-gerard2021pairwise" class="csl-entry">

Gerard, David. 2021. “Pairwise Linkage Disequilibrium Estimation for
Polyploids.” *Molecular Ecology Resources* 21 (4): 1230–42.
<https://doi.org/10.1111/1755-0998.13349>.

</div>

<div id="ref-mcallister2016single" class="csl-entry">

McAllister, Christine A., and Allison J. Miller. 2016. “Single
Nucleotide Polymorphism Discovery via Genotyping by Sequencing to Assess
Population Genetic Structure and Recurrent Polyploidization in
*Andropogon Gerardii*.” *American Journal of Botany* 103 (7): 1314–25.
<https://doi.org/10.3732/ajb.1600146>.

</div>

<div id="ref-uitdewilligen2013next" class="csl-entry">

Uitdewilligen, Anne-Marie A. AND D’hoop, Jan G. A. M. L. AND Wolters.
2013. “A Next-Generation Sequencing Method for Genotyping-by-Sequencing
of Highly Heterozygous Autotetraploid Potato.” *PLOS ONE* 8 (5): 1–14.
<https://doi.org/10.1371/journal.pone.0062355>.

</div>

<div id="ref-voorrips2012simulation" class="csl-entry">

Voorrips, Roeland E., and Chris A. Maliepaard. 2012. “The Simulation of
Meiosis in Diploid and Tetraploid Organisms Using Various Genetic
Models.” *BMC Bioinformatics* 13 (1): 248.
<https://doi.org/10.1186/1471-2105-13-248>.

</div>

</div>
