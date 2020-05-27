
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduce the Results of Gerard (2020)

## Introduction

This repository contains the code and instructions needed to reproduce
all of the results from Gerard (2020).

If you find a bug, please create an
[issue](https://github.com/dcgerard/reproduce_pairwise_ld/issues).

## Instructions

1.  Install the appropriate R packages
    
    ``` r
    install.packages(c("updog",
                       "tidyverse", 
                       "Rcpp",
                       "RcppArmadillo",
                       "roptim",
                       "foreach",
                       "doParallel",
                       "devtools"))
    devtools::install_github("dcgerard/ldsep")
    devtools::install_github("tpbilton/GUSbase")
    devtools::install_github("AgResearch/GUS-LD")
    ```

2.  Get coffee. Running `make sims` should take a few hours. You should
    get some coffee\! Here is a list of some of my favorite places:
    
      - Washington, DC
          - [Colony
            Club](https://www.yelp.com/biz/colony-club-washington)
          - [Grace Street
            Coffee](https://www.yelp.com/biz/grace-street-coffee-georgetown)
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

## References

<div id="refs" class="references">

<div id="ref-gerard2020pairwise">

Gerard, David. 2020. “Pairwise Linkage Disequilibrium Shrinkage
Estimation for Polyploids.”

</div>

</div>
