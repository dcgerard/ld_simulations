
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduce the Results of Gerard (2020)

## Introduction

This repository contains the code and instructions needed to reproduce
all of the results from Gerard (2020).

If you find a bug, please create an
[issue](https://github.com/dcgerard/reproduce_pairwise_ld/issues).

## Instructions

You’ll need to have up-to-date versions of R and python3 to run these
simulatoins.

1.  Install the appropriate R packages
    
    ``` r
    install.packages(c("tidyverse",
                       "ggthemes",
                       "latex2exp",
                       "Rcpp",
                       "RcppArmadillo",
                       "foreach",
                       "doParallel",
                       "corrplot",
                       "devtools"))
    devtools::install_github("dcgerard/updog")
    devtools::install_github("dcgerard/ldsep")
    ```

2.  Install ngsLD in the working directory of the simulation files:
    <https://github.com/fgvieira/ngsLD>
    
    ``` bash
    git clone https://github.com/fgvieira/ngsLD.git
      cd ngsLD
      make
    ```

3.  Install msprime:
    <https://msprime.readthedocs.io/en/stable/installation.html>
    
    I did this via pip:
    
    ``` bash
    pip3 install msprime
    ```
    
    Make sure that you have added `~/.local/bin` to the `PATH` by adding
    the following to your `.bashrc` file:
    
    ``` bash
    export PATH=$PATH:~/.local/bin
    ```

4.  Download the following files from Uitdewilligen et. al. (2013),
    place them in the “data” folder, and extract them.
    
      - <https://doi.org/10.1371/journal.pone.0062355.s004>
      - <https://doi.org/10.1371/journal.pone.0062355.s007>
    
    You can do this with the following bash code
    
    ``` bash
    sudo apt-get install p7zip p7zip-full p7zip-rar
    wget --directory-prefix=data -nc https://doi.org/10.1371/journal.pone.0062355.s004
    mv ./data/journal.pone.0062355.s004 ./data/journal.pone.0062355.s004.gz
    wget --directory-prefix=data -nc https://doi.org/10.1371/journal.pone.0062355.s007
    mv ./data/journal.pone.0062355.s007 ./data/journal.pone.0062355.s007.gz
    7z e ./data/journal.pone.0062355.s004.gz -o./data/
    7z e ./data/journal.pone.0062355.s007.gz -o./data/
    ```
    
    As long as you have 7z installed, `make` will attempt to download
    and extract the files for you.

5.  Run `make`. This will run all of the simulations.
    
    ``` bash
    make
    ```

You may choose to run only part of hte simulations `make mle` or `make
ngsLD`.

5.  Get coffee/sweets. Running `make sims` should take a few hours. You
    should get some coffee\! Here is a list of some of my favorite
    places:
    
      - Washington, DC
          - [Colony
            Club](https://www.yelp.com/biz/colony-club-washington)
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

## References

  - Uitdewilligen, Jan GAML, Anne-Marie A. Wolters, Bjorn B. D’hoop,
    Theo JA Borm, Richard GF Visser, and Herman J. van Eck. “A
    next-generation sequencing method for genotyping-by-sequencing of
    highly heterozygous autotetraploid potato.” *PloS one* 8, no. 5
    (2013).

<div id="refs" class="references">

<div id="ref-gerard2020pairwise">

Gerard, David. 2020. “Pairwise Linkage Disequilibrium Shrinkage
Estimation for Autopolyploids.”

</div>

</div>
