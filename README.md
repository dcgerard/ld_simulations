
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
                       "devtools",
                       "vcfR"))
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

You may choose to run only part of the simulations `make mle` or `make
ngsLD`.

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
