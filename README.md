glmmBoot
================

`glmmBoot` consists of a series of tools that for fitting and
bootstrapping generalized linear mixed models (GLMMs) that takes
advantage of the Template Model Builder (`TMB`), built on`CppAD` and
`Eigen`.

At this stage of development it can handle Mixed Logit and Mixed Poisson
Models, yet it has the potential of being extended intended to handle a
wide range of statistical distributions.

# Mixed Poisson Example

This code to illustrates the implementation of the bootstrapping methods
described in [(Flores-Agreda and Cantoni, Under
Review)](https://www.researchgate.net/publication/315768128_Bootstrapping_Generalized_Linear_Mixed_Models_via_a_Weighted_Laplace_Approximation)
in a real-data context. Data comes from the  example found in
[Molenberghs & Verbeke
(2006)](https://www.springer.com/gp/book/9780387251448) initially
reported by [Faught et. al.
(1996)](https://www.ncbi.nlm.nih.gov/pubmed/8649570).

## Analysis

To start performing the analysis, you need:

1.  the following packages:

<!-- end list -->

``` r
library(lme4)
library(glmmTMB)
library(sas7bdat)
library(TMB)
library(dplyr)
library(ggplot2)
```

2.  to compile `glmm_rirs.cpp` and load the resulting binary into the
    environment:

<!-- end list -->

``` r
compile(file = "glmmBoot/src/glmm_rirs.cpp")
```

    ## [1] 0

``` r
dyn.load(dynlib("glmmBoot/src/glmm_rirs"))
```

3.  the functions in the file `glmmBoot.R`:

<!-- end list -->

``` r
source("glmmBoot/R/glmmBoot.R")
```

4.  the dataset `epilepsy.sas7bdat`

<!-- end list -->

``` r
epilepsy <- read.sas7bdat(file = "glmmBoot/data/epilepsy.sas7bdat")
```
