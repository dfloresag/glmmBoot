glmmBoot
================

`glmmBoot` consists of a series of tools that for fitting and
bootstrapping generalized linear mixed models (GLMMs) that takes
advantage of the Template Model Builder (`TMB`), built on`CppAD` and
`Eigen`.

At this stage of development it can handle Mixed Logit and Mixed Poisson
Models, yet it has the potential of being extended intended to handle a
wide range of statistical distributions. Fixed and random effects models
can be specified for the conditional and zero-inflated components of the
model, as well as fixed effects models for the dispersion parameter.

# Mixed Poisson Example

This code to illustrates the implementation of the bootstrapping methods
described in [(Flores-Agreda and Cantoni, Under
Review)](https://www.researchgate.net/publication/315768128_Bootstrapping_Generalized_Linear_Mixed_Models_via_a_Weighted_Laplace_Approximation)
in a real-data context. Data comes from the  example found in
[Molenberghs & Verbeke
(2006)](https://www.springer.com/gp/book/9780387251448) initially
reported by [Faught et. al.
(1996)](https://www.ncbi.nlm.nih.gov/pubmed/8649570).

The aim of the study was to verify the effects of an new anti-epileptic
drug (AED) compared to a placebo on the number of seizures experienced
by patients during the study. To do this, consider a mixed Poisson model
for the outcome containing two potentially correlated random effects:
one for a random intercept and another one for the visit time i.e.

<img src="img/eq01.png" alt="eq01" height="25">

where:

  - <img src="img/eq02.png" alt="eq02" height="20"> represents the
    effect of the treatment and
  - <img src="img/eq03.png" alt="eq03" height="20"> the visit time.

The variance-covariance structure of the vector of Normal random effects
<img src="img/eq04.png" alt="eq04" height="20"> , comprises a
correlation coefficient <img src="img/eq05.png" alt="eq05" height="14">
, i.e.

<img src="img/eq06.png" alt="eq06" height="50">

To start performing the analysis, you need the following packages in the

``` r
library(lme4)
```

    ## Loading required package: Matrix

``` r
library(glmmTMB)
library(sas7bdat)
library(TMB)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
```

``` r
compile(file = "glmmBoot/src/glmm_rirs.cpp")  
```

    ## [1] 0

``` r
dyn.load(dynlib("glmmBoot/src/glmm_rirs"))

source("glmmBoot/R/glmmBoot.R")
```

``` r
epilepsy <- read.sas7bdat(file = "glmmBoot/data/epilepsy.sas7bdat")
```

``` r
obj.glmerMod <- glmer(nseizw~ trt*studyweek + (studyweek|id), 
                      family = "poisson",
                      data = epilepsy, 
                      control = glmerControl(optimizer = "bobyqa"))

## glmmTMB estimation ####

obj.glmmTMB <- glmmTMB(nseizw~ trt*studyweek + (studyweek|id), 
                       family = "poisson",
                       data = epilepsy)
```

``` r
estimates_1 <-  extract_estimates(obj.glmerMod)
estimates_2 <-  extract_estimates(obj.glmmTMB)
estimates_3 <-  estimate_glmm(obj.glmerMod)
```

``` r
## Bootstrap inference 

rwlb_reps <- bootstrap_glmm(obj.glmerMod,
                            B = 1000,
                            method = "rwlb")
```

    ## Bootstrap Replicates via RWLB

    ## Progress: =                   5%

    ## Progress: ==                  10%

    ## Progress: ===                 15%

    ## Progress: ====                20%

    ## Progress: =====               25%

    ## Progress: ======              30%

    ## Progress: =======             35%

    ## Progress: ========            40%

    ## Progress: =========           45%

    ## Progress: ==========          50%

    ## Progress: ===========         55%

    ## Progress: ============        60%

    ## Progress: =============       65%

    ## Progress: ==============      70%

    ## Progress: ===============     75%

    ## Progress: ================    80%

    ## Progress: =================   85%

    ## Progress: ==================  90%

    ## Progress: =================== 95%

    ## Progress: ====================100%

``` r
rwlb_reps %>% confint(bootstrap_type = "percentile")
```

    ## # A tibble: 7 x 5
    ##   Parameters    Estimate Std_Errors ci_lower ci_upper
    ##   <fct>            <dbl>      <dbl>    <dbl>    <dbl>
    ## 1 (Intercept)     0.896     0.155     0.590   1.18   
    ## 2 rho            -0.324     0.164    -0.634  -0.00546
    ## 3 s_1             1.12      0.129     0.896   1.38   
    ## 4 s_2             0.0473    0.00738   0.0341  0.0632 
    ## 5 studyweek      -0.0273    0.00895  -0.0452 -0.0108 
    ## 6 trt            -0.257     0.247    -0.736   0.218  
    ## 7 trt:studyweek   0.0109    0.0138   -0.0169  0.0391

``` r
rwlb_reps %>% confint(bootstrap_type = "studentized")
```

    ## # A tibble: 7 x 5
    ##   Parameters    Estimate Std_Errors ci_lower ci_upper
    ##   <fct>            <dbl>      <dbl>    <dbl>    <dbl>
    ## 1 (Intercept)     0.896     0.155     0.585   1.19   
    ## 2 rho            -0.324     0.164    -0.742  -0.0566 
    ## 3 s_1             1.12      0.129     0.845   1.33   
    ## 4 s_2             0.0473    0.00738   0.0323  0.0605 
    ## 5 studyweek      -0.0273    0.00895  -0.0445 -0.00987
    ## 6 trt            -0.257     0.247    -0.725   0.216  
    ## 7 trt:studyweek   0.0109    0.0138   -0.0174  0.0395

The Bootstrap distributions were then used to construct (i) Estimates of
the Parameters, by averaging over the replicates (ii) Estimates of the
Standard Errors by computing the standard deviations for each
distribution and (iii) percentile-based Confidence Intervals (CI) with a
level of 95%.
