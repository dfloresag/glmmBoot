Mixed Poisson Implementation
================

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
compile(file = "src/glmm_rirs.cpp")  
```

    ## [1] 0

``` r
dyn.load(dynlib("src/glmm_rirs"))

source("R/glmmBoot.R")
```

``` r
epilepsy <- read.sas7bdat(file = "data/epilepsy.sas7bdat")
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
    ## 1 (Intercept)     0.893     0.152     0.607    1.19  
    ## 2 rho            -0.317     0.168    -0.645    0.0107
    ## 3 s_1             1.11      0.130     0.890    1.40  
    ## 4 s_2             0.0471    0.00756   0.0332   0.0624
    ## 5 studyweek      -0.0272    0.00883  -0.0449  -0.0110
    ## 6 trt            -0.248     0.261    -0.745    0.271 
    ## 7 trt:studyweek   0.0106    0.0139   -0.0161   0.0381

``` r
rwlb_reps %>% confint(bootstrap_type = "studentized")
```

    ## # A tibble: 7 x 5
    ##   Parameters    Estimate Std_Errors ci_lower ci_upper
    ##   <fct>            <dbl>      <dbl>    <dbl>    <dbl>
    ## 1 (Intercept)     0.893     0.152     0.605   1.20   
    ## 2 rho            -0.317     0.168    -0.770  -0.0452 
    ## 3 s_1             1.11      0.130     0.834   1.33   
    ## 4 s_2             0.0471    0.00756   0.0313  0.0600 
    ## 5 studyweek      -0.0272    0.00883  -0.0441 -0.00973
    ## 6 trt            -0.248     0.261    -0.755   0.274  
    ## 7 trt:studyweek   0.0106    0.0139   -0.0157  0.0394

The Bootstrap distributions were then used to construct (i) Estimates of
the Parameters, by averaging over the replicates (ii) Estimates of the
Standard Errors by computing the standard deviations for each
distribution and (iii) percentile-based Confidence Intervals (CI) with a
level of 95%.
