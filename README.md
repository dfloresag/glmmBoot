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

To start performing the analysis, you need:

1.  the following packages:

<!-- end list -->

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
library(kableExtra)
```

    ## 
    ## Attaching package: 'kableExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     group_rows

2.  to compile `glmm_rirs.cpp` and load the resulting binary:

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

The function `estimate_glmm` uses the frame and estimates of objects
obtained with `lme4::glmer()` as starting points

``` r
obj.glmerMod <- lme4::glmer(nseizw~ trt*studyweek + (studyweek|id), 
                            family = "poisson",
                            data = epilepsy, 
                            control = glmerControl(optimizer = "bobyqa"))
```

For comparison,

``` r
est_lmer <- extract_estimates(obj.glmerMod)
est_tmb  <- estimate_glmm(obj.glmerMod)
```

``` r
est_lmer %>% 
  inner_join(est_tmb, by=("Parameters"), suffix=c(" (`glmer`)", " (`TMB`)")) 
```

    ##      Parameters Estimates (`glmer`) Std..Errors (`glmer`) Estimates (`TMB`)
    ## 1   (Intercept)          0.89447741           0.178194833        0.89448617
    ## 2           trt         -0.24422209           0.253813668       -0.24430762
    ## 3     studyweek         -0.02714089           0.009853978       -0.02714543
    ## 4 trt:studyweek          0.01067076           0.013847317        0.01067175
    ## 5           s_1          1.12722970                    NA        1.12726849
    ## 6           s_2          0.04871389                    NA        0.04871498
    ## 7           rho         -0.33394327                    NA       -0.33387552
    ##   Std..Errors (`TMB`)
    ## 1         0.178577425
    ## 2         0.254353981
    ## 3         0.009920940
    ## 4         0.013933245
    ## 5         0.097487567
    ## 6         0.005702092
    ## 7         0.131163679

The generation of replicates using Random Weighted Laplace Bootstrap
[(Flores-Agreda and Cantoni, Under
Review)](https://www.researchgate.net/publication/315768128_Bootstrapping_Generalized_Linear_Mixed_Models_via_a_Weighted_Laplace_Approximation)
is carried by the function `bootstrap_glmm()`, with the option
`method=rwlb`.

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

The resulting Bootstrap replicates are stored in an object of the class
`glmmBoot`, which can be used to construct

  - Estimates of the Parameters: by averaging over the replicates
  - Estimates of the Standard Errors: by computing the standard
    deviations for the replicates and
  - Percentile-based Confidence Intervals (CI) with a level of 95%.

We have implemented the following

  - `confint()` : which produces bootstrap estimates, standard errors
    and confidence intervals with the
      - the `percentile`

<!-- end list -->

``` r
rwlb_reps %>% 
  confint(bootstrap_type = "percentile")
```

    ## # A tibble: 7 x 5
    ##   Parameters    Estimate Std_Errors ci_lower ci_upper
    ##   <fct>            <dbl>      <dbl>    <dbl>    <dbl>
    ## 1 (Intercept)     0.896     0.149     0.610   1.19   
    ## 2 rho            -0.330     0.166    -0.641  -0.00130
    ## 3 s_1             1.12      0.130     0.891   1.41   
    ## 4 s_2             0.0473    0.00743   0.0333  0.0626 
    ## 5 studyweek      -0.0270    0.00867  -0.0446 -0.0104 
    ## 6 trt            -0.259     0.247    -0.747   0.218  
    ## 7 trt:studyweek   0.0115    0.0134   -0.0146  0.0370

    - and the `studentized` methods

``` r
rwlb_reps %>% 
  confint(bootstrap_type = "studentized")
```

    ## # A tibble: 7 x 5
    ##   Parameters    Estimate Std_Errors ci_lower ci_upper
    ##   <fct>            <dbl>      <dbl>    <dbl>    <dbl>
    ## 1 (Intercept)     0.896     0.149     0.618   1.19   
    ## 2 rho            -0.330     0.166    -0.745  -0.0537 
    ## 3 s_1             1.12      0.130     0.835   1.34   
    ## 4 s_2             0.0473    0.00743   0.0313  0.0603 
    ## 5 studyweek      -0.0270    0.00867  -0.0445 -0.00949
    ## 6 trt            -0.259     0.247    -0.752   0.206  
    ## 7 trt:studyweek   0.0115    0.0134   -0.0158  0.0380

  - `plot()` : which provides a visualization of the replicates for all
    the parameters or a subset, e.g.
    
      - fixed effect parameters

<!-- end list -->

``` r
plot(rwlb_reps, 
     parm_subset=c("(Intercept)", "trt", "studyweek", "trt:studyweek")) + 
  ggtitle("Replicates of Fixed Effects") +
  theme_classic()
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

    - variances of random effects 

``` r
plot(rwlb_reps, 
     parm_subset=c("s_1", "s_2", "rho")) + 
  ggtitle("Replicates of Variance Components")+theme_classic()
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->
