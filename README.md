# glmmBoot

`glmmBoot` consists of a series of tools that for fitting and bootstrapping generalized linear mixed models (GLMMs) that takes advantage of the Template Model Builder (`TMB`), built on` CppAD` and `Eigen`. 

At this stage of development it can handle Mixed Logit and Mixed Poisson Models, yet it has the potential of being extended intended to handle a wide range of statistical distributions. Fixed and random effects models can be specified for the conditional and zero-inflated components of the model, as well as fixed effects models for the dispersion parameter. 

