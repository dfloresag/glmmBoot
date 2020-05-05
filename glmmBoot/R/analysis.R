## MFDLD datasets ####

rm(list=ls())

library(lme4)
library(glmmTMB)
library(sas7bdat)
library(TMB)
library(dplyr)

compile(file = "./glmmBoot/src/glmm_rirs.cpp")  
dyn.load(dynlib("./glmmBoot/src/glmm_rirs"))

epilepsy <- read.sas7bdat(file = "./glmmBoot/data/epilepsy.sas7bdat")

source("./glmmBoot/R/glmmBoot.R")

## lme4 estimation ####

obj.glmerMod <- glmer(nseizw~ trt*studyweek + (studyweek|id), 
                      family = "poisson",
                      data = epilepsy, 
                      control = glmerControl(optimizer = "bobyqa"))

## glmmTMB estimation ####

obj.glmmTMB <- glmmTMB(nseizw~ trt*studyweek + (studyweek|id), 
                       family = "poisson",
                       data = epilepsy)

estimates_1 <-  extract_estimates(obj.glmerMod)
estimates_2 <-  extract_estimates(obj.glmerMod)
estimates_3 <-  estimate_glmm(obj.glmerMod)

## Bootstrap inference 

rwlb_reps <- bootstrap_glmm(obj.glmerMod,
                            B = 1000,
                            method = "rwlb")

parb_reps <- bootstrap_glmm(obj.glmerMod,
                            B = 1000,
                            method = "parb")


confint.glmmBoot <- function(obj, 
                             parm, 
                             level = 0.95,
                             boot.method = c("rwlb","parb"),
                             boot.type = c("percentile","studentized","norm"), ...) {
  
  # if(inherits(obj, "glmerMod")) {
  #   boot.method <- match.arg(boot.method)
  #   tmp <- bootstrap_glmm(obj.glmerMod,
  #                         B = 1000,
  #                         method = boot.method)
  # } else if (inherits(obj, "glmerMod")){
    tmp <- obj
  # }
  # extract_estimates(tmp)
  
  boot.type <- match.arg(boot.type)
  
  if(boot.type=="percentile"){
    summary_boot <- tmp$replicates %>% 
      group_by(Parameters) %>% 
      # TODO: add options to select parameters
      summarise(Estimate     = mean(Estimates),
                Std_Errors   =  sd(Estimates), 
                ci_lower     = quantile(Estimates, prob = (1-level)/2), 
                ci_upper     = quantile(Estimates, prob = (1+level)/2)) %>% 
      ungroup()
  } else if (boot.type == "studentized"){
    # TODO
  }
  summary_boot
}

confint(rwlb_reps, boot.method = "rwlb", boot.type = "percentile")
confint(parb_reps, boot.method = "parb", boot.type = "percentile")


parb_reps %>% 
  group_by(Parameters) %>% 
  summarise(Estimate = mean(Estimates),
            Std_Error= sd(Estimates), 
            Quantile_2.5 = quantile(Estimates, prob = 0.025), 
            Quantile_97.5 = quantile(Estimates, prob = 0.975)) %>% 
  ungroup() %>% mutate(Method="parb")


rwlb_reps %>% 
  filter(Parameters == '(Intercept)') %>% 
  select(Estimates) %>% 
  boxplot()
abline(h= est_tmp[est_tmp$Parameters=='(Intercept)','Estimates'], 
       col = 'blue')

parb_reps %>% 
  filter(Parameters == '(Intercept)') %>% 
  select(Estimates) %>% 
  boxplot()
abline(h= est_tmp[est_tmp$Parameters=='(Intercept)','Estimates'], 
       col = 'blue')

rwlb_reps %>% 
  filter(Parameters == 'trt') %>% 
  select(Estimates) %>% 
  boxplot()
abline(h= est_tmp[est_tmp$Parameters=='trt','Estimates'], 
       col = 'blue')

rwlb_reps %>% 
  filter(Parameters == 's_1') %>% 
  select(Estimates) %>% 
  boxplot()
abline(h= est_tmp[est_tmp$Parameters=='s_1','Estimates'], 
       col = 'blue')

rwlb_reps %>% 
  filter(Parameters == 's_2') %>% 
  select(Estimates) %>% 
  boxplot()
abline(h= est_tmp[est_tmp$Parameters=='s_2','Estimates'], 
       col = 'blue')

rwlb_reps %>% 
  filter(Parameters == 'rho') %>% 
  select(Estimates) %>% 
  boxplot()
abline(h= est_tmp[est_tmp$Parameters=='rho','Estimates'], 
       col = 'blue')

