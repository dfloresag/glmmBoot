#' @title extract_info.glmerMod : extract model info and frame from `glmerMod` object.
#'
#' @param obj Object of class `glmerMod` resulting from estimation with `lme4::glmer()`
#' @param data Original dataset. TODO: make it general by extracting it from `obj`
#' @param ... other arguments 
#' @examples
#' @rdname glmmBoot
#' @export
#' @import lme4
extract_info.glmerMod <-  function(obj, data = epilepsy, ...){
  
  if(inherits(obj, "glmerMod")){
    
    bs  = unname(fixef(obj))
    s1  = unname(attr(VarCorr(obj)$id, "stddev")[1])
    s2  = unname(attr(VarCorr(obj)$id, "stddev")[2])
    rh  = unname(attr(VarCorr(obj)$id, "correlation"))[1,2]
    
    se_bs = unname(sqrt(diag(vcov(obj)))) 
    
    fix_nm  = names(fixef(obj))
    
    ng  =  nlevels(factor(obj@frame$id))
    
    group = factor(obj@frame$id, labels = 1:ng)
    y     = obj@resp$y
    X     = as.matrix(obj@pp$X[1:length(y), 1:length(bs)])
    
    list(
      bs  = bs,
      s1  = s1,
      s2  = s2,
      rh  = rh,
      
      se_bs  = se_bs ,
      fix_nm = fix_nm,
      
      ng  =  ng,
      
      group = group,
      y     = y,
      X     = X 
    )  
  } else {
    message('function only works with objects of class glmerMod')
  }
}

#' @title extract_estimates: extracts estimates of `glmmTMB`` and `glmerMod` objects in a tidy way
#'
#' @param obj Object 
#' 
#' @param data Original dataset. TODO: make it general by extracting it from `obj`
#' @param ... other arguments 
#' @examples
#' @rdname glmmBoot
#' @export
#' @import lme4
extract_estimates <- function(obj, ... ){
  if(inherits(obj, "glmerMod")){
    tmp <- extract_info.glmerMod(obj)
    bs  <- tmp$bs
    s1  <- tmp$s1
    s2  <- tmp$s2
    rh  <- tmp$rh
    
    se_bs <- tmp$se_bs
    
    fix_nm <- tmp$fix_nm
    
  } else if(inherits(obj, "glmmTMB")){
    
    bs  <-  unname(fixef(obj)$cond)
    s1  <-  unname(attr(VarCorr(obj)$cond$id, "stddev")[1])
    s2  <-  unname(attr(VarCorr(obj)$cond$id, "stddev")[2])
    rh  <-  unname(attr(VarCorr(obj)$cond[[1]], "correlation"))[1,2]
    
    se_bs <- unname(sqrt(diag(vcov(obj)$cond)))
    
    fix_nm <- names(fixef(obj)[[1]])
    
  } 
  
  data.frame(
    Parameters = c(fix_nm, 
                   "s_1", "s_2", "rho"), 
    Estimates = c(bs, s1, s2, rh),
    'Std. Errors' =  c(se_bs, NA, NA, NA)
  )
}

#' @title estimate_glmm: estimate 
#'
#' @param obj Object of class `glmerMod` resulting from estimation with `lme4::glmer()`
#' @param data Original dataset. TODO: make it general by extracting it from `obj`
#' @param ... other arguments 
#' @examples
#' @rdname glmmBoot
#' @export
#' @import lme4
estimate_glmm <- function(obj.glmerMod, 
                          data = epilepsy, 
                          weights = NULL,
                          predictions = FALSE,
                          newData= NULL, 
                          maxit = 100, ...){
  
  tmp <- extract_info.glmerMod(obj.glmerMod)
  
  bs_0  <- tmp$bs
  s1_0  <- tmp$s1
  s2_0  <- tmp$s2
  rh_0  <- tmp$rh
  
  fix_nm <- tmp$fix_nm
  
  ng <- tmp$ng
  
  group <- tmp$group
  
  if(is.null(newData)){
    y <- tmp$y
    X <- tmp$X
  } else{
    y <- newData$y
    X <- newData$X
  }
  
  if(is.null(weights)) weights = rep(1, times = ng)
  
  cnv = err <- FALSE
  iter <- 1 
  while(!cnv & !err | iter<=maxit){
    
    dta <- list(group = group, 
                y= y, 
                X= X, 
                t= data$studyweek, 
                w= weights,
                family =1 
    )
    
    prm <- list(
      u = cbind(u1 = rep(0, times = ng),
                u2 =rep(0, times = ng)),
      betas = bs_0,
      log_sigma=log(c(s1_0, s2_0)),
      tan_rho = (pi/2)*tan(rh_0)
    )
    
    obj <- MakeADFun(data=dta,
                     parameters=prm,
                     random=c("u"),
                     type = c("ADFun", "ADGrad" , "Fun"),
                     DLL="glmm_rirs",
                     silent = TRUE,
                     hessian = FALSE, 
                     method = "CG",
                     inner.control = list(maxit = 10000),
                     control = list(maxit =10000)
    )
    
    opt <- try(nlminb(obj$par, obj$fn, obj$gr), silent = TRUE)
    
    if(!inherits(opt, "try-error")){
      if(!(opt$convergence==1)){
        
        cnv = err <- TRUE
        
        sdr_fix <- summary(sdreport(obj),select = 'fixed')
        sdr_rep <- summary(sdreport(obj),select = 'report')
        
        bs <- unname(sdr_fix[1:length(bs_0),'Estimate'])
        s1 <- unname(sdr_rep[1,'Estimate'])
        s2 <- unname(sdr_rep[2,'Estimate'])
        rh <- unname(sdr_rep['rho','Estimate'])
        
        se_bs <- unname(sdr_fix[1:length(bs_0),'Std. Error'])
        se_s1 <- unname(sdr_rep[1,'Std. Error'])
        se_s2 <- unname(sdr_rep[2,'Std. Error'])
        se_rh <- unname(sdr_rep['rho','Std. Error'])
        
        res <- data.frame(
          Parameters = c(fix_nm, 
                         "s_1", "s_2", "rho"), 
          Estimates = c(bs, s1, s2, rh),
          'Std. Errors' =  c(se_bs, se_s1, se_s2, se_rh)
        )
        break
      } else { 
        message("Convergence error: Retrying")
        res <- cnv
        iter <- iter + 1
      } 
    } else { 
      message("Estimation error: Retrying")
      res <- err
      iter <- iter + 1
    }
  }
  
  if(iter>=maxit){
    message("No solution found: Review your initial values")
    return(res)
  } else {
    return(res)  
  }
}

#' @title predict_glmm: extracts estimates of `glmmTMB`` and `glmerMod` objects in a tidy way
#'
#' @param obj Object of class `glmerMod` resulting from estimation with `lme4::glmer()`
#' @param data Original dataset. TODO: make it general by extracting it from `obj`
#' @param ... other arguments 
#' @examples
#' @rdname glmmBoot
#' @export
#' @import lme4
# TODO : Finish preparing this function
predict_glmm  <- function(obj,
                          data = epilepsy,
                          newData = NULL, 
                          weights = NULL, ...){
  
  tmp <- extract_info.glmerMod(obj)
  
  bs_0  <- tmp$bs
  s1_0  <- tmp$s1
  s2_0  <- tmp$s2
  rh_0  <- tmp$rh
  
  fix_nm <- tmp$fix_nm
  
  ng <- tmp$ng
  
  group <- tmp$group
  
  if(is.null(newData)){
    y <- tmp$y
    X <- tmp$X
  } else{
    y <- newData$y
    X <- newData$X
  }
  
  if(is.null(weights)) weights = rep(1, times = ng)
  
  dta <- list(group = group,
              y= y,
              X= X,
              t= data$studyweek,
              w= weights,
              family =1
  )
  
  prm <- list(
    u = cbind(u1 = rep(0, times = ng),
              u2 =rep(0, times = ng)),
    betas = bs_0,
    log_sigma=log(c(s1_0, s2_0)),
    tan_rho = (pi/2)*tan(rh_0)
  )
  
  obj <- MakeADFun(data=dta,
                   parameters=prm,
                   random=c("u"),
                   type = c("ADFun", "ADGrad" , "Fun"),
                   DLL="glmm_rirs",
                   silent = TRUE,
                   hessian = FALSE,
                   method = "CG",
                   inner.control = list(maxit = 10000),
                   control = list(maxit =10000)
  )
  
  opt     <- try(nlminb(obj$par, obj$fn, obj$gr), silent = TRUE)
  
  if(!inherits(opt, "try-error")){
    if(!(opt$convergence==1)){
      
      cnv = err <- TRUE
      
      sdr_ref <- try(summary(sdreport(obj),select = 'random'))
      
      if(!inherits(opt, "try-error")){
        res <- data.frame(
          Parameters = c(paste("u1", 1:ng),
                         "s_1", "s_2", "rho"),
          Cond_Mode = c(bs, s1, s2, rh),
          Cond_Var =  c(se_bs, se_s1, se_s2, se_rh)
        )
        break
      } else {
        message("Prediction error: Retrying")
        res <- cnv
        iter <- iter + 1
      }
    } else {
      message("Estimation error: Retrying")
      res <- err
      iter <- iter + 1
    }
  }
}

#' @title simulate_glmm: simulates 
#'
#' @param obj Object of class `glmerMod` resulting from estimation with `lme4::glmer()`
#' @param data Original dataset. TODO: make it general by extracting it from `obj`
#' @param ... other arguments 
#' @examples
#' @rdname glmmBoot
#' @export
#' @import lme4
# TODO : adapt this function to retrieve parameters resulting from *the template* and not the glmerMod object
simulate_glmm <- function(obj, 
                          data = epilepsy, 
                          newData = NULL, ...){
  
  tmp <- extract_info.glmerMod(obj)
  
  bs  <- tmp$bs
  s1  <- tmp$s1
  s2  <- tmp$s2
  rh  <- tmp$rh
  
  
  
  fix_nm <- tmp$fix_nm
  
  ng <- tmp$ng
  
  group <- tmp$group
  
  if(is.null(newData)){
    X <- tmp$X
  } else{
    X <- newData$X
  }
  
  u <- MASS::mvrnorm(
    n=ng, 
    mu = c(0,0), 
    Sigma = matrix(
      c(s1^2, rh*s1*s2, 
        rh*s1*s2 , s2^2),
      nrow = 2, ncol=2)) 
  
  eta <- X%*%bs + u[group, 1] + data$studyweek*u[group, 2]
  
  rpois(n = length(eta), lambda =  exp(eta))
  
}

#' @title simulate_glmm: simulates 
#'
#' @param obj Object of class `glmerMod` resulting from estimation with `lme4::glmer()`
#' @param data Original dataset. TODO: make it general by extracting it from `obj`
#' @param ... other arguments 
#' @examples
#' @rdname glmmBoot
#' @export
#' @import lme4
# TODO : adapt this function to retrieve parameters resulting from *the template* and not the glmerMod object
bootstrap_glmm <- function(obj.glmerMod, 
                           B = 1000, 
                           data = epilepsy,
                           method = c("rwlb", "parb")){
  
  tmp <- extract_info.glmerMod(obj.glmerMod)
  
  bs_0  <- tmp$bs
  s1_0  <- tmp$s1
  s2_0  <- tmp$s2
  rh_0  <- tmp$rh
  
  fix_nm <- tmp$fix_nm
  
  ng <- tmp$ng
  
  group <- tmp$group
  y     <- tmp$y
  X     <- tmp$X
  
  if (method == "rwlb"){
    message("Bootstrap Replicates via RWLB")
    rep_list <- lapply(
      X   = as.list(1:B), 
      FUN = function(x) {
        show_progress(x, B)
        data.frame(
          Replicate = x,
          estimate_glmm(
            obj = obj,
            weights    =  rexp(n = ng, rate = 1))
        )
      }
    )
  } else if (method == "parb"){
    message("Bootstrap Replicates via Parametric Bootstrap")
    rep_list <- lapply(
      X   = as.list(1:B), 
      FUN = function(x) {
        show_progress(x, B)
        y_prb <- simulate_glmm(obj,newData = list(X=tmp$X))
        data.frame(
          Replicate = x,
          estimate_glmm(
            obj, newData = list(y=y_prb, X=tmp$X))
        )
      }
    )
  }
  structure(list(
    replicates = dplyr::bind_rows(rep_list),
    model_object = obj.glmerMod,
    method = method
  ),
  class = "glmmBoot")
}

#' @title simulate_glmm: simulates 
#'
#' @param obj Object of class `glmerMod` resulting from estimation with `lme4::glmer()`
#' @param data Original dataset. TODO: make it general by extracting it from `obj`
#' @param ... other arguments 
#' @examples
#' @rdname glmmBoot
#' @export
#' @import lme4
# TODO : adapt this function to retrieve parameters resulting from *the template* and not the glmerMod object
confint.glmmBoot <- function(obj, 
                             level = 0.95,
                             parm_subset=NULL, 
                             bootstrap_type = c("percentile","studentized","norm"), ...) {
  
  tmp <- obj$replicates
  
  if(!is.null(parm_subset)) tmp <- tmp %>% filter(Parameters %in% parm_subset)
  
  bootstrap_type <- match.arg(bootstrap_type)
  if(bootstrap_type=="percentile"){
    sum_boot <- tmp %>% 
      group_by(Parameters) %>%
      summarise(Estimate     = mean(Estimates),
                Std_Errors   = sd(Estimates), 
                ci_lower     = quantile(Estimates, prob = (1-level)/2), 
                ci_upper     = quantile(Estimates, prob = (1+level)/2)) %>% 
      ungroup()
  } else if (bootstrap_type == "studentized"){
    est <- estimate_glmm(obj$model_object)
    
    tmp <- tmp %>% inner_join(est, by =c("Parameters"),   suffix=c("_star", "_hat"))
    
    sum_boot <- tmp %>% 
      group_by(Parameters) %>%
      summarise(Estimate     = mean(Estimates_star),
                Std_Errors   = sd(Estimates_star), 
                ci_lower     = mean(stud_quant(Estimates_star, Std..Errors_star, 
                                             Estimates_hat , Std..Errors_hat, prob = (1-level)/2)), 
                ci_upper     = mean(stud_quant(Estimates_star, Std..Errors_star, 
                                             Estimates_hat, Std..Errors_hat, prob = (1+level)/2))) %>% 
      ungroup()
  }
  sum_boot
}

#' @title plot.glmmBoot: simulates 
#'
#' @param obj Object of class `glmerMod` resulting from estimation with `lme4::glmer()`
#' @param data Original dataset. TODO: make it general by extracting it from `obj`
#' @param ... other arguments 
#' @examples
#' @rdname glmmBoot
#' @export
#' @import lme4
# TODO : adapt this function to retrieve parameters resulting from *the template* and not the glmerMod object
plot.glmmBoot <- function(obj, parm_subset = NULL, ...){
  
  # TODO: separate between fixed and random effect parameters
  
  tmp <- obj$replicates
  
  if(!is.null(parm_subset)) tmp <- tmp %>% filter(Parameters %in% parm_subset)
  
  library(ggplot2)
  
  est <- estimate_glmm(obj$model_object)
  tmp <-  tmp%>%inner_join(est, by = "Parameters", suffix = c("_boot", "_glmm"))
  
  tmp %>% 
    ggplot() + 
    geom_boxplot(aes(y=Estimates_boot, group=Parameters))+
    geom_hline(aes(yintercept = Estimates_glmm), color = "blue", lty =2)+
    facet_wrap(.~Parameters, scales = "free_y" )+
    ylab(label = "Bootstrap Replicates")
}

#### Internal functions ####
show_progress <- function(x, B) {
  if (x%%(B/20)==0) {
    message("Progress: ", paste(rep("=", times=(x/(B/20)))), paste(rep(" ", times=(20-x/(B/20)))), x/B*100,"%")
  }
}

stud_quant <-  function(x_star, sd_x_star, x_hat, sd_x_hat, prob, ...){
  t_star <-  (x_star-x_hat)/sd_x_star
  z <-  quantile(t_star, prob, ...)
  x_hat + z*sd_x_hat
}




# TODO: 
# print.glmmBoot <- function(){}
# summary.glmmBoot <- function(){}
# plot.glmmBoot <- function(){}
