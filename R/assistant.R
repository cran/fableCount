#' @importFrom tibble as_tibble
split_tibble = function(input_tibble) {
  # Calculate the number of rows needed
  num_rows = nrow(input_tibble) %/% 4

  # Initialize an empty list to store rows
  rows = list()

  # Loop through each row and extract the data for 4 columns
  for (i in 1:num_rows) {
    start_row = (i - 1) * 4 + 1
    end_row = min(i * 4, nrow(input_tibble))
    rows[[i]] = input_tibble[start_row:end_row, 1]
  }

  # Pad rows with NA values if necessary to ensure consistent column lengths
  max_length = max(sapply(rows, length))
  rows = sapply(rows, function(row) c(row, rep(NA, max_length - length(row))))

  # Create a tibble from the list of rows
  result_tibble = as_tibble(do.call(cbind, rows))

  return(result_tibble)
}







matrix_to_list = function(input_matrix) {
  num_rows = nrow(input_matrix)
  row_list = vector("list", length = num_rows)
  for (i in 1:num_rows) {
    row_list[[i]] = input_matrix[i, ]
  }
  row_list = lapply(1:num_rows,
                    function(i){
                      input_matrix[i, ]
                    }
  )
  return(row_list)
}





g = function(x, link=c("identity", "log")){
  link = match.arg(link)
  result = if(is.null(x)) NULL else switch(link,
                                            "identity" = x,
                                            "log" = log(x)
  )
  return(result)
}

#Transformation function:
trafo = function(x, link=c("identity", "log")){
  link = match.arg(link)
  result = if(is.null(x)) NULL else switch(link,
                                            "identity" = x,
                                            "log" = log(x+1)
  )
  return(result)
}

#Inverse of transformation function:
trafo_inv = function(x, link=c("identity", "log")){
  link = match.arg(link)
  result = if(is.null(x)) NULL else switch(link,
                                            "identity" = x,
                                            "log" = exp(x)-1
  )
  return(result)
}

#Inverse of link function:
g_inv = function(x, link=c("identity", "log")){
  link = match.arg(link)
  result = if(is.null(x)) NULL else switch(link,
                                            "identity" = x,
                                            "log" = exp(x)
  )
  return(result)
}

#1st derivative of inverse link function:
g_inv_1st = function(x, link=c("identity", "log")){
  link = match.arg(link)
  result = if(is.null(x)) NULL else switch(link,
                                            "identity" = 1,
                                            "log" = exp(x)
  )
  return(result)
}

#2nd derivative of inverse link function:
g_inv_2nd = function(x, link=c("identity", "log")){
  link = match.arg(link)
  result = if(is.null(x)) NULL else switch(link,
                                            "identity" = 0,
                                            "log" = exp(x)
  )
  return(result)
}

tsglm.parameterlist = function(paramvec, model){
  p = length(model$past_obs)
  P = seq(along=numeric(p)) #sequence 1:p if p>0 and NULL otherwise
  q = length(model$past_mean)
  Q = seq(along=numeric(q)) #sequence 1:p if p>0 and NULL otherwise
  r = length(paramvec) - (1+p+q)
  R = seq(along=numeric(r))
  names(paramvec) = NULL
  result = list(intercept=paramvec[1], past_obs=paramvec[1+P], past_mean=paramvec[1+p+Q], xreg=paramvec[1+p+q+R])
  return(result)
}

tsglm.parametercheck = function(param, link=c("identity", "log"), stopOnError=TRUE, silent=TRUE){
  #Check parameter vector of a count time series following GLMs

  ##############
  #Checks and preparations:
  link = match.arg(link)

  if(link == "identity") parametercheck = function(param){
    stopifnot(
      param$intercept>0,
      param$past_obs>=0,
      param$past_mean>=0,
      param$xreg>=0
    )
    sum_param_past = sum(param$past_obs)+sum(param$past_mean)
    if(sum_param_past>=1) stop(paste("Parameters are outside the stationary region, sum of parameters for regression\non past observations and on past conditional means is", sum_param_past, "> 1"))
    return(TRUE)
  }

  if(link == "log") parametercheck = function(param){
    stopifnot(
      abs(param$past_obs)<1,
      abs(param$past_mean)<1
    )
    sum_param_past = abs(sum(param$past_obs)+sum(param$past_mean))
    if(sum_param_past>=1) stop(paste("Parameters are outside the stationary region, absolute sum of parameters for\nregression on past observations and on past conditional means is", sum_param_past, "> 1"))
    return(TRUE)
  }

  if(stopOnError){
    result = parametercheck(param)
  }else{
    result = try(parametercheck(param), silent=silent)
    if(inherits(result, "try-error")) result = FALSE
  }
  return(result)
}

simcoefs = function(...) UseMethod("simcoefs")
#' @exportS3Method tscount::simcoefs
#' @importFrom tscount tsglm
#' @importFrom parallel parSapply
#' @importFrom stats rnorm
#' @importFrom stats vcov
#' @importFrom stats rnorm
simcoefs.tsglm = function(fit, method=c("bootstrap", "normapprox"), B=1, parallel=FALSE, ...){
  stopifnot(
    length(B)==1,
    B%%1==0,
    B>=1
  )
  method = match.arg(method)
  if(method=="bootstrap"){
    simfit = function(seed, fit, ...){
      set.seed(seed)
      ts_sim = tsglm.sim(fit=fit)$ts
    fit_sim = tsglm(ts=ts_sim, model=fit$model, xreg=fit$xreg, link=fit$link, distr=fit$distr, score=FALSE, info="none", ...)
      if(fit$distr=="nbinom" && fit_sim$distr=="poisson") fit_sim$distrcoefs = c(size=NA)
      result = c(coef(fit_sim), sigmasq=fit_sim$sigmasq, fit_sim$distrcoefs)
      return(result)
    }
    seeds = sample(1e+9, size=B)
    if(parallel){
      Sapply = function(X, FUN, ...) parSapply(cl=NULL, X=X, FUN=FUN, ...)
    }else{
      Sapply = sapply
    }
    coefs = t(Sapply(seeds, simfit, fit=fit, ...))
    result = list(coefs=coefs)
  }
  if(method=="normapprox"){
    rmvnorm_stable = function(n, mean=rep(0, nrow(sigma)), sigma=diag(length(mean))){
      #Function for stable generation of random values from a multivariate normal distribution (is robust against numerical deviations from symmetry of the covariance matrix. Code is taken from function rmvnorm in the package mvtnorm and modified accordingly.
      if(length(mean) != nrow(sigma)) stop("mean and sigma have non-conforming size")
      ev = eigen(sigma, symmetric=TRUE)
      ev$values[ev$values < 0] = 0
      R = t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
      centred = matrix(rnorm(n=n*ncol(sigma)), nrow=n, byrow=TRUE) %*% R
      result = sweep(centred, 2, mean, "+")
      colnames(result) = names(mean)
      return(result)
    }
    f = 1.1 #one could choose this factor according to the probability of a parameter from the multivariate normal distribution to be outside the parameter space
    coefs = rmvnorm_stable(n=ceiling(f*B), mean=coef(fit), sigma=vcov(fit))
    repeat{
      valid_coefs = apply(coefs, 1, function(x) tsglm.parametercheck(tsglm.parameterlist(paramvec=x, model=fit$model), link=fit$link, stopOnError=FALSE))
      if(sum(valid_coefs) >= B) break
      coefs = rbind(coefs, rmvnorm_stable(n=ceiling((B-sum(valid_coefs))*f/mean(valid_coefs)), mean=coef(fit), sigma=vcov(fit)))
    }
    use_coefs = which(valid_coefs)[1:B]
    coefs = coefs[use_coefs, , drop=FALSE]
    n_invalid = max(use_coefs) - B
    distrcoefs_matrix = matrix(rep(c(sigmasq=fit$sigmasq, fit$distrcoefs), B), byrow=TRUE, nrow=B)
    colnames(distrcoefs_matrix) = c("sigmasq", names(fit$distrcoefs))
    coefs = cbind(coefs, distrcoefs_matrix)
    result = list(coefs=coefs, n_invalid=n_invalid)
  }
  return(result)
}


#' @importFrom stats is.ts
#' @importFrom stats start
#' @importFrom stats frequency
#' @importFrom stats coef
#' @importFrom stats window
#' @importFrom stats tsp
#' @importFrom tscount qdistr
#' @importFrom tscount pdistr
#' @importFrom tscount ddistr
#' @importFrom tscount tsglm.sim
predict_call_tsglm = function(object, n.ahead=1, newobs=NULL, newxreg=NULL, level=0.95, global=FALSE, type=c("quantiles", "shortest", "onesided"), method=c("conddistr", "bootstrap"), B=1000, estim=c("ignore", "bootstrap", "normapprox", "given"), B_estim=B, coefs_given=NULL, ...){
  newxreg = if(is.null(newxreg)) matrix(0, nrow=n.ahead, ncol=ncol(object$xreg)) else as.matrix(newxreg)  #if no covariates are provided, these are set to zero
  stopifnot(n.ahead>0,
            n.ahead%%1==0,
            ncol(newxreg)==ncol(object$xreg)
  )
  n = object$n_obs
  link = object$link
  model = object$model
  p = length(model$past_obs)
  P = seq(along=numeric(p)) #sequence 1:p if p>0 and NULL otherwise
  q = length(model$past_mean)
  Q = seq(along=numeric(q)) #sequence 1:p if p>0 and NULL otherwise
  r = ncol(object$xreg)
  R = seq(along=numeric(r))
  xreg = rbind(object$xreg, newxreg)
  if(nrow(xreg) < n+n.ahead) xreg = rbind(xreg, matrix(xreg[nrow(xreg), ], nrow=n+n.ahead-nrow(xreg), ncol=r, byrow=TRUE)) #if not enough future values of the covariates are given, then the values of the last available time point are used, which is usually NOT sensible!
  new = rep(NA, n.ahead)
  new[seq(along=newobs)] = newobs
  ts = c(object$ts, new)
  if(is.ts(object$ts)) ts = ts(ts, start=start(object$ts), frequency=frequency(object$ts))
  nu = c(rep(NA, n-object$n_eff), object$linear.predictors, rep(NA, n.ahead))
  if(is.ts(object$ts)) nu = ts(nu, start=start(object$ts), frequency=frequency(object$ts))
  for(t in n+(1:n.ahead)){
    nu[t] = sum(coef(object)*c(1, trafo(ts[t-model$past_obs], link=link), nu[t-model$past_mean]-if(r>0){sum((as.numeric(model$external)*coef(object)[1+p+q+R])*t(xreg[t-model$past_mean,]))}else{0}, xreg[t,]))
    if(is.na(ts[t])) ts[t] = g_inv(nu[t], link=link) #unobserved future observations are replaced by their prediction (by the conditional mean)
  }
  if(is.ts(object$ts)){
    pred = window(g_inv(nu, link=link), start=tsp(object$ts)[2]+1/frequency(object$ts)) #use time series class if input time series has this class
  }else{
    pred = g_inv(nu, link=link)[n+(1:n.ahead)]
  }
  result = list(pred=pred)

  #Prediction intervals:
  stopifnot(
    length(level)==1,
    !is.na(level),
    level<1,
    level>=0
  )
  if(level>0){ #do not compute prediction intervals for level==0
    warning_messages = character()
    level_local = if(global){1-(1-level)/n.ahead}else{level} #Bonferroni adjustment of the coverage rate of each individual prediction interval such that a global coverage rate is achieved
    method = match.arg(method, several.ok=TRUE)
    conddistr_possible = n.ahead==1 || (length(newobs)>=n.ahead-1 && !any(is.na(newobs[1:(n.ahead-1)]))) #method="conddistr" is only possible if all predictions are 1-step ahead predictions such that the true previous observations are available (this is the case if only one observation is to be predicted or when the first n.ahead-1 observations are given in argument 'newobs')
    if(all(c("conddistr", "bootstrap") %in% method)) method = if(conddistr_possible){"conddistr"}else{"bootstrap"} #automatic choice of the method, prefer method="conddistr" if possible
    type = match.arg(type)
    if(type=="quantiles") a = (1-level_local)/2
    if(type=="onesided") a = 1-level_local
    if(type=="shortest"){
      largestdensityinterval = function(probs, level){ #find shortest interval with given probability, probs are assumed to be the probabilities corresponding to the values seq(along=probs)-1
        values = seq(along=probs)-1
        ord = order(probs, decreasing=TRUE)
        result = range(values[ord][1:which(cumsum(probs[ord])>=level)[1]])
        return(result)
      }
    }
    estim = match.arg(estim)
    if(estim %in% c("bootstrap", "normapprox") && method!="bootstrap") stop("Accounting for the estimation uncertainty is only available if argument 'method'\nis set to\"bootstrap\".")
    if(method=="conddistr"){
      if(!conddistr_possible) stop("Computation of prediction intervals with argument 'method' set to \"conddistr\"\ndoes only work for 1-step-ahead predictions. If argument 'n.ahead' is larger\nthan 1, future observations have to be provided in argument 'newobs'.")
      if(type %in% c("quantiles", "onesided")){
        qdistr_mod = function(p, meanvalue, distr=c("poisson", "nbinom"), distrcoefs, upper=FALSE){ #allows for alternative definition of the quantile
          result = qdistr(p=p, meanvalue=meanvalue, distr=distr, distrcoefs=distrcoefs)
          if(upper) result = result + (pdistr(q=result, meanvalue=meanvalue, distr=distr, distrcoefs=distrcoefs)==p) #alternative definition of the quantile
          return(result)
        }
        lower = if(type=="onesided"){integer(n.ahead)}else{qdistr_mod(a, meanvalue=pred, distr=object$distr, distrcoefs=object$distrcoefs)}
        upper = qdistr_mod(1-a, meanvalue=pred, distr=object$distr, distrcoefs=object$distrcoefs, upper=TRUE)
        predint = cbind(lower=lower, upper=upper)
      }
      if(type=="shortest"){
        cutoff = qdistr(p=1-(1-level_local)/10, meanvalue=max(pred), distr=object$distr, distrcoefs=object$distrcoefs) #very large quantile which is almost certainly larger than the upper bound of the shortest prediction interval
        predint = t(sapply(pred, function(predval) largestdensityinterval(probs=ddistr(x=0:cutoff, meanvalue=predval, distr=object$distr, distrcoefs=object$distrcoefs), level=level_local)))
      }
      predmed = qdistr(p=0.5, meanvalue=pred, distr=object$distr, distrcoefs=object$distrcoefs)
      B = NULL
    }
    if(method=="bootstrap"){
      if(!is.null(newobs)) stop("Computation of prediction intervals with argument 'method' set to \"bootstrap\"\ncan currently not take into account the future observations given in argument\n'newobs'. Either set argument 'method' to \"conddistr\" or do not compute\nprediction intervals at all by setting argument 'level' to zero.")
      stopifnot(
        length(B)==1,
        B%%1==0,
        B>=10
      )
      if(estim %in% c("normapprox", "bootstrap")) stopifnot(length(B_estim)==1, B_estim%%1==0, B_estim>=1)
      if(estim=="normapprox"){
        coefs_complete = simcoefs(object, method="normapprox", B=B_estim)
        coefs = coefs_complete$coefs[, -(1+p+q+r+1)] #remove column 'sigmasq'
        n_invalid = coefs_complete$n_invalid
        if(n_invalid > 0){
          bootstrap_message = paste("In", n_invalid, "cases the bootstrap samples of the regression parameter generated\nfrom the normal approximation was not valid and the respective sample has been\ngenerated again. Note that this might affect the validity of the final result.")
          warning_messages = c(warning_messages, bootstrap_message)
          warning(bootstrap_message)
        }
      }
      if(estim=="bootstrap"){
        coefs_complete = simcoefs(object, method="bootstrap", B=B_estim, ...)
        coefs = coefs_complete$coefs[, -(1+p+q+r+1)] #remove column 'sigmasq'
        n_invalid = sum(coefs_complete$coefs[, 1+p+q+r+1]==1e-9)
        if(n_invalid > 0){
          bootstrap_message = paste("In", n_invalid, "cases of the parametric bootstrap to account for estimation uncertainty\nthe dispersion parameter could not be estimated and a Poisson distribution is\nused instead. Note that this might affect the validity of the final result.")
          warning_messages = c(warning_messages, bootstrap_message)
          warning(bootstrap_message)
        }
      }
      if(estim=="given"){
        if(is.null(coefs_given) || nrow(coefs_given)<1) stop("Argument 'coefs_given' needs to be a matrix with the parameters in the rows.")
        coefs = coefs_given
        B_estim = nrow(coefs)
      }
      if(estim=="ignore"){
        coefs = t(replicate(B, c(coef(object), object$distrcoefs)))
        B_estim = 0
      }
      n_coefs = nrow(coefs)
      use_coefs = sample(n_coefs, size=B, replace=(n_coefs<B))
      coefs = coefs[use_coefs, ]
      simfunc = function(coefs, fit, xreg, n.ahead){
        fit$coefficients = coefs[seq(along=fit$coefficients)]
        fit$distrcoefs = coefs[length(fit$coefficients)+seq(along=fit$distrcoefs)]
        result = tsglm.sim(n=n.ahead, xreg=xreg[-(1:fit$n_obs), , drop=FALSE], fit=fit, n_start=0)$ts
        return(result)
      }
      futureobs = matrix(apply(coefs, 1, simfunc, fit=object, xreg=xreg, n.ahead=n.ahead), nrow=n.ahead, ncol=nrow(coefs))
      if(type %in% c("quantiles", "onesided")){
        quantiles = t(apply(futureobs, 1, quantile, probs=c(a, 1-a), type=1))
        lower = if(type=="onesided"){integer(n.ahead)}else{quantiles[, 1]}
        upper = quantiles[, 2]
        predint = cbind(lower=lower, upper=upper)
      }
      if(type=="shortest") predint = t(apply(futureobs, 1, function(x) largestdensityinterval(tabulate(x+1)/length(x), level=level_local)))
      predmed = apply(futureobs, 1, median)
    }
    colnames(predint) = c("lower", "upper")
    if(is.ts(object$ts)){ #use time series class if input time series has this class
      predint = ts(predint, start=tsp(object$ts)[2]+1/frequency(object$ts), frequency=frequency(object$ts))
      predmed = ts(predmed, start=tsp(object$ts)[2]+1/frequency(object$ts), frequency=frequency(object$ts))
    }
    # result = c(result, list(interval=predint, level=level, global=global, type=type, method=method, B=B, estim=estim, B_estim=B_estim, warning_messages=warning_messages, median=predmed))

  }
  if(n.ahead > 1) ret_obj = futureobs
  else ret_obj = pred
  return(ret_obj)
}





