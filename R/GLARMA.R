globalVariables(c("p", "q","P", "Q"))


arma_to_glarma = function(x, p = 0, q = 0,
                           P = 0, Q = 0,
                           trace = T, ic  = 'aic',
                           link = 'identity', distr = 'poisson',
                           xreg = NULL){
  arma_model =
    tibble::tibble(y_var = x,
                   time_index =
                     lubridate::make_date(year = 1:length(x)) |>
                     tsibble::yearquarter()) |>
    tsibble::as_tsibble(index = time_index) |>
    model(auto = ARIMA(y_var, ic = ic |> tolower()))

  model_report = arma_model$auto[1] |> as.character() |>
    stringr::str_extract("\\[(\\d+)\\]") |>
    stringr::str_extract("\\d+") |>
    as.numeric()

  arma_model = arma_model |>
    tidy() |>
    dplyr::mutate(term = term |> stringr::str_remove_all("[0-9]"),
                  id = 1)


  test_vector = c('ar', 'ma', 'sar', 'sma')
  test_vector = test_vector[!(c('ar', 'ma', 'sar', 'sma') %in% arma_model$term)]

  arma_model = tibble::tibble(term = test_vector,
                              id = rep(0, length(test_vector))) |>
    dplyr::bind_rows(arma_model) |>
    dplyr::group_by(term) |>
    dplyr::summarise(order = sum(id)) |>
    dplyr::bind_rows(
      tibble::tibble(term = 'seas_index',
                     order = model_report)
    )

  params =
    list(pq =
           arma_model |>
           dplyr::filter(term == 'ar' | term == 'ma') |>
           dplyr::pull(order),
         PQ =
           arma_model |>
           dplyr::filter(term == 'sar' | term == 'sma' |  term == 'seas_index') |>
           dplyr::pull(order))

  params =
    c(params$pq[1],
      params$pq[2],
      params$PQ[1],
      params$PQ[2],
      params$PQ[3])
  return(params)
}

#' @importFrom glarma glarma
glarma_call = function(y = y, ic = ic, type = 'Poi',
                       p = 0, q = 0,
                       method = 'NR', residuals = 'Pearson',
                       automatic = automatic, trace = T, xreg = NULL){

  if(is.null(xreg) == F)
      xreg = cbind(rep(1, length(y)), xreg)
     else
       xreg = matrix(rep(1, length(y)), ncol = 1)
  colnames(xreg)[1] = 'Intercept'
  params = c(p, q)
  if(automatic){
    params = arma_to_glarma(x = y, p = 0, q = 0,
                   P = 0, Q = 0,
                   trace = trace, ic = ic,
                   xreg = xreg)
  }
  params_cleaned = clean_params(params)
  gl_model = tryCatch(expr = glarma::glarma(y = y,
                               X = xreg,
                               thetaLags = params_cleaned$p,
                               phiLags = params_cleaned$q,
                               type = type,
                               method = method,
                               residuals = residuals) ,
                      error = function(err){
                        if(automatic) stop('The automatic selection algorithm failed to find stable parameters')
                        else stop('An error has occurred in the model estimation, try changing the model parameters')
                        return(NA)
                      })

return(list(params = params, gl_model = gl_model))
}


#' Estimate a GLARMA model
#'
#' Estimate Generalized Linear Autoregressive Moving Average  model
#' with Poisson or Negative Binomial distribution.
#' Also is provide a automatic parameter algorithm selection for the Autorregressive and Moving Average params
#'
#'
#' @param formula Model specification (see "Specials" section).
#' @param ic Character, can be 'AIC','BIC'. The information criterion used in selecting the model.
#' @param distr Character, can be 'poisson' or 'nbinom'. The probabilty distribution used for the generalized model
#' @param method Character, can be 'FS' (Fisher scoring) or 'NR' (Newton-Raphson). The method of iteration to be used
#' @param residuals Character, can be 'Pearson' or 'Score'. The type of residuals to be used
#' @param trace Logical. If the automatic parameter algorithm is runnig, print the path to the best model estimation
#'
#' @section Specials:
#'
#' \subsection{pq}{
#' pq defines the non-seasonal autoregressive and moving avarages terms,
#' it can be define by the user,
#' or if it's omited, the automatic parameter selection algorithm is trigered
#' The automatic parameter selection algorithm gonna fit the best model based on the information criterion
#' }
#'
#' \subsection{PQ}{
#' PQ defines the seasonal autoregressive and moving avarages terms,
#' it can be define by the user,
#' or if it's omited, the automatic parameter selection algorithm is trigered (only for 'arma_to_GLARMA' algorithm)
#' The automatic parameter selection algorithm gonna fit the best model based on the information criterion
#' }
#'
#' \subsection{xreg}{
#' Exogenous regressors can be included in an GLARMA model without explicitly using the `xreg()` special.
#' Common exogenous regressor specials as specified in [`common_xregs`] can also be used.
#' These regressors are handled using [stats::model.frame()],
#' and so interactions and other functionality behaves similarly to [stats::lm()].
#'
#' The inclusion of a constant in the model follows the similar rules to [`stats::lm()`],
#' where including `1` will add a constant and `0` or `-1` will remove the constant.
#' If left out, the inclusion of a constant will be determined by minimising `ic`.
#'
#' If a xreg is provided, the model forecast is not avaliable
#'
#' \preformatted{
#' xreg(..., fixed = list())
#' }
#'
#' \tabular{ll}{
#'   `...`      \tab Bare expressions for the exogenous regressors (such as `log(x)`)\cr
#'   `fixed`    \tab A named list of fixed parameters for coefficients. The names identify the coefficient, and should match the name of the regressor. For example, `fixed = list(constant = 20)`.
#' }
#' }
#'
#' @examples
#' # Manual GLARMA specification
#' tsibbledata::aus_production |>
#'   fabletools::model(manual_gla = GLARMA(Beer ~ pq(1,0)))
#'
#' # Automatic GLARMA specification
#' tsibbledata::aus_production |>
#'   fabletools::model(auto_gla = GLARMA(Beer, ic = 'aic'))
#'
#'
#' @return A model specification.
#' @importFrom stats fitted
#' @importFrom stats var
#' @export
GLARMA = function(formula,
                   ic = c('aic','bic'),
                   distr = c('Poi', 'NegBin'),
                   method = c('FS', 'NR'),
                   residuals = c('Pearson', 'Score'),
                   trace = FALSE) {

  ic = match.arg(ic)
  distr = match.arg(distr)
  method = match.arg(method)
  residuals = match.arg(residuals)


  model_GLARMA = new_model_class("GLARMA",
                                  train = train_GLARMA,
                                  specials = specials_GLARMA,
                                  check = function(.data) {
                                    if (!tsibble::is_regular(.data)) stop("Data must be regular")
                                  }
  )

  # Return a model definition which stores the user's model specification
  new_model_definition(model_GLARMA, {{formula}},
                       ic = ic,
                       distr = distr,
                       method = method,
                       residuals = residuals,
                       trace = trace)
}

specials_GLARMA = new_specials(
  pq = function(p = 'not choosen', q = 'not choosen',
                p_init = 2, q_init = 2,
                fixed = list()) {

    if(!all(grepl("^(ma|ar)\\d+", names(fixed)))){
      abort("The 'fixed' coefficients for pq() must begin with ar or ma, followed by a lag number.")
    }
    as.list(environment())
  },
  common_xregs,
  xreg = function(..., fixed = list()) {
    dots = enexprs(...)
    env = map(enquos(...), get_env)
    env[map_lgl(env, compose(is_empty, env_parents))] = NULL
    env = if (!is_empty(env)) get_env(env[[1]]) else base_env()

    constants = map_lgl(dots, inherits, "numeric")
    constant_forced = any(map_lgl(dots[constants], `%in%`, 1))

    model_formula = new_formula(
      lhs = NULL,
      rhs = reduce(dots, function(.x, .y) call2("+", .x, .y))
    )

    # Mask user defined lag to retain history when forecasting
    env = env_bury(env, lag = lag)

    xreg = model.frame(model_formula, data = env, na.action = stats::na.pass)
    tm = terms(xreg)
    constant = as.logical(tm %@% "intercept")
    xreg = model.matrix(tm, xreg)

    if (constant) {
      xreg = xreg[, -1, drop = FALSE]
    }

    list(
      constant = if (!constant || constant_forced) constant else c(TRUE, FALSE),
      xreg = if (NCOL(xreg) == 0) NULL else xreg,
      fixed = fixed
    )
  },
  .required_specials = "pq",
  .xreg_specials = names(common_xregs)
)


train_GLARMA = function(.data, specials, ic,
                         distr, residuals,
                         method, trace, ...){
  mv = tsibble::measured_vars(.data)
  if(length(mv) > 1) stop("GLARMA is a univariate model.")
  y = .data[[mv]]


  automatic = F
  if(specials$pq[[1]]$p |>
     is.character()) automatic = T
  xreg = specials$xreg
  glarma_model = glarma_call(y = y,
                               p = specials$pq[[1]]$p,
                               q = specials$pq[[1]]$q,
                               ic = ic,
                               type = distr,
                               method = method,
                               residuals = residuals,
                               automatic = automatic,
                               trace = trace,
                               xreg = xreg[[1]]$xreg)

  # Compute fitted values and residuals
   fit = glarma_model$gl_model |> fitted()
   e = y - fit

  # Create S3 model object
  # It should be small, but contain everything needed for methods below
   structure(
     list(
       coef = glarma_model$params,
       gl_model = glarma_model$gl_model,
       distr = distr,
       residuals_type = residuals,
       method = method,
       n = length(y),
       y_name = mv,
       fitted = fit,
       residuals = e,
       sigma2 = var(e, na.rm = TRUE)
     ),
    class = "GLARMA"
  )
}

#' @export
model_sum.GLARMA = function(x){
  if(x$gl_model$r == 1) out = sprintf("GLARMA(%i, %i)", x$coef[1], x$coef[2] )
  else out = sprintf("GLARMA(%i, %i) w/ covariates", x$coef[1], x$coef[2])
  out
}



#' @export
report.GLARMA = function(object, ...){
  if(object$distr == 'Poi')
    distr_x = 'Poisson'
  else distr_x = 'Negative Binomial'
  x_tidy = object |> tidy()
  cat('\n')
  cat(sprintf("%s GLARMA(%i, %i)", distr_x, object$coef[1], object$coef[2]))
  cat('\n')
  x_tidy |> dplyr::filter(statistic == 'estimate' | statistic == 'std_error') |> print()
  cat('\n')
  cat(paste('log likelihood='), object$gl_model$logLik, sep = '')
  cat('\n')
  cat(paste('AIC='), object$gl_model$aic, sep = '')


}


#' Tidy a fable model
#'
#' Returns the coefficients from the model in a `tibble` format.
#'
#' @inheritParams generics::tidy
#'
#' @return The model's coefficients in a `tibble`.
#'
#' @examples
#' tsibbledata::aus_production |>
#'   fabletools::model(manual_gla = GLARMA(Beer ~ pq(1,0))) |>
#'   dplyr::select(manual_gla) |>
#'   fabletools::tidy()
#'
#' @export
tidy.GLARMA = function(x, ...){
  sum_model = x$gl_model |> summary()
  out = sum_model[14 : (14 + x$gl_model$pq)] |>
    unlist() |>
    tibble::as_tibble()
  if(x$coef[1] == 0) ar_param = NULL
   else ar_param = paste('ar_', rep(1:x$coef[1], 4) |> sort(), sep = '')
  if(x$coef[2] == 0) ma_param = NULL
    else ma_param = paste('ma_', rep(1:x$coef[2], 4) |> sort(), sep = '')
  if(nrow(sum_model$coefficients1) > 1){
  coef1 =
    c(
      rep('intercept', 4),
      paste('xreg_',
            rep(1:(nrow(sum_model$coefficients1) - 1), 4) |>
              sort(), sep = '')
      )
  }
  else coef1 = rep('intercept', 4)
  out = out |>
    dplyr::mutate(term =
                    c(coef1, ar_param, ma_param),
                  statistic =
                    rep(c('estimate', 'std_error', 'z_ratio', 'p_value'),
                        nrow(sum_model$coefficients1)+x$gl_model$pq)
    ) |>
    tidyr::pivot_wider(values_from = value, names_from = term)


}


#' Extract fitted values from a fable model
#'
#' Extracts the fitted values.
#'
#' @inheritParams forecast.GLARMA
#'
#' @return A vector of fitted values.
#'
#' @examples
#' tsibbledata::aus_production |>
#'   fabletools::model(manual_gla = GLARMA(Beer ~ pq(1,0))) |>
#'   dplyr::select(manual_gla) |>
#'   fitted()
#' @export
fitted.GLARMA = function(object, ...){
  object$fitted
}


#' Extract residuals from a fable model
#'
#' Extracts the residuals.
#'
#' @inheritParams forecast.GLARMA
#'
#' @return A vector of fitted residuals.
#'
#' @examples
#' tsibbledata::aus_production |>
#'   fabletools::model(manual_gla = GLARMA(Beer ~ pq(1,0))) |>
#'   dplyr::select(manual_gla) |>
#'   residuals()
#' @export
residuals.GLARMA = function(object, ...){
  object$residuals
}


#' Glance a GLARMA model
#'
#' Construct a single row summary of the GLARMA model.
#'
#' @format A data frame with 1 row, with columns:
#' \describe{
#'   \item{sigma2}{The unbiased variance of residuals. Calculated as `sum(residuals^2) / (num_observations - num_pararameters + 1)`}
#'   \item{log_lik}{The log-likelihood}
#'   \item{AIC}{Akaike information criterion}
#' }
#'
#' @inheritParams generics::glance
#'
#' @return A one row tibble summarising the model's fit.
#'
#' @examples
#' tsibbledata::aus_production |>
#'   fabletools::model(manual_ing = GLARMA(Beer ~ pq(1,1))) |>
#'   dplyr::select(manual_ing) |>
#'   glance()
#' @export
glance.GLARMA = function(x, ...){
  tibble::tibble(sigma2 = x$sigma2,
                 log_lik = x$gl_model$logLik,
                 AIC = x$gl_model$aic,
  )
}




#' Forecast a model from the fable package
#'
#' Produces forecasts from a trained model.
#'
#' Predict future observations based on a fitted GLM-type model for time series of counts.
#' Futher informations about the forecast method can be obtained typing ?glarma::forecast
#'
#' @inheritParams generics::forecast
#' @param new_data Tsibble, it has to contains the time points and exogenous regressors to produce forecasts for.
#'
#' @importFrom stats formula residuals
#'
#' @return A list of forecasts.
#'
#' @examples
#' tsibbledata::aus_production |>
#'   fabletools::model(manual_gla = GLARMA(Beer ~ pq(1,0))) |>
#'   dplyr::select(manual_gla) |>
#'   fabletools::forecast(h = 2)
#' @export
forecast.GLARMA = function(object, new_data,...){
  h = NROW(new_data)
  newdata = matrix(rep(1, h), ncol = 1)
  newoffset  = matrix(rep(0, h), ncol = 1)
  newm = matrix(rep(1, h), ncol = 1)
  values = glarma::forecast(object = object$gl_model,
                            n.ahead = h,
                            newdata = newdata,
                            newoffset = newoffset)
  distributional::dist_poisson(values$Y)
}

clean_params = function(params_vector){
  if(is.na(params_vector[5]) == F){
    if(params_vector[1] == 0) p = NULL
    else p = 1:params_vector[1]
    if(params_vector[2] == 0) q = NULL
    else q = 1:params_vector[2]
    params =
      list(p = c(p, params_vector[5]:(params_vector[5] + params_vector[3] - 1)),
           q = c(q, params_vector[5]:(params_vector[5] + params_vector[4] - 1))
      )
  }
  else{
    if(params_vector[1] == 0) p = NULL
    else p = 1:params_vector[1]
    if(params_vector[2] == 0) q = NULL
    else q = 1:params_vector[2]
    params = list(p = p,
                  q = q)
  }
  return(params)
}



