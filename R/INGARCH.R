globalVariables(c("p", "q","P", "Q"))



naive_search = function(x, p = 0, q = 0, P = 0, Q = 0,
                        trace = T,
                        ic  = 'AIC', link = 'identity', distr = 'poisson', xreg = NULL){
  # search_matrix storage a INGARCH(p,q) model in a matrix[p,q]
  search_matrix = matrix(ncol = 4, nrow = 4)
  for(i in 0:3){
    for(j in 0:3){
      if(i == 0) past_obs_auto = NULL
      else past_obs_auto = 1:i
      if(j == 0) past_mean_auto = NULL
      else past_mean_auto = 1:j
      step_model =
        tryCatch(expr =
                   {tscount::tsglm(x,
                                   model = list(past_obs = past_obs_auto,
                                                past_mean = past_mean_auto),
                                   link = link, distr = distr , xreg = xreg[[1]]$xreg)},
                 error = function(err){
                   return(NA)
                 })
      search_matrix[i + 1, j + 1] = ifelse(is.na(get(ic)(step_model)) == F, get(ic)(step_model), NA)
      if(trace){
        cat(sprintf('INGARCH[%d, %d], %s:%.2f \n', i, j, ic, search_matrix[i+1,j+1]))
      }
    }
  }
  params = which(search_matrix == min(search_matrix), arr.ind = TRUE)
  params = params - 1
  return(params)
}


#' @importFrom fable ARIMA
arma_to_ingarch = function(x, p = 0, q = 0,
                           P = 0, Q = 0,
                           trace = T, ic,
                           xreg = NULL){
  arma_model =
    tibble::tibble(y_var = x,
                 time_index =
                   lubridate::make_date(year = 1:length(x)) |>
                   tsibble::yearquarter()) |>
    tsibble::as_tsibble(index = time_index) |>
    model(auto = ARIMA(y_var, ic = ic |> tolower(),
                       trace = trace))

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



ingarch_tscall = function(x, p = 0, q = 0, P = 0, Q = 0,
                          automatic = T, trace = T,
                          ic  = 'AIC', link = 'identity', distr = 'poisson', xreg = NULL,
                          algorithm = c('naive_search', 'arma_to_ingarch')){

  params = c(p,q)
  if(automatic == T){
    algorithm = match.arg(algorithm)
    params = get(algorithm)(x = x, p = 0, q = 0, P = 0, Q = 0,
                   trace = trace, ic  = 'AIC', link = 'identity',
                   distr = 'poisson', xreg = NULL)
  }
  params = clean_params(params)
  tscount_model = tscount::tsglm(x,
                                 model =
                                   list(past_obs = params[[1]],
                                        past_mean = params[[2]]),
                                 link = link,
                                 distr = distr,
                                 xreg = xreg[[1]]$xreg
                                 )
  params_ret = lapply(params, max) |> unlist()
  return(list(params = params_ret, tscount_model = tscount_model))
}


#' Estimate a INGARCH model
#'
#' Estimate Integer-valued Generalized Autoregressive Conditional Heteroscedasticity model
#' with Poisson or Negative Binomial distribution.
#' Also is provide a automatic parameter algorithm selection for the Autorregressive and Moving Avarege params
#'
#'
#' @param formula Model specification (see "Specials" section).
#' @param ic Character, can be 'aic' 'bic' or 'qic'. The information criterion used in selecting the model.
#' @param link Character, can be 'identity' or 'log' The link function used for the generalized model
#' @param distr Character, can be 'poisson' or 'nbinom'. The probabilty distribution used for the generalized model
#' @param algorithm Character, specifies the automatic parameter selection algorithm. Can be 'naive_search' or 'arma_to_ingarch'.
#' If 'naive_search' is selected, a search in a 4x4 matrix parameter space is performed, where the model to minimize the criterion value is selected.
#' If 'arma_to_ingarch' is selected, uses an auto_arma as the starting point for the selection algorithm.
#' The ‘arma_to_ingarch’ is the only one to perform a seasonal adjustment
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
#' or if it's omited, the automatic parameter selection algorithm is trigered (only for 'arma_to_ingarch' algorithm)
#' The automatic parameter selection algorithm gonna fit the best model based on the information criterion
#' }
#'
#'
#' \subsection{xreg}{
#' Exogenous regressors can be included in a INGARCH model without explicitly using the `xreg()` special.
#' Common exogenous regressor specials as specified in [`common_xregs`] can also be used.
#' These regressors are handled using [stats::model.frame()],
#' and so interactions and other functionality behaves similarly to [stats::lm()].
#'
#' The inclusion of a constant in the model follows the similar rules to [`stats::lm()`],
#' where including `1` will add a constant and `0` or `-1` will remove the constant.
#' If left out, the inclusion of a constant will be determined by minimising `ic`.
#' }
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
#'
#'
#' @return A model specification.
#' @examples
#' \donttest{
#' # Manual INGARCH specification
#' tsibbledata::aus_production |>
#'   fabletools::model(manual_ing = INGARCH(Beer ~ pq(1,1)))
#'
#' # Automatic INGARCH specification
#'  tsibbledata::aus_production |>
#' fabletools::model(auto_ing_naive =
#'                     INGARCH(Beer,
#'                             ic = 'aic',
#'                             trace = TRUE,
#'                            algorithm = 'naive_search'),
#'                   auto_ing_arm_ing =
#'                     INGARCH(Beer,
#'                             ic = 'aic',
#'                             trace = TRUE,
#'                             algorithm = 'arma_to_ingarch'))
#'                             }
#'
#'
#' @importFrom stats fitted
#' @importFrom stats var
#' @export
INGARCH = function(formula,
                    ic = c('aic', 'bic', 'qic'),
                    link = c('identity', 'log'),
                    distr = c('poisson', 'nbinom'),
                    algorithm = c('naive_search', 'arma_to_ingarch'),
                    trace = FALSE) {

  ic = match.arg(ic) |> toupper()
  link = match.arg(link)
  distr = match.arg(distr)

  model_INGARCH = new_model_class("INGARCH",
                                 train = train_INGARCH,
                                 specials = specials_INGARCH,
                                 check = function(.data) {
                                   if (!tsibble::is_regular(.data)) stop("Data must be regular")
                                 }
  )

  # Return a model definition which stores the user's model specification
  new_model_definition(model_INGARCH, {{formula}},
                       ic = ic,
                       link = link,
                       distr = distr,
                       trace = trace)
}

specials_INGARCH = new_specials(
  pq = function(p = 'not choosen', q = 'not choosen',
                 p_init = 2, q_init = 2,
                 fixed = list()) {

    if(!all(grepl("^(ma|ar)\\d+", names(fixed)))){
      abort("The 'fixed' coefficients for pq() must begin with ar or ma, followed by a lag number.")
    }
    as.list(environment())
  },
  PQ = function(P = 'not choosen', Q = 'not choosen',
                P_init = 2, Q_init = 2,
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
  .required_specials = c("pq", "PQ"),
  .xreg_specials = names(common_xregs)
)

train_INGARCH = function(.data, specials, ic,
                         link, distr,
                         trace, ...){
  mv = tsibble::measured_vars(.data)
  if(length(mv) > 1) stop("INGARCH is a univariate model.")
  y = .data[[mv]]

  automatic = F
  if(specials$pq[[1]]$p |>
     is.character()) automatic = T
  xreg = specials$xreg
  tsglm_model = ingarch_tscall(x = y,
                               ic = ic,
                               p = specials$pq[[1]]$p,
                               q = specials$pq[[1]]$q,
                               link = link,
                               distr = distr,
                               automatic = automatic,
                               algorithm = algorithm,
                               trace = trace,
                               xreg = xreg)

  # Compute fitted values and residuals
  fit = tsglm_model$tscount_model |> fitted()
  e = y - fit

  # Create S3 model object
  # It should be small, but contain everything needed for methods below
  structure(
    list(
      coef = tsglm_model$params,
      tsmodel = tsglm_model$tscount_model,
      distr = distr,
      link = link,
      n = length(y),
      y_name = mv,
      fitted = fit,
      residuals = e,
      sigma2 = var(e, na.rm = TRUE)
    ),
    class = "INGARCH"
  )
}

#' @export
model_sum.INGARCH = function(x){
  print(x$coef)
  if(is.na(x$tsmodel$xreg[1]) & length(x$coef) <= 2) out = sprintf("INGARCH(%i, %i)", x$coef[1],x$coef[2])
  else if(is.na(x$tsmodel$xreg[1]) & length(x$coef) > 2) out = sprintf("Seasonal INGARCH(%i, %i)(%i,%i)[%i]", x$coef[1],x$coef[2], x$coef[3], x$coef[4], x$coef[5])
  else if(is.na(x$tsmodel$xreg[1]) == F & length(x$coef) <= 2) out = sprintf("INGARCH(%i, %i) w/ covariates", x$coef[1],x$coef[2])
  else out = sprintf("Seasonal INGARCH(%i, %i)[%i,%i][%i] w/ covariates", x$coef[1],x$coef[2], x$coef[3], x$coef[4], x$coef[5])
  out
}

#' @export
report.INGARCH = function(object, ...){
  model_sum_obj = object$tsmodel |>
    summary()

    rep_model = model_sum_obj$coefficients |>
    as.data.frame() |>
    tibble::rownames_to_column() |>
    dplyr::rename(term = 1) |>
    tidyr::pivot_longer(-term, names_to = 'statistic') |>
    tidyr::pivot_wider(names_from = term) |>
    dplyr::filter(statistic == 'Estimate' | statistic == 'Std.Error')

  cat('\n')
  cat(sprintf("%s INGARCH(%i, %i) w/ %s link", object$distr, object$coef[1], object$coef[2], object$link))
  cat('\n')
  rep_model |> print()
  cat('\n')
  cat(paste('log likelihood='), model_sum_obj$logLik, sep = '')
  cat('\n')
  cat(paste('AIC='), model_sum_obj$AIC, sep = '')
  cat('\n')
  cat(paste('BIC='), model_sum_obj$BIC, sep = '')
  cat('\n')
  cat(paste('QIC='), model_sum_obj$QIC, sep = '')
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
#'   fabletools::model(manual_ing = INGARCH(Beer ~ pq(1,1))) |>
#'   dplyr::select(manual_ing) |>
#'   fabletools::tidy()
#'
#' @import stringr
#' @export
tidy.INGARCH = function(x, ...){
  model_summary = x$tsmodel |>
    summary() |>
    {\(x)x$coefficients}() |>
    tibble::rownames_to_column() |>
    tibble::as_tibble() |>
    dplyr::rename(term = 1,
                  estimate = 2,
                  std.error = 3) |>
    dplyr::mutate(term =
                    stringr::str_replace(term, '(Intercept)', 'constant'),
                  term =
                    stringr::str_replace(term, 'beta', 'ar'),
                  term =
                    stringr::str_replace(term, 'alpha', 'ma'))


  return(model_summary)
}


#' Extract fitted values from a fable model
#'
#' Extracts the fitted values.
#'
#' @inheritParams forecast.INGARCH
#'
#' @return A vector of fitted values.
#'
#' @examples
#' tsibbledata::aus_production |>
#'   fabletools::model(manual_ing = INGARCH(Beer ~ pq(1,1))) |>
#'   dplyr::select(manual_ing) |>
#'   fitted()
#'
#' @import fabletools
#' @export
fitted.INGARCH = function(object, ...){
  object$fitted
}


#' Extract residuals from a fable model
#'
#' Extracts the residuals.
#'
#' @inheritParams forecast.INGARCH
#'
#' @return A vector of fitted residuals.
#'
#' @examples
#' tsibbledata::aus_production |>
#'   fabletools::model(manual_ing = INGARCH(Beer ~ pq(1,1) + PQ(1,1))) |>
#'   dplyr::select(manual_ing) |>
#'   residuals()
#'
#' @import fabletools
#' @export
residuals.INGARCH = function(object, ...){
  object$residuals
}


#' Glance a INGARCH model
#'
#' Construct a single row summary of the INGARCH model.
#'
#' @format A data frame with 1 row, with columns:
#' \describe{
#'   \item{sigma2}{The unbiased variance of residuals. Calculated as `sum(residuals^2) / (num_observations - num_pararameters + 1)`}
#'   \item{log_lik}{The log-likelihood}
#'   \item{AIC}{Akaike information criterion}
#'   \item{BIC}{Bayesian information criterion}
#' }
#'
#' @inheritParams generics::glance
#'
#' @return A one row tibble summarising the model's fit.
#'
#' @examples
#' tsibbledata::aus_production |>
#'   fabletools::model(manual_ing = INGARCH(Beer ~ pq(1,1))) |>
#'   dplyr::select(manual_ing) |>
#'   glance()
#'
#' @importFrom stats AIC BIC
#' @export
glance.INGARCH = function(x, ...){
  tibble::tibble(sigma2 = sum(x$residuals^2) / (x$n - sum(x$coef) + 1),
                 log_lik = x$tsmodel$logLik,
                 AIC = x$tsmodel |> AIC(),
                 BIC = x$tsmodel |> BIC(),
                 )
}

#' Forecast a model from the fable package
#'
#' Produces forecasts from a trained model.
#'
#' Predict future observations based on a fitted GLM-type model for time series of counts.
#' For 1 step ahead, it returns parametric forecast, based on the 'distr' param especified distribution,
#' for multiples steps forecast, the distribution is not know analytically, so it uses a parametric bootstrap
#'
#' @inheritParams generics::forecast
#' @param new_data Tsibble, it has to contains the time points and exogenous regressors to produce forecasts for.
#'
#' @importFrom stats formula residuals
#'
#' @return A list of forecasts.
#'
#' @examples
#'  # 1 step ahead parametric forecast
#' tsibbledata::aus_production |>
#'   fabletools::model(manual_ing = INGARCH(Beer ~ pq(1,1) + PQ(1,1))) |>
#'   dplyr::select(manual_ing) |>
#'   fabletools::forecast(h = 1)
#'
#' # Multiples steap ahead parametric bootstrap forecast
#' tsibbledata::aus_production |>
#'   fabletools::model(manual_ing = INGARCH(Beer ~ pq(1,1) + PQ(1,1))) |>
#'   dplyr::select(manual_ing) |>
#'   fabletools::forecast(h = 4)
#'
#' @import fabletools
#' @export
forecast.INGARCH = function(object, new_data,...){
  h = NROW(new_data)
  values = predict_call_tsglm(object = object$tsmodel,
          n.ahead = h, ...)
  if(h == 1) ret_values = distributional::dist_poisson(values)
  else ret_values = matrix_to_list(values) |> distributional::dist_sample()

  ret_values
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


