
inv_logit <- function(x) 1 / (1 + exp(-x))

#' The Hyper Chain Binomial distribution
#'
#' Probability mass function, expected value, and random generation, for the hyper chain binomial distribution,
#' with parameters s0, prob, i0, and s0_obs for the number of observed infected
#' cases in a population of size s0 with s0_obs observed individuals.
#'
#' @param x numeric vector of the number of infected.
#' @param s0 the number of initial susceptibles.
#' @param prob the transmission, or the per person risk of infection by an infected.
#' @param s0_obs The number of observed individuals.
#' @param i0 the number of primary cases.
#'
#' @export
dcbhyper <- function(x, s0, prob, s0_obs, i0 = 1){

  # Combine input to matrix, to expand/recycle input data.
  inp <- as.matrix(cbind(x, s0, prob, s0_obs, i0))

  stopifnot(any(inp[,'s0_obs'] <= inp[,'s0']))

  res <- numeric(nrow(inp))

  for (ii in 1:nrow(inp)){

    true_escaped <- inp[ii, 's0'] - 0:inp[ii, 's0']

    dd <- dchainbinom(x=0:inp[ii, 's0'], s0 = inp[ii, 's0'], sar = inp[ii, 'prob'], i0 = inp[ii, 'i0'])

    # x: I_obs
    # m: I
    # n = S0 - I
    # k: S0_obs
    res[ii] <- sum(stats::dhyper(x = inp[ii, 'x'], m = 0:inp[ii, 's0'], n = true_escaped, k = inp[ii, 's0_obs']) * dd)

  }

  return(res)
}


#' @rdname dcbhyper
#' @export
ecbhyper <- function(s0, prob, s0_obs, i0 = 1){

  stopifnot(is.numeric(s0) | is.logical(s0),
            is.numeric(prob) | is.logical(prob),
            is.numeric(i0) | is.logical(i0),
            is.numeric(s0_obs) | is.logical(s0_obs),
            all(prob >= 0, na.rm = TRUE),
            all(prob <= 1, na.rm = TRUE),
            all(s0 >= 0, na.rm = TRUE),
            all(i0 >= 1, na.rm = TRUE),
            all(s0_obs >= 1, na.rm = TRUE))

  # Coerce to numericals, in case of logicals being provided.
  s0 <- as.numeric(s0)
  prob <- as.numeric(prob)
  i0 <- as.numeric(i0)
  s0_obs <- as.numeric(s0_obs)

  # Combine input to matrix, to expand/recycle input data.
  inp <- as.matrix(cbind(s0, prob, s0_obs, i0))

  n <- nrow(inp)
  res <- numeric(n)

  for (ii in 1:n){

    if (any(is.na(inp[ii,]))){
      res[ii] <- NA
      next
    }
    xx <- 0:inp[ii, 's0_obs']

    pp <- dcbhyper(x  = xx, s0 = inp[ii, 's0'], i0 = inp[ii, 'i0'], prob = inp[ii, 'prob'], s0_obs = inp[ii, 's0_obs'])
    res[ii] <- sum(xx * pp)
  }

  return(res)

}





# The negative chain binomial log-likelihood, for use with optim.
negloglok_cbhyper <- function(prob, infected, s0, s0_obs, i0, transform_inv_logit = TRUE){

  # The parameter should be transformed when using this function for optimization (estimation),
  # but not when computing the hessian (for standard errors) or when the parameter
  # has already been transformed (in regression modelling).
  if (transform_inv_logit){
    prob <- inv_logit(prob)
  } else {
    if (any(prob < 0) | any(prob > 1)){
      return(Inf)
    }
  }

  nll <- -sum(log(dcbhyper(x = infected, s0 = s0, s0_obs = s0_obs, prob = prob, i0 = i0)), na.rm = TRUE)

  if (is.nan(nll)){
    return(Inf)
  } else {
    return(nll)
  }
}



#' Estimate Secondary Attack Rate of the Hyper Chain Binomial Model
#'
#' Given data on the number of infected after a number of generation, initial
#' number of susceptible, and initial number of infected, estimate the
#' secondary attack rate (SAR) using maximum likelihood.
#'
#' @param infected numeric.
#' @param s0 numeric.
#' @param s0_obs numeric.
#' @param i0 numeric.
#' @param se logical. If TRUE (default), the standard error is computed.
#'
#'
#' @returns A list of class `sar2` with the following components:
#' * `sar_hat` The point estimate of the secondary attack rate.
#' * `se` Standard error of the estimate (if se = TRUE).
#' * `loglikelihood` the log likelihood value at the point estimate.
#' * `data` the input data.
#'
#' @seealso
#' [confint.sar2()] for calculating confidence intervals.
#'
#'@export
estimate_sar_cbhyper <- function(infected, s0, s0_obs, i0 = 1, se = TRUE){

  stopifnot(is.numeric(infected) | is.logical(infected),
            is.numeric(s0) | is.logical(s0),
            is.numeric(i0) | is.logical(i0),
            #is.numeric(generations) | is.logical(generations),
            all(infected <= s0, na.rm = TRUE),
            all(infected >= 0, na.rm = TRUE),
            all(i0 > 0, na.rm = TRUE),
            #all(!is.na(generations)),
            is.logical(se),
            length(se) == 1)


  sar_init <- 0.5
  optim_res <- stats::optim(par = sar_init,
                     fn = negloglok_cbhyper, method = 'L-BFGS-B',
                     hessian = FALSE,
                     infected = infected,
                     s0 = s0,
                     s0_obs = s0_obs,
                     i0 = i0)


  sar_hat <- inv_logit(optim_res$par)

  if (se){
    he <- numDeriv::hessian(negloglok_cbhyper, x = sar_hat,
                            infected = infected,
                            s0 = s0,
                            s0_obs = s0_obs,
                            i0 = i0,
                            #generations = generations,
                            transform_inv_logit = FALSE)

    sar_se <- sqrt(as.numeric(solve(he)))
  } else {
    sar_se <- NA
  }


  inp <- list(infected = infected,
              s0 = s0,
              s0_obs = s0_obs,
              i0 = i0)

  res <- list(sar_hat = sar_hat,
              se = sar_se,
              loglikelihood = -optim_res$value,
              data = inp)

  class(res) <- 'sar2'

  return(res)
}




obj_ci_wilks2 <- function(x, infected, s0, i0, s0_obs, max_loglik, critical_value){

  # Log likelihood at the parameter value x.
  x_loglik <-  -negloglok_cbhyper(prob = x, infected = infected, s0 = s0, i0 = i0, s0_obs = s0_obs,
                                  #generations = generations,
                                  transform_inv_logit = FALSE)

  # -2 log-likelihood ratio.
  nll2 <-  -2 * (x_loglik - max_loglik)

  critical_value - nll2

}


# Find the search interval for use with uniroot.
find_intervall_upr_hcb <- function(sh){

  if (sh >= 0.99){
    search_int <- c(0, 1)
  } else {
    search_int <- c(sh, 1)
  }

  return(search_int)

}

find_intervall_lwr_hcb <- function(sh){

  if (sh <= 0.01){
    search_int <- c(0, 1)
  } else {
    search_int <- c(0, sh)
  }

  return(search_int)

}


#' Confidence intervals for sar2 Object.
#'
#' @param object a sar2 object.
#' @param parm Ignored.
#' @param level Default is 0.95, for 95% confidence intervals.
#' @param method Either 'chisq'(default) or 'normal'.
#' @param ... other arguments. Ignored.
#'
#' @returns A numeric of length 2 with the lower and upper end of the confidence interval.
#'
#'@export
confint.sar2 <- function(object, parm = NULL, level = 0.95, method = 'chisq', ...){

  stopifnot(method %in% c('chisq', 'normal'))

  if (method == 'normal'){

    # Compute Standard error if needed.
    if (is.na(object$se)){

      he <- numDeriv::hessian(negloglok_cbhyper, x = object$sar_hat,
                              infected = object$data$infected,
                              s0 = object$data$s0,
                              i0 = object$data$i0,
                              s0_obs = object$data$s0_obs,
                              #generations = object$data$generations,
                              transform_inv_logit = FALSE)
      sar_se <- sqrt(as.numeric(solve(he)))

    } else{
      sar_se <- object$se
    }

    # compute CI using normal approximation.
    ci_lwr <- max(object$sar_hat + (sar_se * stats::qnorm((1-level)/2, lower.tail = TRUE)), 0)
    ci_upr <- min(object$sar_hat + (sar_se * stats::qnorm((1-level)/2, lower.tail = FALSE)) , 1)

  } else if (method == 'chisq'){

    plwr <- (1-level)/2
    pupr <- 1 - plwr


    critical_value_lower <- stats::qchisq(plwr, df = 1, lower.tail = FALSE)

    uniroot_res_lwr <- stats::uniroot(f = obj_ci_wilks2, interval = find_intervall_lwr_hcb(object$sar_hat),
                                     infected = object$data$infected, s0 = object$data$s0,
                                     i0 = object$data$i0,
                                     s0_obs = object$data$s0_obs,
                                     #generations = object$data$generations,
                                     max_loglik = object$loglikelihood,
                                     critical_value = critical_value_lower,
                                     tol = 0.0000001)

    critical_value_upper <- stats::qchisq(pupr, df = 1, lower.tail = TRUE)

    uniroot_res_upr <- stats::uniroot(f = obj_ci_wilks2, interval = find_intervall_upr_hcb(object$sar_hat),
                               infected = object$data$infected, s0 = object$data$s0,
                               i0 = object$data$i0,
                               s0_obs = object$data$s0_obs,
                               #generations = object$data$generations,
                               max_loglik = object$loglikelihood,
                               critical_value = critical_value_upper,
                               tol = 0.0000001)

    ci_lwr <- min(object$sar_hat, uniroot_res_lwr$root)
    ci_upr <- max(object$sar_hat, uniroot_res_upr$root)

  }


  cn <- sprintf('%.1f %%', c(100 * (1-level)/2, 100 * (1-((1-level)/2)) ))

  res <- c(ci_lwr, ci_upr)
  names(res) <- cn
  return(res)
}


