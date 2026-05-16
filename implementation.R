library(survival)
library(icenReg)

# Helper functions -------------------------------------------------------------

#' Generate lower and upper limits for interval-censored data
#' 
#' @param exit Survival time at exit in days
#' @param width Interval size in days for each individual
#' @return A data frame with two columns: lower and upper censoring limit
generate_intervals <- function(exit, width) {
  G <- sample(1:width, length(exit), replace = TRUE)
  L <- pmax(0, exit - G)
  U <- exit + width - G
  df <- data.frame("lower" = L, "upper" = U, "obs" = U, "mid" = U - (U - L)/2)
  return(df)
}


#' Generate the data set
#' 
#' @param n Number of observations
#' @param pars A named list of dist. parameter values
#' @param width Vector of interval widths
#' @return A list of generated data frames with interval-censored times
generate_df <- function(n, pars, width) {
  exit <- ceiling(rweibull(n, pars[["shape"]], pars[["scale"]]))
  df.intervals <- lapply(width, function(x)
    cbind(exit, generate_intervals(exit, x)))
  return(df.intervals)
}


# Transform parameters from survreg to pweibull parametrisation
transform_pars <- function(icoefs) {
  return(c(1/exp(icoefs[[2]]), exp(icoefs[[1]])))
}


#' Estimate survival for different time points and distribution parameters
#' 
#' @param D Time points for which to estimate survival
#' @param df Data frame containing interval-censored survival times
#' @return A list of arrays of survival and parameter estimates
get_estimates <- function(D, df) {
  status <- rep(1, nrow(df))
  
  # Fit non-parametric models
  KM.mid <- survfit(Surv(mid, status) ~ 1, data = df)
  KM.obs <- survfit(Surv(obs, status) ~ 1, data = df)
  NPMLE <- ic_np(df[, c("lower", "upper")])
  scurv <- getSCurves(NPMLE)
  
  # Fit parametric models
  mid <- survreg(Surv(mid, status) ~ 1, data = df)
  ign <- survreg(Surv(obs, status) ~ 1, data = df)
  ic <- survreg(Surv(lower, upper, type = "interval2") ~ 1, data = df)
  
  # Estimates from model fit differ from base R parametrisation, 
  # thus we need to transform them.
  mid.pars <- transform_pars(mid$icoef)
  ign.pars <- transform_pars(ign$icoef)
  ic.pars <- transform_pars(ic$icoef)
  
  # Estimate and save probabilities into an array
  pD <- array(c(
    summary(KM.mid, times = D)[[6]],
    summary(KM.obs, times = D)[[6]],
    sapply(D, function(d) 
      scurv$S_curves$baseline[scurv$Tbull_ints[, 2] >= d][1]),
    1 - pweibull(D, mid.pars[[1]], mid.pars[[2]]),
    1 - pweibull(D, ign.pars[[1]], ign.pars[[2]]),
    1 - pweibull(D, ic.pars[[1]], ic.pars[[2]])
  ),
  dim = c(length(D), 6), 
  dimnames = list(D, c("KM.MID", "KM.IGN", "NPMLE", "MID", "IGN", "IC"))
  )
  
  # Save parameters into an array.
  pars <- array(c(mid.pars, ign.pars, ic.pars), 
                dim = c(2, 3),
                dimnames = list(c("shape", "scale"), c("MID", "IGN", "IC")))
  
  return(list(pD, pars))
}


#' Repeated simulation and parameter estimation
#' 
#' @param iter Number of simulated data sets
#' @param days Vector of survival times in days to estimate
#' @param df.pars Named list of parameters to generate the data sets
#' @return List of arrays containing est. survival and dist. parameter values
repeat_estimation <- function(iter, days, df.pars) {
  # Empty array to save survival probability estimates
  model.names <- c("KM.MID", "KM.IGN", "NPMLE", "MID", "IGN", "IC")
  est.prob <- array(
    NA, 
    dim = c(iter, length(days), 6, 3),
    dimnames = list(1:iter, days, model.names, 1:3)
  )
  
  # Empty array to save distribution parameter estimates
  est.pars <- array(
    NA, 
    dim = c(iter, 2, 3, 3),
    dimnames = list(1:iter, c("shape", "scale"), c("MID", "IGN", "IC"), 1:3)
  )
  
  # Repeat simulation and estimation `iter` times
  for(i in 1:iter) {
    df.lst <- do.call(generate_df, df.pars)
    for(j in 1:3) {
      est <- get_estimates(days, df.lst[[j]])
      est.prob[i,,,j] <- est[[1]]
      est.pars[i,,,j] <- est[[2]]
    }
  }
  
  return(list("prob" = est.prob, "pars" = est.pars))
}


#' Calculate summary statistics for the bias
#' 
#' @param arr Array of estimates
#' @param true.vals Vector of true parameter values
#' @return Array of summery statistics for the bias estimates
summary_bias <- function(arr, true.vals) {
  bias <- apply(arr, c(1, 3, 4), function(v) v - true.vals)
  stats.bias <- apply(bias, c(1, 3, 4), function(v)
    c(
      "bias" = mean(v),
      "rel. bias" = 0,
      "sd" = sd(v),
      "rmse" = sqrt(mean(v^2))
    ))
  stats.bias[2,,,] <- stats.bias[1,,,] / true.vals
  return(stats.bias)
}


# Global parameters ------------------------------------------------------------
set.seed(1984)
widths <- c("1" = 1, "2" = 30, "3" = 365)  # interval widths
days <- c(181, 365, 546, 730, 911, 1095, 1276, 1461, 1643, 1826)

# Short-term survival ----------------------------------------------------------
pars1 <- c("shape" = 0.5, "scale" = 200)
S1 <- 1 - pweibull(days, pars1[["shape"]], pars1[["scale"]])
df.pars1 <- list("n" = 2000, "pars" = pars1, "width" = widths)

# Estimation and derivation of statistics
est1 <- repeat_estimation(1000, days, df.pars1)
par.bias1 <- summary_bias(est1[["pars"]], pars1)
prob.bias1 <- summary_bias(est1[["prob"]], S1)

# Long-term survival -----------------------------------------------------------
pars5 <- c("shape" = 0.6, "scale" = 2000)
S5 <- 1 - pweibull(days, pars5[["shape"]], pars5[["scale"]])
df.pars5 <- list("n" = 2000, "pars" = pars5, "width" = widths)

# Estimation and derivation of statistics
est5 <- repeat_estimation(1000, days, df.pars5)
par.bias5 <- summary_bias(est5[["pars"]], pars5)
prob.bias5 <- summary_bias(est5[["prob"]], S5)

# Figures ----------------------------------------------------------------------
source("./plots.R")

# Boxplots for estimated 1-year and 5-year survival probabilities
plot_boxplots(est1[["prob"]][,2,,], S1[2])  # Short survival, 1-year
plot_boxplots(est1[["prob"]][,10,,], S1[10])  # Short survival, 5-year
plot_boxplots(est5[["prob"]][,2,,], S5[2])  # Long survival, 1-year
plot_boxplots(est5[["prob"]][,10,,], S5[10])  # Long survival, 5-year

# True survival curves and average estimated survival probabilities
plot_means(prob.bias1[1,,,1] + S1, prob.bias5[1,,,1] + S5, pars1, pars5)  # 1
plot_means(prob.bias1[1,,,2] + S1, prob.bias5[1,,,2] + S5, pars1, pars5)  # 30
plot_means(prob.bias1[1,,,3] + S1, prob.bias5[1,,,3] + S5, pars1, pars5)  # 365

# Probability bias tables
table_prob(prob.bias1[,2,,], prob.bias5[,2,,], "prob.1year.tex")  # 1-year
table_prob(prob.bias1[,10,,], prob.bias5[,10,,], "prob.5year.tex")  # 5-year

# Parameter bias table
table_par(par.bias1, par.bias5, pars1, pars5, "par.table.tex")

