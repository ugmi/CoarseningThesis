library(dplyr)
library(gt)

# Figures ----------------------------------------------------------------------

#' Plot survival curves of the Weibull distribution
#' 
#' @param pars1 A named list of dist. parameters for short survival setting
#' @param pars5 A named list of dist. parameters for long survival setting
#' @return NULL
plot_curves <- function(pars1, pars5) {
  # Helper variables
  times <- 1:2000
  surv1 <- 1 - pweibull(times, pars1[["shape"]], pars1[["scale"]])
  surv5 <- 1 - pweibull(times, pars5[["shape"]], pars5[["scale"]])
  legendtxt <- c(
    paste0("shape=", pars5[["shape"]], ", scale=", pars5[["scale"]]),
    paste0("shape=", pars1[["shape"]], ", scale=", pars1[["scale"]]))
  
  # Plot
  par(mar = c(2, 4, 1, 1))
  plot(times, surv1, xlim = c(0, followup5), pch = 20, cex = 0.5, xaxt = "n",
       xlab = "", ylab = "Survival probability", las = 2, ylim = c(0, 1))
  points(times, surv5, pch = 20, cex = 0.5)
  axis(1, at = c(0, 365, 1826), labels = c("0", "1 year",  "5 years"))
  # Mark 1-year and 5-year survival probabilities
  points(c(365, 1826, 365, 1826), c(surv1[c(365, 1826)], surv5[c(365, 1826)]), 
         pch = "|")
  text(c(365, 1826, 365, 1826), c(surv1[c(365, 1826)], surv5[c(365, 1826)]), 
       labels = round(c(surv1[c(365, 1826)], surv5[c(365, 1826)]), 2), 
       pos = 3)
  # Add legend
  legend("topright", bty = "n", pch = c(20, 20), legend = legendtxt)
}


#' Plot the bias of survival probabilities
#' 
#' @param bias Array of bias estimates
#' @return NULL
plot_bias <- function(bias) {
  # Extract needed info from given objects
  labs <- dimnames(bias)
  days <- as.numeric(labs[[1]])
  ymax <- max(max(bias), max(-bias)) + 0.0005
  pchs <- c(20, 8, 21, 22, 24, 25)
  cols <- c("darkorange", "red", "brown", "forestgreen", "royalblue", "cyan3")
  
  # Make the plot
  plot(days, bias[, 1], pch = 20, col = cols[1], ylim = c(-ymax, ymax),
       xaxt = "n", xlab = "t, in days", ylab = "Bias")
  for(i in 2:6) 
    points(days, bias[, i], pch = pchs[i], col = cols[i])
  axis(1, at = days, labels = TRUE)
  lines(c(0, 2000), c(0, 0), col = "grey", lwd = 0.5)
  legend("bottom", pch = pchs, bty = "n", col = cols, legend = labs[[2]], ncol = 3)
}


#' Plot boxplots for probability estimates
#' 
#' @param ests Array of probability estimates
#' @param true.surv True value of the survival probability
#' @return NULL
plot_boxplots <- function(ests, true.surv) {
  # For setting y-axis limits and label positioning
  ymax <- max(ests)
  ymin <- min(ests)
  offset <- (ymax - ymin)/10
  
  # Plot
  par(mar = c(2, 4, 1, 1))
  boxplot(ests[,,1], boxwex = 0.15, at = 1:6 - 0.2, xaxt = "n", lty = 1,
          col = "#00000055", ylim = c(ymin - offset, ymax), xlim = c(0.7, 6.3),
          range = 0, las = 1, ylab = "Survival probability")
  boxplot(ests[,,2], boxwex = 0.15, at = 1:6, add = TRUE, lty = 1, 
          col = "#00000055", range = 0, las = 1)
  boxplot(ests[,,3], boxwex = 0.15, at = 1:6 + 0.2, add = TRUE, xaxt = "n",
          lty = 1, col = "#00000055", range = 0, las = 1)
  # Line for true survival probability
  lines(c(0, 7), rep(true.surv, 2), col = "#99999999", las = 1)
  # Additional x-axis for censoring settings
  axis(1, at = c(sapply(1:6, function(x) x - c(0.2, 0, -0.2))),
       labels = FALSE, las = 2, tck = 0.01)
  text(c(sapply(1:6, function(x) x - c(0.2, 0, -0.2))), 
       rep(ymin - offset - offset/5, 18), labels = rep(c(1, 30, 365), 6), 
       adj = 0, srt = 90, cex = 0.9)
}


#' Plot the means of estimated survival probabilities
#' 
#' @param means1 Array of means for short survival setting
#' @param means5 Array of means for long survival setting
#' @param pars1 A named list of dist. parameters for short survival setting
#' @param pars5 A named list of dist. parameters for long survival setting
#' @return NULL
plot_means <- function(means1, means5, pars1, pars5) {
  surv1 <- 1 - pweibull(0:2000, pars1[["shape"]], pars1[["scale"]])
  surv5 <- 1 - pweibull(0:2000, pars5[["shape"]], pars5[["scale"]])
  ymax <- max(max(means1), max(means5))
  labs <- dimnames(means1)
  days <- as.numeric(labs[[1]])
  pchs <- c(20, 8, 21, 22, 24, 25)
  cols <- c("darkorange", "red", "brown", "forestgreen", "royalblue", "cyan3")
  
  # Make the plot
  par(mar = c(4, 4, 1, 1), xpd = FALSE)
  plot(days, means1[, 1], pch = 20, col = cols[1], ylim = c(0, 1), xaxt = "n",
       xlim = c(0, 1900), xlab = "Time in days", ylab = "Survival probability", 
       cex = 0.7, las = 2)
  axis(1, at = days, labels = TRUE)
  for(i in 2:6) 
    points(days, means1[, i], pch = pchs[i], col = cols[i], cex = 0.7)
  for(i in 1:6) 
    points(days, means5[, i], pch = pchs[i], col = cols[i], cex = 0.7)
  points(0:2000, surv1, col = "grey", pch = 20, cex = 0.01)
  points(0:2000, surv5, col = "grey", pch = 20, cex = 0.01)
  legend("topright", pch = pchs, bty = "n", col = cols, 
         legend = labs[[2]], ncol = 2)
}

# Tables -----------------------------------------------------------------------

#' Make table for statistics of survival probability bias
#' 
#' @param bias1 Array of statistics for short survival
#' @param bia5 Array of statistics for long survival
#' @param filename Name of file (with extension) to save the table to
#' @return NULL
table_prob <- function(bias1, bias5, filename) {
  table1 <- as.data.frame(apply(bias1, 2, rbind)) %>%
    mutate(width = rep(c(1, 30, 365), each = 4),
           rown = rep(c("Bias", "Rel. bias", "SD", "RMSE"), 3),
           surv = "Short survival")
  
  table5 <- as.data.frame(apply(bias5, 2, rbind)) %>%
    mutate(width = rep(c(1, 30, 365), each = 4),
           rown = rep(c("Bias", "Rel. bias", "SD", "RMSE"), 3),
           surv = "Long survival")
  
  rbind(table1, table5) %>% group_by(surv, width) %>%
    gt(rowname_col = "rown") %>% 
    tab_spanner(label = "Non-parametric", columns = 1:3) %>%
    tab_spanner(label = "Parametric (Weibull)", columns = 4:6) %>%
    fmt_number(decimals = 5) %>% cols_align("center") %>%
    gtsave(filename)
}


#' Make data frame for statistics of parameter bias
#' 
#' @param arr Array of statistics for the parameter
#' @param partxt Label for the parameter
#' @param survtxt Label for the survival setting
#' @return A formatted data frame
make_subtable <- function(arr, partxt, survtxt) {
  sub <- as.data.frame(t(apply(arr, 1, cbind))) %>%
    mutate(rown = c("Bias", "Rel. bias", "SD", "RMSE"),
           par = partxt,
           surv = survtxt)
  colnames(sub) <- c(paste0(c("MID", "IGN", "IC"), rep(c(1, 30, 365), each = 3)), 
                     "rown", "par", "surv")
  return(sub)
}


#' Make table for statistics of parameter bias
#' 
#' @param bias1 Array of statistics for short survival
#' @param bias5 Array of statistics for long survival
#' @param pars1 A named list of dist. parameters for short survival setting
#' @param pars5 A named list of dist. parameters for long survival setting
#' @param filename Name of file (with extension) to save the table to
#' @return NULL
table_par <- function(bias1, bias5, pars1, pars5, filename) {
  # Format subtables for each parameter and survival setting
  partxt1 <- paste(c("Shape =", "Scale ="), pars1)
  sub11 <- make_subtable(bias1[,1,,], partxt1[1], "Short survival")
  sub12 <- make_subtable(bias1[,2,,], partxt1[2], "Short survival")
  partxt5 <- paste(c("Shape =", "Scale ="), pars5)
  sub51 <- make_subtable(bias5[,1,,], partxt1[1], "Short survival")
  sub52 <- make_subtable(bias5[,2,,], partxt1[2], "Short survival")
  
  table <- rbind(sub11, sub12, sub51, sub52) %>% group_by(surv, par) %>%
    gt(rowname_col = "rown") %>%
    tab_spanner(label = "1", columns = 1:3) %>%
    tab_spanner(label = "30", columns = 4:6) %>%
    tab_spanner(label = "365", columns = 7:9) %>%
    fmt_number(decimals = 5) %>% cols_align("center") %>%
    gtsave(filename)
}


