library(robustbase)
library(investr)
library(ggplot2)
library(gridExtra)
data("delivery")


dist.lm <- lm(delTime ~ distance, data = delivery)
prod.lm <- lm(delTime ~ n.prod, data = delivery)


par(mfrow=c(2,2))
#
plotFit(dist.lm, interval = 'confidence', adjust = 'Scheffe', main = 'Working-Hotelling DelTime ~ Distance')
plotFit(prod.lm, interval = 'confidence', adjust = 'Scheffe', main = 'Working-Hotelling DelTime ~ Products')
#
plotFit(dist.lm, interval = 'confidence', k = 0.95, adjust = 'Bonferroni', main = 'Bonferroni DelTime ~ Distance')
plotFit(prod.lm, interval = 'confidence', k = 0.95, adjust = 'Bonferroni', main = 'Bonferroni DelTime ~ Products')

working.hotelling.bonferroni.intervals <- function(x, y) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  n <- length(y)

  # Get the fitted values of the linear model
  fit <- lm(y ~ x)
  fit <- fit$fitted.values
  
  # Find standard error as defined above
  se <- sqrt(sum((y - fit)^2) / (n - 2)) * 
    sqrt(1 / n + (x - mean(x))^2 / 
           sum((x - mean(x))^2))

  # Calculate B and W statistics for both procedures.
  W <- sqrt(2 * qf(p = 0.95, df1 = 2, df2 = n - 2))
  B <- 1-qt(.95/(2 * 3), n - 1)

  # Compute the simultaneous confidence intervals
  
  # Working-Hotelling
  wh.upper <- fit + W * se
  wh.lower <- fit - W * se
  
  # Bonferroni
  bon.upper <- fit + B * se
  bon.lower <- fit - B * se
  
  xy <- data.frame(cbind(x,y))
  
  # Plot the Working-Hotelling intervals
  wh <- ggplot(xy, aes(x=x, y=y)) + 
    geom_point(size=2.5) + 
    geom_line(aes(y=fit, x=x), size=1) + 
    geom_line(aes(x=x, y=wh.upper), colour='blue', linetype='dashed', size=1) + 
    geom_line(aes(x=x, wh.lower), colour='blue', linetype='dashed', size=1) +
    labs(title='Working-Hotelling')
  
  # Plot the Bonferroni intervals
  bonn <- ggplot(xy, aes(x=x, y=y)) + 
    geom_point(size=2.5) + 
    geom_line(aes(y=fit, x=x), size=1) + 
    geom_line(aes(x=x, y=bon.upper), colour='blue', linetype='dashed', size=1) + 
    geom_line(aes(x=x, bon.lower), colour='blue', linetype='dashed', size=1) +
    labs(title='Bonferroni')
  
  grid.arrange(wh, bonn, ncol = 2)
  
  # Collect results of procedures into a data.frame and return
  res <- data.frame(round(cbind(W, B), 3), row.names = c('Result'))
  colnames(res) <- c('W', 'B')
  
  return(res)
}

working.hotelling.bonferroni.intervals(delivery$n.prod, delivery$delTime)

working.hotelling.bonferroni.intervals(delivery$distance, delivery$delTime)
