# Import libraries --------------------------------------------------------
library(ggplot2)
library(GGally)
library(reshape)
library(kableExtra)
library(knitr)
library(dplyr)
library(plotly)
library(RColorBrewer)
# Assignment 1 ------------------------------------------------------------

dt = read.delim("T1-9.dat", header=FALSE)

colnames(dt) = c('country', '100m', '200m', '400m', '800m', '1500m', '3000m','marathon')

head(dt)

col_means = sapply(dt[, -1], mean)
print(col_means, round(3))

col_sd = sapply(dt[, -1], sd)
print(col_sd)


ggpairs(dt[, -1])

# QQ-plots
par(mfrow=c(4,2), bg='whitesmoke')
for (i in 2:8){
par(mar = c(2, 2, 2, 2))
qqnorm(dt[, i], main=paste('Q-Q plot for', colnames(dt)[i]),  panel.first = grid(25,25))
qqline(dt[, i], col='red')
}


# Histograms
# Values for the normal distribution.
x = seq(-5, 5, 0.1)
y = dnorm(x)
par(mfrow=c(4,2), bg='whitesmoke')
for (i in 2:8){
hist(scale(dt[, i]),
freq=FALSE,
breaks=10,
main=paste('Histogram of variable', colnames(dt)[i]),
col='gray', 
border='blue', panel.first = grid(25,25))
lines(x, y, col='tomato4')
}

# Boxplots
par(mfrow=c(4,2), bg='whitesmoke')
for(i in 2:9){
  if(i!=9){
  boxplot(dt[, i], horizontal = TRUE, main = paste('Boxplot for variable', colnames(dt)[i]))
  # Add mean line
  segments(x0 = mean(dt[, i]), y0 = 0.8,
           x1 = mean(dt[, i]), y1 = 1.2,
           col = "red", lwd = 2)
  # Add mean point
  # points(mean(dt[, i]), 1, col = 3, pch = 19, cex=2)
  stripchart(dt[, i], method = "jitter", pch = 19, add = TRUE, col = "blue", cex =0.5)}else{
    par(mai=c(0,0,0,0))
    plot.new()
    legend('center',legend=c('points','mean'), col=c('blue', 'red'), pch=c(19, NA), lwd=c(NA, 2), cex=0.7)
  }
  
}


# Assignment 2 ------------------------------------------------------------

# a) ---------------------------------------------
# calculate matrices
corr_mat=cor(dt[, 2:8]) ; cov_mat=cov(dt[, 2:8]) 
# print correlation mat
print(corr_mat)
kable(corr_mat)
# print covariance mat
print(cov_mat)
kable(cov_mat)

# b) ---------------------------------------------
par(mfrow=c(3,2), bg='whitesmoke')
for(i in 2:7){
  name1=colnames(dt)[i+1]
  name0=colnames(dt)[i]
  title=paste0(name1," vs ",name0) 
  print(title)
  plot(dt[, i], dt[, i+1], 
       xlab=colnames(dt)[i], ylab=colnames(dt)[i+1], 
       col='red', pch =19,
       main=paste("Scatterplot ", title))
  lm_model=lm(dt[,i+1]~dt[,i], data=dt)
  abline(lm_model,lty=2, lwd=2)
  }


# c) ---------------------------------------------

my_cols= colorRampPalette(brewer.pal(8, "PiYG"))(25)
heatmap(as.matrix(dt[, 2:8]), labRow=dt$country, scale='column', col = my_cols)


# Assignment 3 ------------------------------------------------------------

euclidean_dist=function(X){
  X_centered=sweep(X, 2, colMeans(X))
  X_dist=sqrt(diag(X_centered %*% t(X_centered)))
return(X_dist)
}

distances_ed = euclidean_dist(as.matrix(dt[, 2:8])); distances_ed
idxs = sort(distances_ed, decreasing=TRUE, index.return=TRUE)$ix;idxs
countries = dt$country[idxs[1:5]]
print(countries)














