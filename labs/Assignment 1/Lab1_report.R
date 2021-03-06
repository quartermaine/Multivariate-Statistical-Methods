## ----message=FALSE,echo=FALSE-------------------------------------------------
# Import libraries --------------------------------------------------------
library(ggplot2)
library(GGally)
library(reshape)
# library(kableExtra)
library(knitr)
library(dplyr)
library(plotly)
library(RColorBrewer)


## ---- echo=FALSE--------------------------------------------------------------
dt = read.delim("T1-9.dat", header=FALSE)

colnames(dt) = c('country', '100m', '200m', '400m', '800m', '1500m', '3000m','marathon')

kable(dt[1:3,],
      caption = "First 3 rows of the data")



## ----echo=F-------------------------------------------------------------------
col_means = sapply(dt[, -1], mean)
kable(col_means,
       caption = "Column means")



## ----echo=F-------------------------------------------------------------------
col_sd = sapply(dt[, -1], sd)
kable(col_sd,
      caption = "Column standard deviations")



## ----echo=F-------------------------------------------------------------------
# Histograms
# Values for the normal distribution.

x = seq(-5, 5, 0.1)
y = dnorm(x)
par(mar=rep(2,4))
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


## ----echo=F-------------------------------------------------------------------
# Boxplots
par(mar=rep(2,4))
par(mfrow=c(4,2), bg='whitesmoke')
for(i in 2:9){
  if(i!=9){
  boxplot(dt[, i], horizontal = TRUE, 
          main = paste('Boxplot for variable', colnames(dt)[i]))
  # Add mean line
  segments(x0 = mean(dt[, i]), y0 = 0.8,
           x1 = mean(dt[, i]), y1 = 1.2,
           col = "red", lwd = 2)
  # Add mean point
  # points(mean(dt[, i]), 1, col = 3, pch = 19, cex=2)
  stripchart(dt[, i], method = "jitter", 
             pch = 19, add = TRUE, 
             col = "blue", cex =0.5)}else{
    par(mai=c(0,0,0,0))
    plot.new()
    legend('center',legend=c('points','mean'), 
           col=c('blue', 'red'), pch=c(19, NA), 
           lwd=c(NA, 2), cex=0.7)
  }

}



## ----echo=FALSE---------------------------------------------------------------

# a) ---------------------------------------------
# calculate matrices
corr_mat=cor(dt[, 2:8]) ; cov_mat=cov(dt[, 2:8]) 
# print correlation mat
# print(corr_mat)
kable(corr_mat,
      caption = "Correlation matrix")
# print covariance mat
# print(cov_mat)
kable(cov_mat,
      caption = "Covariance matrix")



## ----echo=FALSE---------------------------------------------------------------

# b) ---------------------------------------------
par(mfrow=c(3,2), bg='whitesmoke',  mar=c(2,2,2,2))
for(i in 2:7){
  name1=colnames(dt)[i+1]
  name0=colnames(dt)[i]
  title=paste0(name1," vs ",name0) 
  # print(title)
  plot(dt[, i], dt[, i+1], 
       xlab=colnames(dt)[i], ylab=colnames(dt)[i+1], 
       col='purple', bg ='green' ,pch =21, lwd = 1,
       main=paste("Scatterplot ", title), cex = 1.2)
  lm_model=lm(dt[,i+1]~dt[,i], data=dt)
  abline(lm_model,lty=2, lwd=2)
  }



## ----echo=FALSE---------------------------------------------------------------

# c) ---------------------------------------------
my_cols= colorRampPalette(brewer.pal(8, "PiYG"))(25)
heatmap(as.matrix(dt[, 2:8]), labRow=dt$country, scale='column', col = my_cols)




## ----echo=FALSE---------------------------------------------------------------

euclidean_dist=function(X){
  X_centered=sweep(X, 2, colMeans(X))
  X_dist=sqrt(diag(X_centered %*% t(X_centered)))
return(X_dist)
}

distances_ed = euclidean_dist(as.matrix(dt[, 2:8]));
idxs = sort(distances_ed, decreasing=TRUE, index.return=TRUE)$ix;
countries = dt$country[idxs[1:5]]
kable(as.data.frame(countries))



## ----code=readLines(knitr::purl("/home/quartermaine/Desktop/multivariate_statistical_methods-732A97/labs/Assignment 1/Lab1_report.Rmd",documentation = 1)), echo = T, eval = F----
## NA

