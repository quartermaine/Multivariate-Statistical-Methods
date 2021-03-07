## ----message=FALSE, echo=F----------------------------------------------------
# import libraries --------------------------------------------------------
library(knitr)
library(corrplot)
library(ggbiplot)
library(tidyverse)
library(gridExtra)
library(psych)
knitr::opts_chunk$set(echo = F)


## -----------------------------------------------------------------------------
# Problem 1 ---------------------------------------------------------------

# a) ---------------------------------------------------

dt = read.table("T1-9.dat")
colnames(dt) = c("Country", "100m", "200m", "400m", "800m", "1500m", "3000m", "Marathon")
kable(head(dt), caption = "First lines of national track data")



## -----------------------------------------------------------------------------
# sample correlation matrix
dt_corr = cor(dt[, 2:8])
dt_eigen = eigen(dt_corr)


## -----------------------------------------------------------------------------
# correlation matrix
kable(dt_corr, caption="Correlation matrix")


## -----------------------------------------------------------------------------
# eigen values 
cat("Eigenvalues: \n", dt_eigen$values)


## -----------------------------------------------------------------------------
# eigen vectors
dt_eigen_vectors = dt_eigen$vectors
row.names(dt_eigen_vectors) = colnames(dt[,2:8])
colnames(dt_eigen_vectors) = c("COMP1", "COMP2", "COMP3","COMP4","COMP5","COMP6","COMP7")

kable(dt_eigen_vectors, caption = "Eigenvectors")



## -----------------------------------------------------------------------------
# variance importance 
pca_obj <- prcomp(dt[,2:8], scale. = F)
var_explained_dt <- data.frame(PC= paste0("PC",1:7),
                               var_explained=(pca_obj$sdev)^2/sum((pca_obj$sdev)^2))

kable(head(var_explained_dt))



## -----------------------------------------------------------------------------
cat("Cumulative percentage of the total variance explained by the first two components: ", 
    sum(var_explained_dt$var_explained[1:2])*100,"%")


## ----fig.height=6.5, fig.width=12, position="h"-------------------------------
# scree plots
g1 = var_explained_dt %>%
  ggplot(aes(x=PC,y=var_explained, group=1))+
  geom_point(size=4, col='darkturquoise')+
  geom_line(col='deeppink4')+
  labs(title="Scree plot: PCA on unscaled data")

g2 = var_explained_dt %>%
  ggplot(aes(x=PC,y=var_explained))+
  geom_col(fill="cornflowerblue")+
  labs(title="Scree plot: PCA on unscaled data")+
   geom_text(aes(label = round(var_explained,4), vjust = -1))

grid.arrange(g1, g2, ncol = 2)   



## ---- fig.align="center"------------------------------------------------------
## PCA computation
dt.pca = dt %>% 
  select(2:8) %>%
  prcomp(scale. = F, center = TRUE)

dt.pca %>% 
  ggbiplot::ggbiplot(scale = 1,
                     groups=dt$Country[which(unique(dt$Country)%in%c("SWE","USA"))],
                     ellipse = T)+labs(title = "PCA plot on unscaled data for Sweden and USA")



## -----------------------------------------------------------------------------
dt_scaled = scale(dt[,2:8])
dt_scaled_corr = cor(dt_scaled)
dt_scaled_eigen = eigen(dt_scaled_corr)
cat("===========================================================\n")
cat("First principal component for scaled variables: \n", dt_scaled_eigen$vectors[,1],"\n")
cat("===========================================================\n")
cat("Second  principal component for scaled variables: \n", dt_scaled_eigen$vectors[,2])


## -----------------------------------------------------------------------------
# ------------------------------------------------
# calculate correlation of std. variables with components
e = matrix(dt_eigen$vectors[,1:2], ncol=2) 
l = diag(dt_eigen$values[1:2]%>%sqrt())
cor_mat = as.matrix(e%*%l) 
cor_dt = as.data.frame(cor_mat)
colnames(cor_dt) = c('Comp1', 'Comp2')

kable(cor_dt, caption ='Correlation of standardized variables with components')



## -----------------------------------------------------------------------------
kable(dt_scaled_corr, caption="Correlation matrix scaled variables")



## -----------------------------------------------------------------------------
# ------------------------------------------------
pca_obj.scaled <- prcomp(dt[,2:8], scale. = T)
var_explained_dt.scaled <- data.frame(PC= paste0("PC",1:7),
                               var_explained=(pca_obj.scaled$sdev)^2/sum((pca_obj.scaled$sdev)^2))

kable(head(var_explained_dt.scaled))



## -----------------------------------------------------------------------------
cat("Cumulative percentage of the total variance explained by the first two components: ", 
    sum(var_explained_dt.scaled$var_explained[1:2])*100,"%")


## ---- fig.height=6.5, fig.width=12,position="h"-------------------------------
g1 = var_explained_dt.scaled %>%
  ggplot(aes(x=PC,y=var_explained, group=1))+
  geom_point(size=4, col='darkturquoise')+
  geom_line(col='deeppink4')+
  labs(title="Scree plot: PCA on scaled data")

g2 = var_explained_dt.scaled %>%
  ggplot(aes(x=PC,y=var_explained))+
  geom_col(fill="cornflowerblue")+
  labs(title="Scree plot: PCA on scaled data")+
   geom_text(aes(label = round(var_explained,4), vjust = -1))


grid.arrange(g1, g2, ncol = 2)   



## ---- fig.align="center"------------------------------------------------------
## PCA computation
dt_scaled.pca = dt %>% 
  select(2:8) %>%
  prcomp(scale. = TRUE, center = TRUE)


dt_scaled.pca %>% 
  ggbiplot::ggbiplot(scale = 1,
                     groups=dt$Country[which(unique(dt$Country)%in%c("SWE","USA"))],
                     ellipse = T)+labs(title = "PCA plot on scaled data for Sweden and USA")



## -----------------------------------------------------------------------------
Y_1 = as.matrix(dt_scaled) %*% dt_scaled_eigen$vectors[, 1]
rank = list(Country=dt$Country, Score=Y_1)
rank = data.frame(rank)
ordered_idxs = order(rank$Score, decreasing=TRUE)
ordered_rank = rank[ordered_idxs, ]
# top 10 countries
#  cat("Top 10 countries: \n")
kable(ordered_rank[1:10,], caption = "Top 10 countries")



## -----------------------------------------------------------------------------
# last 10 countries
# cat("Last 10 countries: \n")
kable(ordered_rank[44:54,], caption = "Last 10 countries")



## ---- fig.width=12,fig.height=8-----------------------------------------------
ordered_rank$dummy<-ifelse(ordered_rank$Score>0,"positive","negative")

ggplot(ordered_rank,aes(Country,Score,fill=dummy))+
  geom_bar(stat="identity",col="black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(title = "Ranking score plot")



## -----------------------------------------------------------------------------

# Maximum likelihood estimation and PCA for S
fit_fact_S = factanal(dt[,2:8], 
                      factors=2, 
                      covmat = cov(dt[, 2:8]))

fit_pca_S = principal(cov(dt[, 2:8]), nfactors=2, covar=T)

# Maximum likelihood estimation and PCA for R 
fit_fact_R = factanal(dt[,2:8],
                      factor=2,
                      covmat = cor(dt[,2:8]))

fit_pca_R = principal(cor(dt[,2:8]),nfactors = 2, covar = FALSE, scores = TRUE)



## -----------------------------------------------------------------------------
fit_fact_S

# cat ("Degres of freedom: ",fit_fact_S$dof, 
#      "and fit for model is: ", fit_fact_S$criteria[1])

# scores 
fit_fact_S_scores = factanal(dt[,2:8],
                             factors=2,
                             scores="Bartlett")$scores


## -----------------------------------------------------------------------------
fit_pca_S

fit_pca_S$loadings

fit_pca_S_scores = factor.scores(dt[,2:8], 
                                 f=fit_pca_S)$scores


## -----------------------------------------------------------------------------
cov_mat_loadings = cbind(fit_fact_S$loadings,fit_pca_S$loadings) 

kable(cov_mat_loadings, caption = "Loadings for covariance matrix S from factor analysis ML and PCA")



## ---- fig.width=16, fig.height=12---------------------------------------------

par(bg="whitesmoke",mfrow=c(1,2),mar=c(2,2,2,2))
# factor scores plot ML
plot(fit_fact_S_scores[,1], fit_fact_S_scores[,2],
     pch=3,
     panel.first = grid(25,25),
     main="Factor Scores (ML with Covariance)",
     xlab = "Component 1", ylab = "Component 2")
text(x = fit_fact_S_scores[,1], 
     y = fit_fact_S_scores[,2],
     labels = dt[,1], adj = c(0,1.5),cex=0.7)
# text(3, 2, "#SAM is extreme\n#KORN and COK are outiers\n#USE and KEN are on the boarder")
# We can see that the outliers are “SAM”, “KORN” and “COK”

# factors scores plot PCA
plot(fit_pca_S_scores[,1], fit_pca_S_scores[,2],
     pch=3,
     panel.first = grid(25,25),
     main="Factor Scores (PCA with Covariance)",
     xlab = "Component 1", ylab = "Component 2")
text(x = fit_pca_S_scores[,1], 
     y = fit_pca_S_scores[,2],
     labels = dt[,1], adj = c(0,1.5),cex=0.7)
# text(3, 2, "#SAM is extreme\n#KORN and COK are outiers\n#USE and KEN are on the boarder")
# The outliers this time are “MEX”, “COK” and “PNG”.


## -----------------------------------------------------------------------------

fit_fact_R
#  
# cat ("Degres of freedom: ",fit_fact_R$dof, 
#      "and fit for model is: ", fit_fact_R$criteria[1])

# scores 
fit_fact_R_scores = factanal(dt[,2:8],
                             factors=2,
                             scores="Bartlett")$scores



## -----------------------------------------------------------------------------
fit_pca_R

fit_pca_R$loadings

fit_pca_R_scores = factor.scores(dt[,2:8], 
                                 f=fit_pca_R)$scores



## -----------------------------------------------------------------------------

cov_mat_loadings = cbind(fit_fact_R$loadings,fit_pca_R$loadings) 

kable(cov_mat_loadings, caption = "Loadings for correlation matrix R from factor analysis ML and PCA")



## ----fig.width=16, fig.height=12----------------------------------------------

par(bg="whitesmoke",mfrow=c(1,2))
# factor scores plot ML
plot(fit_fact_R_scores[,1], fit_fact_R_scores[,2],
     pch=3,
     panel.first = grid(25,25),
     main="Factor Scores (ML with Correlation)",
     xlab = "Component", ylab = "Component 2")
text(x = fit_fact_R_scores[,1], 
     y = fit_fact_R_scores[,2],
     labels = dt[,1], adj = c(0,1.5),cex=0.7)
# We can see that the outliers are “SAM”, “KORN” and “COK”

# factors scores plot PCA
plot(fit_pca_R_scores[,1], fit_pca_R_scores[,2],
     pch=3,
     panel.first = grid(25,25),
     main="Factor Scores (PCA with Correlation)",
     xlab = "Component 1", ylab = "Component 2")
text(x = fit_pca_R_scores[,1], 
     y = fit_pca_R_scores[,2],
     labels = dt[,1], adj = c(0,1.5),cex=0.7)
# The outliers this time are “MEX”, “COK” and “PNG”.



## ----code=readLines(knitr::purl("/home/quartermaine/Courses/Multivariate-Statistical-Methods/labs/Assignment 3/Lab3_report.Rmd",documentation = 1)), echo = T, eval = F----
## NA

