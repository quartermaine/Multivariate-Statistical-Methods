# import libraries --------------------------------------------------------
library(car)
library(latticeExtra)
library(heplots)
library(GGally)
library(knitr)
library(CARS)
# Problem 1 --------------------------------------------------------------------------------------------

# # a) ----------------------------------------------------------------------
# 
# 
# dt = read.table("T5-12.DAT")
# 
# colnames(dt) = c('x1(Tail length)', 'x2(Wing length)')
# 
# kable(head(dt))
# 
# par(bg='whitesmoke')
# plot(dt$`x1(Tail length)`, 
#      dt$`x2(Wing length)` ,
#      xlim = c(160,220),ylim=c(240,320),
#      xlab='Tail length',
#      ylab='Wing length',
#      main='Ellipse plot for 95% confidence',
#      panel.first = grid(25,25))
# 
# ellipse(dt,center = colMeans(dt), shape = cov(dt),radius = sqrt(qchisq(0.95, df=2)),col = 'purple')
# 
# 
# # -------------------------------------------------------------------------
# # https://biotoolbox.binghamton.edu/Multivariate%20Methods/Multivariate%20Tools%20and%20Background/pdf%20files/MTB%20070.pdf 
# 
# xyplot(dt$`x2(Wing length)`~dt$`x1(Tail length)`,
#        scales="free", aspect='yx',
#        xlim=c(165,220), ylim = c(240,320),
#        xlab='Tail length',ylab='Wing length',
#        par.seƫngs=list(plot.symbol=list(col='blue',cex=0.5,pch=19)),
#        panel=function(x,y,...){
#          panel.xyplot(x,y,...)
#          panel.ellipse(x,y,lwd=2,lty = 2,level=0.68,col='green')
#          panel.ellipse(x,y,lwd=2,lty=2,level=0.90,col='blue')
#          panel.ellipse(x,y,lwd=2,lty=2,level=0.95,col='red')
#          },
#        # auto.key=list(x=.1,y=.8,corner=c(0,0),
#        # lines=list(col=c('green','brown','red'),lty=c(3,3,3),lwd=6),
#        # text=list(c('68 CI','90 CI','95 CI')))
#        key=list(space="right",
#          lines=list(col=c("green","blue","red"), lty=c(2,2,2), lwd=2),
#          text=list(c("68% CI","90% CI","95% CI")))
#          )
# 
# 
# 
# # -------------------------------------------------------------------------
# 
# # b) ----------------------------------------------------------------------
# 
# # https://rpubs.com/aaronsc32/simultaneous-confidence-intervals
# 
# working.hotelling.bonferroni.intervals <- function(x, y) {
#   y <- as.matrix(y)
#   x <- as.matrix(x)
#   n <- length(y)
# 
#   # Get the fitted values of the linear model
#   fit <- lm(y ~ x)
#   fit <- fit$fitted.values
#   
#   # Find standard error as defined above
#   se <- sqrt(sum((y - fit)^2) / (n - 2)) * 
#     sqrt(1 / n + (x - mean(x))^2 / 
#            sum((x - mean(x))^2))
# 
#   # Calculate B and W statistics for both procedures.
#   W <- sqrt(2 * qf(p = 0.95, df1 = 2, df2 = n - 2))
#   B <- 1-qt(.95/(2 * 3), n - 1)
# 
#   # Compute the simultaneous confidence intervals
#   
#   # Working-Hotelling
#   wh.upper <- fit + W * se
#   wh.lower <- fit - W * se
#   
#   # Bonferroni
#   bon.upper <- fit + B * se
#   bon.lower <- fit - B * se
#   
#   xy <- data.frame(cbind(x,y))
#   
#   # Plot the Working-Hotelling intervals
#   wh <- ggplot(xy, aes(x=x, y=y)) + 
#     geom_point(size=2.5) + 
#     geom_line(aes(y=fit, x=x), size=1) + 
#     geom_line(aes(x=x, y=wh.upper), colour='blue', linetype='dashed', size=1) + 
#     geom_line(aes(x=x, wh.lower), colour='blue', linetype='dashed', size=1) +
#     labs(title='Working-Hotelling')
#   
#   # Plot the Bonferroni intervals
#   bonn <- ggplot(xy, aes(x=x, y=y)) + 
#     geom_point(size=2.5) + 
#     geom_line(aes(y=fit, x=x), size=1) + 
#     geom_line(aes(x=x, y=bon.upper), colour='blue', linetype='dashed', size=1) + 
#     geom_line(aes(x=x, bon.lower), colour='blue', linetype='dashed', size=1) +
#     labs(title='Bonferroni')
#   
#   grid.arrange(wh, bonn, ncol = 2)
#   
#   # Collect results of procedures into a data.frame and return
#   res <- data.frame(round(cbind(W, B), 3), row.names = c('Result'))
#   colnames(res) <- c('W', 'B')
#   
#   return(res)
# }
# 
# # plots
# working.hotelling.bonferroni.intervals(dt$`x2(Wing length)`, dt$`x1(Tail length)`)


# a) ----------------------------------------------------------------------

birds=read.table("T5-12.DAT")
colnames(birds)=c("Tail Length","Wing Length")
kable(head(birds))


x1=birds$`Tail Length` ; x2=birds$`Wing Length`
n=dim(birds)[1] ; p=dim(birds)[2]

a=0.05
mu0=c(190,275)
xbar=colMeans(birds)
S=cov(birds)

dist=xbar-mu0
# crit_val=sqrt(p*(n-1)/(n*(n-p)))*qf(1-a,p,n-p)
crit_val=sqrt(p*(n-1)/(n*(n-p))*qf(1-a,p,n-p))


angles=seq(0,2*pi,length.out=200)
# eigen values and eigen vectors of covariance matrix
eigVal <-eigen(S)$values
eigVec <- eigen(S)$vectors
eigScl <- eigVec%*%diag(sqrt(eigVal))
xMat <- rbind(xbar[1] + eigScl[1,]**crit_val, xbar[1]- eigScl[1,]*crit_val)
yMat <- rbind(xbar[2] + eigScl[2,]**crit_val, xbar[2]- eigScl[2,]*crit_val)
ellBase <- cbind(sqrt(eigVal[1])*crit_val*cos(angles), sqrt(eigVal[2])* crit_val*sin(angles))
ellRot <- eigVec%*%t(ellBase)

par(bg='whitesmoke')
plot(birds$`Tail Length`,
     birds$`Wing Length` ,
     xlim = c(170,220),ylim=c(240,320),
     xlab='Tail length',
     ylab='Wing length',
     main='Ellipse plot for 95% confidence',
     panel.first = grid(20,20),pch=20)
lines( (ellRot+xbar)[1,],(ellRot+xbar)[2,],
       asp=1,type="l",lwd=2,col="orange") 
points(xbar[1],xbar[2],pch=4,col="blue",lwd=3)
points(mu0[1],mu0[2],pch=3,col="red",lwd=3)
legend('topleft',legend=c('95%ellipse','sample means','test means'), 
           col=c('orange','blue', 'red'), pch=c(NA,4,3), 
           lwd=c(3,NA, NA), cex=0.7)

# b) ----------------------------------------------------------------------
# Simultaneous Intervals
f <- sqrt(((n-1)*p/(n-p))*qf(1-a, p, n-p))
sim_low <- round((t(xbar) - f * sqrt(diag(S)/n)),2)
sim_up  <- round((t(xbar) + f * sqrt(diag(S)/n)),2)
sim_interval=rbind(sim_low, sim_up) ; sim_interval
rownames(sim_interval)=c("lower band", "upper band")

kable(sim_interval)



#Bonferroni Intervals
t <- qt((1-a/(2)), df = (n-1))
bon_low <- round((t(xbar) - t * sqrt(diag(S)/n)),2)
bon_up  <- round((t(xbar) + t * sqrt(diag(S)/n)),2)
bon_interval=rbind(bon_low, bon_up)
rownames(bon_interval) <- c("lower band", "upper band")

kable(bon_interval)



# Hector's solution -------------------------------------------------------

# simultaneous intervals
# fsqrt = function(i){
# return(sqrt(p*(n-1)/ ((n-p)) * qf(0.95,df1=p, df=n-p)) * sqrt(S[i,i]/n))
# }
# 
# sim_CI = NULL
# sim_CI$lower = (xbar[1]) - fsqrt(1)
# sim_CI$upper =(xbar[1]) + fsqrt(1)
# sim_CI = as.data.frame(sim_CI)
# temp = c ((xbar[2]) - fsqrt(2), (xbar[2]) + fsqrt(2))
# sim_CI = rbind(sim_CI, temp)
# rownames(sim_CI) = c( "Mu_1_tail","Mu_2_wing")
# 
# kable(t(sim_CI))
# 

# Bonferroni CI
# set.seed(12345)
# ben_alpha = 0.95
# Ber_CI = NULL
# Ber_CI$Lower = xbar[1] - abs(qt(0.05/2*p, df=n-1)) * sqrt(S[1,1]/n)
# Ber_CI$Upper= xbar[1] + abs(qt(0.05/2*p, df=n-1) ) * sqrt(S[1,1]/n)
# Ber_CI = as.data.frame(Ber_CI)
# temp = c(xbar[2] - abs(qt(0.05/2*p, df=n-1) ) * sqrt(S[2,2]/n),
# xbar[2] + abs(qt(0.05/2*p, df=n-1) ) * sqrt(S[2,2]/n))
# Ber_CI = rbind(Ber_CI, temp)
# rownames(Ber_CI) = c("Mu_1_tail", "Mu_2_wing")
# 
# kable(Ber_CI)


# c) ----------------------------------------------------------------------

# qqnorm(birds[,1], main = "Q-Q plot for x1")
# qqline(birds[,1])
# qqnorm(birds[,2], main = "Q-Q plot for x2")
# qqline(birds[,2])
# plot(birds[,1], birds[,2])

qqPlot(birds$`Tail Length`, main ="QQ plot for X1: Tail length",id=F)

qqPlot(birds$`Wing Length`, main="QQ plot for X2: Wing Length",id=F )


# Problem 2 --------------------------------------------------------------------------------------------

# a) ----------------------------------------------------------------------
kable(head(Skulls, n=4))

ggpairs(Skulls, mapping = aes(color = epoch)) + theme_bw() 

# b) ----------------------------------------------------------------------

res = manova(cbind(mb,bh,bl,nh)~epoch, Skulls)
res
summary(res)

# The mean vectors do defer except nasal height. All other means are different between the epochs with
# a significante level between of 5 percent.
# Also, we reject the hypothesis that all the in-group means (per-epoch means) are equal.

summary.aov(res)

# c) ----------------------------------------------------------------------

w_mb =sum(res$residuals[,1]^2)
w_bh=sum(res$residuals[,2]^2)
w_bl=sum(res$residuals[,3]^2)
w_nh=sum(res$residuals[,4]^2)
w=c(w_mb, w_bh, w_bl, w_nh)
epoch =as.character(unique(Skulls$epoch))
g =length(unique(Skulls$epoch))
p =ncol(Skulls)
n = 150
a = 0.05
C = -qt(a/((p-1)*g*(g-1)), (n-g))* sqrt(2*w/(30*(n-g)))
C_mat =matrix(c(1,1,1,1),4)%*%C


#Calculating the mean values of the samples and the differences between them.
xbar=matrix(0, nrow = 5, ncol = 4, dimnames =list(epoch,names(Skulls[,-1])))
dist=matrix(0, 4,4)

for(i in 2:p){
        for(j in 1:g){
                xbar[j,(i-1)] =mean(Skulls[which(Skulls$epoch==epoch[j]),i])
        }
        for(k in 1:4) {dist[k, i-1] <- xbar[1, i-1]-xbar[k+1, i-1]
        }
}

SI_lower = dist-C_mat
SI_upper = dist+C_mat
e1 =round(t(SI_lower),3)
e2 =round(t(SI_upper),3)
interval =matrix(0, 4,4)
for(i in 1:4){
        for(j in 1:4) {
                interval[i,j] =paste("(",e1[i,j],e2[i,j],")" ,sep = " ")
                }
        }
colnames(interval) =names(Skulls[,-1])
rownames(interval) =c("(epoch1 - epoch2)", "(epoch1 - epoch3)","(epoch1 - epoch4)", "(epoch1 - epoch5)")
knitr::kable(interval)



