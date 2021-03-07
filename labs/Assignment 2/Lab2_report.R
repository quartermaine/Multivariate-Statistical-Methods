## ----message=FALSE,echo=FALSE-------------------------------------------------
# Import libraries --------------------------------------------------------
# import libraries --------------------------------------------------------
library(car)
library(latticeExtra)
library(heplots)
library(GGally)
library(knitr)
library(CARS)


## ---- echo=FALSE--------------------------------------------------------------
# a) ----------------------------------------------------------------------

birds=read.table("T5-12.DAT")
colnames(birds)=c("Tail Length","Wing Length")
kable(head(birds))



## ----echo=F-------------------------------------------------------------------

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



## ---- echo=F------------------------------------------------------------------
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
legend('topleft',legend=c('95%ellipse','sample means',expression(mu[0])), 
           col=c('orange','blue', 'red'), pch=c(NA,4,3), 
           lwd=c(2,NA, NA), cex=0.7)



## ---- echo=F------------------------------------------------------------------
# Simultaneous Intervals

f <- sqrt(((n-1)*p/(n-p))*qf(1-a, p, n-p))
sim_low <- round((t(xbar) - f * sqrt(diag(S)/n)),2)
sim_up  <- round((t(xbar) + f * sqrt(diag(S)/n)),2)
sim_interval=rbind(sim_low, sim_up)
rownames(sim_interval)=c("lower band", "upper band")
kable(sim_interval, caption = "Simultaneous Intervals Table")



## ---- echo=F------------------------------------------------------------------
#Bonferroni Intervals
t <- qt((1-a/(2)), df = (n-1))
bon_low <- round((t(xbar) - t * sqrt(diag(S)/n)),2)
bon_up  <- round((t(xbar) + t * sqrt(diag(S)/n)),2)
bon_interval=rbind(bon_low, bon_up)
rownames(bon_interval) <- c("lower band", "upper band")

kable(bon_interval, caption = "Bonferroni Intervals Table")



## ---- echo=F------------------------------------------------------------------
# c) ----------------------------------------------------------------------

# qqnorm(birds[,1], main = "Q-Q plot for x1", 
#        col="purple", pch=19, panel.first=grid(25, 25))
# qqline(birds[,1], col="orange", lwd=2)

qqPlot(birds$`Tail Length`, main ="QQ plot for X1: Tail length",id=F)



## ---- echo=F------------------------------------------------------------------
# qqnorm(birds[,2], main = "Q-Q plot for x2",
#        col="mediumaquamarine", pch =19,
#        panel.first=grid(25, 25))
# qqline(birds[,2], col="mediumslateblue", lwd=2)

qqPlot(birds$`Wing Length`, main="QQ plot for X2: Wing Length",id=F )



## ---- echo=F------------------------------------------------------------------
par(bg='whitesmoke')
plot(birds[,1], birds[,2], 
     xlab=colnames(birds)[1], ylab = colnames(birds)[2],
     col="tomato",pch=19, panel.first = grid(25,25),
     main="Scatter plot wing length~tail length")



## ---- echo=F------------------------------------------------------------------
# a) ----------------------------------------------------------------------
kable(head(Skulls, n=4))



## ---- echo=F, messenge=F, fig.width=12, fig.height=12-------------------------

ggpairs(Skulls, mapping = aes(color = epoch)) + theme_bw() 


## ---- echo=F------------------------------------------------------------------
res = manova(cbind(mb,bh,bl,nh)~epoch, Skulls)
print("=================== Result ========================")
res
cat("\n")
cat("\n")
print("=================== Summary =======================")
summary(res)


## ---- echo=F------------------------------------------------------------------

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


## ----code=readLines(knitr::purl("/home/quartermaine/Courses/Multivariate-Statistical-Methods/labs/Assignment 2/Lab2_report.Rmd",documentation = 1)), echo = T, eval = F----
## NA

