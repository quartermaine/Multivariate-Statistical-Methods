## 732A97 Multivariate Statistical Methods

data<-read.table("T1-9.dat",row.names=1)
## column names are for some reason not in the dat file
colnames(data)<-c("100m (s)","200m (s)","400m (s)","800m (s)","1500m (min)","3000m (min)","Marathon (min)")

#Problem 1

#a)

R<-cor(data)
eigen(R)
lam<-eigen(R)$values
e<-eigen(R)$vectors

#b)

y1<-e[,1]
y2<-e[,2]
prop1<-lam[1]/7
prop2<-lam[2]/7
prop<-(lam[1]+lam[2])/7

r11<-y1[1]*sqrt(lam[1])
r12<-y1[2]*sqrt(lam[1])
r13<-y1[3]*sqrt(lam[1])

r21<-y2[1]*sqrt(lam[2])
r22<-y2[2]*sqrt(lam[2])
r23<-y2[3]*sqrt(lam[2])

#...


#d) 

#Standardization of data

x<-data.matrix(data)
mean<-apply(x,2,mean)
sd<-apply(data,2,sd)
D12<-diag(sd)
u<-c(rep(1,54))
v<-kronecker(u,t(mean))
z<-(x-v)%*%solve(D12)

score<-z%*%y1
a<-cbind(score,rank(score))

a

#Problem 2

S<-cov(data)
eigen(S)
lams<-eigen(S)$values
es<-eigen(S)$vectors

#m=2 

ys1<-es[,1]
ys2<-es[,2]
props1<-lams[1]/sum(diag(S))
props2<-lams[2]/sum(diag(S))
props<-(lams[1]+lams[2])/sum(diag(S))

L1<-sqrt(lams[1])*ys1
L2<-sqrt(lams[2])*ys2
L<-cbind(L1,L2)
LL<-L%*%t(L)

hs1<-LL[1,1]
hs2<-LL[2,2]
#...

psis1<-S[1,1]-LL[1,1]
psis2<-S[2,2]-LL[2,2]
#...

fs1<-solve(t(L)%*%L)%*%t(L)%*%(x[1,]-mean)

FS0<-c(rep(0,54))
FS<-cbind(FS0,FS0)

for (i in 1: 54) {fs<-solve(t(L)%*%L)%*%t(L)%*%(x[i,]-mean)
FS[i,]<-fs}

#For R

#m=2

LR1<-sqrt(lam[1])*y1
LR2<-sqrt(lam[2])*y2
LR<-cbind(LR1,LR2)
LLR<-LR%*%t(LR)

h1<-LLR[1,1]
h2<-LLR[2,2]
psi1<-R[1,1]-LLR[1,1]
psi2<-R[2,2]-LLR[2,2]

f11<-1/sqrt(lam[1])*t(y1)%*%solve(D12)%*%(x[1,]-mean)
f12<-1/sqrt(lam[2])*t(y2)%*%solve(D12)%*%(x[1,]-mean)

f1<-solve(t(LR)%*%LR)%*%t(LR)%*%solve(D12)%*%(x[1,]-mean)

F0<-c(rep(0,54))
F<-cbind(F0,F0)

for (i in 1: 54) {f<-solve(t(LR)%*%LR)%*%t(LR)%*%solve(D12)%*%(x[i,]-mean)
F[i,]<-f}



