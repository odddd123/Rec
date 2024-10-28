#########################################################
#parameter---------------------------------------------
#########################################################
library(splines)
library(truncnorm)
load('exp1data.RData')
time=100 #number of simulation

#########################################################
#result--------------------------------------------------
#########################################################
n=800 #number of data set we chose, we could change to be 400, 600 to see result of experiment I and II
load("exp1nnn.RData") #upload estimation result, we could change to be exp1n.RData, exp1nn.RData to see result of experiment I and II

rectheta<-matrix(0,((q)+degreefreedom*p),time) #estimation of varying coefficient model 
recvaralpha<-array(0,dim = c((degreefreedom*p),(degreefreedom*p),time)) #covariance of varying coefficient 
recvarbeta<-array(0,dim = c(q,q,time)) #covariance of constant coefficient
recztheta<-matrix(0,((q)+p),time) #estimation of constant coefficient model 
reczvarbeta<-array(0,dim = c(((q)+p),((q)+p),time)) #covariance of constant coefficient model 
for (i in 1:time) {
  rectheta[,i]=t(results[[i]][[1]][[1]])
  recvaralpha[,,i]=results[[i]][[1]][[3]]
  recvarbeta[,,i]=results[[i]][[1]][[4]]
  
  recztheta[,i]=t(results[[i]][[2]][[1]])
  reczvarbeta[,,i]=results[[i]][[2]][[3]]
}

#########################################################
#RASE----------------------------------------------------
#########################################################
hgrid=0.01 
ngrid=seq(hgrid,1,hgrid) #select point

ASE1=c() #ASE 
ktill=inner
for (i in 1:time) {
  degreefreedom=1+ktill+ddgree
  sknots=seq(1/(1+ktill),1-1/(1+ktill),1/(1+ktill))
  bgrid=bs(ngrid,df=degreefreedom,degree = ddgree,knots = sknots,Boundary.knots = c(0,1),intercept = T)

  ASE1[i]=mean(((2+2*sin(ngrid*pi*2))-bgrid%*%rectheta[1:(degreefreedom),i])^2)+mean((2+2*cos(pi*ngrid*2)-bgrid%*%rectheta[(1+degreefreedom):(2*degreefreedom),i])^2)+mean((4+2*sin(pi*ngrid*2)+2*cos(2*pi*ngrid)-bgrid%*%rectheta[(1+2*degreefreedom):(3*degreefreedom),i])^2)
}

mean((ASE1));sd((ASE1)) #mean and sd of ASE

#########################################################
#censort----------------------------------------------------
#########################################################
#rate of censor
median(unlist(lapply(1:time, function(i,datas){sum(datas[[i]][,'censorstatus']==0)/datas[[i]][dim(datas[[i]])[1],1]},datas=datas)))
#average number of recurrent event for each individual
median(unlist(lapply(1:time, function(i,datas){dim(datas[[i]])[1]/datas[[i]][dim(datas[[i]])[1],1]},datas=datas)))


#########################################################
#theta----------------------------------------------------
#########################################################
selectk=ktill
ualpha=runif(10000,0,1)
sknots=seq(1/(1+selectk),1-1/(1+selectk),1/(1+selectk))
degreefreedom=1+selectk+ddgree
Balpha=bs(ualpha,df=degreefreedom,degree = ddgree,knots = sknots,Boundary.knots = c(0,1),intercept = T)
balpha1=coef(lm((2+2*sin(2*pi*ualpha))~Balpha-1)) #verify spline function form
balpha2=coef(lm((2+2*cos(2*pi*ualpha))~Balpha-1)) #verify spline function form
balpha3=coef(lm((4+2*sin(2*pi*ualpha)+2*cos(2*pi*ualpha))~Balpha-1)) #verify spline function form
balpha=c(balpha1,balpha2,balpha3)

##rectheta##########################################################r
theta0=c(balpha,beta)
means=rbind(apply(rectheta, 1, mean),theta0);means #estimation of theta

#########################################################
#rectheta##########################################################r
#########################################################
alphaasymvar=apply(recvaralpha,c(1,2),mean)
alphaasymstds=apply(recvaralpha,c(1,2),var)
betaasymvar=apply(recvarbeta,c(1,2),mean)
betaasymstds=apply(recvarbeta,c(1,2),var)
recvars=rbind(apply(rectheta, 1, var),c(diag(alphaasymvar),diag(betaasymvar)))

#########################################################
#estimation result summary----------------------------------------------------
#########################################################
tp=seq(0.1,0.9,0.1)# varying coefficient point
#beta
recbeta=apply(rectheta, 1, mean)[(1+(degreefreedom)*p):(q+(degreefreedom)*p)];round(recbeta,3)# beta estimation
sdbeta=sqrt(recvars[1,(1+(degreefreedom)*p):(q+(degreefreedom)*p)]);round(sdbeta,3)#sd of beta
sebeta=sqrt(recvars[2,(1+(degreefreedom)*p):(q+(degreefreedom)*p)]);round(sebeta,3)#se of beta
balphap=bs(tp,df=degreefreedom,degree = ddgree,knots = sknots,Boundary.knots = c(0,1),intercept = T)
balphap=cbind(balphap,balphap,balphap)
balphap# bsplines at each point
sdalpha1=apply(balphap[,1:degreefreedom]%*%rectheta[1:(degreefreedom),],1,var)#sd of alpha 1
sdalpha2=apply(balphap[,(degreefreedom+1):(degreefreedom*2)]%*%rectheta[(degreefreedom+1):(degreefreedom*2),],1,var)#se of alpha 1
sdalpha3=apply(balphap[,(degreefreedom*2+1):(degreefreedom*p)]%*%rectheta[(degreefreedom*2+1):(degreefreedom*p),],1,var) #sd of alpha 2
sealpha1=diag(balphap[,1:degreefreedom]%*%alphaasymvar[1:degreefreedom,1:degreefreedom]%*%t(balphap[,1:degreefreedom]))#se of alpha 2
sealpha2=diag(balphap[,(degreefreedom+1):(degreefreedom*2)]%*%alphaasymvar[(degreefreedom+1):(degreefreedom*2),(degreefreedom+1):(degreefreedom*2)]%*%t(balphap[,(degreefreedom+1):(degreefreedom*2)]))#sd of alpha 3
sealpha3=diag(balphap[,(2*degreefreedom+1):(degreefreedom*p)]%*%alphaasymvar[(2*degreefreedom+1):(degreefreedom*p),(2*degreefreedom+1):(degreefreedom*p)]%*%t(balphap[,(2*degreefreedom+1):(degreefreedom*p)]))#se of alpha 3
round(sqrt(sdalpha1),3);round(sqrt(sealpha1),3) #round sd and se of alpha 1
round(sqrt(sdalpha2),3);round(sqrt(sealpha2),3) #round sd and se of alpha 2
round(sqrt(sdalpha3),3);round(sqrt(sealpha3),3) #round sd and se of alpha 3

alphamean1=apply(balphap[,1:degreefreedom]%*%rectheta[(1:degreefreedom),],1,median);t(round(alphamean1,3))#median of estimation alpha 1
alphamean2=apply(balphap[,(degreefreedom+1):(2*degreefreedom)]%*%rectheta[((degreefreedom+1):(2*degreefreedom)),],1,median);t(round(alphamean2,3))#median of estimation alpha 2
alphamean3=apply(balphap[,(2*degreefreedom+1):(3*degreefreedom)]%*%rectheta[((2*degreefreedom+1):(3*degreefreedom)),],1,median);t(round(alphamean3,3))#median of estimation alpha 3

#plot of varying coefficient
par(mfrow = c(1,3))
plot(tp,alphamean1,type = 'l',ylim = c(-15,15),xlab =expression(u[i]),ylab=expression(alpha[1]))
lines(tp,2+2*sin(2*pi*tp),col="red",lty=1,lwd=1.5)
lines(tp,alphamean1+1.96*sqrt(sealpha1),col="green",lty=1,lwd=1.5)
lines(tp,alphamean1-1.96*sqrt(sealpha1),col="green",lty=1,lwd=1.5)
legend('top', legend=c("Estimated varying coefficient function",
                       "True mean functions",
                       "95% asymptotic confidence intervals"),
       col=c("black","red","green"),
       lty=c(1,1,2,2,3,3),
       cex=0.8,box.col="black",horiz = F
       ,text.font = 1)
abline(h=0,col='blue')

plot(tp,alphamean2,type = 'l',ylim = c(-10,10),xlab =expression(u[i]),ylab=expression(alpha[2]))
lines(tp,2+2*cos(2*pi*tp),col="red",lty=1,lwd=1.5)
lines(tp,alphamean2+1.96*sqrt(sealpha2),col="green",lty=1,lwd=1.5)
lines(tp,alphamean2-1.96*sqrt(sealpha2),col="green",lty=1,lwd=1.5)
legend('top', legend=c("Estimated varying coefficient function",
                       "True mean functions",
                       "95% asymptotic confidence intervals"),
       col=c("black","red","green"),
       lty=c(1,1,2,2,3,3),
       cex=0.8,box.col="black",horiz = F
       ,text.font = 1)
abline(h=0,col='blue')

plot(tp,alphamean3,type = 'l',ylim = c(-5,15),xlab =expression(u[i]),ylab=expression(alpha[3]))
lines(tp,4+2*cos(2*pi*tp)+2*sin(2*pi*tp),col="red",lty=1,lwd=1.5)
lines(tp,alphamean3+1.96*sqrt(sealpha3),col="green",lty=1,lwd=1.5)
lines(tp,alphamean3-1.96*sqrt(sealpha3),col="green",lty=1,lwd=1.5)
legend('top', legend=c("Estimated varying coefficient function",
                       "True mean functions",
                       "95% asymptotic confidence intervals"),
       col=c("black","red","green"),
       lty=c(1,1,2,2,3,3),
       cex=0.8,box.col="black",horiz = F
       ,text.font = 1)
abline(h=0,col='blue')

#########################################################
#coverage probability------------------------------------
#########################################################
qtn=-qnorm(0.025)
coverage1=(alphamean1-qtn*sqrt(sealpha1))<(balphap[,1:degreefreedom]%*%rectheta[(1:degreefreedom),1])&(alphamean1+qtn*sqrt(sealpha1))>(balphap[,1:degreefreedom]%*%rectheta[(1:degreefreedom),1])
coverage2=(alphamean2-qtn*sqrt(sealpha2))<(balphap[,(degreefreedom+1):(2*degreefreedom)]%*%rectheta[((degreefreedom+1):(2*degreefreedom)),1])&(alphamean2+qtn*sqrt(sealpha2))>(balphap[,(degreefreedom+1):(2*degreefreedom)]%*%rectheta[((degreefreedom+1):(2*degreefreedom)),1])
coverage3=(alphamean3-qtn*sqrt(sealpha3))<(balphap[,(2*degreefreedom+1):(3*degreefreedom)]%*%rectheta[((2*degreefreedom+1):(3*degreefreedom)),1])&(alphamean3+qtn*sqrt(sealpha3))>(balphap[,(2*degreefreedom+1):(3*degreefreedom)]%*%rectheta[((2*degreefreedom+1):(3*degreefreedom)),1])
for (j in 2:time) {
  coverage1=cbind(coverage1,(alphamean1-qtn*sqrt(sealpha1))<(balphap[,1:degreefreedom]%*%rectheta[(1:degreefreedom),j])&(alphamean1+qtn*sqrt(sealpha1))>(balphap[,1:degreefreedom]%*%rectheta[(1:degreefreedom),j]))
  coverage2=cbind(coverage2,(alphamean2-qtn*sqrt(sealpha2))<(balphap[,(degreefreedom+1):(2*degreefreedom)]%*%rectheta[((degreefreedom+1):(2*degreefreedom)),j])&(alphamean2+qtn*sqrt(sealpha2))>(balphap[,(degreefreedom+1):(2*degreefreedom)]%*%rectheta[((degreefreedom+1):(2*degreefreedom)),j]))
  coverage3=cbind(coverage3,(alphamean3-qtn*sqrt(sealpha3))<(balphap[,(2*degreefreedom+1):(3*degreefreedom)]%*%rectheta[((2*degreefreedom+1):(3*degreefreedom)),j])&(alphamean3+qtn*sqrt(sealpha3))>(balphap[,(2*degreefreedom+1):(3*degreefreedom)]%*%rectheta[((2*degreefreedom+1):(3*degreefreedom)),j]))
}
cp1=apply(coverage1, 1, mean);cp1 #coverage probability of alpha 1
cp2=apply(coverage2, 1, mean);cp2 #coverage probability of alpha 2
cp3=apply(coverage3, 1, mean);cp3 #coverage probability of alpha 3

coverage=(recbeta-qtn*sqrt(sebeta))<rectheta[(degreefreedom*p+1):(degreefreedom*p+q),1]&(recbeta+qtn*sqrt(sebeta))>rectheta[(degreefreedom*p+1):(degreefreedom*p+q),1]
for (j in 2:time) {
  coverage=cbind(coverage,(recbeta-qtn*sqrt(sebeta))<rectheta[(degreefreedom*p+1):(degreefreedom*p+q),j]&(recbeta+qtn*sqrt(sebeta))>rectheta[(degreefreedom*p+1):(degreefreedom*p+q),j])
}
cp=apply(coverage, 1, mean);round(cp,3) #coverage probability of beta


