library(np)
library(MASS)

##Data
table=read.table('C:/Users/Peter-Jack/Dropbox/work/Bootstrap/BootyStrapy/mcycle.txt',header=TRUE)
attach(table)

##number of boostrap samples
B=1000

##spline smoothing
plot(mcycle$times,mcycle$accel)
accel.spline=smooth.spline(mcycle$times,mcycle$accel,df=12)
lines(accel.spline,col="blue",type="l")

##We need a smoother function + spline smoother "smoother"
find_b=function(b){
  #  if (b<1/2){return(2*b)}
  #    else{return(b+(1-b)/2)}
  if (b<0.96){return(b+0.04)}
  else{return(b+(1-b)/2)}
}


accel.spline=smooth.spline(mcycle$times,mcycle$accel,df=12)
b=accel.spline$spar

b=b+0.1
accel.spline2=smooth.spline(mcycle$times,mcycle$accel,df=12,spar=b)
b=b-0.2
accel.spline3=smooth.spline(mcycle$times,mcycle$accel,df=12,spar=b)

##Residual
par(mfrow=c(1, 3))
plot(times,accel-accel.spline$y,main="Cross-validation",xlab="time",ylab="accel")
plot(times,accel-accel.spline2$y,main="Smoother",xlab="time",ylab="accel")
plot(times,accel-accel.spline3$y,main="Less-smoother",xlab="time",ylab="accel")
par(mfrow=c(1, 1))

##Getting the different values of mu

mu_B_sm=function(B,data,values_x){
  
  n_x=length(values_x)
  n=length(data$times)
  output_hat=matrix(0,ncol=n_x,nrow=B)
  
  curve_smoothing_normal=smooth.spline(data$times,data$accel,df=12)
  b=curve_smoothing_normal$spar
  
  curve_smoothing=smooth.spline(data$times,data$accel,df=12,spar=b+0.1)
  
  curve_smoothing2=smooth.spline(data$times,data$accel,df=12,spar=b-0.1)
  
  epsilon_star=sample(residuals(curve_smoothing2),n,replace=T)
  y_star=fitted(curve_smoothing)+epsilon_star
  
  mu_hat_x   = predict(curve_smoothing  , values_x)$y
  mu_tilde_x = predict(curve_smoothing2 , values_x)$y
  
  for (i in 1:B){
    epsilon_star=sample(residuals(curve_smoothing2),n,replace=T)
    y_star=fitted(curve_smoothing)+epsilon_star
    
    curve_smoothing_star=smooth.spline(data$times,y_star,df=12)
    
    output_hat[i,]=predict(curve_smoothing_star,values_x)$y
  }
  return(list(hat=mu_hat_x,tilde=mu_tilde_x,hat_star=output_hat))
}

CI=function(output_mu_B){
  n_x=dim(output_mu_B$hat_star)[2]
  B=dim(output_mu_B$hat_star)[1]
  res=matrix(0,nrow=2,ncol=n_x)
  for (i in 1:n_x){
    mu.hat.x=output_mu_B$hat[i]
    mu.til.x=output_mu_B$tilde[i]
    
    mu.x.star=sort(na.omit(output_mu_B$hat_star[,i]))
    
    res[1,i]=mu.hat.x-(mu.x.star[as.integer(length(mu.x.star)*0.95)]-mu.til.x)
                       
    res[2,i]=mu.hat.x-(mu.x.star[as.integer(length(mu.x.star)*0.05)]-mu.til.x)
                       
  }
  return(res)
}

##Plotting confidence points for smothing spline

x_point=c(10,20,25,30,35,45,50)
B=1000

mu_sm=mu_B_sm(B,mcycle,x_point)

CI_mu_sm=CI(mu_sm)
CI_mu_sm

accel.spline=smooth.spline(mcycle$times,mcycle$accel)
plot(mcycle$times,mcycle$accel,main="Smoothing spline",xlab="time",ylab="accel")
lines(accel.spline,col="green",type="l")

i=0
for (x in x_point){
  i=i+1
  arrows(x, CI_mu_sm[1,i], x, CI_mu_sm[2,i], col= 'orange',length = 0.1, angle = 90)
  arrows(x, CI_mu_sm[2,i], x, CI_mu_sm[1,i], col= 'orange',length = 0.1, angle = 90)
}

## Finding the right h / cross validation

gamma_i=function(x,X,h,i){
  gam=dnorm((x-X)/h)
  return((gam/sum(gam))[i])
}

quantity=function(y,y.hat,X,h){
  num=(y-y.hat)^2
  n=length(y)
  dem=rep(0,n)
  for (i in 1:n){
    dem[i]=(1-gamma_i(X[i],X,h,i))
  }
  dem=dem*dem
  val=sum(num/dem)
  return(val)
}
cross_validation=function(x,y){
  h=1
  best_h=h
  
  curve_smooth=ksmooth(x,y, kernel="normal",bandwidth=best_h,
                       range.x = range(x),
                       n.points = max(100L, length(x)),x)
  best_val=quantity(curve_smooth$y,y,x,best_h)
  for (i in range(10000)){
    h=h+0.0001
    test_curve_smooth=ksmooth(x,y, kernel="normal", 
                         bandwidth=h,
                         range.x = range(x),
                         n.points = max(100L, length(x)),x)
    test_best_val=quantity(curve_smooth$y,y,x,h)
    if (test_best_val<best_val){
      curve_smooth=test_curve_smooth
      best_h=h
    }
  }
  return(best_h)
  
}
## Nadaraya et Watson

mu_B_ks=function(B,data,values_x){
  
  n_x=length(values_x)
  n=length(data$times)
  output_hat=matrix(0,ncol=n_x,nrow=B)
  
  h=npregbw(xdat=data$times, ydat=data$accel,regtype="lc")$bw
  
  k_smooth=ksmooth(data$times,data$accel, kernel="normal", 
                   bandwidth=2*h,range.x = range(times_x),
                   n.points = max(100L, length(times_x)),data$times)
  

  k_smooth2=ksmooth(data$times,data$accel, kernel="normal", 
             bandwidth=h/2,range.x = range(times_x),
             n.points = max(100L, length(times_x)),data$times)
  
  epsilon_s=data$accel-k_smooth2$y
  
  mu_hat_x   = ksmooth(data$times,data$accel, kernel="normal", 
                       bandwidth=h,range.x = range(times_x),
                       n.points = max(100L, length(times_x)),values_x)$y
  mu_tilde_x = ksmooth(data$times,data$accel, kernel="normal", 
                       bandwidth=2*h,range.x = range(times_x),
                       n.points = max(100L, length(times_x)),values_x)$y
    
  for (i in 1:B){
    epsilon_star=sample(epsilon_s,n,replace=T)
    y_star=k_smooth$y+epsilon_star
    
    h=npregbw(xdat=data$times, ydat=y_star,regtype="lc")$bw
    accel.ksm_x_hat=ksmooth(data$times,y_star, kernel="normal", 
                            bandwidth=h,range.x = range(times_x),
                            n.points = max(100L, length(times_x)),values_x)
    
    output_hat[i,]=accel.ksm_x_hat$y
  }
  return(list(hat=mu_hat_x,tilde=mu_tilde_x,hat_star=output_hat))
}

##Plotting confidence intervals for certain points

B=1000

x_point=c(10,20,25,30,35,45,50)

mu_ks=mu_B_ks(B,mcycle,x_point)

CI_mu_ks=CI(mu_ks)
CI_mu_ks


h=npregbw(xdat=mcycle$times, ydat=mcycle$accel,regtype="lc")$bw

plot(mcycle$times,mcycle$accel,main="kernel smoothing",xlab="time",ylab="accel")
kernel_smooth=ksmooth(mcycle$times,mcycle$accel, kernel="normal", 
                      bandwidth=2*h)
lines(kernel_smooth,col="blue",type="l")

i=0
for (x in x_point){
  i=i+1
  arrows(x, CI_mu_ks[1,i], x, CI_mu_ks[2,i], col= 'red',length = 0.1, angle = 90)
  arrows(x, CI_mu_ks[2,i], x, CI_mu_ks[1,i], col= 'red',length = 0.1, angle = 90)
}




##Plotting both at the same time:

plot(mcycle$times,mcycle$accel,main="Non parametric estimation",xlab="time",ylab="accel")

h=npregbw(xdat=mcylce$times, ydat=mcycle$accel,regtype="lc")$bw
h=2
kernel_smooth=ksmooth(mcycle$times,mcycle$accel, kernel="normal", 
                      bandwidth=h)
lines(kernel_smooth,col="blue",type="l")

i=0
for (x in x_point){
  i=i+1
  arrows(x, CI_mu_ks[1,i], x, CI_mu_ks[2,i], col= 'red',length = 0.1, angle = 90)
  arrows(x, CI_mu_ks[2,i], x, CI_mu_ks[1,i], col= 'red',length = 0.1, angle = 90)
}

accel.spline=smooth.spline(mcycle$times,mcycle$accel)
lines(accel.spline,col="green",type="l")

i=0
for (x in x_point){
  i=i+1
  arrows(x, CI_mu_sm[1,i], x, CI_mu_sm[2,i], col= 'orange',length = 0.1, angle = 90)
  arrows(x, CI_mu_sm[2,i], x, CI_mu_sm[1,i], col= 'orange',length = 0.1, angle = 90)
}

legend(40,-50,legend =c("kernel smoothing", "spline smoothing", "CI for k_s","CI for s_s"), col = c("blue","green","red","orange"),
       pch = 15, bty = "n", pt.cex = 2, cex = 0.8,, inset = c(0.1, 0.1))


##############################################################

# WILD BOOTSTRAP 

rm(list=ls())

library(bbemkr)
library(fANCOVA)
table=read.table("mcycle.txt",header=TRUE)
attach(table)

x=times
y=accel
NW=NadarayaWatsonkernel(x, y, h = 1, gridpoint = x)
mat <- cbind(y,NW$mh)
matplot(x, mat,pch=3)
# plot(NW$gridpoint,NW$mh,col='blue')
reshat<- y-NW$mh
reshat
res.boot <- wild.boot(reshat, nboot=1)
res.boot
NW2<-NadarayaWatsonkernel(x, y, h = 3, gridpoint = x)
ystar<-NW2$mh+res.boot
result<-NadarayaWatsonkernel(x, ystar, h = 2, gridpoint = x)
mat2<- cbind(y,result$mh)
matplot(x, mat2,pch=3)
