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

b=accel.spline$spar
b=find_b(b)
accel.spline2=smooth.spline(mcycle$times,mcycle$accel,df=12,spar=b)
lines(accel.spline2,col="red",type="l")

##Residual
plot(times,accel-accel.spline$y)

#bootstrap address in table : (x,y)
bootystrapy=function(data,B){
  n=dim(data)[1]
  matrice_x=c()
  for ( i in 1:B ){
    matrice_x=cbind(matrice_x,sample(c(1:n),n,replace=T))
  }
  return(matrice_x)
}

##Getting the different values of mu
mu_B_sm=function(matrice_x,data,values_x){
  n_x=length(values_x)
  B=dim(matrice_x)[2]
  output_hat=matrix(0,ncol=n_x,nrow=B)
  output_tilde=matrix(0,ncol=n_x,nrow=B)
  for (i in 1:B){
    times_x=data$times[matrice_x[,i]]
    accel_x=data$accel[matrice_x[,i]]
    accel.spline_x_hat=smooth.spline(times_x,accel_x,df=12)
    output_hat[i,]=predict(accel.spline_x_hat,values_x)$y
    b=find_b(accel.spline_x_hat$spar)
    accel.spline_x_tilde=smooth.spline(times_x,accel_x,df=12,spar=b)
    output_tilde[i,]=predict(accel.spline_x_tilde,values_x)$y
  }
  return(list(hat=output_hat,tilde=output_tilde))
}

mu_B_sm=function(B,data,values_x){
  
  n_x=length(values_x)
  n=length(data$times)
  output_hat=matrix(0,ncol=n_x,nrow=B)
  
  curve_smoothing=smooth.spline(data$times,data$accel,df=12)
  b=curve_smoothing$spar
  b=find_b(b)
  curve_smoothing2=smooth.spline(data$times,data$accel,df=12,spar=b)
  
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

CI(output_mu_B)
##Plotting confidence points for smothing spline

x_point=c(10,20,25,30,35,45,50)
B=1000

mu_sm=mu_B(B,mcycle,x_point)

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

## Finding the right h
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
cross_validation=function(list_x,data){
  h=1
  best_h=h
  times_b=data$times[list_x]
  accel_b=data$accel[list_x]
  
  curve_smooth=ksmooth(times_b,accel_b, kernel="normal",bandwidth=best_h,
                       range.x = range(times_b),
                       n.points = max(100L, length(times_b)),times_b)
  best_val=quantity(curve_smooth$y,accel_b,times_b,best_h)
  for (i in range(10000)){
    h=h+0.0001
    test_curve_smooth=ksmooth(times_b,accel_b, kernel="normal", 
                         bandwidth=h,
                         range.x = range(times_b),
                         n.points = max(100L, length(times_b)),times_b)
    test_best_val=quantity(curve_smooth$y,accel_b,times_b,h)
    if (test_best_val<best_val){
      curve_smooth=test_curve_smooth
      best_h=h
    }
  }
  return(best_h)
  
}

mu_B_ks=function(B,data,values_x){
  
  n_x=length(values_x)
  n=length(data$times)
  output_hat=matrix(0,ncol=n_x,nrow=B)
  h=3
  k_smooth=ksmooth(data$times,data$accel, kernel="normal", 
                   bandwidth=h,range.x = range(times_x),
                   n.points = max(100L, length(times_x)),data$times)
  
  h=h*2
  k_smooth2=ksmooth(data$times,data$accel, kernel="normal", 
             bandwidth=h,range.x = range(times_x),
             n.points = max(100L, length(times_x)),data$times)
  
  output_tilde=matrix(0,ncol=n_x,nrow=B)
  
  for (i in 1:B){
    times_x=data$times[matrice_x[,i]]
    accel_x=data$accel[matrice_x[,i]]
    ##h=cross_validation(matrice_x[,i],mcycle)
    h=3
    accel.ksm_x_hat=ksmooth(times_x,accel_x, kernel="normal", 
                            bandwidth=h,range.x = range(times_x),
                            n.points = max(100L, length(times_x)),values_x)
    output_hat[i,]=accel.ksm_x_hat$y
    b=h*2
    accel.ksm_x_tilde=ksmooth(times_x,accel_x, kernel="normal", 
                                bandwidth=b,range.x = range(times_x),
                                n.points = max(100L, length(times_x)),values_x)
    output_tilde[i,]=accel.ksm_x_hat$y
  }
  return(list(hat=output_hat,tilde=output_tilde))
}
mu_B_sm=function(B,data,values_x){
  
  n_x=length(values_x)
  n=length(data$times)
  output_hat=matrix(0,ncol=n_x,nrow=B)
  
  curve_smoothing=smooth.spline(data$times,data$accel,df=12)
  b=curve_smoothing$spar
  b=find_b(b)
  curve_smoothing2=smooth.spline(data$times,data$accel,df=12,spar=b)
  
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

##Plotting confidence intervals for certain points

B=1000
mcycle_B=bootystrapy(mcycle,B)


x_point=c(10,20,25,30,35,45,50)

mu_ks=mu_B_ks(mcycle_B,mcycle,x_point)

CI_mu_ks=CI(mu_ks)
CI_mu_ks

plot(mcycle$times,mcycle$accel,main="kernel smoothing",xlab="time",ylab="accel")
kernel_smooth=ksmooth(mcycle$times,mcycle$accel, kernel="normal", 
                      bandwidth=3)
lines(kernel_smooth,col="blue",type="l")

i=0
for (x in x_point){
  i=i+1
  arrows(x, CI_mu_ks[1,i], x, CI_mu_ks[2,i], col= 'red',length = 0.1, angle = 90)
  arrows(x, CI_mu_ks[2,i], x, CI_mu_ks[1,i], col= 'red',length = 0.1, angle = 90)
}




##Plotting both at the same time:

plot(mcycle$times,mcycle$accel,main="Comparaison",xlab="time",ylab="accel")

kernel_smooth=ksmooth(mcycle$times,mcycle$accel, kernel="normal", 
                      bandwidth=3)
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

legend(3,0.6,legend = c("(a,b)=(1,0.5)","(a,b)=(1,2)","(a,b)=(1,10)","(a,b)=(2,2)","Jeffrey"), col = colours, pch = 15, bty = "n", pt.cex = 2, cex = 0.8,, inset = c(0.1, 0.1))
