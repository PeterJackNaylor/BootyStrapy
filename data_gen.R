library(np)
library(MASS)

##Data
table=read.table('C:/Users/Peter-Jack/Dropbox/work/Bootstrap et ré echantillonage/BootyStrapy/mcycle.txt',header=TRUE)
attach(table)

##number of boostrap samples
B=1000

##spline smoothing
plot(times,accel)
accel.spline=smooth.spline(times,accel,df=12)
lines(accel.spline,col="blue",type="l")

##We need a smoother function + spline smoother "smoother"
find_b=function(b){
  #  if (b<1/2){return(2*b)}
  #    else{return(b+(1-b)/2)}
  if (b<0.94){return(b+0.06)}
  else{return(b+(1-b)/2)}
}

b=accel.spline$spar
b=find_b(b)
accel.spline2=smooth.spline(times,accel,df=12,spar=b)
lines(accel.spline2,col="red",type="l")

##Residual
plot(times,accel-accel.spline$y)

##N-W
res=npreg(accel ~ times)
res$mean
length(c(res$eval))
res$eval
lines(times,res$mean,col="red")


#bootstrap address in table : (x,y)
bootystrapy=function(data,B){
  n=dim(data)[1]
  matrice_x=c()
  for ( i in 1:B ){
    matrice_x=cbind(matrice_x,sample(c(1:n),n,replace=T))
  }
  return(matrice_x)
}

mcycle_B=bootystrapy(mcycle,B)


mu_B=function(matrice_x,data,values_x){
  n_x=length(values_x)
  B=dim(matrice_x)[2]
  output_hat=matrix(0,ncol=n_x,nrow=B)
  output_tilde=matrix(0,ncol=n_x,nrow=B)
  for (i in 1:B){
    times_x=data$times[matrice_x[,i]]
    accel_x=data$accel[matrice_x[,i]]
    accel.spline_x_hat=smooth.spline(times_x,accel_x,df=12)
    output_hat[i,]=predict(accel.spline_x_hat,values_x)$y
    print(b)
    b=find_b(accel.spline_x_hat$spar)
    accel.spline_x_tilde=smooth.spline(times_x,accel_x,df=12,spar=b)
    output_tilde[i,]=predict(accel.spline_x_tilde,values_x)$y
  }
  return(list(hat=output_hat,tilde=output_tilde))
}

x_point=c(10,20,25,30,35,45,50)

mu=mu_B(mcycle_B,mcycle,x_point)

mean(values_x$hat[,2])
mean(values_x$tilde[,2])
i=2
output_mu_B=mu

CI=function(output_mu_B){
  n_x=dim(output_mu_B$hat)[2]
  res=matrix(0,nrow=2,ncol=n_x)
  for (i in 1:n_x){
    mu_tilde=output_mu_B$tilde[,i]
    mu_hat=output_mu_B$hat[,i]
    mu_hat=sort(mu_hat)
    mu_tilde_mean=mean(mu_tilde)
    mu_hat_mean=mean(mu_hat)
    res[1,i]=mu_hat_mean-(mu_hat[950]-mu_tilde_mean)
    res[2,i]=mu_hat_mean-(mu_hat[50]-mu_tilde_mean)
  }
  return(res)
}

CI_mu=CI(mu)

plot(times,accel)
lines(accel.spline,col="blue",type="l")

i=0
for (x in x_point){
  i=i+1
  arrows(x, CI_mu[1,i], x, CI_mu[2,i], col= 'red',length = 0.1, angle = 90)
  arrows(x, CI_mu[2,i], x, CI_mu[1,i], col= 'red',length = 0.1, angle = 90)
}