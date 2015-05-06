
mcycle {MASS}
library(MASS)
table=data(mcycle)
print(table)
attach(table)
table=read.table('mcycle.txt',header=TRUE)
attach(table)
plot(times,accel)

accel.spline=smooth.spline(times,accel,df=12)
plot(times,accel)
lines(accel.spline,col="blue",type="l")

plot(times,accel-accel.spline$y)
library(MASS)
library(np)
res$mean
length(c(res$eval))
res$eval$type
res=npreg(accel ~ times)
lines(times,res$mean,col="red")

accel.spline
smooth.spline$