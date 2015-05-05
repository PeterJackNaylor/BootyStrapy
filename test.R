n.s=c(1000,5000,10000)
B.s=c(10,50,100,250,500)

n=n.s[1]
B=B.s[1]


##Observation

mu=0
sigma=1
X=rnorm(n,mu,sigma)

curve(dnorm(x,mu,sigma),from=mu-2*sigma^2,to=mu+2*sigma^2)

##Noyau
K=function(x,X,h){return(dnorm((x-X)/h))}


##AMISE et Sheather Jones

## Theorical

##with normal mu and sigma and gaussian kernel

R_K=1
pi=3.14
R_f_2_prime=3/(8*sqrt(pi)*sigma^5)
sigma_K=1
n

h_star.s= (R_K/(n.s*R_f_2_prime*sigma_K^4))^1/5

#With bootstrap

X_10=sample(X,B.s[1])
X_50=sample(X,B.s[2])
X_100=sample(X,B.s[3])
X_250=sample(X,B.s[4])
X_500=sample(X,B.s[5])


## Sheather and Jones
L_4=function(x){
  
  value_1=exp(-x*x/2)/sqrt(2*pi)
  value_2=x^4-6*x*x+3
  value=value_1*value_2
  return(value)
  
}


L_6=function(x){
  
  value_1=exp(-x*x/2)/sqrt(2*pi)
  value_2=x^6-11*x^4-4*x^3+25*x*x+20*x-15
  value=value_1*value_2
  return(value)
  
}


X_val=X
n=length(X_val)
lambda.hat=IQR(X_val)
a=0.920*lambda.hat*n^(-1/7)
b=0.912*lambda.hat*n^(-1/9)


diff_table=t(outer(X_val,X_val, `-`)) 
T_b=-sum(L_6(diff_table/a))/(n*(n-1)*b^7)
S=function(a){return(sum(L_4(diff_table/a))/(n*(n-1)*a^5))}
S_a=S(a)
alpha_func=function(h){return(1.357*(S_a/T_b)^(1/7)*h^(5/7))}

equation=function(h){
  value=R_K/(sigma_K*S(alpha_func(h)))
  
  return((value/n)^(1/5)-h)
}

equation(0.9)
ans=rep(0,1000)
h=rep(0,1000)
for (i in 1:1000){
  if (i==1){
    h[i]=0.001
  }else{
    h[i]=0.001+h[i-1]
  }
  ans[i]=equation(h[i])
  
} 
plot(h,ans)
h_opti_hat <- uniroot(equation,c(0.001,1))$root


###
