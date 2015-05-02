n.s=c(1000,5000,10000)
B.s=c(10,50,100,250,500)

n=n.s[1]
B=B.s[1]

alpha=2
beta=2
X=rbeta(n,alpha,beta)

curve(dbeta(x,alpha,beta))

K=function(x,X,h){return(dnorm((x-X)/h))}

R_K=1
R_f_2_prime=4
sigma_k=1
n

h_star= (R_K/(n*R_f_2_prime*sigma_k^4))^1/5
?sample
