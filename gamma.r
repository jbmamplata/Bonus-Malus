#setwd("~/Documents/Diane/r")

library(plotly)

maxage <- 130 # maximum age
age <- 20
B = 100 # maximum death benefit
L = 10 # maximum hospitalization benefit
N <- 10      # number of states
rate <- 0.05
v <- 1/(1+rate)
BB=rep(0, N)
H = rep(0, N)
for (i in 1:N){
  BB[i] = (B/N)*(N-i+1)
  H[i] = (L/N)*(N-i+1)
}
# use Makeham formula for survival functions, set A=0 for Gompertz assumption
A = 0.0001
B = 0.00035
c = 1.075

p <- function(t){
  p <- exp(-A*t)*exp(-(B/log(c))*(c^age)*(c^t-1))
  return(p)
}

# q_x+t
qq <- function(t){
  qq <- p(t+age)-p(t+age+1)
  return(qq)
}

# parameters of Gamma
alph = 10
bet = 10

alpha=seq(1,50)
beta=seq(1,50)
Prem=matrix(0,nrow=50,ncol=50)

for (alph in 1:50){
  for (bet in 1:50){
  
    f <- function(x){
      f <- x^(alph-1)*exp(-x/bet)/(bet^alph * gamma(alph))
      return(f)
    }
  
    T=array(0, dim=c(N,N))
  for (i in 1:(N-1)){
    for (j in 1:(N-1)){
      if (i==j){
        T[i,j] = integrate(f,((i-1)*L/N),(i*L/N))$value
      }
      else if (i<j) {
        T[i,j] = integrate(f,((j-1)*L/N),(j*L/N))$value
      }
      else if (i == (j+1)) {
        T[i,j] = integrate(f,0,(j*L/N))$value
      }
    }
    T[i,N] = 1-sum(T[i,])
  }
  T[N,N-1] = integrate(f, 0, ((N-1)*L/N))$value
  T[N,N] = 1-T[N,N-1]
  
  TPM = array(0, dim=c(maxage,N,N))
  TPM[1,,]=T
  
  for (i in 2:maxage){
    TPM[i,,]=TPM[i-1,,]%*%T
  }
 

    #premium
    P=BB[1]*v*qq(0) + min(H[1],alph/bet)*v
  for (i in 1:N){
    for (j in 1:maxage-age){
      P = P + prod(BB[i],v^(j+1),TPM[j,1,i],p(j),qq(j)) + prod(min(H[i],alph/bet),v^(j+1),TPM[j,1,i],p(j))
    }
  }
   
  ##premium
  #P=BB[1]*p(0)*qq(0) + H[1]*p(0)*qq(0)
  
  #for (i in 1:N){
  #  for (j in 1:maxage-age){
  #    P = P + prod(BB[i],v^(j+1),TPM[j,1,i],p(j),qq(j)) + prod(min(alph/bet,H[i]),v^j,TPM[j,1,i],p(j),qq(j))
  #  }
  #}
  
  ann=1
  for (k in 1:(maxage-age)){
    ann = ann + v^k*p(k)
  }
  
  P=P/ann
  #print(P)
  Prem[alph,bet] = P
  }
}

persp(alpha,beta,Prem, main="premium price", zlab="premium", theta=30, phi=15, col="springgreen", shade=0.5)

#write.csv(Prem, "premium.csv")

#premium=read.csv("premium.csv")

#data = read.csv("premium_2.csv")
#plot(data$a_b, data$prem, type="l", xlab="a/b", ylab="premium")
