maxage <- 130 # maximum age
age <- 20
daily <- 10000
N <- 10      # number of states
L = 10000 # maximum hospitalization benefit
rate <- 0.05
v <- 1/(1+rate)
ben=1000000
BB=rep(0, N)
H = rep(0, N)
for (i in 1:N){
  BB[i] = (ben/N)*(N-i+1)
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


#lam <- 10   # parameter of the Poisson distribution
f <- function(lam,x){
  f <- exp(-lam)*lam^x/factorial(x)
  return(f)
}

lambda=seq(1,20)
Prem=rep(0,20)

for (lam in 1:20){
T=array(0, dim=c(N,N))
for (i in 1:(N-1)){
  for (j in 1:(N-1)){
    if (i==j){
      T[i,j] = f(lam,i-1)
    }
    else if (i<j) {
      T[i,j] = f(lam,j-1)
    }
    else if (i == (j+1)) {
      T[i,j] = sum(f(lam,0:(j-1)))
    }
  }
  T[i,N] = 1-sum(T[i,])
}
T[N,N-1] = sum(f(lam,0:(N-2)))
T[N,N] = 1-T[N,N-1]

TPM = array(0, dim=c(maxage,N,N))
TPM[1,,]=T

for (i in 2:maxage){
  TPM[i,,]=TPM[i-1,,]%*%T
}

#premium
P=BB[1]*v*qq(0) + H[1]*lam*v
for (i in 1:N){
  for (j in 1:maxage-age){
    P = P + prod(BB[i],v^(j+1),TPM[j,1,i],p(j),qq(j)) + prod(H[i],lam,v^(j+1),TPM[j,1,i],p(j))
  }
}
ann=1
for (k in 1:maxage-age){
  ann = ann + (v^k)*p(k)
}
P = P/ann
#print(P)
Prem[lam] = P

##premium
#P=BB[1]*p(0)*qq(0)

#for (i in 1:N){
#  for (j in 1:maxage-age){
#    P = P + prod(BB[i],v^j,TPM[j,1,i],p(j),qq(j))
#  }
#}

#P <- (P+daily*lam*sum(p(1:maxage-age)))/(v^c(0:maxage-age)%*%p(0:maxage-age))
##print(P)
#Prem[lam] = P
}
plot(lambda,Prem, type="l", xlab="lambda", ylab="premium")
