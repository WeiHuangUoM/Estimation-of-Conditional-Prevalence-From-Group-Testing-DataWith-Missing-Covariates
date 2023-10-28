

EqG=function(h,X,delta,Z_q,ker,sigk2,Rk,p0hat,q){
  sigk4=sigk2^2
  n=length(X)

  kerd4=function(x){(x^4-6*x^2+3)*dnorm(x)}
  kerd6=function(x){(x^6-15*x^4+45*x^2-15)*dnorm(x)}
  k4=function(x,y){sum(delta*Z_q*kerd4((x-X)/y),na.rm = T)}
  k6=function(x,y){sum(delta*Z_q*kerd6((x-X)/y),na.rm = T)}
  
  ShatD=function(y){
    s=0
    for (j in 1:n)
    { if (delta[j]==1)
      {s = s+Z_q[j]*k4(X[j],y)}
    }
    (n*(n-1))^(-1)*y^(-5)*s
  }
  ThatD=function(y){
    s=0
    for (j in 1:n)
    { if (delta[j]==1){
      s=s+Z_q[j]*k6(X[j],y)
    }
    }
    
    -(n*(n-1))^(-1)*y^(-7)*s
  }
  lambda=IQR(X[delta==1&Z_q>0])
  a=0.92*lambda*n^(-1/7)
  b=0.912*lambda*n^(-1/9)
  c=ShatD(a)/ThatD(b)
  alpha2h=1.357*c^(1/7)*h^(5/7)
  out=(Rk*p0hat*q/(sigk4*ShatD(alpha2h)))^(1/5)*n^(-1/5)-h
  out
}

