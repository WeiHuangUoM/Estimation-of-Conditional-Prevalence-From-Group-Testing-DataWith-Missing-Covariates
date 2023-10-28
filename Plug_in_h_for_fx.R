

Eqf=function(h,X,delta,ker,sigk2,Rk,rhohat){
  sigk4=sigk2^2
  Xobs=X[delta==1]
  n=length(X)
  
  kerd4=function(x){(x^4-6*x^2+3)*dnorm(x)}
  kerd6=function(x){(x^6-15*x^4+45*x^2-15)*dnorm(x)}
  k4=function(x,y){sum(kerd4((x-Xobs)/y))}
  k6=function(x,y){sum(kerd6((x-Xobs)/y))}
  ShatD=function(y){
    s=sapply(Xobs,k4,y)
    (n*(n-1))^(-1)*y^(-5)*sum(s)
  }
  ThatD=function(y){
    s=sapply(Xobs,k6,y)
    -(n*(n-1))^(-1)*y^(-7)*sum(s)
  }
  lambda=IQR(Xobs)
  a=0.92*lambda*n^(-1/7)
  b=0.912*lambda*n^(-1/9)
  c=ShatD(a)/ThatD(b)
  alpha2h=1.357*c^(1/7)*h^(5/7)
  out=(Rk*rhohat/(sigk4*ShatD(alpha2h)))^(1/5)*n^(-1/5)-h
  out
}

newton=function(f,tol=1e-5,x0=1,N=100){
  h=1e-7
  i=1
  x1=x0
  p=numeric(N)
  while(i <= N){
    df.dx=(f(x0+h)-f(x0))/h
    x1=(x0-(f(x0)/df.dx))
    p[i]=x1
    i=i+1
    if(abs(x1-x0)<tol)break
    x0=x1
    
  }
  return(p[(i-1)])
}

