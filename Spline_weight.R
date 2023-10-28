#--------------------------------------------------------------------------------------------------
#Spline estimator with weight
#--------------------------------------------------------------------------------------------------

#Function to set up the matrix D_l (See Wand and Ormerod (2007), Theorem on page 11)
formOmega = function(a,b,intKnots)
{
  allKnots = c(rep(a,4),intKnots,rep(b,4))
  K = length(intKnots) ; L = 3*(K+8)
  xtilde = (rep(allKnots,each=3)[-c(1,(L-1),L)]+
              rep(allKnots,each=3)[-c(1,2,L)])/2
  wts = rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
  Bdd = spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),
                   outer.ok=TRUE)$design
  Omega = t(Bdd*wts)%*%Bdd
  return(Omega)
}


#betahat of the penalized spline estimator of E(Z|X=x)
betahat=function(K,lambda,a,b,X,delta,Z_q,w)
{
  m=length(X)
  X[is.na(X)|X<a|X>b]=a
  varphi=m/sum(w)
  d = diag(delta)
  intKnots = quantile(unique(X[delta==1]),seq(0,1,length=
                                      (K+2))[-c(1,(K+2))])
  allknots = c(rep(a,4),intKnots,rep(b,4))
  B = bs(X,knots=intKnots,degree=3,
         Boundary.knots=c(a,b),intercept=TRUE)
  Bd = crossprod(B,d)
  Bwd = tcrossprod(Bd,w)
  BTB = tcrossprod(Bwd,t(B))
  Dl = formOmega(a,b,intKnots)
  Inverse = solve(varphi*BTB+lambda*Dl,tol=1e-30)
  
  Hat=varphi*t(Bd)%*%Inverse%*%Bwd
  coef=varphi*Inverse%*%Bwd%*%Z_q
  estimator=t(Bd)%*%coef
  
  result=new.env()
  result$knots=intKnots
  result$hatmatrix=Hat
  result$Inv = Inverse
  result$BTB=BTB*varphi
  result$Bwd=Bwd
  result$beta=coef
  result$fhat=estimator
  result
}

#Leave-one-group-out cross validation 
formcv=function(K,lambda,a,b,n,X,delta,Z_q,J,nj,w)
{
  
    est = betahat(K,lambda,a,b,X,delta,Z_q,w)
    A=est$hatmatrix
    nc=length(delta)
    Lev= diag(nc)-A
    fhat=est$fhat
    dif_ym=(Z_q*delta-fhat)
    dif_leave=NULL
    index1=c(0,cumsum(nj[1:(J-1)]))+1
    index2=cumsum(nj)
    for (j in 1:J)
    {
      Lev_tmp = Lev[index1[j]:index2[j],index1[j]:index2[j]]
      dif_tmp = dif_ym[index1[j]:index2[j]]
      dif_leave_tmp = solve(Lev_tmp, dif_tmp)
      dif_leave = c(dif_leave,dif_leave_tmp)
    }
    options(digits=15)
    cv  = crossprod(dif_leave,dif_leave)
  
  cv
}


