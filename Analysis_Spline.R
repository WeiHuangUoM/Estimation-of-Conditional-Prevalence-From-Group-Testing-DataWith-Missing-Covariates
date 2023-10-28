#---------------------------------------------------------------------------------------------------------
# Wei Huang's R code for calculating the spline estimators (unweighted and weighted), for the paper by Delaigle, Huang and Lei
# "Estimation of conditional prevalence from group testing data with missing covariates", 2019.
#---------------------------------------------------------------------------------------------------------

library(splines)
source("Spline_weight.R")

#----------------------------------------------------------------
#Read data
#----------------------------------------------------------------
nj <- 			#vector of group sizes
J=length(nj) 		#number of group
n=sum(nj) 		#sample size
X <-  			#vector of individual covariate data (n*1)
Ystar <-  		#vector of the grouped response data (J*1)
delta <- 		#vector of indicator of missingness (n*1, 1:observed; 0:missing)

a=quantile(X[!is.na(X)],probs = 0)
b=quantile(X[!is.na(X)],probs = 1)  #Interval of covariate that is of interest
  

x <-                    #vector of covariates for prediction

#-----------------------------------------
#Point estimations
#-----------------------------------------
    Zstar=1-Ystar
    ZZ=rep(Zstar,nj)

  #calculate the ML estimator of q
  phi.q=function(q)
  {
    
    sumj=rep(0,J)
    for(j in 1:J)
      sumj[j]=sum(q^(0:(nj[j]-1)))
    out=sum(nj/sum(sumj)*(Zstar-q^nj))
    out
  }
  q=uniroot(f=phi.q,lower=0,upper=1)$root
  
  Z_q=q^(1-Nj)*ZZ
  
  rhohat=mean(delta)
  p0hat=sum(ZZ[delta==1])/(sum(ZZ))
  if(p0hat>1)
  {p0hat=1}
  if(is.na(p0hat)){p0hat=0.01}
  p1hat=(rhohat-p0hat*q)/(1-q)
  if(q==1)
  {p1hat=1}
  if(p1hat>1)
  {p1hat=1}
  if(p1hat<=0)
  {p1hat=0.01}
  
  muZ=mean(ZZ)
  
  #----------------------------------------------------------------------------------
  #Calculate the unweighted spline estimators
  #----------------------------------------------------------------------------------
  #Spline estimator for nonmissing
  p_spl_w=function(x,X,delta,Z_q,a,b,K,lambda,w)
  {
    est=betahat(K,lambda,a,b,X,delta,Z_q,w)
    intKnots=est$knots
    Bx= bs(x,knots=intKnots,degree=3,
           Boundary.knots=c(a,b),intercept=TRUE)
    mhat=Bx%*%est$beta
    phat=1-mhat
    as.double(phat)
  }
  
  #Corrected spline for missing
  p_cor_w=function(x,X,delta,Z_q,a,b,Kobs,lambdaobs,w)
  {
    pobs=p_spl_w(x,X,delta,Z_q,a,b,Kobs,lambdaobs,w)
    den=1+(((p1hat/p0hat)-1)*(1-pobs))
    phat=pobs/den
    as.double(phat)
  }
  
  #CV parameters for unweighted spline estimator of p 
  qantX = quantile(X[delta==1],c(0.1,0.9))
  Wxi=dunif(X,qantX[1],qantX[2])*(qantX[2]-qantX[1])
  Wxi[is.na(Wxi)]=0
  ind = Wxi>0
  njq=rep(0,J)
  index1=c(0,cumsum(nj[1:(J-1)]))+1
  index2=cumsum(nj)
  for(j in 1:J)
  {
    Wxi_j=Wxi[index1[j]:index2[j]]
    njq[j]=sum(Wxi_j)
  }
  njq=njq[njq>0]
  Jq=length(njq)
  nq=sum(njq)
  
  Vcri=function(x1,x2)
  {     formcv(x1,x2,a,b,nq,X[ind],delta[ind],Z_q[ind],Jq,njq,diag(length(X[ind])))
  }
  
  Kp = c(3,4,5,6,7,8)                        #Candidates for number of knots K
  lambdap = exp(c(0.05,0.1,0.5,0,1.5))       #Candidates for smoothing parameter lambda
  
  Kpp = rep(Kp,each=length(lambdap))
  lambdapp = rep(lambdap,length(Kp))
  
  CV = mapply(Vcri,Kpp,lambdapp)
  cvind = which(CV==min(CV))
  
  Kunw=Kpp[cvind]
  lambdaunw=lambdapp[cvind]
 
  #Predict p(x) using unweighted spline estimator
  sp_unw_X=p_cor_w(x,X,delta,Z_q,a,b,Kunw,lambdaunw,diag(length(X)))
  sp_unw_X[sp_unw_X<0] = 0
  sp_unw_X[sp_unw_X>1] = 1


#-------------------------------------------------------------------------------------
#Calculate the weighted spline estimator if data are of unequal group size
#------------------------------------------------------------------------------------- 
  source("Plug_in_h_for_fx.R", local=TRUE)  
#Estimate f_X
  ker=function(x){dnorm(x)}
  Rk=integrate(function(x){ker(x)^2},-Inf,Inf)$value
  sigk2=integrate(function(x){x^2*ker(x)},-Inf,Inf)$value
  hf=newton(f=function(h){Eqf(h,X,delta,ker,Rk,sigk2,rhohat)}) #Bandwidth for f_X
  
  fxobs=function(x){
    kd=ker((x-X[delta==1])/hf)
    out=sum(kd)/(n*hf)
    out
  }
  
  
  intervaly = (b-a)/1000
  y_int = seq((a+intervaly/2),(b-intervaly/2),intervaly)
  
  
  gpilot = 1-p_spl_w(y_int,X,delta,Z_q,a,b,Kunw,lambdaunw,diag(length(X)))
  gpilot[gpilot<0]=0
  gpilot[gpilot>1]=1
  
  gpilotx = 1-p_spl_w(x_int,X,delta,Z_q,a,b,Kunw,lambdaunw,diag(length(X)))
  gpilotx[gpilotx<0]=0
  gpilotx[gpilotx>1]=1
  
  pihat = p1hat/(1+(p1hat/p0hat-1)*gpilotx)
    
  wl = function(K,lambda)
  {
    
    
    intKnots = quantile(unique(X[delta==1]),seq(0,1,length=(K+2))[-c(1,(K+2))])
    
    B= function(x)
    {
      bs(x,knots=intKnots,degree=3,
         Boundary.knots=c(a,b),intercept=TRUE)}
    
    Npif = function(x)
    {
      t(B(x))%*%B(x)*fxobs(x)
    }
    
    
    Npifint = 0
    for (i in 1:(length(y_int)))
    {
      Npifint = Npifint+Npif(y_int[i])
    }
    Npifint = Npifint*intervaly
    
    H = lambda*formOmega(a,b,intKnots)/n+Npifint
    Hinv = solve(H)
    
    Npifgint = 0
    for (i in 1:(length(y_int)))
    {
      Npifgint = Npifgint+Npif(y_int[i])*gpilot[i]
    }
    Npifgint = Npifgint*intervaly
    
    Npifg2int = 0
    for (i in 1:(length(y_int)))
    {
      Npifg2int = Npifg2int+Npif(y_int[i])*((gpilot[i])^2)
    }
    Npifg2int = Npifg2int*intervaly
    
    V1 =function(x)
    {
      A1 = crossprod(Npifgint,Hinv)
      A2 = crossprod(A1,Hinv)
      B(x)%*%A2%*%t(B(x))
    }
    
    V2 = function(x)
    {
      A1 = crossprod(Npifg2int,Hinv)
      A2 = crossprod(A1,Hinv)
      B(x)%*%A2%*%t(B(x))
    }
    
    V1int = sum(sapply(x_int,V1)*pihat^4)*interval
    V2int = sum(sapply(x_int,V2)*pihat^4)*interval
    
    w = 1/(q^(1-Nj)*V1int-V2int)
    
    
  }
  
  Vcri=function(x1,x2)
  {
    w=wl(x1,x2)
    cv=formcv(x1,x2,a,b,nq,X[ind],delta[ind],Z_q[ind],Jq,njq,diag(w[ind]))
    c(cv,w)
  }
  
  
  
  Kp = c(3,4,5,6,7,8)                        #Candidates for number of knots K
  lambdap = exp(c(0.05,0.1,0.5,0,1.5))       #Candidates for smoothing parameter lambda
  

  Kpp = rep(Kp,each=length(lambdap))
  lambdapp = rep(lambdap,length(Kp))
  
  CV = mapply(Vcri,Kpp,lambdapp)
  cvind = which(CV[1,]==min(CV[1,]))
  
  Kw=Kpp[cvind]
  lambdaw=lambdapp[cvind]
  w = CV[2:(n+1),cvind]
  
#Predict p(x) using weighted spline estimator
  sp_w_X=p_cor_w(x,X,delta,Z_q,a,b,Kw,lambdaw,diag(w))
  sp_w_X[sp_w_X<0] = 0
  sp_w_X[sp_w_X>1] = 1