#---------------------------------------------------------------------------------------------------------
# Wei Huang's R code for calculating the parametric estimator, for the paper by Delaigle, Huang and Lei
# "Estimation of conditional prevalence from group testing data with missing covariates", 2019.
# Important note: All bandwidth selectors are valid only for standard normal kenel but can be easily adapted to another 
# kernel. The local linear estimators are also all calculated with a standard normal kernel (dnorm function) but can be  
# easily adapted to another kernel.
#---------------------------------------------------------------------------------------------------------

library(DEoptim)
source("Plug_in_h_for_Gn.R", local=TRUE)
source("Plug_in_h_for_fx.R", local=TRUE)
source("Likelihood_intn.R",local=TRUE)

#---------------------------------
#Logistic parametric model
#--------------------------------
p_phat=function(x,theta)
{
  theta1=theta[1]
  theta2=theta[2]
  num=exp(theta1+theta2*x)
  den=(1+num)
  out = num/den
  out
}  

#----------------------------------------------------------------
#Read data
#----------------------------------------------------------------
nj <- 			#vector of group sizes
J=length(nj) 		#number of group
n=sum(nj) 		#sample size
X <-  			#vector of individual covariate data (n*1)
Ystar <-  		#vector of the grouped response data (J*1)
delta <- 		#vector of indicator of missingness (n*1, 1:observed; 0:missing)

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
 
  #Estimte p using local constant estimator
  X[delta==0]=1
  ker=function(x){dnorm(x)}
  Rk=integrate(function(x){ker(x)^2},-Inf,Inf)$value
  sigk2=integrate(function(x){x^2*ker(x)},-Inf,Inf)$value
  hG=newton(f=function(h){EqG(h,X,delta,Z_q,ker,Rk,sigk2,p0hat,q)}) #Bandwidth for G_n
  
  G_n=function(x){
    k=ker((x-X)/hG)
    kd=k*delta*Z_q
    out=sum(kd)/(n*hG)
    out
  }
  
  #Estimate f_X
  sigk2=integrate(function(x){x^2*ker(x)},-Inf,Inf)$value
  hf=newton(f=function(h){Eqf(h,X,delta,ker,Rk,sigk2,rhohat)}) #Bandwidth for f_X
  
  fxobs=function(x){
    k=ker((x-X)/hf)
    kd=k*delta
    out=sum(kd)/(n*hf)
    out
  }
  
  fxn = function(x){
    out1 = fxobs(x)/p1hat
    out2 = (p1hat/p0hat-1)/p1hat*G_n(x)
    out1+out2
    
  }
  
  f_int <- sapply(x_int,fxn)
  f_X=sapply(X,fxn)
  
  #Full likelihood with integration and non-parametric f_X
  op_intn=DEoptim(sl_intn,lower=c(-10,-2),upper=c(10,0)
                  ,control = DEoptim.control(itermax=60,trace=F))
  
  #Predict p(x)
  p_MLEX=sapply(x,p_phat,op_intn$optim$bestmem)
  