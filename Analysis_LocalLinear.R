#---------------------------------------------------------------------------------------------------------
# Shaoke Lei's R code for calculating the local linear estimator, for the paper by Delaigle, Huang and Lei
# "Estimation of conditional prevalence from group testing data with missing covariates", 2019.
# Important note: All bandwidth selectors are valid only for standard normal kenel but can be easily adapted to another 
# kernel. The local linear estimators are also all calculated with a standard normal kernel (dnorm function) but can be  
# easily adapted to another kernel.
#---------------------------------------------------------------------------------------------------------

source("DM.R",local=TRUE)
  
  
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
  
  #----------------------------------------------------------------------------------
  #Local linear stimators
  #----------------------------------------------------------------------------------
  
 
  
  
  #Calculate the local linear estimator  of g based on full data using DM method
  g_ful_DM=function(x,X,ZZ,h){
    out=locpolvec(x,X,ZZ,h)
    out
  }
  
  X_obs=X[delta==1]
  ZZ_obs=ZZ[delta==1]
  #Calculates the corrected local linear estimator of g 
  px_1=function(x,X_obs,ZZ_obs,h){
    EZX=q*g_ful_DM(x,X_obs,ZZ_obs,h)/muZ
    num=1-EZX
    den=1+(p1hat/p0hat-1)*EZX
    out=num/den
    out
  }
  
  
#----------------------------------
  #Calulate the CV bandwidth
#----------------------------------

  njobs=rep(0,J)
  index1=c(0,cumsum(nj[1:(J-1)]))+1
  index2=cumsum(nj)
  for(j in 1:J)
  {
    w_j=delta[index1[j]:index2[j]]
    njobs[j]=sum(w_j)
  }
  njobs=njobs[njobs>0]
  Jobs=length(njobs)
  nobs=sum(njobs)
  Njobs = rep(njobs,njobs)
 
  CV=function(X,ZZ,h,n,J,nj,q,muZ){
    
    index1=c(0,cumsum(nj[1:(J-1)]))+1
    index2=cumsum(nj)
    
    index1 = rep(index1,nj)
    index2 = rep(index2,nj)
    
    qantX=quantile(X,c(0.1,0.9))
    WXi=dunif(X,qantX[1],qantX[2])*(qantX[2]-qantX[1])
    
    
    RSS=0
    for(i in 1:n){
      lgo=index1[i]:index2[i]
      
      if( i %in% which(WXi>0))
        RSS=RSS+(ZZ[i]-q^Nj[i]*g_ful_DM(X[i],X[-lgo],ZZ[-lgo],h)/muZ)^2
    }
    if(is.na(RSS)){RSS<-10000}
    RSS
  }	
  
  CV_obs_DM=function(h){CV(X_obs,ZZ_obs,h,nobs,Jobs,njobs,q,muZ)}
  hseq = seq(0.001,10,length=50)
  cv = sapply(hseq,CV_obs_DM)
  hind = which(cv==min(cv))
  h.obs_DM = hseq[hind]

#-------------------------------
  #Predict p(x)
#-------------------------------
  p_LLX=sapply(x,px_1,X_obs,ZZ_obs,h.obs_DM)
  p_LLX[p_LLX<0] = 0
  p_LLX[p_LLX>1] = 1
