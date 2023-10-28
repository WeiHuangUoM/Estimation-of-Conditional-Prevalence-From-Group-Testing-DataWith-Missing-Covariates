#Likelihood function for group testing data with missing covariate 

sl_intn=function(theta)
{
  X[delta==0]=1
  p_phatX=sapply(X, p_phat,theta)
  num=(1-p_phatX)*p0hat
  den=num+p_phatX*p1hat
  p1f=function(x){p_phat(x,theta)}
  p1f_int <- sapply(x_int,p1f)
  inte=(1-p1f_int)*f_int
  pdxint=(1-p1f_int*p1hat-p0hat+p1f_int*p0hat)*f_int
  
  p0d=num
  py0=sum(inte*interval)
  p0star.vec=(((1-p0hat)*py0)^(1-delta))*((p0d*f_X)^delta)
  p0star=rep(0,J)
  for (j in 1:J)
  {
    p0star[j]=prod(p0star.vec[index1[j]:index2[j]])
  }
  lp0=p0star^(1-Ystar)
  lp0[lp0==0]=min(lp0[lp0>0])/(1e10)
  sl1=sum(log(lp0))
  
  pd0xint=sum(pdxint*interval)
  if (pd0xint<((1-p0hat)*py0))
  {pd0xint=(1-p0hat)*py0}
  
  pdx=(den*f_X)^delta*(pd0xint^(1-delta))
  pdxj=rep(0,J)
  for (j in 1:J)
  {
    pdxj[j]=prod(pdx[index1[j]:index2[j]])
  }
  p1star=pdxj-p0star
  lp1=p1star^Ystar
  lp1[lp1==0]=min(lp1[lp1>0])/(1e10)
  sl2=sum(log(lp1))
  
  out=sl1+sl2
  ifelse(is.na(out),1000000000,-out)
}