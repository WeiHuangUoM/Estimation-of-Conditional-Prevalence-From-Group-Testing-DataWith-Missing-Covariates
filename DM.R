#---------------------------------------------------------------------------------------------------------
# Aurore Delaigle's R code for calculating the estimator and ROT and PI bandwidths, for the paper by Delaigle and Meister
# "Nonparametric regression analysis for group testing data", 2011.
# HERE ALL GROUPS ARE OF EQUAL SIZE
# Important note: All bandwidth selectors are valid only for standard normal kenel but can be easily adapted to another 
# kernel. The local poly estimators are also all calculated with a standard normal kernel (dnorm function) but can be  
# easily adapted to another kernel.
#---------------------------------------------------------------------------------------------------------




#This routine is needed to calculate the ML estimator of q
fobject=function(qq)
{

	sumj=rep(0,J)
	for(j in 1:J)
		sumj[j]=sum(qq^(0:(nj[j]-1)))

	som=sum(nj/sum(sumj)*(Zstar-qq^nj))
	som
}


#This routine calculates the local linear estimator  of g at a single point x
fnorm=function(x,X,ZZ,h)
{
	n=length(X)
	a0=dnorm((x-X)/h)
	a1=a0*(X-x)/h
	a2=a1*(X-x)/h
	
	S0=sum(a0)/(n*h)
	S1=sum(a1)/(n*h)
	S2=sum(a2)/(n*h)
	T0=sum(ZZ*a0)/(n*h)
	T1=sum(ZZ*a1)/(n*h)
	
	(T0*S2-T1*S1)/(S0*S2-S1^2)
}



#This routine calculates the local linear estimator of g at a vector x
locpolvec=function(x,X,ZZ,h)
{

	sapply(x,fnorm,X,ZZ,h)

}


#--------------------------------------------------------------------------------------------------
#DM plugin/ROT bandwidth
#--------------------------------------------------------------------------------------------------
#Local poly estimator of order p of the dth derivatiave g^{(d)}
g.lpder=function(x,X,ZZ,h,d,p){
	      n=length(X)
	      a0=dnorm((X-x)/h)
	      S0=sum(a0)/(n*h)
	      T0=sum(ZZ*a0)/(n*h)

	      if(p>0){
	   	    for(i in 1:(2*p)){
			   eval(parse(text=paste("a",i,"=a",i-1,"*(X-x)/h", sep = "")))
			   eval(parse(text=paste("S",i,"=sum(a",i,")/(n*h)", sep = "")))
			   eval(parse(text=paste("T",i,"=sum(ZZ*a",i,")/(n*h)", sep = "")))
		    }
		  }

	      TT=c()
	      SS=matrix(0,nrow=p+1,ncol=p+1)
	      for(i in 0:p){
		     eval(parse(text=paste("TT=c(TT,T",i,")", sep = "")))
		     for(j in 0:p)
		     eval(parse(text=paste("SS[",i+1,",",j+1,"]=S",i+j, sep = "")))
		  }

	      out=solve(SS)%*%TT
	      out=factorial(d)*h^(-d)*out[d+1]
	      out
}

#Fitted values of local poly estimator of order p of the dth derivatiave g^{(d)}
g.lpderfit=function(x,X,ZZ,h,d,p){
	sapply(x,g.lpder,X,ZZ,h,d,p)
}

#Calculates a global polynomial estimator or order p of the d derivative of g, where d can be 0, 2 or 4 and p can be d, d+1 or d+2. Can easily be extended to other derivatives d and other values of p

globpoly=function(x,X,ZZ,d,p,q){
	sdX=sqrt(var(X))
	x=x/sdX
	X=X/sdX
	#find the index of the first component of each group
	cumsom=cumsum(nj)+1
	cumsom=cumsom[-J]
	cumsom=c(1,cumsom)
	hatbeta=rep(0,p+1)
	#These are the T_j*'s
	ZZ=ZZ*q^(-Nj)*mean(ZZ)
	#Design matrix
	XD=0*ZZ+1
	for(j in 1:p)
	   XD=cbind(XD,X^j)
	#Estimator of beta
	hatbeta=hatbeta+solve(t(XD)%*%XD)%*%t(XD)%*%ZZ
	rm(XD)
	#Estimator of g
	gofx=hatbeta[1]
	for(j in 1:p)
	   gofx=gofx+hatbeta[j+1]*x^j
	#Take the dth derivative
	if((d==2)&&(p==4))
	    gofxpp=2*hatbeta[3]+6*hatbeta[4]*x+12*hatbeta[5]*x^2
	if((d==2)&&(p==3))
	    gofxpp=2*hatbeta[3]+6*hatbeta[4]*x
	if((d==2)&&(p==2))
	    gofxpp=2*hatbeta[3]
	if((d==4)&&(p==6))
		gofxpp=24*hatbeta[5]+120*hatbeta[6]*x+720*hatbeta[7]*x^2
	if((d==4)&&(p==5))
		gofxpp=24*hatbeta[5]+120*hatbeta[6]*x
	if((d==4)&&(p==4))
		gofxpp=24*hatbeta[5]
	if(d==0)
		gofxpp=gofx
	gofxpp
}
#-----------------------------------------------------------------------------------------------
# calculates the weighted PI bandwidth for a standard normal kernel K, with weight weight=c(a,b)
#------------------------------------------------------------------------------------------------
h.plugin=function(X,ZZ,weight,q,th24spline){
	N=length(X)
	#Second moment of a standard normal kernel, mu2 in the paper
	muK2=1
	#\int K^2 for a standard normal kernel
	RK=1/(2*sqrt(pi))

	qantX=range(X)
	x=seq(qantX[1],qantX[2],(qantX[2]-qantX[1])/100)
	qantX=quantile(X,weight)
	WXi=sapply(X,dunif,qantX[1],qantX[2])*(qantX[2]-qantX[1])
    #Calculation of the bandwidth h2 for nonparam estimation of the bias  term b.
    #Global polynomial estimation of th24spline
	theta24=sum(th24spline*WXi)/N
	#Constants related to the kernel K
	C2K=(3/(8*sqrt(pi)))^(1/7)
	if(theta24>0)
		C2K=(15/(16*sqrt(pi)))^(1/7)

	#Estimation of the variance term v
	#Find the index of the first component of each group
	cumsom=cumsum(nj)+1
	cumsom=cumsom[-J]
	cumsom=c(1,cumsom)
	mult=q^(-nj)*mean(ZZ)
	intvar=0

	sumwi=0
	sumwi2=0
	for(i in 1:max(nj)){
		#keep only groups of size nj>=i
		condnj=(nj>=i)
		Ji=sum(condnj)
		wi=sqrt(Ji)
		sumwi=sumwi+wi
		sumwi2=sumwi2+1/wi
		cumsomi=cumsom[condnj]
		indices=cumsomi+i-1
		indord=order(X[indices])
		Xord=X[indices[indord]]
		Tord=mult[condnj]*ZZ[indices[indord]]
		intvar=intvar+wi*sum(Tord[1:(Ji-1)]*(1-Tord[2:Ji])*(Xord[2:Ji]-Xord[1:(Ji-1)]))
	}
	intvar=intvar/sumwi

	#Bandwidth h2
	h2=C2K*(intvar/abs(theta24))^(1/7)*(sumwi/sumwi2)^(-1/7)
	#Calculation of the nonpar estimator of the bias term b.
	b.np=0
	ZZ=ZZ*q^(-Nj)*mean(ZZ)

	for(i in 1:max(nj)){
		#keep only groups of size nj>=i
		condnj=(nj>=i)
		Ji=sum(condnj)
		wi=sqrt(Ji)
		cumsomi=cumsom[condnj]
		indices=cumsomi+i-1

		Xi=X[indices]
		Zi=ZZ[indices]
		qantX=quantile(Xi,weight)
		WXi=dunif(Xi,qantX[1],qantX[2])*(qantX[2]-qantX[1])
		g2=g.lpderfit(Xi,Xi,Zi,h2,2,3)
		b.np=b.np+sum(g2^2*WXi)/wi	
	}
	b.np=b.np/sumwi
	#Calculation of the PI bandwidth.
	Num=RK*intvar
	Den=muK2^2*b.np
	h=(Num/Den)^(1/5)*N^(-1/5)
	h
}
#---------------------------------------------------------------
#calculates the unweighted ROT
#--------------------------------------------------------------
h.ROT=function(X,ZZ,q){

      N=length(X)
      muK2=1
      RK=1/(2*sqrt(pi))

      #Gobal polynomial estimation of the bias term
      g2.globpoly=globpoly(X,X,ZZ,2,3,q)
      b.globpoly=sum(g2.globpoly^2)

	#Estimation of the variance term

	cumsom=cumsum(nj)+1
	cumsom=cumsom[-J]
	cumsom=c(1,cumsom)
	mult=q^(-nj)*mean(ZZ)
	intvar=0

	for(i in 1:max(nj))
	{
		indices=cumsom+i-1
		indord=order(X[indices])
		Xord=X[indices[indord]]
		Tord=mult*ZZ[indices[indord]]
		intvar=intvar+sum(Tord[1:(J-1)]*(1-Tord[2:J])*(Xord[2:J]-Xord[1:(J-1)]))
	}

	intvar=intvar/max(nj)

	#Bandwidth
	Num=RK*intvar
	Den=muK2^2*b.globpoly
	h=(Num/Den)^(1/5)
	h
}

#---------------------------------------------------------------
#calculates the weighted ROT
#--------------------------------------------------------------
h.ROT2=function(X,Y,q){

      N=length(X)
      muK2=1
      RK=1/(2*sqrt(pi))

      #Gobal polynomial estimation of the bias term
      g2.globpoly=globpoly(X,X,Y,2,3,q)
      qantX=quantile(X,c(0.1,0.9))
      wofXi=dunif(X,qantX[1],qantX[2])*(qantX[2]-qantX[1])
      b.globpoly=sum(g2.globpoly^2*wofXi)

	#Estimation of the variance term

	cumsom=cumsum(nj)+1
	cumsom=cumsom[-J]
	cumsom=c(1,cumsom)
	mult=q^(-nj)*mean(Y)
	intvar=0

	for(i in 1:max(nj))
	{
		indices=cumsom+i-1
		indord=order(X[indices])
		Xord=X[indices[indord]]
		Tord=mult*Y[indices[indord]]
		intvar=intvar+sum(Tord[1:(J-1)]*(1-Tord[2:J])*(Xord[2:J]-Xord[1:(J-1)]))
	}

	intvar=intvar/max(nj)

	#Bandwidth
	Num=RK*intvar
	Den=muK2^2*b.globpoly
	h=(Num/Den)^(1/5)
	h

}

