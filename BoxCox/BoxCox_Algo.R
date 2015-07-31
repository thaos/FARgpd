library(MASS)
library(tseries)
library(evmix)
library(extRemes)
source("../Oridea/Oridea_Algo.R", chdir=TRUE)
g=function(x,lambda){
  y=NA
  if(lambda!=0){y=(x^lambda-1)/lambda}
  if(lambda==0){y=log(x)}
  return(y)
} 

gauss_lik <- function(mle,covariate,y){
	#         print("*****************************************")
	require(mvtnorm)
	mu0=mle[1]
	mu1=mle[2]
	sig0=mle[3]
	sig1=mle[4]
	mu=mu1*covariate+mu0
	#         print(system.time({
	n=length(mu)
	sig2=sig1*covariate+sig0
	negll=-0.5*(-n*log(2*pi)-sum(log(sig2))-sum((y-mu)^2/sig2))
	if(is.na(negll)) negll=10^6
	negll
		#         }))
	#         print(system.time({
	#                 sig2=diag(sig1*covariate+sig0)
	#                 negdmv=-dmvnorm(y,mu,sig2,log=TRUE)
	#         }))
	#         print(negdmv)
		#         print(system.time({
		#                 sig=sqrt(sig1*covariate+sig0)
		#                 negloglik=-sum(mapply(dnorm,x=y,mean=mu,sd=sig,MoreArgs=list(log=TRUE)))
		#                 if(is.na(negloglik)) negloglik=10^6
		#                 print(negloglik)
		#         }))
		#         print(negll)
		#         print(negloglik)
		#         print(all.equal(negloglik,negll))
		#         negloglik
	#                         sig=diag(sig)
	#                         negloglik=-sum(mapply(dnorm,x=y,mean=mu,sd=sig,MoreArgs=list(log=TRUE)))
	#                         negloglik
}

gauss_lik_sansvar <- function(mle,covariate,y){
	require(mvtnorm)
	mu0=mle[1]
	mu1=mle[2]
	sig2=diag(mle[3],length(covariate))
	mu=mu1*covariate+mu0
	negdmv=-dmvnorm(y,mu,sig2,log=TRUE)
	negdmv
	#                         sig=diag(sig)
	#                         negloglik=-sum(mapply(dnorm,x=y,mean=mu,sd=sig,MoreArgs=list(log=TRUE)))
	#                         negloglik
}

lambda_lik_sansvar=function(lambda,covariate,y){
	#         print("*******************************")
	#         print("lambda")
	#         print(lambda)
	#         print("*******************************")
	ylambda=g(y,lambda)
	lm.fit=lm(ylambda~covariate)
	lm.fit.s=summary(lm.fit)
	mle <- c(coefficients(lm.fit),lm.fit.s$sigma^2)
	#         print(lm.fit)
	#         print(logLik(lm.fit))
	#         print(mle)
	#         print("FIT2")
	#         print(gauss_lik_sansvar(mle,covariate,ylambda))
	y.fit2=nlminb(start=mle,gauss_lik_sansvar,covariate=covariate,y=ylambda)
	y.fit2=within(y.fit2,{objective=objective+(1-lambda)*sum(log(y))})
	#         if(y.fit2$convergence!=0 )
	#                 y.fit2$objective=10^6
	#         y.fit2$objective=logLik(lm.fit)
	       #         e=lm.fit$residuals
	       #         elY=(exp(mean(log(y))))
	       #         loglik=-n/2*log(lm.fit.s$sigma^2)+(lambda-1)*sum(log(y))
	       #         print(loglik)
	       #         loglik=-n/2*log(sum((e/elY^lambda)^2))
	       #         print(loglik)
	#         y.fit2$objective=logLik(lm.fit)
	       #         y.fit2$objective=loglik
	#         print(y.fit2)
	y.fit2
}

lambda_lik=function(lambda,covariate,y){
	#         print("*******************************")
	#         print("lambda")
	#         print(lambda)
	#         print("*******************************")
	ylambda=g(y,lambda)
	lm.fit=lm(ylambda~covariate)
	lm.fit.s=summary(lm.fit)
	fit_res=residuals(lm.fit)
	var.fit=lm(fit_res^2~covariate)
	#                 mle=c(coef(lm.fit),lm.fit.s$sigma^2,1)
	mle <- c(coefficients(lm.fit),coefficients(var.fit))
	#         print(mle)
	#         print("FIT2")
	#         print(gauss_lik(mle,covariate,ylambda))
	y.fit2=nlminb(start=mle,gauss_lik,covariate=covariate,y=ylambda)
	y.fit2=within(y.fit2,{objective=objective+(1-lambda)*sum(log(y))})
	#         print(y.fit2)
	y.fit2
}

bc.fit_sansvar=function(covariate,y,range=seq(0.5,1.5,by=0.01)){
	y <- y / exp(mean(log(y)))
	#         browser()
	L=range
	#         LV=Vectorize(lambda_lik)(L)
	LV=-unlist(lapply(L,function(L)lambda_lik_sansvar(L,covariate,y)$objective))
	overallmax <- max(LV)
	crit <- overallmax - qchisq(0.999, 1)/2
	cond <- LV > crit
	Lc <- L[cond]
	LVc <- LV[cond]
	plot(L,LV,type="l",ylab="")
	boxcox(y~covariate,lambda=seq(-10,10,.1))
	#         browser()
	lambda=(maxf=optimize(function(L)lambda_lik_sansvar(L,covariate,y)$objective,range(L),maximum=FALSE))$minimum
	print(lambda)
	res=lambda_lik_sansvar(lambda,covariate,y)
	res$lambda=lambda
	res
}


bc.fit_old=function(covariate,y,range=seq(0.5,1.5,by=0.01)){
	#         browser()
	L=range
	#         LV=Vectorize(lambda_lik)(L)
	LV=-unlist(lapply(L,function(L)lambda_lik(L,covariate,y)$objective))
	overallmax <- max(LV)
	crit <- overallmax - qchisq(0.999, 1)/2
	cond <- LV > crit
	Lc <- L[cond]
	LVc <- LV[cond]
	plot(Lc,LVc,type="l",ylab="")
	#         browser()
	lambda=(maxf=optimize(function(L)lambda_lik(L,covariate,y)$objective,range(L),maximum=FALSE))$minimum
	#         print(lambda)
	res=lambda_lik(lambda,covariate,y)
	res$lambda=lambda
	res
}

bc.fit=function(covariate,y,range=seq(0.5,1.5,by=0.01),to.plot=FALSE){
	#for the stability of the boxcox trannsformation
	y <- y / exp(mean(log(y)))
	#         browser()
	L=range
	if(to.plot){
		#         LV=Vectorize(lambda_lik)(L)
		LV=-unlist(lapply(L,function(L)lambda_lik(L,covariate,y)$objective))
		overallmax <- max(LV)
		crit <- overallmax - qchisq(0.999, 1)/2
		cond <- LV > crit
		Lc <- L[cond]
		LVc <- LV[cond]
		par(mfrow=c(1,2))
		plot(Lc,LVc,type="l",ylab="")
		boxcox(y~covariate,lambda=L)
	}
	lambda=(maxf=optimize(function(l)lambda_lik(l,covariate,y)$objective,range(L),maximum=FALSE))$minimum
	#         print(lambda)
	res=lambda_lik(lambda,covariate,y)
	res$lambda=lambda
	res
}


zp2yp <- function(lambda,mu,sigma,zp) 
	(lambda*(mu+sigma*zp)+1)^(1/lambda)


bc.fit.original=function(t,y,range=seq(-1,1,by=0.01)){
	logv.original=function(lambda){
		ylambda=g(y,lambda)
		lm.fit=lm(ylambda~t)
		e=lm.fit$residuals
		sigma=summary(lm.fit)$sigma
		res=-n/2*log(2 * pi)-n*(log(sigma))-.5/sigma^2*(sum(e^2))+(lambda-1)*sum(log(y))
		res
	}
	L=range
	#         LV=Vectorize(logv)(L)
	LV=unlist(lapply(L,logv.original))
	plot(L,LV,type="l",ylab="")
	lambda=(maxf=optimize(logv.original,range(L),maximum=TRUE))$maximum
	lambda
}

getFAR.boxcox <- function(pt0,pt1,xp,ydat,L=seq(-10,10,0.1),to.plot=FALSE){
	if(any(ydat$y<= 0)){
		offs=+ceiling(abs(min(ydat$y)))
		ydat$y=ydat$y+offs
		xp=xp+offs
	}
	covariate=ydat$mua
	sf=exp(mean(log(ydat$y)))
	#         res=with(ydat,bc.fit(mua,y,range=prange,to.plot=TRUE))
	res=with(ydat,bc.fit(mua,y,range=L))
	lambda=res$lambda
	print(lambda)
	#         res=with(ydat,lambda_lik(lambda,covariate,y/sf))
	#         res$lambda=lambda
	ylambda=g(ydat$y/sf,lambda)
	xplambda=g(xp/sf,lambda)
	mu.pred=with(res,par[1]+par[2]*covariate)
	sigma2.pred=with(res,par[3]+par[4]*covariate)
	mu.pred0=with(res,par[1]+par[2]*pt0[2])
	mu.pred1=with(res,par[1]+par[2]*pt1[2])
	residus=(ylambda-mu.pred)/sqrt(sigma2.pred)
	sigma2.pred0=with(res,par[3]+par[4]*pt0[2])
	sigma2.pred1=with(res,par[3]+par[4]*pt1[2])
	r0=(xplambda-mu.pred0)/sqrt(sigma2.pred0)
	r1=(xplambda-mu.pred1)/sqrt(sigma2.pred1)
	#         print(shapiro.test(residus))
	#         print(Box.test(residus))
	#         print(PP.test(c(residus)))
	#         print(adf.test(residus))
	tc=tcplot_sthao(residus)
	# Easton avec threshold fixe 
	threshold=select.mu(tc)
	if(is.finite(threshold)){
		print("Path 1")
		rate=mean(residus>=threshold)
		residus.fit=fevd(residus,threshold=threshold,type="GP",time.units="years")
		mle=residus.fit$results$par
		p0=pevd(r0, scale=mle[1], shape=mle[2], threshold=threshold, type="GP",lower.tail=FALSE)
		p1=pevd(r1, scale=mle[1], shape=mle[2], threshold=threshold, type="GP",lower.tail=FALSE)
	}
	#         print(c(r0,r1))
	#         if(r0<threshold) p0=mean(r0<residus) # ou pnorm(r0,lower.tail=FALSE
	#         if(r0<threshold) p0=mean(r0<residus) # ou pnorm(r0,lower.tail=FALSE
		print("Direct 2 Path 2")
	if(r0<threshold) p0=pnorm(r0,lower.tail=FALSE)
	if(r1<threshold) p1=pnorm(r1,lower.tail=FALSE)
	if(to.plot){
		par(mfrow=c(1,3))
		plot(ydat$year,ydat$y)
		plot(ydat$year,ylambda)
		plot(ydat$year,residus)
	}
	res=FAR(p0,p1)
	res
}
