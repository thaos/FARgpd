library(MASS)
library(evmix)
library(extRemes)

g=function(x,lambda){
  y=NA
  if(lambda!=0){y=(x^lambda-1)/lambda}
  if(lambda==0){y=log(x)}
  return(y)
} 


n=350
t=seq(1,n)
mu=(t>100)*.01*(t-100)+100
sigma2=(t>100)*.01*(t-100)+1
sigma2=1
y=mapply(rnorm,1,mu,sqrt(sigma2))
covariate=c(scale(mu))

gauss_lik <- function(mle,covariate,y){
	require(mvtnorm)
	mu0=mle[1]
	mu1=mle[2]
	sig0=mle[3]
	sig1=mle[4]
	mu=mu1*covariate+mu0
	sig2=diag(sig1*covariate+sig0)
	negdmv=-dmvnorm(y,mu,sig2,log=TRUE)
	negdmv
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
	print("*******************************")
	print("lambda")
	print(lambda)
	print("*******************************")
	ylambda=g(y,lambda)
	lm.fit=lm(ylambda~covariate)
	lm.fit.s=summary(lm.fit)
	mle <- c(coefficients(lm.fit),lm.fit.s$sigma^2)
	print(lm.fit)
	print(logLik(lm.fit))
	print(mle)
	print("FIT2")
	print(gauss_lik_sansvar(mle,covariate,ylambda))
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
	print("*******************************")
	print("lambda")
	print(lambda)
	print("*******************************")
	ylambda=g(y,lambda)
	lm.fit=lm(ylambda~covariate)
	lm.fit.s=summary(lm.fit)
	fit_res=residuals(lm.fit)
	var.fit=lm(fit_res^2~covariate)
	#                 mle=c(coef(lm.fit),lm.fit.s$sigma^2,1)
	mle <- c(coefficients(lm.fit),coefficients(var.fit))
	print(mle)
	print("FIT2")
	print(gauss_lik(mle,covariate,ylambda))
	y.fit2=nlminb(start=mle,gauss_lik,covariate=covariate,y=ylambda)
	y.fit2=within(y.fit2,{objective=objective+(1-lambda)*sum(log(y))})
	print(y.fit2)
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
	browser()
	lambda=(maxf=optimize(function(L)lambda_lik_sansvar(L,covariate,y)$objective,range(L),maximum=FALSE))$minimum
	print(lambda)
	res=lambda_lik_sansvar(lambda,covariate,y)
	res$lambda=lambda
	res
}

ylambda=y
m.fit=lm(ylambda~covariate)
lm.fit.s=summary(lm.fit)
fit_res=residuals(lm.fit)
var.fit=lm(fit_res^2~covariate)
#                 mle=c(coef(lm.fit),lm.fit.s$sigma^2,1)
mle <- c(coefficients(lm.fit),coefficients(var.fit))
print(mle)
print("FIT2")
print(gauss_lik(mle,covariate,ylambda))
y.fit2=nlminb(start=mle,gauss_lik,covariate=covariate,ylambda=ylambda)
print(y.fit2)
y.fit2



bc.fit=function(covariate,y,range=seq(0.5,1.5,by=0.01)){
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
	browser()
	lambda=(maxf=optimize(function(L)lambda_lik(L,covariate,y)$objective,range(L),maximum=FALSE))$minimum
	print(lambda)
	res=lambda_lik(lambda,covariate,y)
	res$lambda=lambda
	res
}
bc.fit=function(covariate,y,range=seq(0.5,1.5,by=0.01)){
	#for the stability of the boxcox trannsformation
	y <- y / exp(mean(log(y)))
	#         browser()
	L=range
	#         LV=Vectorize(lambda_lik)(L)
	LV=-unlist(lapply(L,function(L)lambda_lik(L,covariate,y)$objective))
	overallmax <- max(LV)
	crit <- overallmax - qchisq(0.999, 1)/2
	cond <- LV > crit
	Lc <- L[cond]
	LVc <- LV[cond]
	plot(L,LV,type="l",ylab="")
	boxcox(y~covariate,lambda=seq(-10,10,.1))
	browser()
	lambda=(maxf=optimize(function(L)lambda_lik(L,covariate,y)$objective,range(L),maximum=FALSE))$minimum
	print(lambda)
	res=lambda_lik(lambda,covariate,y)
	res$lambda=lambda
	res
}


zp2yp <- function(lambda,mu,sigma,zp) 
	(lambda*(mu+sigma*zp)+1)^(1/lambda)

sf=exp(mean(log(y)))
res=bc.fit(covariate,y,range=seq(-10,10,.1))
lambda=res$lambda
res=lambda_lik(lambda,covariate,y/sf)
res$lambda=lambda
res
ylambda=g(y/sf,lambda)
mu.pred=with(res,par[1]+par[2]*covariate)
sigma2.pred=with(res,par[3]+par[4]*covariate)
residus=(ylambda-mu.pred)/sqrt(sigma.pred)
print(shapiro.test(residus))
print(Box.test(residus))
print(PP.test(c(residus)))
print(adf.test(residus))
tc=tcplot(residus,c(1,max(residus)))
tc=tcplot(residus,c(quantile(residus,0.5),max(residus)))
tc=threshrange.plot(residus,c(quantile(residus,0.5),max(residus)))
# Easton avec threshold fixe 
threshold=0.009
rate=mean(residus>=threshold)
ydat=cbind(rep(1,n),rep(1,n))
residus.fit=fevd(residus,threshold=threshold,type="GP",time.units="years")
muscsh=findpars(residus.fit)
muscsh$mu=rep(threshold,length(muscsh$scale))
z95f=evmix::qgpd(rep(0.95,length(muscsh$mu)),u=muscsh$mu,sigmau=muscsh$scale,xi=muscsh$shape,phiu=rate)
y95f=zp2yp(lambda,mu=mu.pred,sigma=sigma.pred,zp=z95f)*sf
q.theo95=mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=0.95))


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
getparam <- function(residus.fit,ydat,sigl){
