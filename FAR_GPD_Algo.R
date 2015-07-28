packages <- c("boot", "quantreg", "evmix", "ismev", "parallel","foreach","doMC","snow")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
	  install.packages(setdiff(packages, rownames(installed.packages())),repos="http://cran.us.r-project.org")
}
library(boot)
library(quantreg)
library(evmix)
library(ismev)
library(extRemes)
library(parallel)
library(foreach)
library(doMC)
library(snow)

getP <- function(p2,y.fit,ydat,to.plot=FALSE){
	#         print(p2)
	tcord=p2[1]
	ccord=p2[2]
	ycord=p2[3]
	mcord=p2[4]
	mle <- y.fit$results$par
	mu=mcord
	sigma=mle[1]+ccord*mle[2]
	shape=mle[3]
	phi=y.fit$rate
	if (ycord>mu)
		res=pevd(ycord,threshold=mu,scale=sigma,shape=shape,type="GP",lower.tail=FALSE)*phi
	else {
		findP <- function(par,ccord,ycord){
			rq.fitted=rq(y~mua,data=ydat,tau=par)
			ndata=data.frame("mua"=ccord)
			predicted=predict.rq(rq.fitted,newdata=ndata)
			#                         print((ycord-predicted)^2)
			#                         print(par)
			abs(ycord-predicted)
		}
		res.optim=optimize(findP,interval=c(0,1-phi),ccord=ccord,ycord=ycord,tol=0.001)
		to.plot=TRUE
		res=res.optim$minimum
		if(to.plot){
			with(ydat,plot(year,y,col="grey"))
			lines(unique(cbind(ydat$year,predict(rq(y~mua,tau=res,data=ydat)))),col="red")
			abline(h=mu)
			points(tcord,ycord,lwd=2,col="blue")
		}
		res=1-res
	}		
	res=c(res,mu,sigma,shape)
	names(res)=c("p","mu","sigma","shape")
	res
}
# getP(c(1920,110),y.fit,ydat,to.plot=TRUE)


getFAR <- function(p1,p2,y.fit,ydat){
	prob1=getP(p1,y.fit,ydat)
	prob2=getP(p2,y.fit,ydat)
	if(prob1[1]==0 & prob2[1]==0)
		FAR=1
	else
	FAR=1-(prob2[1]/prob1[1])
	res=c(FAR,prob2,prob1)
	names(res)=c("FAR","p0","mu0","sc0","sh0","p1","mu1","sc1","sh1")
	#         print(res)
	res
}
# getFAR(c(1920,110),2102,y.fit,ydat)

getFAR.theo=function(xp,i0,i1,mu,sigma){
	p0=1-pnorm(xp,mean=mu[i0],sd=sigma[i0])
	p1=1-pnorm(xp,mean=mu[i1],sd=sigma[i1])
	if(p0==0 & p1==0)
		p0p1=1
	else
		p0p1=p0/p1
	res=1-p0p1
	attr(res,"p0")=p0
	attr(res,"p1")=p1
	res
}

FARBoot <- function(ydat,indice,p1,x2){
	y.fit.dat=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,method="MLE")
	init=as.list(y.fit.dat$results$par)
	data.b=ydat[indice,]
	#         y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga)
	y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga,method="MLE",initial=init)
	getFAR(p1,x2,y.fit,ydat)
}

FARBoot <- function(ydat,indice,qthreshold,pt1,pt0){
	rq.fitted=rq(y~mua,data=ydat,tau=qthreshold)
	threshold=predict(rq.fitted)
	y.fit.dat=fevd(y,ydat,threshold=threshold,scale.fun=~siga,type="GP",method="MLE")
	init=as.list(y.fit.dat$results$par)
	data.b=ydat[indice,]
	data.b=data.b[order(data.b$year),]
	rq.fitted.b=rq(y~mua,data=data.b,tau=qthreshold)
	threshold.b=predict(rq.fitted.b)
	ndata=data.frame("mua"=c(pt0[2],pt1[2]))
	m=predict(rq.fitted.b,ndata)
	m0=m[1]
	m1=m[2]
	pt0=c(pt0,m0)
	pt1=c(pt1,m1)
	y.fit=fevd(y,data.b,threshold=threshold,scale.fun=~siga,method="MLE",type="GP",initial=init)
	getFAR(pt1,pt0,y.fit,data.b)
}

FARBoot.bc <- function(ydat,indice,qthreshold,pt1,pt0){
	data.b=ydat[indice,]
	data.b=data.b[order(data.b$year),]
	res=getFAR.boxcox(pt0,pt1,xp,data.b,to.plot=FALSE)
	res
}

FARBoot.or <- function(ydat,indice,qthreshold,pt1,pt0){
	data.b=ydat[indice,]
	data.b=data.b[order(data.b$year),]
	far.or=getFAR.oridea(pt0,pt1,xp,ydat,to.plot=FALSE)
}


FARBoot_explore <- function(ydat,indice,p1,x2){
	data.b=ydat[indice,]
	y.fit.dat=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,method="MLE")
	init=as.list(y.fit.dat$results$par)
	y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga,method="MLE",initial=init)
	res=getFAR(p1,x2,y.fit,ydat)
	i0=min(which((abs(ydat$year-t0))==min(abs(ydat$year-t0))))
	i1=min(which((abs(ydat$year-t1))==min(abs(ydat$year-t1))))
	r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu=ydat$mu,sigma=ydat$sigma,xi=ydat$shape)
	#         with(ydat,plot(year,y,yaxt="n",col="red",ylim=c(100,108)))
	#         par(new=TRUE)
	#         with(data.b,hist(year,breaks=250))
	#         par(new=TRUE)
	#         with(data.b,plot(year,y,yaxt="n",col="green", ylim=c(100,108)))
	if (abs(res[1] - (r.theo[1])) > 0.2){
		print("chelou")
		print(res)
		#                 browser()
	}
	if (abs(res[1] - (r.theo[1])) <= 0.2) {
		print("normal")
		print(res)
		#                 browser()
	}
	res
}

FARBoot.Spline <- function(ydat,indice,p1,x2){
	data.b=ydat[indice,]
	covariate=with(data.b,predict(rq( y~ bs(year, df=3,degree=3), tau=0.5)))
	data.b$mua=covariate
	data.b$siga=covariate
	y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga)
	getFAR(p1,x2,y.fit,ydat)
}

function2draw <- function(condition){
	function(ydat,nrepet,year){
		cond.e=eval(condition)
		sb=which(cond.e)
		sample(sb,size=nrepet,replace=TRUE)
	}
}


rbind_list <- function(liste){
	if(is.data.frame(liste)|length(liste)==1)
	   return(liste)
	if(length(liste)>1){
		liste[[2]]= rbind(liste[[1]],liste[[2]])
		return(rbind_list(liste[-1]))
	}
}

cond=quote(ydat$year==year)
draw_same_year <- function2draw(cond)
		
cond2=function(tol){
	substitute(ydat$year <= year+tol &ydat$year >= year-tol,list(tol=tol))
}
draw_around_year <- function2draw(cond2(5))

year2mu <- function(year,ydat){
	i=which(year==ydat$year)
	ydat[i[1],"mua"]
}
cond3 <- function(tol){
	substitute(ydat$mua <= year2mu(year,ydat)+tol &ydat$mua >= year2mu(year,ydat)-tol,list(tol=tol))
}
draw_around_mua <- function2draw(cond3(0.1))

FARBoot_gen  <- function(to_draw){
	function(ydat,indice,p1,x2){
		y.fit.dat=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,method="MLE")
		init=as.list(y.fit.dat$results$par)
		years=aggregate(y~year,data=ydat,FUN=length)
		names(years)=c("year","eff")
		indice=c(with(years,mapply(to_draw,nrepet=eff,year=year,MoreArgs=list(ydat=ydat),SIMPLIFY=TRUE)))
		data.b=ydat[indice,]
		years.v=c(with(years,mapply(rep,x=year,times=eff)))
		data.b$year=years.v
		y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga,method="MLE",initial=init)
		getFAR(p1,x2,y.fit,ydat)
	}
}
FARBoot_ax=FARBoot_gen(draw_around_mua) 
FARBoot_sy=FARBoot_gen(draw_same_year) 
FARBoot_ay=FARBoot_gen(draw_around_year) 
