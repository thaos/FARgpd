source("FAR_GPD_Algo.R")

Coverage <- function(Code2Execute,n=100){
	function(env){
		Code2Execute=get("Code2Execute",parent.env(environment()))
		assign("Code2Execute",Code2Execute,envir=env)
		res=t(sapply(1:n,function(x){
			   if(x==n)
				   Code2Execute(env,to.plot=TRUE)
			   else{
				   print(x)
				   Code2Execute(env)
			   }
			}
	))
		res=as.data.frame(res)
		invisible(res)
	}
}

make_Boot2Execute <- function(f_boot){
	function(env,to.plot=FALSE){
		list2env(as.list(env),envir=environment())
		y=mapply(rnorm,1,mean=mu,sd=sigma)
		ydat=data.frame(year,y,rep(1,n*repet),covariate,rep(1,n*repet),covariate,rep(1,n*repet))
		names(ydat)=c("year","y","mub","mua","sigb","siga","xi")
		i0=min(which.min((abs(ydat$year-t0))))
		i1=min(which.min((abs(ydat$year-t1))))
		c0=ydat$mua[i0]
		c1=ydat$mua[i1]
		pt0=c(t0,c0,xp)
		pt1=c(t1,c1,xp)
		boot.res=boot(data=ydat,statistic=f_boot,R=250,pt1=pt1,pt0=pt0,qthreshold=qthreshold, parallel="snow")
		theta.boot=colMeans(boot.res$t)
		r.boot=mean(boot.res$t)
		r.boot=theta.boot[1]
		if(is.na(r.boot))
			browser()
		alpha=0.05
		ic.boot=quantile(boot.res$t,p=c(alpha/2,1-alpha/2))
		ic.boot=apply(boot.res$t,2,quantile,p=c(alpha/2,1-alpha/2))
		boot.ic=c(ic.boot[1,1],r.boot,ic.boot[2,1])
		r.theo=getFAR.theo(xp=xp,i0=i0,i1=i1,mu,sigma)
		names(boot.ic) <- c("LowerCI", "Estimate", "UpperCI")
		res=c(r.theo,boot.ic)
		names(res)[1]="Theoric"
		res
	}
} 
Boot2Execute <- make_Boot2Execute(FARBoot)
Boot2Execute.bc <- make_Boot2Execute(FARBoot.bc)
Boot2Execute.or <- make_Boot2Execute(FARBoot.or)
#Boot2Execute_sy <- make_Boot2Execute(FARBoot_sy)
#Boot2Execute_ay <- make_Boot2Execute(FARBoot_ay)
#Boot2Execute_ax <- make_Boot2Execute(FARBoot_ax)


Prof2Execute <- function(env,to.plot=FALSE){
	list2env(as.list(env),envir=environment())
	y=mapply(revd,1,mu,sigma,shape)
	covariate=covariate-mean(covariate)
	ydat=data.frame(year,y,rep(1,n*repet),covariate,rep(1,n*repet),covariate,rep(1,n*repet))
	names(ydat)=c("year","y","mub","mua","sigb","siga","xi")
	y.fit=fevd(y,ydat,location.fun=~mua,scale.fun=~siga)
	i0=min(which((abs(ydat$year-t0))==min(abs(ydat$year-t0))))
	i1=min(which((abs(ydat$year-t1))==min(abs(ydat$year-t1))))
	r.ic=gpd.ratio.ic(xp=xp,t0=i0,t1=i1,y.fit=y.fit,qthreshold=qthreshold,ydat=ydat)
	#         r.ic=r.ic[1:3]
	#         r.ic=r.ic[c(3,2,1)]
	r.ic=r.ic[3:1]
	r.ic=1-r.ic
	names(r.ic) <- c("LowerCI", "Estimate", "UpperCI")
	r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu,sigma)
	# print(r.theo)
	# print(r.ic)
	res=c(r.theo,r.ic)
	names(res)[1]="Theoric"
	print(res)
	res
}


