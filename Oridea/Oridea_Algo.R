library(extRemes)
library(evmix)
library(quantregGrowth)
source("../tcplot_sthao.R") 

getP.or <- function(p,m){
	x=p[1]
	y=p[2]
	coeff=m$coefficients
	pred=coeff[1,]+x*coeff[2,]
	dif=abs(y-pred)
	i=which.min(dif)
	#         if(dif[i]>10^-1){
	#                 print(dif[i])
	#                 stop("Cant find the point : diff too big")
	#         }
	as.numeric(colnames(coeff)[i])
}
	
getPs <- function(x,y,m){
	mapply(function(x,y)getP.or(c(x,y),m),x=x,y=y)
}

send2x0 <- function(x0,m,alpha){
	coeff=m$coefficients
	extract  <- function(alpha,coeff){
		ec=coeff[,paste(alpha)]
		ec[1]+x0*ec[2]
	}
	unlist(lapply(alpha,extract,coeff=coeff))
}

select.mu <- function(tc){
	i.na=which(apply(tc,1,function(x)any(is.na(x))))
	ifelse(length(i.na)!= 0,{tcc=tc[-c(min(i.na):nrow(tc)),]},{tcc=tc})
	#         ifelse(length(i.na)!= 0,{tcc=tc[-i.na,]},{tcc=tc})
	mu2rm=which(apply(tcc,1,function(x)x["nu"]<11))
	if(length(mu2rm)!=0) tcc=tcc[-mu2rm,]
	list2env(tcc,env=environment())
	found=FALSE
	i=1
	while(!found){
		if(i> length(u)){
			par(mfrow=c(2,1))
			with(tc,plot(u,xi))
			with(tc,lines(u,cil.xi,col="red"))
			with(tc,lines(u,ciu.xi,col="red"))
			with(tc,plot(u,mod.sigmau))
			with(tc,lines(u,mod.cil.sigmau,col="red"))
			with(tc,lines(u,mod.ciu.sigmau,col="red"))
			browser()
			stop("no threshold found")
		}
		mu.cur=u[i]
		sig.cur=mod.sigmau[i]
		xi.cur=xi[i]
		i2r=seq(1,i,1)
		tsig=all(sig.cur>mod.cil.sigmau[-i2r] & sig.cur < mod.ciu.sigmau[-i2r])
		txi=all(xi.cur>cil.xi[-i2r] & xi.cur < ciu.xi[-i2r])
		found=tsig&txi
		if(found == FALSE) i=i+1
	}
	mu.cur
}

getRP1 <- function(x,y,m,x0,threshold,gpd.fit){
	p=getP.or(c(x,y),m)
	s=send2x0(x0,m,p)
	if(s < threshold){
		message("value inferior to threshold")
		return(1-p)
	}
	param=gpd.fit$results$par
	sc=param[1]
	sh=param[2]
	phi=gpd.fit$rate
	res=pevd(s, scale=sc, shape=sh, threshold=threshold, type="GP",lower.tail=FALSE)
	res=res*phi
	res
}
getRP0 <- function(y,m,x0,threshold,gpd.fit){
	param=gpd.fit$results$par
	sc=param[1]
	sh=param[2]
	phi=gpd.fit$rate
	res=pevd(y, scale=sc, shape=sh, threshold=threshold, type="GP",lower.tail=FALSE)
	res=res*phi
	res
}

getFAR.oridea <- function(pt0,pt1,xp,ydat,to.plot=FALSE){
	pas=.01
	tauss<-seq(pas,1-pas,by=pas) #fix the percentiles of interest
	# Marche seulement si il y a une corrélation postive entre la variance des données et la covarriable ?
	my<-gcrq(y~mua, tau=tauss, data=ydat) #linear 
	alpha=with(ydat,getPs(mua,y,my))
	sent=send2x0(pt0[2],my,alpha) 
	tc=tcplot_sthao(sent)
	threshold=select.mu(tc)
	#         print(threshold)
	gpd.fit=fevd(sent,threshold=threshold,type="GP")
	p1=with(ydat,getRP1(pt1[2],xp,my,pt0[2],threshold=threshold,gpd.fit))
	p0=with(ydat,getRP0(xp,my,pt0[2],threshold=threshold,gpd.fit))
	if(to.plot){
		par(mfrow=c(2,2))
		with(ydat,plot(year,y))
		apply(my$coefficients,2,function(x)with(ydat,lines(year,x[1]+x[2]*mua)))
		with(ydat,plot(mua,y))
		apply(my$coefficients,2,function(x)with(ydat,lines(mua,x[1]+x[2]*mua)))
		with(my,apply(coefficients,1,function(y)plot(taus,y)))
	}
	FAR(p0,p1)
}
	
FAR <- function(p0,p1){
	if(p0==0 & p1==0)
		FAR=1
	else
	FAR=1-(p0/p1)
	attr(FAR,"p0")=p0
	attr(FAR,"p1")=p1
	FAR
}


