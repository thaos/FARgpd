library(extRemes)
library(evmix)
library(quantregGrowth)
source("../FARgpd/FAR_GPD_Algo.R")
source("../BoxCox/BoxCox_Algo.R")
source("proftable.R")
data(growthData) #load data
pas=.001
tauss<-seq(pas,1-pas,by=pas) #fix the percentiles of interest
m2<-gcrq(rev(y)~ps(x, mon=0, lambda=1000, pdiff=2), tau=tauss, data=growthData) #linear 
# Marche seulement si il y a une corrélation postive entre la variance des données et la covarriable ?
m2<-gcrq(y~x, tau=tauss, data=growthData) #linear 
with(growthData,plot(x,y))
apply(m2$coefficients,2,function(x)abline(b=x[2],a=x[1]))

par(mfrow=c(nrow(m2$coefficients),1))
with(m2,apply(coefficients,1,function(y)plot(taus,y)))

getP <- function(p,m){
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
	
getP(with(growthData[1,],c(x,y)),m2)

getPs <- function(x,y,m){
	mapply(function(x,y)getP(c(x,y),m),x=x,y=y)
}
alpha=apply(growthData,1,function(l){getP(c(l["x"],l["y"]),m2)})

send2x0 <- function(x0,m,alpha){
	coeff=m$coefficients
	extract  <- function(alpha,coeff){
		ec=coeff[,paste(alpha)]
		ec[1]+x0*ec[2]
	}
	unlist(lapply(alpha,extract,coeff=coeff))
}

sent=send2x0(0.5,m2,alpha) 
hist(sent)

select.mu <- function(tc){
	i.na=which(apply(tc,1,function(x)any(is.na(x))))
	ifelse(length(i.na)!= 0,{tcc=tc[-i.na,]},{tcc=tc})
	mu2rm=which(apply(tcc,1,function(x)x["nu"]<11))
	if(length(mu2rm)!=0) tcc=tcc[-mu2rm,]
	list2env(tcc,env=environment())
	found=FALSE
	i=1
	while(!found){
		if(i> length(u))
			stop("no threshold found")
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

changement=1950
years=1850:2200
n=length(unique(years))
changement=100
repet=4
t=seq(1,n)
mu=(t>changement)*.020*(t-changement)+100
sigma=(t>changement)*.001*(t-changement)+1
sigma=(t>changement)*.005*(t-changement)+1
sigma=sqrt(sigma)
shape=rep(-0.10,n)
mu=rep(mu,repet)
sigma=rep(sigma,repet)
t=rep(t,repet)
shape=rep(shape,repet)
year=rep(years,repet)
covariate=(mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=0.5)))
#On ajoute un bruit insignifiant pour éviter que la covariable soit exactementy identique sur la première période, sans quoi, ce la fait bugguer le profile de vraisemblance: division par zéro pour une différence en covariable nulle.
# covariate=covariate+rnorm(n*repet,sd=0.000001)
covariate=scale(covariate,scale=FALSE)
y=mapply(rnorm,1,mean=mu,sd=sigma)
ydat=data.frame(year,y,mu,sigma,shape,rep(1,n*repet),covariate,rep(1,n*repet),covariate,rep(1,n*repet))
names(ydat)=c("year","y","mu","sigma","shape","mub","mua","sigb","siga","xi")
ydat=ydat[order(year),]
		
	

pas=.001
tauss<-seq(pas,1-pas,by=pas) #fix the percentiles of interest
# Marche seulement si il y a une corrélation postive entre la variance des données et la covarriable ?
my<-gcrq(y~mua, tau=tauss, data=ydat) #linear 
with(ydat,plot(year,y))
apply(my$coefficients,2,function(x)with(ydat,lines(year,x[1]+x[2]*mua)))

with(ydat,plot(mua,y))
apply(my$coefficients,2,function(x)with(ydat,lines(mua,x[1]+x[2]*mua)))

par(mfrow=c(nrow(my$coefficients),1))
with(my,apply(coefficients,1,function(y)plot(taus,y)))

alpha=with(ydat,getPs(mua,y,my))
t0=1875
i0=min(which((abs(ydat$year-t0))==min(abs(ydat$year-t0))))
sent=send2x0(ydat$mua[i0],my,alpha) 
hist(sent)
ks.test(sent+rnorm(length(sent),sd=.00000001),"pnorm",mean=mean(sent),sd=sd(sent))
# tc=threshrange.plot(sent,c(quantile(sent,0.5),max(sent)))
tc=tcplot(sent,c(quantile(sent,0.5),max(sent)))
threshold=select.mu(tc)
gpd.fit=fevd(sent,threshold=threshold,type="GP")
lapply(findpars(gpd.fit),unique)

getRP1 <- function(x,y,m,x0,threshold,gpd.fit){
	p=getP(c(x,y),my)
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

t0=1900
i0=min(which((abs(ydat$year-t0))==min(abs(ydat$year-t0))))
t1=2150
xp=106
i1=min(which.min(abs(ydat$year-t1)))
p1=with(ydat,getRP1(mua[i1],xp,my,mua[i0],threshold=threshold,gpd.fit))
p0=with(ydat,getRP0(xp,my,mua[i0],threshold=threshold,gpd.fit))
r.oridea=FAR(p0,p1)
r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu=ydat$mu,sigma=ydat$sigma)
print(p0)
print(p1)
print(r.oridea)
print(r.theo)

getFAR.oridea <- function(i0,i1,xp,ydat,to.plot=FALSE){
	pas=.001
	tauss<-seq(pas,1-pas,by=pas) #fix the percentiles of interest
	# Marche seulement si il y a une corrélation postive entre la variance des données et la covarriable ?
	my<-gcrq(y~mua, tau=tauss, data=ydat) #linear 
	alpha=with(ydat,getPs(mua,y,my))
	sent=send2x0(ydat$mua[i0],my,alpha) 
	tc=tcplot(sent,c(quantile(sent,0.5),max(sent)))
	threshold=select.mu(tc)
	print(threshold)
	gpd.fit=fevd(sent,threshold=threshold,type="GP")
	p1=with(ydat,getRP1(mua[i1],xp,my,mua[i0],threshold=threshold,gpd.fit))
	p0=with(ydat,getRP0(xp,my,mua[i0],threshold=threshold,gpd.fit))
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
far.or=getFAR.oridea(i0,i1,xp,ydat,to.plot=TRUE)
	

covariate=ydat$mua
sf=exp(mean(log(ydat$y)))
res=with(ydat,bc.fit(mua,y,range=seq(-10,10,.1)))
lambda=res$lambda
res=with(ydat,lambda_lik(lambda,covariate,y/sf))
res$lambda=lambda
res
ylambda=g(ydat$y/sf,lambda)
xplambda=g(xp/sf,lambda)
mu.pred=with(res,par[1]+par[2]*covariate)
sigma2.pred=with(res,par[3]+par[4]*covariate)
residus=(ylambda-mu.pred)/sqrt(sigma2.pred)
print(shapiro.test(residus))
print(Box.test(residus))
print(PP.test(c(residus)))
print(adf.test(residus))
# tc=tcplot(residus,c(1,max(residus)))
tc=tcplot(residus,c(quantile(residus,0.5),max(residus)))
# tc=threshrange.plot(residus,c(quantile(residus,0.5),max(residus)))
# Easton avec threshold fixe 
threshold=select.mu(tc)
rate=mean(residus>=threshold)
residus.fit=fevd(residus,threshold=threshold,type="GP",time.units="years")
muscsh=lapply(findpars(residus.fit),unique)
muscsh$mu= threshold

	
	
FAR <- function(p0,p1){
	if(p0==0 & p1==0)
		FAR=1
	else
	FAR=1-(p0/p1)
	attr(FAR,"p0")=p0
	attr(FAR,"p1")=p1
	FAR
}

r0=(xplambda-mu.pred[i0])/sqrt(sigma2.pred[i0])
r1=(xplambda-mu.pred[i1])/sqrt(sigma2.pred[i1])
p1=with(muscsh,pevd(r1, scale=scale, shape=shape, threshold=mu, type="GP",lower.tail=FALSE))
p0=with(muscsh,pevd(r0, scale=scale, shape=shape, threshold=mu, type="GP",lower.tail=FALSE))
if(r0<threshold) p0=mean(r0>residus)
if(r1<threshold) p1=mean(r1>residus)
r.boxcox=FAR(p0,p1)
r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu=ydat$mu,sigma=ydat$sigma)
print(r.boxcox)
print(r.oridea)
print(r.theo)

mr<-gcrq(residus~ydat$mua, tau=tauss) #linear 
plot(ydat$year,residus)
apply(mr$coefficients,2,function(x)abline(b=x[2],a=x[1]))
apply(mr$coefficients,2,function(x)with(ydat,lines(mua,x[1]+x[2]*mua)))

par(mfrow=c(nrow(mr$coefficients),1))
with(mr,apply(coefficients,1,function(y)plot(taus,y)))

#Marche seulement pour variance évoluant linéairement
Rprof(tmp <- tempfile(), line.profiling=TRUE)
far.or=getFAR.oridea(pt0,pt1,xp,ydat,to.plot=FALSE)
Rprof()
summaryRprof(tmp)
proftable(tmp)

#Marche seulement pour variance évoluant linéairement
Rprof(tmp <- tempfile(), line.profiling=TRUE)
far.bc=getFAR.boxcox(pt0,pt1,xp,ydat,to.plot=FALSE)
Rprof()
summaryRprof(tmp)
proftable(tmp)

system.time({far.bc=getFAR.boxcox(pt0,pt1,xp,ydat,to.plot=TRUE)})
system.time({far.or=getFAR.oridea(pt0,pt1,xp,ydat,to.plot=TRUE)})
