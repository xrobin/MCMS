# This is Jesper's variance model (globalBayesVar) for MS analysis and related functions

#' Variance model for MS data
#' This function creates a model of the variance based on the global proteomics data.
#' The data is cut in bins of similar |log(ratio)| where the variance is asseased. The binning is controlled by the \code{n.cuts}, \code{counts.pr.cut}, \code{cutVector} and \code{cut.type} parameters
#' @param data the normalized data
#' @param mean the name of the column with the mean ratio
#' @param q the name of the column with q = (n-1) * sd^2
#' @param n the name of the column with the number of effective observations
#' @param n.cuts is the number of cuts along the x=|log(ratio)| axis to estimate conditional shape and rate parameters for the inverse variances
#' @param plot whether to show diagnostic plots
#' @param fixed.shape,fixed.rate whether to keep shape and rate fixed
#' @param min.n use only observations with a number of effective observations greater or equal to that
#' @param alpha proportion of the tails of the quantiles to ignore
#' @param counts.pr.cut the target number of ratios in each bin
#' @param cutVector the direct vector of cut points. Overrides the other binning arguments (\code{n.cuts}, \code{counts.pr.cut} and \code{cut.type})
#' @param cut.type whether to cut by \dQuote{quantile} or cuts of \dQuote{equal} size
#' @param weight.type type of weighting
#' @param n.quantiles.pr.cut number of quantiles for the gamma function
#' @param counts.pr.quantile.pr.cut the target number of observation in each quantile for the gamma function
#' @param alpha.pr.cut proportion of the tails to ignore in the gamma function
#' @param min.variance floor levell for the variance
#' @param probs.pr.cut probabilities for the gamma function
#' @param moment.fit whether to fit the moment. Values of 0, 1 or 2 are accepted
#' @param rate.from.mean whether to use the rate in the mean fit
#' @param weighted.mean.shape whether to use weighting to calculate the mean shape
#' @param max.rate maximum allowed value for the rate
#' @author Jesper Ferkinghoff-Borg
#' @export
variance.model <- function(data, mean="mean", q = "q", n = "n", n.cuts=NULL, plot = FALSE, fixed.shape=FALSE, fixed.rate=FALSE,
					  min.n=5, alpha=c(0,0),
					  #n.quantiles=NULL,
					  counts.pr.cut=400, cutVector=NULL, cut.type=c("quantile", "equal"),
					  weight.type=c("normal", "none", "log"),
					  n.quantiles.pr.cut=NULL, counts.pr.quantile.pr.cut=20, alpha.pr.cut=0, min.variance=0.005,
					  probs.pr.cut=NULL, moment.fit=0, rate.from.mean=FALSE, weighted.mean.shape=TRUE, max.rate=25, ...) {
	data=data[data[[n]]>=min.n & !is.na(data[[mean]]),]
	data$variance=data[[q]]/(data[[n]]-1)
	data=data[data$variance>min.variance,]
	filt=data.frame(mean=abs(data[[mean]]),variance=data$variance)
	filt$n=data[[n]]
	init.multiple.plots(1)
	plot(filt$mean,filt$variance,cex=0.3,...)

	cut.type <- match.arg(cut.type)
	weight.type <- match.arg(weight.type)

	if (length(alpha)==1) alpha=rep(alpha,2)
	q.a=quantile(filt$mean,probs=c(alpha[1],1-alpha[2]))
	if (!is.null(cutVector))
	{
		q.a[1]=max(q.a[1],cutVector[1])
		q.a[2]=min(q.a[2],cutVector[length(cutVector)])
	}
	filt=filt[filt$mean>=q.a[1] & filt$mean<=q.a[2],]

	if (is.null(cutVector))
	{
		if (is.null(n.cuts)) n.cuts=max(2,floor(length(filt$variance)/counts.pr.cut))

		if (cut.type=="quantile")
		{
			print("cutting by quantiles")
			cut.probs=seq(n.cuts)/n.cuts
			cutVector=c(0,quantile(filt$mean,probs=cut.probs))
			# Merge duplicate cuts (eg if many 0 ratios)
			cutVector <- cutVector[!duplicated(cutVector)]
		} else {
			print("cutting equal")
			#q.c=quantile(filt$mean,probs=(n.cuts-1)/n.cuts)
			#cutVector=c(seq(0,q.c,q.c/(n.cuts-1)),max(filt$mean))
			m=max(filt$mean)
			cutVector=seq(0,m,m/n.cuts)
		}
	}	else {n.cuts=length(cutVector)-1}
	#if (is.null(nquantiles)) nquantiles=(nrow(filt$n)/(20*ncuts)) #20 counts pr. quantile of the variance-distribution
	#quantiles=seq(nquantiles)/nquantiles
	filt$groups=cut(filt$mean,cutVector, include.lowest = TRUE)
	tg<-table(filt$groups)
	if (plot) print(tg)
	if (min(tg)<0.8*counts.pr.cut)
	{
		if (!plot) print(tg)
		warning(paste("Some of the bins have counts less than ",0.8*counts.pr.cut," -> make cut strategy more coarse-grained"))
	}

	#fit=nls(y~pgamma(x,shape,1./scale),start=list(scale=0.1,shape=1.1),data=his,lower=c(0.01,0.5),algorithm="port")
	if (plot) init.multiple.plots(n.cuts)

	result <- (filt %>% group_by(groups) %>%
			   	do(
			   		gamma.fit.dplyr(., plot=plot, as.dataframe=TRUE,alpha=alpha.pr.cut,n.quantiles=n.quantiles.pr.cut,
			   						counts.pr.quantile=counts.pr.quantile.pr.cut,
			   						probs=probs.pr.cut,moment.fit=moment.fit,min.val=min.variance,...)
			   	))
	direct.rate=result$rate
	direct.shape=result$shape

	if (weight.type=="normal")
	{
		result$shape.w=1/(result$shape.std)^2
		result$rate.w=1/(result$rate.std)^2
	}
	if (weight.type=="none")
	{
		result$shape.w=rep(1,length(result$shape.std))
		result$rate.w=rep(1,length(result$rate.std))
	}
	if (weight.type=="log")
	{
		result$shape.w=result$shape^2/(result$shape.std)^2
		result$rate.w=result$rate^2/(result$rate.std)^2
	}

	direct.rate.w=result$rate.w
	direct.shape.w=result$shape.w

	#weighted shape average
	mean=list()
	if (weighted.mean.shape) mean$shape=sum(result$shape.w*result$shape)/sum(result$shape.w) else 	mean$shape=mean(result$shape)
	mean$rate=sum(result$rate.w*result$rate)/sum(result$rate.w)

	if (is.numeric(fixed.shape)) f.shape=fixed.shape
	else if (fixed.shape) f.shape=mean$shape else f.shape=NULL
	if (is.numeric(fixed.rate)) f.rate=fixed.rate
	else if (fixed.rate) f.rate=mean$rate else f.rate=NULL

	if (plot) init.multiple.plots(n.cuts)
	if (fixed.rate || fixed.shape)
	{
		print("redo fit with fixed shape")
		result <- (filt %>% group_by(groups) %>%
				   	do(
				   		gamma.fit.dplyr(.,plot=plot, as.dataframe=TRUE,alpha=alpha.pr.cut,n.quantiles=n.quantiles.pr.cut,
				   						counts.pr.quantile=counts.pr.quantile.pr.cut,
				   						probs=probs.pr.cut,moment.fit=moment.fit,fixed.shape=f.shape,fixed.rate=f.rate,min.val=min.variance,...)

				   	))
	}
	result$direct.rate=direct.rate
	result$direct.rate.w=direct.rate.w
	result$direct.shape=direct.shape
	result$direct.shape.w=direct.shape.w
	if (fixed.shape || weight.type=="none") result$shape.w=rep(1,length(result$shape)) else
	{
		if (weight.type=="normal") result$shape.w=1/(result$shape.std)^2
		if (weight.type=="log") result$shape.w=1/(result$shape.std)^2*result$shape^2
	}
	if (fixed.rate || weight.type=="none") result$rate.w=rep(1,length(result$rate)) else
	{
		if (weight.type=="normal") result$rate.w=1/(result$rate.std)^2
		if (weight.type=="log") result$rate.w=result$rate^2/(result$rate.std)^2
	}
	result2<- (filt %>% group_by(groups) %>%
			   	dplyr::summarise(
			   		#mean.x = mean(mean),
			   		#mean.var=mean(variance),
			   		mean.x=median(mean),
			   		mean.var=mean(variance),
			   		mean.n=mean(n),
			   		#var.var=0.5*mean(variance)^2/(mean(n)-1),
			   		var.var=0.5*median(variance)^2/(mean(n)-1),
			   		n = n()))

	result<-left_join(result,result2)
	curve=list()
	if (rate.from.mean)
		curve$rate=data.frame(x=result$mean.x,y=(result$shape-1)*result$mean.var, w=1/((result$shape-1)^2*result$var.var+result$mean.var^2*(result$shape.std^2)))
	else curve$rate=data.frame(x=result$mean.x,y=result$rate,w=1/(result$rate.std^2+(fixed.rate)))
	curve$shape=data.frame(x=result$mean.x,y=result$shape,w=1/(result$shape.std^2+(fixed.shape)))

	m.offset=0
	#curve$mode=data.frame(x=result$mean.x,y=result$rate/(result$shape+m.offset),
	#					  w=1/(result$rate.std^2/(result$shape+m.offset)^2+(result$rate/(result$shape+m.offset)^2)^2*result$shape.std^2))

	curve$mode=data.frame(x=result$mean.x,y=result$rate/(result$shape+m.offset),w=1/(result$rate.std^2/(result$shape+m.offset)^2+result$shape.std^2*result$rate^2/(result$shape+m.offset)^4))
	#curve$mean=data.frame(x=result$mean.x,y=log(result$rate/(result$shape)),w=1/(1/curve$rate$w+result$shape.std^2/(result$shape^2)))
	#not really mean, actually 1/<lambda>=rate/shape (as opposed to <1/lambda>=rate/(shape-1))
	offset=1
	mean.var=result$rate/(result$shape-offset)
	#mean.var=result$mean.variance
	ind=(result$mean.x<30 & mean.var>0)
	#curve$mean=data.frame(x=result$mean.x[ind],y=log(mean.var[ind]),w=1/(1/curve$rate$w[ind]+result$shape.std[ind]^2/(result$shape[ind]-offset)^2))
	#curve$mean=data.frame(x=result$mean.x[ind],y=log(result$mean.var[ind]),w=1/(result$var.var/result$mean.var[ind]^2))
	#curve$mean=data.frame(x=result$mean.x[ind],y=log(mean.var[ind]),w=1/(1/curve$rate$w[ind]+result$shape.std[ind]^2/(result$shape[ind]-offset)^2))
	curve$mean=data.frame(x=result$mean.x[ind],y=mean.var[ind],w=1/(result$rate.std[ind]^2/(result$shape[ind]-offset)^2+result$shape.std[ind]^2*result$rate[ind]^2/(result$shape[ind]-offset)^4))

	#do linear fit first, in all cases
	fit=list()
	mdl=list()
	curve.types=c("rate","shape","mode")

	if (plot) init.multiple.plots(1)

	for (ct in curve.types)
	{
		if (plot) plot(curve[[ct]]$x,curve[[ct]]$y,main=ct)
		if (nrow(curve[[ct]])>0)
		{
			print(paste("doing curvetype",ct))
			if (ct=="mean") {
				#fit[[ct]]<-nls(y~a+(x/b)^(1+nu*(b/(x+b))),data =curve[[ct]], weights = curve[[ct]]$w,start=list(nu=0.5,b=0.5,a=-0.2))
				#fit[[ct]]<-nls(y~a-log(exp(b)+1)+log(exp(b)+exp(c*x)),data =curve[[ct]], weights = curve[[ct]]$w,start=list(a=-0.5,b=1,c=1))
				attempt<-try(fit[[ct]]<-nls(y~a+b*x^nu,data =curve[[ct]], weights = curve[[ct]]$w,start=list(nu=1,a=1,b=1)))
			}	else if (ct=="shape")
			{
				#l=length(curve[[ct]]$w)
				#m=max(curve[[ct]]$w)
				#curve[[ct]]$w[1]=5*m
				i.max=which.max(result$shape)
				b=result$mean.x[i.max]
				#curve[[ct]]$w[i.max]=m
				#curve[[ct]]$w[l]=0.1*min(curve[[ct]]$w)
				#fit[[ct]]<-nls(y~A*(1+((x-b)/g)^2)^(-nu)+2,data=curve[[ct]],start=list(nu=0.7,b=b,g=0.2,A=25),lower=c(0,0,0,0),
				#			   weights=curve[[ct]]$w,algorithm="port")
				#fit[[ct]]<-nls(y~A*(1+(x/g)^2)^(-nu)+1,data=curve[[ct]],start=list(nu=0.7,g=0.2,A=25),lower=c(0,0,0),
				#			   weights=curve[[ct]]$w,algorithm="port")
				attempt<-try(fit[[ct]]<-nls(y~stretched.gamma(x,a,s,b=1,n,k),data=curve[[ct]],start=list(s=2,k=0.9,a=0.7,n=0.1),
							   lower=rep(0.01,5),algorithm="port")) #weights=curve[[ct]]$w,
				#fit[[ct]]<-lm(y~x,data =curve[[ct]],weights=curve[[ct]]$w)
			}
			else if (ct=="mode"){
				attempt<-try(fit[[ct]] <- lm(y ~ x, data =curve[[ct]], weights = curve[[ct]]$w))
				#fit[[ct]]<-nls(y~a*(1+g*x)^nu,data =curve[[ct]], weights = curve[[ct]]$w,start=list(g=1,nu=1,a=1),lower=c(0,0,0),algorithm="port")
				#fit[[ct]]<-nls(y~a+c*tanh(b*(x-d)),data =curve[[ct]], weights = curve[[ct]]$w,start=list(c=0.5,b=1,d=0,a=-0.2),lower=c(0,0.001,-10,-1),algorithm="port")
				#fit[[ct]]<-nls(y~a+(x/b)^(r-nu*(x/(x+b))),data =curve[[ct]], weights = curve[[ct]]$w,start=list(r=1,nu=0.5,b=0.5,a=-0.2))
				#fit[[ct]]<-nls(y~a+(b*x)^nu/((b*x)^nu+1),data =curve[[ct]], weights = curve[[ct]]$w,start=list(nu=1,a=-0.2,b=1/0.5),lower=c(0,-5,0.1),algorithm="port")
				#fit[[ct]]<-nls(y~a+(x/b)^(c+exp(-x/b)/(1+c)),data =curve[[ct]], weights = curve[[ct]]$w,start=list(b=1,a=-0.2,c=0.5),lower=c(0.1,-3,0),algorithm="port")
				#fit[[ct]]<-nls(y~a+(x/b)^(1-nu*(x/(x+b))),data =curve[[ct]], weights = curve[[ct]]$w,start=list(nu=0.5,b=0.5,a=-0.2))
			}else if (ct=="rate")
			{
				attempt<-try(fit[[ct]]<-nls(y~rate0+(x/rate1)^rate2,data =curve[[ct]], weights = curve[[ct]]$w,
							   start=list(rate2=0.5,rate0=curve[[ct]]$y[1],rate1=3),lower=c(0.1,0,0.01),algorithm="port"))
				#fit[[ct]]<-nls(y~k+0.5*(1+tanh(s*(x/x0-1)))*(x/x0)^b,data=curve[[ct]],weights=curve[[ct]]$w,start=list(s=3,k=0.01,b=1,x0=0.7),algorithm="port",
				#		 lower=rep(0.001,5))
			}
			if (is(attempt,"try-error")) {
				warning(paste("Could not perform model fit for",ct,". Reverting to constant fit..."))
				fit[[ct]]<-lm(y~1,data=curve[[ct]],weights=curve[[ct]]$w)
			}
			curve[[ct]]$predict=predict(fit[[ct]])
			mdl[[ct]]=summary(fit[[ct]])$coefficients
		}
	}
	print("prediction ok ")
	#if (fit.type=="stretched")
	#{
	#	a0=b_fit$coefficients["(Intercept)"]
	#	a1=b_fit$coefficients["x"]
	#	b_fit=nls(y~a0+a1*x^nu, data=fit,weights=result$w,start=list(nu=1,a0=a0,a1=a1))
	#}

	ret=list(data=result,mdl=mdl,
			 get.direct.params = function(x) {
			 	data.frame(shape=predict(fit$shape,data.frame(x=abs(x))),rate=predict(fit$rate, data.frame(x=abs(x))))},
			 #get.mparams=function(x) {
			 #	data.frame(shape=1/(exp(predict(fit$mean,data.frame(x=abs(x)))-predict(fit$mode,data.frame(x=abs(x))))-1),
			 #			   rate=1/(exp(-predict(fit$mode,data.frame(x=abs(x))))-exp(-predict(fit$mean,data.frame(x=abs(x))))))}
			 get.params=function(x)
			 {
			 	data.frame(shape=predict(fit$shape,data.frame(x=abs(x))),
			 			   #rate=(predict(fit$shape,data.frame(x=abs(x)))+m.offset)*exp(predict(fit$mode,data.frame(x=abs(x)))))
			 			   rate=pmin(predict(fit$rate, data.frame(x=abs(x))),max.rate))
			 	#rate=(predict(fit$shape,data.frame(x=abs(x)))-offset)*predict(fit$mean,data.frame(x=abs(x))))
			 	#rate=regular.at.zero(predict(fit$shape,data.frame(x=abs(x)))-offset,0.1,0.1)*exp(predict(fit$mean,data.frame(x=abs(x)))))
			 	#1/(exp(predict(fit$mean,data.frame(x=abs(x)))-predict(fit$mode,data.frame(x=abs(x))))-1),
			 	#			   rate=1/(exp(-predict(fit$mode,data.frame(x=abs(x))))-exp(-predict(fit$mean,data.frame(x=abs(x))))))}
			 }
	)

	print("return ok")
	if (plot) {

		init.multiple.plots(n.cuts)
		for (i in seq(1,nrow(result)))
		{

			m=result$mean.x[i]
			param=ret$get.params(m)
			g=result$groups[i]
			print(paste("doing group",g," with m=",m,"shape=",param$shape,"rate=",param$rate))
			l=filt[filt$groups==g,]
			dummy=gamma.fit.dplyr(l,
								  plot=TRUE, as.dataframe=TRUE,alpha=alpha.pr.cut,n.quantiles=n.quantiles.pr.cut,
								  counts.pr.quantile=counts.pr.quantile.pr.cut,
								  probs=probs.pr.cut,moment.fit=moment.fit,fixed.shape=param$shape,fixed.rate=param$rate,min.val=min.variance,...)
		}

		#init.multiple.plots(length(curve.types))
		init.multiple.plots(1)
		print("start plotting")
		for (ct in curve.types)
		{
			if (nrow(curve[[ct]])>0) {
				xl=c(0,max(curve[[ct]]$x))
				#if (ct=="mode") xl[2]=4
				errbar(curve[[ct]]$x,curve[[ct]]$y,curve[[ct]]$y+1/sqrt(curve[[ct]]$w),curve[[ct]]$y-1/sqrt(curve[[ct]]$w),xlab="|log(ratio)|",ylab=paste0('log(',ct,')'),xlim=xl)
				lines(curve[[ct]]$x,curve[[ct]]$predict,col="blue",type='b')
				if (ct=="rate")
				{
					m.predict=(ret$get.params(curve[[ct]]$x))[[ct]]
					lines(curve[[ct]]$x,m.predict,col="red",type='b')
				}
				if (ct=="mode")
				{
					params=ret$get.params(curve[[ct]]$x)
					predict=params$rate/(params$shape+m.offset)
					lines(curve[[ct]]$x,predict,col="red")
				}
			}
		}
	}
	return(ret)
}


init.multiple.plots<-function(n.plots,max.n=4)
{
	if (n.plots>max.n) n.plots=max.n
	if (n.plots==1) par(mfrow=c(1,1))
	if (n.plots==2) par(mfrow=c(2,1))
	if (n.plots>2 && n.plots<=4 ) par(mfrow=c(2,2))
	if (n.plots>4 && n.plots<=6) par(mfrow=c(3,2))
	if (n.plots>6 && n.plots<=8) par(mfrow=c(4,2))
	if (n.plots>8) par(mfrow=c(5,2))
}

stretched.gamma<-function(x,a=1,s=2.7795,b=0.07945,n=4.93786,k=0.947836)
{
	return((x/a)^s*exp(-(x/b)^n)+k)
}


gamma.fit.dplyr <- function(x,...) {
	gamma.fit(variances=x$variance,interval=c(min(x$mean),max(x$mean)),...)
}


gamma.fit<-function(variances,n.quantiles=NULL,probs=NULL,counts.pr.quantile=20,plot=FALSE,alpha=0.0,
					fixed.shape=NULL,fixed.rate=NULL,moment.fit=0,as.dataframe=FALSE,min.val=0.005,min.shape=0.25,
					interval=NULL,...)

{

	do.fit=(is.null(fixed.rate) || is.null(fixed.shape))
	#vars <- variances %>% filter(!is.na(.)) %>% filter(variances>0)
	if (length(alpha)==1) alpha=rep(alpha,2)
	vars=variances[!is.na(variances) & variances>min.val]
	lambdas=1./vars
	if (is.null(probs))
	{
		if (!is.null(n.quantiles)) nbins=n.quantiles+1 else
			nbins=max(6, floor(length(lambdas)/counts.pr.quantile+1))
		probs=seq(1,nbins-1)/nbins
	}
	probs=probs[probs>=alpha[1] & probs<=1-alpha[2]]

	qs=quantile(lambdas,probs=probs)
	p.first=probs[1]
	first=qs[1]
	p.last=probs[length(probs)-1]
	last=qs[length(qs)-1]
	data=data.frame(y=probs,x=qs)
	shape=fixed.shape
	rate=fixed.rate
	lambda.restrict=lambdas[lambdas>=first & lambdas<=last]
	if (do.fit)
	{
		if (moment.fit>0)
		{
			m=mean(lambda.restrict)
			n.r=length(lambda.restrict)
			m.std=1/sqrt(n.r)*sd(lambda.restrict)
			v=sd(lambda.restrict)^2
			v.std=sqrt(2*v^4/(n.r-1))
			if (is.null(shape)) {
				#get initial guesses from moments
				shape=m^2/v
				shape.std=sqrt((2*m/v)^2*m.std^2+(m/v)^4*v.std^2)
			} else
			{
				shape.std=0
			}
			if (is.null(rate)) {
				rate=(shape/m)
				rate.std=sqrt(m.std^2/v^2+(m/v^2)^2*v.std^2)
			}
			else {
				shape=rate*m;
				shape.std=rate*m.std;
				rate.std=0
			}
		} else
		{
			q.s=quantile(lambdas,probs=c(0.16,0.5,0.84))
			m=q.s[2]
			if (is.null(rate))
			{
				V=(0.5*(q.s[3]-q.s[1]))^2
				rate=1/(2*V)*(m+sqrt(m^2+4*V))
			}
			if (is.null(shape)) shape=rate*m+1 else rate=(shape-1)/m
		}
	} else {shape.std=0;rate.std=0}

	if (shape<min.shape) shape=min.shape

	shape0=shape
	rate0=rate

	if (moment.fit<2)
	{
		if (do.fit)
		{
			s.list=list()
			lows=NULL
			if (is.null(fixed.shape)) {s.list$shape=max(1E-6, shape+0.001);lows=max(1E-6,min.shape);}
			if (is.null(fixed.rate))  {s.list$rate=max(1E-6,rate);lows=c(lows,1E-6);}
			fit <- try(nls(y ~ pgamma(x, shape, rate), data=data,
						   start=s.list,lower=lows, algorithm="port", control = nls.control(warnOnly = TRUE)))
			if (is(fit, "try-error")) stop(fit)
			s=summary(fit)$parameters
		}
		if (is.null(fixed.shape))
		{
			shape=s["shape","Estimate"]
			shape.std=s["shape","Std. Error"]
		} else shape.std=0
		if (is.null(fixed.rate))
		{
			rate=s["rate","Estimate"]
			rate.std=s["rate","Std. Error"]
		} else rate.std=0
	}

	if (plot)
	{
		s=sprintf("%.2f",shape)
		r=sprintf("%.2f",rate)
		int.s=""
		if (!is.null(interval)) int.s<-paste(",",sprintf("%.3g",interval[1]),"< x <",sprintf("%.3g",interval[2]))
		txt=paste("s=",s,", r=",r,", n=",length(vars),int.s)


		plot(data$x,data$y,main=txt,xlab='1/V',ylab='P(X<1/V)')
		if (do.fit) {
			if (moment.fit<2) fity=predict(fit,data$x) else fity=pgamma(data$x,shape,rate)
			lines(data$x,fity,col="blue")
		}
		lines(data$x,pgamma(data$x,shape0,rate0),col=c("blue","red")[as.numeric(do.fit)+1])
		#qlim=quantile(lambdas,c(min(alpha[1],1-alpha[2]))
		#lambda.restrict=lambdas[lambdas>=first && lambdas<=last]
		br=max(15,length(lambda.restrict)/50)
		#hist(lambda.restrict,freq=FALSE,breaks=data$x,main=txt)
		hist(lambda.restrict,freq=FALSE,breaks="FD",main=txt,xlab='1/V',ylab='P(1/V)')
		#h=hist(lambda.restrict,freq=FALSE,breaks=data$x,main=txt)
		x=seq(min(lambda.restrict),max(lambda.restrict),(max(lambda.restrict)-min(lambda.restrict))/50)
		lines(x,1/(p.last-p.first)*dgamma(x,shape=shape,rate=rate),col="blue")
		#lines(x,dgamma(x,shape=shape,rate=rate),col="blue")
	}
	median=quantile(lambdas,probs=0.5)
	mean=shape/rate


	if (as.dataframe) return(data.frame(shape=shape,shape.std=shape.std,rate=rate,rate.std=rate.std,median=median,mean=mean))
	return(c(shape=shape,shape.std=shape.std,rate=rate,rate.std=rate.std,median=median,mean=mean))
}


