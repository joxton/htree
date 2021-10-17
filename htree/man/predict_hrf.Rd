\name{predict_hrf}
\alias{predict_htb}
\alias{predict_hrf}
\title{Prediction}
\description{
  Prediction functions for \code{hrf} and \code{htb}.
}
\usage{

predict_htb(object,x=NULL,time=NULL,id=NULL,
	all.trees=FALSE,type="response",ntrees=NULL,se=FALSE)

predict_hrf(object,x=NULL,time=NULL,id=NULL,
	all.trees=FALSE,se=FALSE)
}
\arguments{
  \item{object}{An object of class \code{htree}.}
  \item{x}{A data frame or matrix containing new data.  If \code{NULL} then training data in \code{object} is used. }
   \item{time}{A vector of observation times.}
   \item{id}{A vector of length \code{nrow(x)} identifying the subjects.}
    \item{type}{If \code{response} then predictions on same scale as response are given.}
  \item{ntrees}{The number of trees to use in prediction.}
   \item{all.trees}{If \code{TRUE} then predictions associated with each tree are returned.}
  \item{se}{If \code{TRUE} then standard errors of predictions are returned.}

}

\value{
Returns predictions.
 }


\seealso{\code{\link{hrf}},\code{\link{htb}}}

\examples{

\dontrun{

# ------------------------------- #
# Simulated data example          #
# ------------------------------- #
# library(htree)
p=5;sigma_e=.5;sigma_a=.5;n=500;pnoise=2
random_intercept=as.numeric(mapply(rep,rnorm(n,sd=sigma_a),times=p))
dat=data.frame(time=rep(1:p,n),
	x=(random_intercept+rnorm(n*p,sd=sigma_e)),
	znoise=matrix(rnorm(n*p*pnoise),ncol=pnoise))
id=sort(rep(1:n,p))

# fit historical random forest
hb=hrf(x=dat,time=dat$time,id=id,yindx=2,se=TRUE)

# get predictions with standard errors
pred=predict_hrf(hb,se=TRUE)


# ------------------------------------------------------------------ #
# Comparison of SE-estimates with actual standard errors for 'hrf'
# ------------------------------------------------------------------ #  


## -- evaluation points
n=200
datp=data.frame(y=rep(0,n),w=seq(-2,2,length=n),z=rep(0,n))

## -- estimate model on 50 simulated data sets
pred=NULL
pred_se=NULL
nsim=20

## -- B=100 bootstrap samples, ensemble size of R=10 on each
control=list(ntrees=500,B=100,R=10,se=TRUE,nodesize=5)
for(k in 1:nsim){
	
	if(is.element(k,seq(1,nsim,by=10)))
		cat(paste("simulation: ",k," of ",nsim," \n",sep=""))
	# -- simulation model -- #
	dat=data.frame(y=(4*datp$w+rnorm(n)),x=datp$w,z=rnorm(n))
	# ---------------------- #
	h=hrf(x=dat,yindx="y",control=control)
	mm=predict_hrf(object=h,x=datp,se=TRUE)
	pred=cbind(pred,mm[,1])
	pred_se=cbind(pred_se,mm[,2])
}


# --- Actual Standard errors at datp 
pred_se_true=apply(pred,1,sd)

# --- Mean of estimated standard errors
pred_se_est=apply(pred_se,1,mean)
pred_se_lower=apply(pred_se,1,quantile,prob=.1)
pred_se_upper=apply(pred_se,1,quantile,prob=.9) 


# -- Plot estimated SE and true SE (+smooth)
z=c(pred_se_true,pred_se_est,pred_se_lower,pred_se_upper)
ylim=c(min(z),max(z))
plot(datp$w,pred_se_est,ylim=ylim,col="blue",xlab="w",
	ylab="Standard error",type="l",main=" SE-true (red) SE-est (blue)")
points(datp$w,pred_se_lower,col="blue",type="l",lty=2)
points(datp$w,pred_se_upper,col="blue",type="l",lty=2)

points(datp$w,pred_se_true,col="red",type="l")

 



# ------------------------------------------------------------------ #
# Comparison of SE-estimates with actual standard errors for 'htb'
# ------------------------------------------------------------------ #  


## -- evaluation points
n=200
datp=data.frame(y=rep(0,n),w=seq(-2,2,length=n),z=rep(0,n))

## -- estimate model on 50 simulated data sets
pred=NULL
pred_se=NULL
nsim=20
for(k in 1:nsim){
	
	if(is.element(k,seq(1,nsim,by=10)))
		cat(paste("simulation: ",k," of ",nsim," \n",sep=""))
	# -- simulation model -- #
	dat=data.frame(y=(4*datp$w+rnorm(n)),x=datp$w,z=rnorm(n))
	# ---------------------- #
	h=htb(x=dat,yindx="y",ntrees=200,cv.fold=10)
	mm=predict_htb(object=h,x=datp,se=TRUE)
	pred=cbind(pred,mm[,1])
	pred_se=cbind(pred_se,mm[,2])
}


# --- Actual Standard errors at datp 
pred_se_true=apply(pred,1,sd)

# --- Mean of estimated standard errors
pred_se_est=apply(pred_se,1,mean)
pred_se_lower=apply(pred_se,1,quantile,prob=.1)
pred_se_upper=apply(pred_se,1,quantile,prob=.9) 


# -- Plot estimated SE and true SE (+smooth)
z=c(pred_se_true,pred_se_est,pred_se_lower,pred_se_upper)
ylim=c(min(z),max(z))
plot(datp$w,pred_se_est,ylim=ylim,col="blue",xlab="w",
	ylab="Standard error",type="l",main=" SE-true (red) SE-est (blue)")
points(datp$w,pred_se_lower,col="blue",type="l",lty=2)
points(datp$w,pred_se_upper,col="blue",type="l",lty=2)

points(datp$w,pred_se_true,col="red",type="l")




}
}
\keyword{classif}% at least one, from doc/KEYWORDS
\keyword{regression}
