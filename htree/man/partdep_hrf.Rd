\name{partdep_hrf}
\alias{partdep_hrf}
\alias{partdep_htb}
\title{Partial dependence}
\description{
  Marginal effects for historical tree ensembles.
}
\usage{
partdep_htb(object,xindx,xlim=NULL,ngrid=10,subsample=.5,
	plot.it=TRUE,cond=NULL)
partdep_hrf(object,xindx,xlim=NULL,ngrid=10,subsample=.5,
	plot.it=TRUE,which.class=1,cond=NULL)



}
\arguments{
  \item{object}{An object of class \code{htree}.}
  \item{xindx}{Name of predictor to compute marginal effect, alternatively the column corresponding this predictor in \code{object$x}.}
  \item{xlim}{Optional range for partial dependence.}
  \item{ngrid}{Number of values in grid for partial dependence.}
   \item{subsample}{Fraction of training data sampled for calculation.}
   \item{plot.it}{If \code{TRUE} then a plot is produced.}
  \item{which.class}{Only used with \code{partdep_hrf} and when \code{hrf} is run with \code{classify=TRUE}. Determines for which class to show the partial dependence.}
   \item{cond}{A logical condition for restricting the partial dependence, see below for details}

}


\value{
Returns \code{x} and \code{y}, the grid and partial dependence values. If \code{htb} run includes cross-validation then 
approximate standard errors 
are also output (based on the leave-out-m jack-knife estimate, where leave-out refers to subjects). 
If \code{hrf} estimation is run with \code{se=TRUE} then approximate bootstrap 
standard errors are returned (resampling subjects with replacement). 
 }
\references{
 J.H. Friedman (2001). \dQuote{Greedy Function Approximation: A Gradient Boosting Machine,} \emph{Annals of Statistics} 29(5):1189-1232.

}

\seealso{\code{\link{hrf}}, \code{\link{htb}}}

\examples{

\dontrun{


# library(htree)
data(mscm)
x=as.data.frame(na.omit(mscm))


# historical predictors 
vh=c("stress","illness")

# concurrent predictors
vc=names(x)[-which(is.element(names(x),vh))]
control=list(vc=vc,vh=vh,nodesize=20)


### -- hrf 
ff=hrf(x=x,yindx="illness",id=x$id,time=x$day,control=control)

## -- baseline illness
pp=partdep_hrf(ff,xindx="billness")


## -- with bootstrap standard errors
control$se=TRUE
ff=hrf(x=x,yindx="illness",id=x$id,time=x$day,control=control)

## -- baseline illness
pp=partdep_hrf(ff,xindx="billness")

## -- mothers stress
pp=partdep_hrf(ff,xindx="stress")


## -- partial dependence for a subset is done using the 'cond' argument

## ... only first week 
pp=partdep_hrf(ff,xindx="billness",cond="day<=7")

## ... last week
pp=partdep_hrf(ff,xindx="billness",cond="day>23",plot.it=FALSE)
points(pp$x,pp$y,type="l",lwd=2,col="blue")


## ... first week and employed mothers
pp=partdep_hrf(ff,xindx="billness",cond="day<=7&emp==1")



### -- hbt -----
# library(htree)
data(mscm)
x=as.data.frame(na.omit(mscm))

# historic predictors
vh=c("stress","illness")
# concurrent predictors
vc=names(x)[-which(is.element(names(x),vh))]

control=list(vc=vc,vh=vh,cvfold=10,family="bernoulli",ntrees=200,lambda=.1)
ff=htb(x=x,yindx="illness",id=x$id,time=x$day,control=control)

vn=c("illness","billness","day","stress")
par(mfrow=c(2,2))
for(k in vn)
	pp=partdep_htb(ff,xindx=k)
par(mfrow=c(1,1))


### --- standard error bands and model estimates 
# library(htree)
sim_data=function(n=100,p=5,pnoise=2,sigma_e=.5,sigma_a=.5){
	## -- random intercept data simulator
	random_intercept=as.numeric(mapply(rep,rnorm(n,sd=sigma_a),times=p))
	dat=data.frame(time=rep(1:p,n),x=(random_intercept+rnorm(n*p,sd=sigma_e)),
		znoise=matrix(rnorm(n*p*pnoise),ncol=pnoise),id=sort(rep(1:n,p)))
	dat
}

## simulate data and estimate model and partial-dependence with standard errors
sdat=sim_data()
control=list(se=TRUE)
h=hrf(x=sdat,yindx="x",id=sdat$id,time=sdat$time,control=control)
pp=partdep_hrf(h,xindx="x",xlim=c(-2,2),ngrid=20)


## estimate and plot partial dependence on simulated data sets 
p_est=NULL
nsim=10
for(s in 1:nsim)
{
	sdat=sim_data()
	hs=hrf(x=sdat,yindx="x",id=sdat$id,time=sdat$time)
	ps=partdep_hrf(hs,xindx="x",plot.it=FALSE,xlim=c(-2,2),ngrid=20)
	p_est=cbind(p_est,ps$y)
	points(ps$x,ps$y,type="l",lty=2,col="blue")
}

## plot of estimated and true standard errors 
res=data.frame(se_est=pp$se,se_obs=apply(p_est,1,sd))
plot(pp$x,res$se_est,ylim=c(min(res),max(res)),type="l",col="blue",
	main="SE-est(blue) SE-obs(red)",xlab="x",ylab="standard error")
points(pp$x,res$se_obs,ylim=c(min(res),max(res)),type="l",col="red")



}

}

