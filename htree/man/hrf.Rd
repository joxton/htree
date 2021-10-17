\name{hrf}
\alias{hrf}
\title{Random forest for longitudinal data}
\description{Fits a random forest of historical regression trees to longitudinal data.}
\usage{



hrf(x,
    time=NULL,
    id=NULL,
    yindx,
    ntrees = 100,
    method="freq",
    mtry=NULL,
    se=FALSE,
    B=100,
    R=10,
    nsamp=5,
    historical=TRUE,
    vh=NULL,
    vc=NULL,
    delta=NULL,
    classify=FALSE,
    control=list())
}
\arguments{
\item{x}{A data frame containing response and predictors}
\item{time}{A vector of observation times associated with rows of \code{x}.}
\item{id}{Subject identifier, if \code{NULL} observations are assumed independent}
\item{yindx}{Name of response variable, alt. its column number in \code{x}}
\item{ntrees}{Number of trees in ensemble}
\item{method}{Historical summary method, can be \code{freq}, \code{frac}, \code{mean0}, \code{freqw}, \code{fracw} and \code{mean0w} (see below).}
\item{mtry}{Number of predictors sampled at each split}
\item{se}{If \code{TRUE} then bootstrap standard errors are computed. Total number of trees fit for bootstrapping is \code{B*R}. }
\item{B}{Only used if \code{se=TRUE}, number of bootstrap samples, defaults to \code{100}.}
\item{R}{Only used if \code{se=TRUE}, forest size for each bootstrap sample, defaults to \code{R=10}.}
\item{nsamp}{Number of sampled \code{delta} values, see below}
\item{historical}{If \code{TRUE} then historical splitting is done, else only standard (ie concurrent predictor) splitting.}
\item{vh}{Optional vector of variable names to be used as historical predictors, can be a numeric vector giving column numbers of \code{x}.}
\item{vc}{Optional vector of variable names to be used as concurrent predictors, can be a numeric vector giving column numbers of \code{x}.}
\item{delta}{An optional vector of time-lags to be used (see below).}
\item{classify}{If \code{TRUE} then a classification tree is built (using gini-impurity node splitting).}
\item{control}{A list of control parameters (see below). All arguments, except those describing the data, can be set in \code{control}. Arguments in \code{control} 
are used if both are given.}

}

\details{

The \code{hrf} function fits a random forest model to longitudinal data. Data is assumed to be of form: 
\eqn{z_{ij}=(y_{ij},t_{ij},x_{ij})} for \eqn{i=1,..,n} and \eqn{j=1,..,n_i}, with \eqn{y_{ij}} being the response for the \eqn{i}-th subject 
at the \eqn{j}-th observation time \eqn{t_{ij}}. The vector of predictors at time \eqn{t_{ij}} are \eqn{x_{ij}}. The number of observations can vary across 
subjects, and the sampling in time can be irregular.


\code{hrf} estimates a model for the response \eqn{y_{ij}} using both \eqn{(t_{ij},x_{ij})} (the observations concurrent with \eqn{y_{ij}}) and all preceeding observations of the \eqn{i}-th subject up to (but not including) time \eqn{t_{ij}}. The model is fit using \code{historical} regression (alt. classification) trees. Here a predictor is one 
of two types, either \code{concurrent} or \code{historical}. The concurrent predictors for \eqn{y_{ij}} are the elements of the vector (\eqn{(t_{ij},x_{ij})}), while a 
historic predictor is the set of all preceeding values (i.e. prior to time \eqn{t_{ij}}) of a given element of \eqn{(y_{ij},t_{ij},x_{ij})}, for subject \eqn{i}. In 
a historic regression tree, node splitting on a \code{concurrent} predictor follows the approach in standard regression (classification) trees. For \code{historical} predictors 
the splitting is modified since, associated with each observed response \eqn{y_{ij}}, the number (and observation times) of observations of a \code{historical} predictor will vary according to 
\eqn{i} and \eqn{j}. For these, the splitting is done by first transforming the preceeding values of a predictor using a \code{summary function}. This summary is 
invertible, in the sense that knowledge of it is equivalent to knowledge of the covariate history. Letting \eqn{\bar{z}_{ijk}} denote the set of historical values of the 
\eqn{k}-th element of \eqn{z_{ij}}, the summary function is denoted \eqn{s(\eta;\bar{z}_{ijk})} where \eqn{\eta} is the argument vector of the summary 
function. Node splitting based on a historical predictor is done by solving \deqn{\mbox{argmin}\sum_{(ij)\in Node} (y_{ij}-\mu_L I(s(\eta;\bar{z}_{ijk})<c)-\mu_R I(s(\eta;\bar{z}_{ijk})\geq c))^2,} where the minimization is over the vector \eqn{(k,\mu,c,\eta)}. Each node of \code{historical} regression tree is split using the best split 
among all splits of concurrent and historical predictors. 




Different \code{summary functions} are available in \code{hrf}, specified by the argument \code{method}. If \code{method="freq"} the summary function summarizing a covariate history is: 
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}} I(z_{ihk}<\eta_2);}  \code{method="frac"}:  
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}} I(z_{ihk}<\eta_2)/n_{ij}(\eta);}  \code{method="mean0"}: 
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}} z_{ihk}/n_{ij}(\eta);}  \code{method="freqw"}: 
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}-\eta_2} I(z_{ihk}<\eta_3);}  \code{method="fracw"}: 
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}-\eta_2} I(z_{ihk}<\eta_3)/n_{ij}(\eta);}  \code{method="meanw0"}:  
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}-\eta_2} z_{ihk}/n_{ij}(\eta).} Here \eqn{n_{ij}(\eta)} denotes the number of observations of subject 
\eqn{i} in the 
time window \eqn{[t_{ij}-\eta_1,t_{ij}-\eta_2)}. In the case  \eqn{n_{ij}(\eta)=0}, the summary function is set to zero, i.e \eqn{s(\eta;\bar{z}_{ijk})=0}.   The default is \code{method="freq"}. The possible values of \eqn{\eta_1} in the summary function can be set by the argument \code{delta}. If not supplied, the set of possible values of 
\eqn{\eta_1} is 
determined by the difference in time between within-subject successive observations in the data. When a split is attempted on a historical predictor, a sample of this set is 
taken from which the best split is selected. The size of this set equals that of the \code{nsamp} argument. See below on \code{control} for futher arguments governing the 
historical splitting. See below on \code{control} for futher arguments governing the 
historical splitting.



Setting \code{se=TRUE} performs standard error estimation. The number of bootstrap samples (sampling subjects with replacement) is determined by 
\code{B}. For each bootstrap sample a random forest with \code{R} trees is built, which defaults to \code{R=10}. The bias induced by 
using   smaller bootstrap ensemble sizes is corrected for in the estimate. Using \code{se=TRUE} will influence summaries from 
the fitted model, such as providing approximate confidence intervals for partial dependence plots (when running \code{partdep_hrf}), and give 
standard errors for predictions when \code{predict_hrf} is used.  



All arguments (except those decribing the data \code{x}, \code{yindx},\code{time} and \code{id}) can be set in the \code{control} list. The arguments supplied 
in \code{control} are used if both are supplied. So if \code{ntrees=300} and \code{control=list(ntrees=500)} then \code{500} trees are fit. Besides the arguments 
described above, a number of other parameters can be set in control. These are: \code{nodesize} giving the minimum number of training observations in a terminal node; 
\code{sample_fraction} giving the fraction of data sample to train each tree; \code{dtry} the number of sampled \code{delta} values used when splitting based on 
a historical variable (note this is alternatively controlled by the argument \code{nsamp} above ; \code{ndelta} the number of \code{delta} values to use if \code{delta} is not supplied, these are taken as the quantiles from the 
distribution of observed \code{delta} values in the data; \code{qtry} the number of sampled values of \eqn{\eta_2} for \code{method=freq/frac}, \eqn{\eta_3} for 
\code{method=freqw/fracw}; \code{quantiles} is a vector of probabilities, and  is used when \code{method=frac} or \code{method=fracw}, ie when covariate histories are split on their windowed empirical distributions. Splitting is restricted to these quantiles. 


 
}

\value{ Returns a \code{list} with elements, the most important being: \code{h} is a list returned by \code{parLapply} for fitting trees in parallel, \code{error} gives the OOB error (mse for regression, misclassification rate for classification), \code{control} a list containing control arguments, \code{boot} 
gives bootstrapped model estimates (if \code{se=TRUE}), \code{x,id,time} give the training data. }  

\references{
L. Breiman (2001). \dQuote{Random Forests,} \emph{Machine Learning} 45(1):5-32.

Zhang and Singer (2010) \dQuote{Recursive Partioning and Applications} \emph{Springer}.

Sexton and Laake (2009) \dQuote{Standard errors for bagged and random forest estimators,} \emph{Computational Statistics and Data Analysis}.
}

\author{Joe Sexton \email{joesexton0@gmail.com}}

\seealso{ \code{\link{predict_hrf}}, \code{\link{partdep_hrf}},
          \code{\link{varimp_hrf}}. }

\examples{

\dontrun{

# ----------------------------------------------------------------------------------------- ##
# Mother's stress on child illness:
# 	Investigate whether mother's stress is (Granger) causal for child illness 
#	'hrf' model is fit using previous observations of mother's stress to predict 
#	child's illness at given time point, but not mother's stress at that time point
#
#	Predictor variables are classified into "historical" and "concurrent"
#
#	A predictor is "historical" if its prior realizations can be used to predict 
#	the outcome.   
#
#	A predictor is "concurrent" if its realization at the same timepoint as the outcome
#	can be used to predict the outcome at that timepoint
#
#	A predictor can be both "concurrent" and "historical", the status of the predictors 
#	can be set by the 'vh' and 'vc' arguments of 'hrf'. 
#	(if not set these are automatically determined) 
#	 
# ------------------------------------------------------------------------------------------- ##

data(mscm) 
mscm=as.data.frame(na.omit(mscm))


# -- set concurrent and historical predictors 
historical_predictors=match(c("stress","illness"),names(mscm))
concurrent_predictors=which(names(mscm)!="stress")

control=list(vh=historical_predictors,vc=concurrent_predictors)


# -- fit historical random forest
#	(NOTE: response is 0/1 so a regression tree is
# 	 the same as  a classification tree with Gini-index splitting) 
ff=hrf(x=mscm,id=mscm$id,time=mscm$day,yindx="illness",control=control)

# out-of-bag error 
plot(1:length(ff$error),ff$error,type="l",main="OOB error",xlab="forest size",ylab="mse")

# .. larger nodesize works slightly better
control$nodesize=20
ff=hrf(x=mscm,id=mscm$id,time=mscm$day,yindx="illness",control=control)
points(1:length(ff$error),ff$error,type="l",col="blue")


# -- variable importance table 
vi=varimp_hrf(ff)
vi



# -- fit historical random forest with 'se=TRUE'
control$se=TRUE
ff=hrf(x=mscm,id=mscm$id,time=mscm$day,yindx="illness",control=control)

# -- partial dependence for top 4 predictors (with +/-2 SE estimates) 
par(mfrow=c(2,2))
for(k in 1:4)
	pd=partdep_hrf(ff,xindx=as.character(vi$Predictor[k]))

par(mfrow=c(1,1))




## -- Classification trees 

## setting classify=TRUE builds classification tree (gini-impurity node splitting)
control$classify=TRUE
## ... standard error estimation not implemented .. turn off bootstrapping 
control$se=FALSE

ff=hrf(x=mscm,id=mscm$id,time=mscm$day,yindx="illness",control=control)

# -- plot oob classification error 
plot(1:length(ff$error),ff$error,type="l",xlab="forest size",ylab="oob classification error")
abline(mean(mscm$illness),0,lty=2)  ## error of constant model

p=predict_hrf(ff)

## variable importance table (model error measured by gini-impurity)
vi=varimp_hrf(ff)
vi


# -------------------------------- #
# Data w/irregular observation times 
# ------------------------------- #
data(cd4)

control=list(se=TRUE)
ff=hrf(x=cd4,id=cd4$id,time=cd4$time,yindx="count",control=control)

vi=varimp_hrf(ff)

# -- partial dependence for top 4 predictors (with +/-2 SE estimates) 
par(mfrow=c(2,2))
for(k in 1:4)
	pd=partdep_hrf(ff,xindx=as.character(vi$Predictor[k]))
par(mfrow=c(1,1))

plot(1:length(ff$error),ff$error,xlab="forest size",ylab="oob mse",type="l")

## by default, the number of delta values (parameter 'eta_1' above) is 20
## can set this using 'ndelta'  
control$ndelta=50

control$se=FALSE # -- turning off bootstrapping .. 
ff=hrf(x=cd4,id=cd4$id,time=cd4$time,yindx="count",control=control)
points(1:length(ff$error),ff$error,type="l",lty=2)

# the grid of delta values 
ff$control$delta


# --------------------------------------- ##
# Boston Housing data (not longitudinal)
# --------------------------------------- ##
# library(htree)
library(mlbench)
library(randomForest)

data(BostonHousing)
dat=as.data.frame(na.omit(BostonHousing))


## omitting arguments time/id assumes rows are iid 
control=list(ntrees=500,sample_fraction=.5,nodesize=1)
h=hrf(x=dat,yindx="medv",control=control)

## randomForest comparison 
##	(by default, randomForest samples with replacement, while hrf samples without) 
r=randomForest(medv~.,data=dat,replace=F,sampsize=ceiling(.5*nrow(dat)),nodesize=1)

## plot oob-error for both
plot(1:length(r$mse),r$mse,type="l",ylim=c(min(r$mse,h$error),max(r$mse,h$error)),
	main="BostonHousing",xlab="forest size",ylab="out-of-bag mse")
points(1:length(h$error),h$error,type="l",col="blue")

## -- variable importance table
vi=varimp_hrf(h)
vi

## -- partial dependence plots with approximate 95% C.I 
control$se=TRUE
h=hrf(x=dat,yindx="medv",control=control)

par(mfrow=c(2,2))
for(k in 1:4)
	pd=partdep_hrf(h,xindx=as.character(vi$Predictor[k]))

par(mfrow=c(1,1))



}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{longitudinal}


