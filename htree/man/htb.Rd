\name{htb}
\alias{htb}
\title{Tree boosting for longitudinal data}
\description{Fits a boosted ensemble of historical regression trees to longitudinal data.}
\usage{



htb(x,
    time=NULL,
    id=NULL,
    yindx,
    ntrees = 100,
    method="freq",
    nsplit=1,
    lambda=.05,
    family="gaussian",
    cv.fold=0,
    cv.rep=NULL,
    nsamp=5,
    historical=TRUE,
    vh=NULL,
    vc=NULL,
    delta=NULL,
    control=list())


}
\arguments{
\item{x}{A data frame containing response and predictors}
\item{time}{A vector of length \code{nrow(x)} of observation times}
\item{id}{A vector of subject identifiers (length equal to \code{nrow(x)}), if \code{NULL} observations are assumed independent}
\item{yindx}{Name of response variable, alt. its column number in \code{x}}
\item{ntrees}{Number of trees in ensemble}
\item{method}{Historical summary method, can be \code{freq}, \code{frac}, \code{mean0}, \code{freqw}, \code{fracw} and \code{mean0w}}
\item{nsplit}{Number of splits in each regression tree.}
\item{lambda}{Shrinkage parameter applied to each tree.}
\item{family}{Either "gaussian" (default) or "bernoulli".}
\item{cv.fold}{Number of cross-validation folds, if \code{cv.fold<=1} no cross-validation is run. }
\item{cv.rep}{Number of times to repeat the cross-validation. If not given set to \code{cv.fold}.}
\item{historical}{If \code{TRUE} then historical splitting is done, else standard splitting.}
\item{nsamp}{Number of sampled \code{delta} values, see below}
\item{vh}{Optional vector of indexes giving the historical predictors.}
\item{vc}{Optional vector of indexes giving strictly concurrent effects.}
\item{delta}{A vector of history lags to be used (see below), defaults to \code{NULL} in which case all possible observed lags are available for splitting.}
\item{control}{A list of control parameters (see below). All arguments, except those describing the data, can be set in \code{control}. Arguments in \code{control} 
are used if both are given.}
}

\details{

The \code{htb} function fits a boosted tree ensemble to longitudinal data. Data are assumed to be of form: 
\eqn{z_{ij}=(y_{ij},t_{ij},x_{ij})} for \eqn{i=1,..,n} and \eqn{j=1,..,n_i}, with \eqn{y_{ij}} being the response for the \eqn{i}-th subject 
at the \eqn{j}-th observation time \eqn{t_{ij}}. The predictors at time \eqn{t_{ij}} are \eqn{x_{ij}}. The number of observations can vary across 
subjects, and the sampling in time can be irregular.


\code{htb} estimates a model for the response \eqn{y_{ij}} using both \eqn{(t_{ij},x_{ij})} (the observations concurrent with \eqn{y_{ij}}) and all preceeding observations of the \eqn{i}-th subject up to (but not including) time \eqn{t_{ij}}. The model is fit using \code{historical} regression (alt. classification) trees. Here a predictor is one 
of two types, either \code{concurrent} or \code{historical}. The concurrent predictors for \eqn{y_{ij}} are the elements of the vector (\eqn{(t_{ij},x_{ij})}), while a 
historic predictor is the set of all preceeding values (preceeding time \eqn{t_{ij}}) of a given element of \eqn{(y_{ij},t_{ij},x_{ij})}, for subject \eqn{i}. In 
a historic regression tree, node splitting on a \code{concurrent} predictor follows the approach in standard regression (classification) trees. For \code{historical} predictors 
the splitting is modified since, associated with each observed response \eqn{y_{ij}}, the set of observations of a \code{historical} predictor will vary according to 
\eqn{i} and \eqn{j}. For these, the splitting is done by first transforming the preceeding values of a predictor using a \code{summary function}. This summary is 
invertible, in the sense that knowledge of it is equivalent to knowledge of the covariate history. Letting \eqn{\bar{z}_{ijk}} denote the set of historical values of the 
\eqn{k}-th element of \eqn{z_{ij}}, the summary function is denoted \eqn{s(\eta;\bar{z}_{ijk})} where \eqn{\eta} is the argument vector of the summary 
function. Node splitting based on a historical predictor is done by solving \deqn{\mbox{argmin}\sum_{(ij)\in Node} (y_{ij}-\mu_L I(s(\eta;\bar{z}_{ijk})<c)-\mu_R I(s(\eta;\bar{z}_{ijk})\geq c))^2,} where the minimization is over the vector \eqn{(k,\mu,c,\eta)}. Each node of \code{historical} regression tree is split using the best split 
among all splits of concurrent and historical predictors. 



Different \code{summary functions} are available in \code{htb}, specified by the argument \code{method}. Setting \code{method="freq"} corresponds the summary 
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}} I(z_{ihk}<\eta_2);}  \code{method="frac"}:  
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}} I(z_{ihk}<\eta_2)/n_{ij}(\eta);}  \code{method="mean0"}: 
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}} z_{ihk}/n_{ij}(\eta);}  \code{method="freqw"}: 
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}-\eta_2} I(z_{ihk}<\eta_3);}  \code{method="fracw"}: 
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}-\eta_2} I(z_{ihk}<\eta_3)/n_{ij}(\eta);}  \code{method="meanw0"}:  
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}-\eta_2} z_{ihk}/n_{ij}(\eta).} Here \eqn{n_{ij}(\eta)} denotes the number of observations of subject 
\eqn{i} in the 
time window \eqn{[t_{ij}-\eta_1,t_{ij}-\eta_2)}. In the case  \eqn{n_{ij}(\eta)=0}, the summary function is set to zero, i.e \eqn{s(\eta;\bar{z}_{ijk})=0}.   The default is \code{method="freqw"}. The possible values of \eqn{\eta_1} in the summary function can be set by the argument \code{delta}. If not supplied, the set of possible values of 
\eqn{\eta_1} is 
determined by the difference in time between successive observations in the data. When a split is attempted on a historical predictor, a sample of this set is 
taken from which the best split is selected. The size of this set equals that of the \code{nsamp} argument.  See below on \code{control} for futher arguments governing the 
historical splitting.



When \code{cv.fold>1} then cross-validation is performed. In subsequent summaries of the model, say the partial dependence plots from \code{partdep_htb}, 
these cross-validation model fits are used to estimate the standard error. This is done using the delete-d jackknife estimator (where deletion refers to 
subjects, instead of rows of the training data). Each cross-validation model fit is performed by removing a random sample of \code{1/cv.fold} of the  subjects. The number of cross-validation runs is determined by \code{cv.rep} which defaults to the value of \code{cv.fold}. 



All arguments (except those decribing the data \code{x}, \code{yindx},\code{time} and \code{id}) can be set in the \code{control} list. The arguments supplied 
in \code{control} are used if both are supplied, so if \code{ntrees=300} and \code{control=list(ntrees=500)} then \code{500} trees are fit. Besides the arguments 
described above, a number of other parameters can be set in control. These are: \code{nodesize} giving the minimum number of training observations in a terminal node; 
\code{sample_fraction} giving the fraction of data sample to train each tree; \code{dtry} the number of sampled \code{delta} values used when splitting a variable 
based on a covariate history (note this is alternatively controlled by the argument \code{nsamp} above ;
 \code{ndelta} the number of \code{delta} values to use if \code{delta} is not supplied, these are taken as the quantiles from the 
distribution of observed \code{delta} values in the data; \code{qtry} the number of sampled values of \eqn{\eta_2} for \code{method=freq/frac}, \eqn{\eta_3} for 
\code{method=freqw/fracw}; \code{quantiles} is a vector of probabilities, and  is used when \code{method=frac} or \code{method=fracw}, ie when covariate histories are split on their windowed empirical distributions. Splitting is restricted to these quantiles. The fold-size for cross-validation can be set by \code{control$cvfold}, and the number 
of repetitions by \code{control$nfold}.  


 
}

\value{ Returns a \code{list} whose more important elements are: \code{trees} giving the tree ensemble, \code{cv} (if \code{cv.fold>1}) the cross-validation model estimates, \code{cv_error} cross-validated error (mse when \code{family="gaussian"} and negative log-likelihood when \code{family="bernoulli"}, \code{control} a list controlling the fit, \code{x,id,time} give the training data.}

\references{

J.H. Friedman (2001). \dQuote{Greedy Function Approximation: A Gradient Boosting Machine,} \emph{Annals of Statistics} 29(5):1189-1232.

J.H. Friedman (2002). \dQuote{Stochastic Gradient Boosting,} \emph{Computational Statistics and Data Analysis} 38(4):367-378.
}

\author{Joe Sexton \email{joesexton0@gmail.com}}
\seealso{ \code{\link{predict_htb}}, \code{\link{partdep_htb}},
          \code{\link{varimp_htb}}. }
\examples{

\dontrun{


# ----------------------------------------------------------------------------------------- ##
# Mother's stress on child illness:
# 	Investigate whether mother's stress is (Granger) causal of child illness 
#	'htb' model is fit using previous observations of mother's stress to predict 
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
historical_predictors=match(c("illness","stress"),names(mscm))
concurrent_predictors=which(names(mscm)!="stress")
control=list(vh=historical_predictors,vc=concurrent_predictors,
	ntrees=200,family="bernoulli",cvfold=10,lambda=.1)

# logistic regression
ff=htb(x=mscm,id=mscm$id,time=mscm$day,yindx="illness",control=control)

# cross-validated negative log-likelihood
plot(1:ff$ntrees,ff$cv_error,type="l",#col="blue",
	xlab="iterations",ylab="cross-validation error")

# -- variable importance table 
vi=varimp_htb(ff)
vi

# -- plot partial dependence (with +/-2 approximate standard errors)
par(mfrow=c(2,2))
for(k in 1:4)
	partdep_htb(ff,xindx=as.character(vi$Predictor[k]))

par(mfrow=c(1,1))


# -- Standard errors are based on cross-validation runs (using delete-d 
#	(subjects) jackknife estimator, d="number-of-subjects"/cvfold)
# -- increasing nfold (which defaults to equal cvfold) gives less 
#	noisy standard error estimates ...    
control$nfold=50  
ff=htb(x=mscm,id=mscm$id,time=mscm$day,yindx="illness",control=control)

par(mfrow=c(2,2))
for(k in 1:4)
	partdep_htb(ff,xindx=as.character(vi$Predictor[k]))

par(mfrow=c(1,1))


# -------------------------------- #
# Data w/irregular observation times 
# ------------------------------- #
data(cd4)

control=list(cvfold=10,lambda=.1,nsplit=3,ntrees=200)
ff=htb(x=cd4,id=cd4$id,time=cd4$time,yindx="count",control=control)

vi=varimp_htb(ff)

# -- partial dependence for top 4 predictors (with +/-2 SE estimates) 
par(mfrow=c(2,2))
for(k in 1:4)
	pd=partdep_htb(ff,xindx=as.character(vi$Predictor[k]))
par(mfrow=c(1,1))


# -- cv error
plot(1:ff$control$ntrees,ff$cv_error,xlab="boosting iteration",
	ylab="cross-validated mse",type="l")
# estimated boosting iteration
abline(v=ff$nboost_cv,lty=2)

## by default, the number of delta values (parameter 'eta_1' above) is 20
## can set this using 'ndelta'  
control$ndelta=50

ff=htb(x=cd4,id=cd4$id,time=cd4$time,yindx="count",control=control)
points(1:ff$control$ntrees,ff$cv_error,type="l",lty=2)  
## note: differences in cv_error can be due (in part) to randomness
##	 in out-of-sample selection

# the grid of delta values 
ff$control$delta


# ------------------------------------------ #
# Boston Housing data (not longitudinal)
# ------------------------------------------ #
# library(htree)
library(mlbench)
data(BostonHousing)
dat=as.data.frame(na.omit(BostonHousing))


# omitting arguments 'time' and 'id' assumes rows are iid 
control=list(cvfold=10,ntrees=500)
h=htb(x=dat,yindx="medv",control=control)


# -- plot cross-validated Mean-squared error --- #
plot(1:h$ntrees,h$cv_error,type="l",xlab="Boosting iterations",
	ylab="MSE",main="Cross-validated Mean-squared error")

# -- variable importance table
vi=varimp_htb(h,nperm=20)
vi 

# -- partial dependence of top 4 predictors with +/-2 S.E.
par(mfrow=c(2,2))
for(k in 1:4)
	partdep_htb(h,xindx=as.character(vi$Predictor[k]))


par(mfrow=c(1,1))

}

}
\keyword{nonparametric}
\keyword{tree}
\keyword{longitudinal}


