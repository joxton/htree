\name{varimp_hrf}
\alias{varimp_htb}
\alias{varimp_hrf}
\title{Variable importance}
\description{Z-score variable importance for \code{hrf} and \code{htb}}
\usage{
varimp_hrf(object,nperm=20,parallel=TRUE)
varimp_htb(object,nperm=20)

}
\arguments{
\item{object}{Return list from \code{hrf} or \code{htb}}
\item{nperm}{Number of permutations.}
\item{parallel}{If \code{TRUE}, run in parallel.}

}

\details{
To measure the importance of a predictor, \code{varimp_hrf} and \code{varimp_htb} compare the prediction errors of the estimated model with the 
prediction errors obtained after integrating the predictor out of the model. If \eqn{F} denotes the estimated model, the model obtained by 
integrating out predictor \code{k} is \eqn{F_k(x)=\int F(x) dP(x_k)}, where \eqn{P(x_k)} is the marginal distribution of \eqn{x_k}. In 
practice, the integration is done by averaging over multiple predictions from \eqn{F},  each obtained using a random permutation of the observed values of \eqn{x_k}. The number of 
permutations is set by \code{nperm}. Letting \eqn{L(y,y_{hat})}) be the loss of predicting \eqn{y} with \eqn{y_{hat}}, 
the vector \eqn{w_i=L(y_i,F_k(x_i))-L(y_i,F(x_i))} for \eqn{i=1,..,n} gives the difference in the prediction error between the original and marginalized model. 
The corresponding z-score 
\eqn{z=mean(w_i)/se(w_i)} corresponds a paired test for the equality of the prediction errors, in which case it is approximately distributed as N(0,1). Larger z-score values indicate that the prediction error increases if \eqn{x_k} is marginalized out, and thus that \eqn{x_k} is useful. On the other hand, large negative values of the z-score indicate that the integrated model is more accurate.  For longitudinal data, the \code{w_i} are computed by averaging across all observations from the \code{i}-th subject. For \code{htb} the prediction error is calculated based on the cross-validation model estimates, for \code{hrf} out-of-bag predictions 
are used. 
}
\value{A \code{data.frame} with columns: \code{Predictor} giving predictor being marginalized; \code{Marginalized error} gives the prediction error of model with \code{Predictor} marginalized out;  \code{Model error} the prediction error with original model;  \code{Relative change} gives relative change in prediction error due to marginalization; \code{Z-value}: \code{Z} value from test comparing prediction errors of original and marginalized models.   }

\references{
L. Breiman (2001). \dQuote{Random Forests,} \emph{Machine Learning} 45(1):5-32.
}

\seealso{ \code{\link{hrf}}, \code{\link{htb}}}

\examples{
\dontrun{

data(mscm) 
mscm=as.data.frame(na.omit(mscm))


# -- set concurrent and historical predictors 
historical_predictors=match(c("stress","illness"),names(mscm))
concurrent_predictors=which(names(mscm)!="stress")
control=list(vh=historical_predictors,vc=concurrent_predictors,nodesize=20)

## -- fit model
ff=hrf(x=mscm,id=mscm$id,time=mscm$day,yindx="illness",control=control)

# -- variable importance table
vi=varimp_hrf(ff)
vi


## same with htb

control=list(vh=historical_predictors,vc=concurrent_predictors,
	lambda=.1,ntrees=200,nsplit=3,family="bernoulli")
control$cvfold=10 ## need cross-validation runs to run varimp_htb
ff=htb(x=mscm,id=mscm$id,time=mscm$day,yindx="illness",control=control)

# -- variable importance table
vi=varimp_htb(ff)
vi




# --------------------------------------------------------------------------------------------- ##
# Boston Housing data 
#	Comparison of Z-score variable importance with coefficient Z-scores from linear model
# --------------------------------------------------------------------------------------------- ##

# Boston Housing data 
library(mlbench)
data(BostonHousing)
dat=as.data.frame(na.omit(BostonHousing))
dat$chas=as.numeric(dat$chas)

# -- random forest 
h=hrf(x=dat,yindx="medv")


# -- tree boosting
hb=htb(x=dat,yindx="medv",ntrees=1000,cv.fold=10,nsplit=3)


# -- Comparison of variable importance Z-scores and Z-scores from linear model 
vi=varimp_hrf(h)
vb=varimp_htb(hb)
dvi=data.frame(var=as.character(vi$Predictor),Z_hrf=vi$Z)
dvb=data.frame(var=as.character(vb$Predictor),Z_htb=vb$Z)

dlm=summary(lm(medv~.,dat))$coeffi
dlm=data.frame(var=rownames(dlm),Z_lm=round(abs(dlm[,3]),3))
dlm=merge(dlm[-1,],dvi,by="var",all.x=TRUE)

# -- Z-scores of hrf and lm for predictor variables 
merge(dlm,dvb,by="var",all.x=TRUE)



}

}

