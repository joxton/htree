\name{misc}
\alias{hrf_boot}
\alias{partdep}
\alias{predict_viaux}
\title{Internal helper functions}
\description{
  Helper functions
}
\usage{

hrf_boot(object,B=50,R=10)
predict_viaux(trees,x,oob,id,time,yindx,concurrent,nconcurrent,
	historic,nhistoric,delta,delta0,type,quantiles,nperm,rf=1,ncat=-10)
partdep(object,xindx,xlim=NULL,ngrid=100,ntrees=NULL,
	subsample=1,which.class=1,se=FALSE,cond=NULL)

}
\arguments{
\item{object}{An object of class \code{htree}.}
\item{B}{Number of bootstrap samples.}
\item{R}{Ensemble size for each bootstrap sample.}
\item{trees}{Matrix with trees.}
\item{x}{Data matrix.}
\item{oob}{OOB matrix.}
\item{id}{Vector of subject id's.}
\item{time}{Vector of observation times.}
\item{yindx}{Column index of response variable in \code{x}.}
\item{concurrent}{Concurrent predictors.}
\item{nconcurrent}{Number of concurrent predictors.}
\item{historic}{Historic predictors.}
\item{nhistoric}{Number of historic predictors.}
\item{delta}{Delta values.}
\item{delta0}{Delta0 values.}
\item{type}{Type of historical summary function.}
\item{quantiles}{Quantiles to be used in summary function, if applicable.}
\item{nperm}{Number of permutations for variable importance.}
\item{rf}{Indicator for random forest model or not.}
\item{ncat}{Number of categories.}
\item{xindx}{Column index for partial dependence variable.}
\item{xlim}{Range of values in plot.}
\item{ngrid}{Number of grid points.}
\item{ntrees}{Number of trees.}
\item{subsample}{Sub-sample fraction.}
\item{which.class}{Which class to plot partial dependence for.}
\item{se}{Logical for whether to return standard errors.}
\item{cond}{Logical subsetting.}
}
\details{
These functions are used internally by \code{htree}.

}
