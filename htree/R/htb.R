


htb=function(x,time=NULL,id=NULL,yindx,ntrees=100,method="freq",nsplit=1,
		lambda=.05,family="gaussian",cv.fold=0,cv.rep=NULL,nsamp=5,historical=TRUE,vh=NULL,vc=NULL,delta=NULL,control=list())
{


	## arguments can be passed 'directly' or via control-list, control arguments are used if both supplied
	dtry=nsamp
	cvfold=cv.fold
	nfold=cv.rep
	mtry=ncol(x)		

	argnames=c("ntrees","family","method","vh","vc","delta","lambda","family","dtry","cvfold","nfold","historical","nsplit","mtry")
	for(a in argnames)
	{

		if(is.null(control[[a]]))
			eval(parse(text=paste("control[[a]]=",a,sep="")))
		eval(parse(text=paste(a,"=control[[a]]",sep="")))
	}


	if(control$cvfold>0&is.null(control$nfold))
		control$nfold=control$cvfold        ## number of cv runs 


	if(control$family=="bernoulli"&(sum(is.element(x[,yindx],c(0,1)))!=nrow(x)))
		stop(" For 'family=bernoulli' response must be 0-1 vector.")


	if(is.null(control$type))
		control$type=method_convert(method)


	if((is.null(id)|is.null(time)))
	{
		## if time/id then data assumed iid at same time point (ie standard random forest)
		id=c(1:nrow(x))
		time=rep(1,nrow(x))
	}
	

	if(!is.element(class(yindx),c("numeric","character")))
	{
		stop(" 'yindx' must give column number, or name, of response in 'x'.") 
	}else{
		if(class(yindx)=="character")
			yindx=which(names(x)==yindx)

		if(length(yindx)==0)
			stop(" 'yindx' not found ")

	}



	 ## Check missing values
	  if (any(is.na(x))) {
	    stop("Missing data in 'x'.", call. = FALSE)
	  }
	  if (any(is.na(id))) {
	    stop("Missing data in 'id'.", call. = FALSE)
	  }
	  if (any(is.na(time))) {
	    stop("Missing data in 'time'.", call. = FALSE)
	  }


 	 if(yindx>ncol(x)|yindx<1)
		stop(" 'yindx' out of range")

 	if(!is.numeric(time)){
		stop("'time' variable must be numeric. ")
	}  


	vtype=unlist(lapply(x,class))

	if(is.character(id))
		id=as.factor(id)

	if(is.factor(id))
		id=as.numeric(id)

	## -- convert factors/strings if any
	rd=reformat_data(x=x)
	x=rd$x

	ii=order(id,time)
	x=x[ii,]
	id=id[ii]
	time=time[ii]
	flist=rd$flist


	
	if(length(unique(c(nrow(x),length(id),length(time))))!=1) 
		stop(" length of 'id', 'time' and number of rows of 'x' must be the same.")



	ntrees_orig=control$ntrees

	if(is.null(control$nsplit))
		control$nsplit=1  ## stumps 
	pc=proc_control(x=x,yindx=yindx,id=id,time=time,control=control)


	NCOL_TREE=16
	nsplit=nrow(x)


	findx=1
	if(control$family=="bernoulli")
		findx=2

	nodesize=pc$nodesize
	concurrent=pc$concurrent
	nconcurrent=pc$nconcurrent
	historic=pc$historic
	nhistoric=pc$nhistoric
	delta=pc$delta
	delta0=pc$delta0
	quantiles=pc$quantiles
	dtry=pc$dtry
	nsplit=pc$nsplit
	qtry=pc$qtry
	quantiles=pc$quantiles
	type=pc$type
	rf=0
	

	h=cfaux(X=1,x,yindx,id,time,ntrees,mtry,nsplit,
		nodesize=nodesize,concurrent=concurrent,nconcurrent=nconcurrent,historic,nhistoric,
		delta,delta0,type,quantiles,dtry,qtry,rf,lambda,findx)

	hcv=NULL
	cv_error=NULL
	if(control$cvfold>0)
	{

	if(is.null(control$ncores))
	{
		control$ncores=detectCores()
	}
	cluster=makeCluster(control$ncores)
	on.exit(stopCluster(cluster))

	
	hcv=parLapply(cl=cluster, X=1:control$nfold,cfaux,x,yindx,id,time,ntrees,mtry,nsplit,
		nodesize=nodesize,concurrent=concurrent,nconcurrent=nconcurrent,historic,nhistoric,
		delta,delta0,type,quantiles,dtry,qtry,rf,control$lambda,findx,boot=FALSE,cvfold=control$cvfold,ncat=-10)

	
	
	cverror=function(z){
			
		z$yindx=yindx
		z$rf=0
		z$ncat=-10
		x=x[-z$row_train,]
		id=id[-z$row_train]
		time=time[-z$row_train]
		z$control=pc  ### CHANGED FROM CONTROL
		pp=predict_cfn(object=z,x=x,id=id,time=time,all.trees=TRUE)
		y=x[[z$yindx]]

		ac=t(apply(pp,1,cumsum))
		if(findx==1)
			h=lapply(as.data.frame(ac),function(w){ mean((y-w)^2)})
		if(findx==2)
			h=lapply(as.data.frame(ac),function(w){ pr=1/(1+exp(-w));mean(-(y*log(pr)+(1-y)*log(1-pr)))})   # mean((y-pr)^2)

		h=unlist(h)
		h

	}

	m=lapply(hcv,cverror)
	m=matrix(unlist(m),ncol=control$ntrees,byrow=TRUE)
	cve=colMeans(m)
	cv_error=cve

	}
	

	nboost_cv=NULL
	if(!is.null(cv_error))
		nboost_cv=order(cv_error)[1]
	
	h$yindx=yindx
	h$x=deformat_data(x=x,flist=flist,vtype=vtype)
	h$id=id
	h$time=time
	h$x[ii,]=h$x  ## place in original order
	h$time[ii]=h$time  
	h$id[ii]=h$id 

	h$rf=0
	h$cv_error=cv_error
	h$cvfold=cvfold
	h$nfold=nfold
	h$cv=hcv
	h$nboost_cv=nboost_cv
	h$control=pc
	h$dimx=dim(x)
	h$ntrees=ntrees
	h$family=family
	h$findx=h$findx
	h$ncat=-10
	h$flist=flist
	h$vtype=vtype
	class(h)="htree"
	h
}






predict_htb=function(object,x=NULL,time=NULL,id=NULL,all.trees=FALSE,type="response",ntrees=NULL,se=FALSE)
{

	#se=FALSE
	if(is.null(x))
	{
		x=object$x
		id=object$id
		time=object$time

	}else{

		if(is.null(time)|is.null(id))
		{
			id=c(1:nrow(x))
			time=rep(1,nrow(x))
		}
	}

	if(class(object)!="htree")
		stop(" 'object' not of class 'htree' ")

	if(object$rf!=0)
		stop(" 'object' not fit by 'htb' ")

	yindx=object$yindx
	
	ii=NULL 
	if(!is.null(x))
	{
	
		  if (any(is.na(x))) {
		    	stop("Missing data in 'x'.", call. = FALSE)
  			}
		  if (any(is.na(id))) {
 			   stop("Missing data in 'id'.", call. = FALSE)
			  }
		  if (any(is.na(time))) {
 			   stop("Missing data in 'time'.", call. = FALSE)
 			 }


		if(sum(unlist(lapply(x,class))!=object$vtype&(object$vtype=="character"|object$vtype=="factor"))>0)
			stop("Variable class mismatch with training data.")

		if(ncol(x)!=object$dimx[2])
			stop("Number of columns in 'x' differs from training data.")

		if(nrow(x)==0)
			stop("No rows in 'x'.") 
		
	
		if(length(object$flist)>0)
		{

			# map strings/factors into integers 
			x=format_data(x=x,flist=object$flist)
		}
		if(is.null(id)|is.null(time))
			stop(" Arguments 'id' and 'time' cannot be empty.")
		if(is.character(id))
			id=as.factor(id)

		if(is.factor(id))
			id=as.numeric(id)

		if(!is.numeric(time))
			stop(" 'time' must be numeric.")

		ii=order(id,time)
		x=x[ii,,drop=FALSE]
		id=id[ii]
		time=time[ii]

		  if(yindx>ncol(x)|yindx<0)
			stop(" 'yindx' out of range")
		
	
	}else{
		id=time=NULL
	}

		

	if(!is.null(ntrees)) 
		ntrees=ntrees
	
	if(is.null(ntrees)&(!is.null(object$nboost_cv)))
		ntrees=object$nboost_cv

	pred=predict_cfn(object=object,x=x,id=id,time=time,all.trees=all.trees,se=se,ntrees=ntrees)

	if(is.null(ncol(pred)))
		pred=matrix(pred,ncol=1)

	if(type=="response"&object$control$family=="bernoulli")
	{

		if(all.trees)
		{
			cat(" For 'all.trees=TRUE' predictions on logit scale. \n")
		}else{
			pred=1/(1+exp(-pred))

		}

	}


	if(se&(!is.null(object$cv))){
		ha=NULL
		for(k in 1:length(object$cv))
		{
			object$trees=object$cv[[k]]$trees
			hh=predict_htb_aux(object,x=x,time=time,id=id,all.trees=FALSE,ntrees=NULL)	
			ha=cbind(ha,hh)
		}



		#aa=aggregate(ha,by=list(tt),mean)
		pm=ha #(as.matrix(aa[,-1]))
		B=object$control$nfold
		v=apply(pm,1,var)*((B-1)/B)*(object$control$cvfold-1)
		stderr=sqrt(v)

		pred=cbind(pred,stderr)

	}


	retlist=pred
	if(!is.null(ii))
		retlist[ii,]=pred

	if(is.null(x))  # -- put predictions on training data in original order 
		retlist[object$indx_original,]=pred


	
	
	retlist
}








partdep_htb=function(object,xindx,xlim=NULL,ngrid=10,subsample=.5,plot.it=TRUE,cond=NULL)
{
	main="Marginal effect"



	if((subsample>1)|(subsample<=0))
		stop(" Invalid 'subsample' value.")

	if(class(xindx)=="character")
	{
		xindx=which(names(object$x)==xindx)
	
		if(length(xindx)==0)
			stop(" 'xindx' not found.")
	}

	if(xindx>object$dimx[2])
		stop(" Invalid 'xindx' ")

	pd=partdep_aux(object=object,xindx=xindx,xlim=xlim,ngrid=ngrid,ntrees=NULL,subsample=subsample,se=(!is.null(object$cv)),cond=cond)


	if(plot.it)
	{
		nx=names(object$x)[xindx]
		if((!is.element(object$vtype[xindx],c("factor","character")))&(length(unique(object$x[,xindx]))>5))
		{
		if(!is.null(pd$se))
		{
				upper=pd$y+2*pd$se;lower=pd$y-2*pd$se
				ylim=c(min(lower),max(upper))
				plot(pd$x,pd$y,type="l",lwd=3,xlab=nx,ylab="",main=main,ylim=ylim,cex.axis=1.5,cex.lab=2,cex.main=2)
				points(pd$x,upper,type="l",lwd=3,lty=2)
				points(pd$x,lower,type="l",lwd=3,lty=2)
	
		}else{
			plot(pd$x,pd$y,type="l",lwd=3,xlab=nx,ylab="",main=main,cex.axis=1.5,cex.lab=2,cex.main=2)
		}
		}else{

			plot.means(pd$y,categories=pd$x,SEs=pd$se,x.axis.label=nx)
		}
	}
pd
}






partdep_aux=function(object,xindx,xlim=NULL,ngrid=100,ntrees=NULL,subsample=1,se=FALSE,cond=NULL)
{

	# -- subsample subjects
	uid=unique(object$id)
	ns=round(subsample*length(uid))
	sid=sample(uid,size=ns,replace=F)


	uv=unique(object$x[,xindx])
	nuv=length(uv)
	if(nuv<ngrid)
	{
		tt=sort(uv)
		ngrid=nuv
	}else{
		if(is.null(xlim)){
			tt=seq(min(object$x[,xindx]),max(object$x[,xindx]),length=ngrid)
		}else{
			tt=seq(xlim[1],xlim[2],length=ngrid)
		}
	}
	

	ii=which(is.element(object$id,sid))
	nii=length(ii)
	mid=max(object$id)+1
	idaux=object$id[ii]
	id=NULL
	for(k in 1:ngrid)
		id=c(id,idaux+mid*k)

	ii=rep(ii,ngrid)
	x=object$x[ii,]
	time=object$time[ii]



	if(is.null(ntrees)&(!is.null(object$nboost_cv))){
		ntrees=object$nboost_cv
	}else{
		ntrees=object$control$ntrees
	}

	fit=object  

	tt=sort(rep(tt,nii)) 

	x[,xindx]=tt

	vv=NULL
	xx=as.matrix(x)
	yindx=fit$yindx


	hh=predict_htb(object,x=x,time=time,id=id)


	if(!is.null(cond))
	{
		txt=paste("ii=with(x,which(",cond,"))",sep="")
		eval(parse(text=txt))
		
		if(length(ii)>0)
		{
			dd=data.frame(y=hh[ii],x[ii,xindx])
			tt=tt[ii]
		}else{
			cat("Condition 'cond' not satisfied by any observations, ignoring ....")
			dd=data.frame(y=hh,x[,xindx])
		
		}

	}else{
		dd=data.frame(y=hh,x[,xindx])

	}		
	ds=aggregate(dd$y,list(tt),mean)

	if(se==TRUE&length(object$cv)==0)
	{
		cat("Need run 'htb' with 'cvfold>0' to get standard errors.\n")
		se=FALSE

	}

	stderr=NULL
	zscore=NULL
	if(se&!is.null(object$cv)){
		ha=NULL
		for(k in 1:length(object$cv))
		{
			object$trees=object$cv[[k]]$trees
			hh=predict_htb(object,x=x,time=time,id=id,all.trees=FALSE,ntrees=NULL)
			if(!is.null(cond)&(length(ii)>0))
				hh=hh[ii,]
	
			ha=cbind(ha,hh)
		}



		aa=aggregate(ha,by=list(tt),mean)
		pm=(as.matrix(aa[,-1]))
		B=object$control$nfold
		v=apply(pm,1,var)*((B-1)/B)*(object$control$cvfold-1)
		stderr=sqrt(v)



	}

	list(y=ds[,2],x=ds[,1],se=stderr)
}
















varimp_htb=function(object,nperm=20)
{

	if(is.null(object[["cv"]]))
		stop(" No 'cross-validation' runs, cannot compute variable importance...")

	x=object$x
	id=object$id
	time=object$time
	yindx=object$yindx
	control=object$control
	#nti=ceiling(ntrees/control$ncores)
	xorig=x
	x=format_data(x=x,flist=object$flist)
	## -- order 
	ii=order(id,time)
	x=x[ii,]
	id=id[ii]
	time=time[ii]
	xorig=xorig[ii,]	

	pc=proc_control(x=x,yindx=yindx,id=id,time=time,control=control)

	nodesize=pc$nodesize
	concurrent=pc$concurrent
	nconcurrent=pc$nconcurrent
	historic=pc$historic
	nhistoric=pc$nhistoric
	delta=pc$delta
	delta0=pc$delta0
	quantiles=pc$quantiles
	dtry=pc$dtry
	nsplit=pc$nsplit
	qtry=pc$qtry
	quantiles=pc$quantiles
	type=pc$type
	rf=0

	## --- ... 
	
	gg=function(X,object,pc,nperm){
		trees=object$cv[[X]]$trees
		trees=trees[which(trees[,1]<object$nboost_cv),]
		ntrees=object$nboost_cv
		x_train=x[-object$cv[[X]]$row_train,]
		id_train=id[-object$cv[[X]]$row_train]
		time_train=time[-object$cv[[X]]$row_train]
		yindx=object$yindx 
		rf=0
		oob=array(1,dim=c(nrow(x),ntrees))
		
		h_vi=htree::predict_viaux(trees,x_train,oob,id_train,time_train,yindx,
			pc$concurrent,pc$nconcurrent,pc$historic,pc$nhistoric,
			pc$delta,pc$delta0,pc$type,pc$quantiles,nperm,rf)
		h_vi
	
	}

	if(is.null(control$ncores))
	{
		control$ncores=detectCores()
	}
	cluster=makeCluster(control$ncores)
	on.exit(stopCluster(cluster))


	h=parLapply(cluster,X=1:length(object$cv),gg,object=object,pc=pc,nperm=nperm)	

	pp=array(0,dim=dim(object$x))
	n=rep(0,nrow(object$x))
	pr=rep(0,nrow(object$x))

	for(k in 1:length(object$cv))
	{
		ii=c(1:nrow(object$x))[-object$cv[[k]]$row_train]
		pp[ii,]=pp[ii,]+h[[k]]
		n[ii]=n[ii]+1
		
		oc=object
		oc$trees=oc$cv[[k]]$trees
		pr[ii]=pr[ii]+predict_htb(oc,x=xorig[ii,],id=id[ii],time=time[ii],ntrees=object$nboost_cv,type="link")

	}	

	ii=which(n>0)
	y=x[[object$yindx]][ii]
	pp=pp[ii,]/n[ii]
	pr=pr[ii]/n[ii]



	zscore=function(z){

		if(object$family=="gaussian")
		{
			error_perm=mean((y-z)^2)
			error_oob=mean((y-pr)^2)
			d=abs(y-z)-abs(y-pr)
		}

		if(object$family=="bernoulli")
		{
			zr=1/(1+exp(-z))
			prr=1/(1+exp(-pr))
			error_perm=mean((-y*log(zr)-(1-y)*log(1-zr)))
			error_oob=mean((-y*log(prr)-(1-y)*log(1-prr)))
			d=sqrt(-y*log(zr)-(1-y)*log(1-zr))-sqrt(-y*log(prr)-(1-y)*log(1-prr))

		}
		
		#d=(object$x[[object$yindx]]-z)^2-(object$x[[object$yindx]]-pred_oob)^2
		da=aggregate(d,by=list(id=id[ii]),mean)[,2]
		zs=mean(da)/(sd(da)/sqrt(length(da)))
		cret=c(error_oob,error_perm,(error_perm-error_oob)/error_oob,zs)
		return(cret)
	}
	h=matrix(unlist(lapply(as.data.frame(pp),zscore)),ncol=4,byrow=TRUE)
	dat=data.frame(vn=names(x),error_perm=round(h[,2],4),error_oob=round(h[,1],4),rchange=round(h[,3],3),zscore=round(h[,4],3))
	

	names(dat)=c("Predictor","Marginalized error","Model error","Relative change","Z-value") 
	dat=dat[order(dat[,2],decreasing=TRUE),]
	dat
}














predict_htb_aux=function(object,x=NULL,time=NULL,id=NULL,all.trees=FALSE,type="response",ntrees=NULL)
{

	## same as predict_htb, but without se

	se=FALSE
	if(is.null(x))
	{
		x=object$x
		id=object$id
		time=object$time

	}

	if(class(object)!="htree")
		stop(" 'object' not of class 'htree' ")

	if(object$rf!=0)
		stop(" 'object' not fit by 'htb' ")

	yindx=object$yindx
	
	ii=NULL 
	if(!is.null(x))
	{
	
		  if (any(is.na(x))) {
		    	stop("Missing data in 'x'.", call. = FALSE)
  			}
		  if (any(is.na(id))) {
 			   stop("Missing data in 'id'.", call. = FALSE)
			  }
		  if (any(is.na(time))) {
 			   stop("Missing data in 'time'.", call. = FALSE)
 			 }


		if(sum(unlist(lapply(x,class))!=object$vtype&(object$vtype=="character"|object$vtype=="factor"))>0)
			stop("Variable class mismatch with training data.")

		if(ncol(x)!=object$dimx[2])
			stop("Number of columns in 'x' differs from training data.")

		if(nrow(x)==0)
			stop("No rows in 'x'.") 
		
	
		if(length(object$flist)>0)
		{

			# map strings/factors into integers 
			x=format_data(x=x,flist=object$flist)
		}
		if(is.null(id)|is.null(time))
			stop(" Arguments 'id' and 'time' cannot be empty.")
		if(is.character(id))
			id=as.factor(id)

		if(is.factor(id))
			id=as.numeric(id)

		if(!is.numeric(time))
			stop(" 'time' must be numeric.")

		ii=order(id,time)
		x=x[ii,,drop=FALSE]
		id=id[ii]
		time=time[ii]

		  if(yindx>ncol(x)|yindx<0)
			stop(" 'yindx' out of range")
		
	
	}else{
		id=time=NULL
	}

		

	if(!is.null(ntrees)) 
		ntrees=ntrees
	
	if(is.null(ntrees)&(!is.null(object$nboost_cv)))
		ntrees=object$nboost_cv

	pred=predict_cfn(object=object,x=x,id=id,time=time,all.trees=all.trees,se=se,ntrees=ntrees)

	if(is.null(ncol(pred)))
		pred=matrix(pred,ncol=1)

	if(type=="response"&object$control$family=="bernoulli")
	{

		if(all.trees)
		{
	##		cat(" For 'all.trees=TRUE' predictions on logit scale. \n")
		}else{
			pred=1/(1+exp(-pred))

		}

	}



	retlist=pred
	if(!is.null(ii))
		retlist[ii,]=pred

	if(is.null(x))  # -- put predictions on training data in original order 
		retlist[object$indx_original,]=pred


	
	
	retlist
}



