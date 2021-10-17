

hrf=function(x,time=NULL,id=NULL,yindx,ntrees=100,method="freq",mtry=NULL,se=FALSE,B=100,R=10,nsamp=5,historical=TRUE,vh=NULL,vc=NULL,delta=NULL,classify=FALSE,control=list())
{

	## arguments can be passed 'directly' or via control-list, control arguments are used if both supplied
	dtry=nsamp
	argnames=c("ntrees","mtry","se","B","R","method","vh","vc","delta","classify","dtry","historical")
	for(a in argnames)
	{

		if(is.null(control[[a]]))
			eval(parse(text=paste("control[[a]]=",a,sep="")))
		eval(parse(text=paste(a,"=control[[a]]",sep="")))
	}


	if(is.null(control$type))
		control$type=method_convert(method)


	if(is.null(mtry)&is.null(control$vc))
		mtry=floor(sqrt(ncol(x)-1))

	if(is.null(mtry)&(!is.null(control$vc)))
		mtry=floor(sqrt(length(control$vc)))

	if(!is.null(mtry)){
		mtry=min(c(mtry,ncol(x)))
		if(!is.null(control$vc))
			mtry=min(c(mtry,length(control$vc)))
	}
	control$mtry=mtry


	if((is.null(id)|is.null(time)))
	{
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

	if(is.null(control$ncores))
	{
		control$ncores=detectCores()
	}
	cluster=makeCluster(control$ncores)
	on.exit(stopCluster(cluster))
	nti=ceiling(control$ntrees/control$ncores)

	
	pc=proc_control(x=x,yindx=yindx,id=id,time=time,control=control)

	NCOL_TREE=16
	nsplit=nrow(x)


	ntrees=nti
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
	sample_fraction=pc$sample_fraction

	


	ncat=-10
	if(classify)
	{
		ncat=length(unique(x[,yindx]))
		if(ncat>50)
			stop("Max number of classes is 50.")
	}

	h=parLapply(cl=cluster, X=1:control$ncores,cfaux,x,yindx,id,time,ntrees,mtry,nsplit,
		nodesize=nodesize,concurrent=concurrent,nconcurrent=nconcurrent,historic,nhistoric,
		delta,delta0,type,quantiles,dtry,qtry,ncat=ncat,sample_fraction=sample_fraction)


	pred_oob=NULL
	ss=NULL
	if(ncat<0){
		## -- oob error 
		oob<-NULL
		hh=paste("oob=cbind(",paste("h[[",1:length(h),"]]$oob",collapse=","),")",sep="")
		eval(parse(text=hh))
		pred<-NULL
		hh=paste("pred=cbind(",paste("h[[",1:length(h),"]]$pred",collapse=","),")",sep="")
		eval(parse(text=hh))

		cr=apply(oob,1,cumsum)
		pr=apply(oob*pred,1,cumsum)
		mr=pr/cr
		yy=kronecker(t(rep(1,ncol(oob))),x[,yindx])
		ss=apply((yy-t(mr))^2,2,mean,na.rm=T)
		#pred_oob=t(mr)[,ntrees]
		# above line is an error, should be following (18/9/2019)
		pred_oob=t(mr)[,ncol(oob)]

	}

	#pc$vh=pc$historic
	#pc$vc=pc$concurrent
	retlist=list(h=h,error=ss,control=pc,yindx=yindx,x=x,id=id,time=time,pred_oob=pred_oob)
	
	
	class(retlist)="htree"
	retlist$dimx=dim(x)
	retlist$rf=1	
	retlist$ncat=ncat
	retlist$vtype=vtype
	retlist$flist=flist
	x=deformat_data(x=x,flist=flist,vtype=vtype)
	retlist$x=x	
	# -- put in original order 
	retlist$id[ii]=id
	retlist$time[ii]=time
	retlist$x[ii,]=x
	retlist$mtry=mtry
	

	if(ncat>0)
		retlist$error=class_error_hrf(object=retlist)
	

	if(se)
		retlist$boot=hrf_boot(object=retlist,B=control$B,R=control$R)

	
	retlist
}



class_error_hrf=function(object)
{

	# Helper function: compute classification oob error for all trees 
	pp=predict_hrf(object=object,x=object$x,id=object$id,time=object$time,all.trees=T)

	oobm<-NULL
	## -- combine oob matrices
	hh=paste("oobm=cbind(",paste("object$h[[",1:length(object$h),"]]$oob",collapse=","),")",sep="")
	eval(parse(text=hh))

	## dup rows
	oobm_rep=oobm[sort(rep(1:nrow(object$x),object$ncat)),]
	pm=pp*oobm_rep
	cp=t(apply(pm,1,cumsum))
	coob=t(apply(oobm_rep,1,cumsum))
	coob[which(coob>0)]=1
	coob[coob==0]=NA
	cp=cp*coob

	y=object$x[[object$yindx]]
	ncat=object$ncat
	gg=function(x){
		xm=matrix(x,ncol=ncat,byrow=TRUE)
		majority_vote=apply(xm,1,which.max)
		majority_vote=as.numeric(majority_vote)-1
		mean(majority_vote!=y,na.rm=T)
	}

	error=as.numeric(unlist(lapply(as.data.frame(cp),gg)))


	error

}



hrf_boot=function(object,B=50,R=10)
{
	x=object$x
	yindx=object$yindx
	time=object$time
	id=object$id
	mtry=object$control$mtry
	control=object$control
	x=format_data(x=x,flist=object$flist)

	if(is.null(control$ncores))
	{
		control$ncores=detectCores()
	}
	cluster=makeCluster(control$ncores)
	on.exit(stopCluster(cluster))

	pc=proc_control(x=x,yindx=yindx,id=id,time=time,control=control)

	NCOL_TREE=16
	nsplit=nrow(x)


	ntrees=R
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
	boot=TRUE
	sample_fraction=pc$sample_fraction

	h=parLapply(cl=cluster, X=1:B,cfaux,x,yindx,id,time,ntrees,mtry,nsplit,
		nodesize=nodesize,concurrent=concurrent,nconcurrent=nconcurrent,historic,nhistoric,
		delta,delta0,type,quantiles,dtry,qtry,boot=boot,sample_fraction=sample_fraction)


	h
}









factor_conversion=function(y)
{
	# convert factor-levels to integer: 0:(length(levels(y))-1)
	# 


	flist=list()
	flist[levels(y)]=c(0:(length(levels(y))-1))
	ilist=list()
	ilist[as.character(c(0:(length(levels(y))-1)))]=levels(y) 

	yint=as.numeric(unlist(flist[y]))
	
	classInfo=list()
	classInfo=list(levels=levels(y),flist=flist,ilist=ilist)
	rlist=list(classInfo=classInfo,y=yint)

	rlist

}


reformat_data=function(x)
{
	# Convert factors/character-variables to numeric 

	fvar=NULL
	flist=list()
	nx=names(x)
	for(k in nx)
	{
		if(is.character(x[[k]]))
			x[[k]]=as.factor(x[[k]])

		if(is.factor(x[[k]]))
		{
			h=factor_conversion(x[[k]])
			x[[k]]=h$y
			flist[[k]]=h$classInfo
		}
		
	} 

	list(x=x,flist=flist)
}


deformat_data=function(x,flist,vtype)
{

	if(!is.null(flist))
	{

		for(k in names(flist)){
			z=as.character(unlist(flist[[k]]$ilist[as.character(x[[k]])]))
			vtypek=vtype[which(names(vtype)==k)]
			txt=paste("x[[k]]=as.",vtypek,"(z)",sep="")
			eval(parse(text=txt))
			#x[[k]]=as.factor()

		}
		
	}
	

	x
}


format_data=function(x,flist)
{
	# Apply data-formatting used in model fit 
	# .. convert strings/factors to integers 

	if(!is.null(flist))
	{

		for(k in names(flist))
		{
			uc=unique(as.character(x[[k]]))
			if(sum(!is.element(uc,flist[[k]]$levels))>0)
				stop(paste(" Factor levels in ",k," not in training data, exiting...",sep="")) 
			
			x[[k]]=(as.numeric(unlist(flist[[k]]$flist[as.character(x[[k]])])))
			
		}


	}


	x
}






predict_hrf=function(object,x=NULL,time=NULL,id=NULL,all.trees=FALSE,se=FALSE)
{

	yindx=NULL

	if((is.null(x)))
	{
		id=object$id
		time=object$time
		x=object$x
	}


	if((is.null(id)|is.null(time)))
	{
		id=c(1:nrow(x))
		time=rep(1,nrow(x))
	}
	


	if(class(object)!="htree")
		stop(" 'object' not of class 'htree' ")

	if(object$rf!=1)
		stop(" 'object' not fit by 'hrf' ")



	if(is.null(yindx))
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

		
		xclass=unlist(lapply(x,class))
		## print(xclass)
		if(sum(xclass!=object$vtype&(is.element(xclass,c("factor","character"))|is.element(object$vtype,c("factor","character"))))>0)
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

		

	pred=predict_cfn(object=object,x=x,id=id,time=time,all.trees=all.trees,se=se)

	if(is.null(ncol(pred)))
		pred=matrix(pred,ncol=1)

	retlist=pred
	if((!is.null(ii))&object$ncat<0){
		retlist[ii,]=pred
	}

	if((object$ncat>0)){

		if(all.trees){
		fa=function(x,ncat) return(((x-1)*ncat+1):(x*ncat))
		jj=as.numeric(mapply(fa,ii,ncat=rep(object$ncat,length(ii))))
		retlist=pred
		retlist[jj,]=pred
		}else{
			retlist[ii,]=pred
			colnames(retlist)=paste("Prob_",0:(object$ncat-1),sep="")
		}

	}
	
	retlist
}



partdep_hrf=function(object,xindx,xlim=NULL,ngrid=10,subsample=.5,plot.it=TRUE,which.class=1,cond=NULL)
{
	main="Marginal effect"

	if(object$ncat>0)
	{
		if(class(which.class)=="character")
			which.class=which(object$classInfo$levels==which.class)

		if(length(which.class)==0)
			stop(" 'which.class' not found. ")

	}

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


	pd=partdep(object=object,xindx=xindx,xlim=xlim,ngrid=ngrid,ntrees=NULL,subsample=subsample,which.class=which.class,se=(!is.null(object$boot)),cond=cond)


	if(plot.it)
	{
		nx=names(object$x)[xindx]
		if((!is.element(object$vtype[xindx],c("factor","character")))&(length(unique(object$x[,xindx]))>5)){
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
















partdep=function(object,xindx,xlim=NULL,ngrid=100,ntrees=NULL,subsample=1,which.class=1,se=FALSE,cond=NULL)
{
	# -- subsample subjects
	# subsample=.5;ngrid=20;object=ff;xindx=20
	uid=unique(object$id)
	ns=round(subsample*length(uid))
	sid=sample(uid,size=ns,replace=F)


	uv=unique(object$x[,xindx])
	nuv=length(uv)
	if((nuv<ngrid)|is.element(object$vtype[xindx],c("character","factor")))
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



	if(is.null(ntrees))
		ntrees=object$ntrees


	tt=sort(rep(tt,nii)) 

	x[,xindx]=tt

	vv=NULL
	xx=as.matrix(x)
	yindx=object$yindx

	hh=predict_hrf(object,x=x,time=time,id=id)



	if(object$ncat>0)
		hh=hh[,which.class]
	

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


	if(se==TRUE&length(object$boot)==0)
	{
		cat("Need to run 'hrf' with 'se=TRUE' to get standard errors.\n")
		se=FALSE

	}

	stderr=NULL
	boot_tree=NULL
	if(se){
		## -- Get standard errors 
		object$h=object$boot
		hh=predict_hrf(object,x=x,time=time,id=id,all.trees=TRUE)
		Mlist=list()
		if(!is.null(cond)&(length(ii)>0))
			hh=hh[ii,]

		
		Mlist$M=hh
		Mlist$tt=tt
	
		M=aggregate(hh,by=list(tt),mean)
		M=(as.matrix(M[,-1]))
		B=length(object$boot);
		R=length(unique(object$boot[[1]]$trees[,1]))

		start=seq(1,((B-1)*R+1),by=R)
		end=start-1+R
		
		mean_v=mapply(function(x,y){ rowMeans(M[,c((x):(y))]) },x=start,y=end)
		var_v=mapply(function(x,y){ rowSums(M[,c((x):(y))]^2) },x=start,y=end)
		var_v=rowSums((var_v-R*mean_v^2)/(R-1))
		bias=var_v/B*(1/R) 
		var_biased=apply(mean_v,1,var)

		## -- bias-correct 
		var_hat=var_biased-bias
		mean_var=mean(var_hat[which(var_hat>=0)],na.rm=T)
		var_hat[var_hat<0]=mean_var
		stderr=sqrt(var_hat)
		
		## -- data frame of partial dependence 
		boot_tree=data.frame(boot=sort(rep(1:B,R*nrow(ds))),tree=rep(sort(rep(1:R,nrow(ds))),B),x=rep(ds[,1],B),y=as.numeric(M))



	}

	#list(y=ds[,2],x=ds[,1],se=stderr,chisq=zzu,chisq_biased=zz,df=df,boot_raw=mean_v,tree_raw=M)
	list(y=ds[,2],x=ds[,1],se=stderr,boot_tree=boot_tree)
}























predict_viaux=function(trees,x,oob,id,time,yindx,concurrent,nconcurrent,historic,nhistoric,delta,delta0,type,quantiles,nperm,rf=1,ncat=-10)
{

#	h_vi=predict_viaux(combt,x,oob,id,time,yindx,
#			concurrent,nconcurrent,historic,nhistoric,
#			delta,delta0,type,quantiles,nperm,ncat=object$ncat)



	NCOL_TREES=16
	ntrees=sum(trees[,2]==(-99))


	npred=nrow(x)*ncol(x)
	if(ncat>0)
		npred=npred*ncat

	h=.C("ctree_varimp",as.double(as.matrix(x)),as.integer(id),as.double(time),as.integer(nrow(x)),
		as.integer(ncol(x)),as.integer(yindx-1),as.integer(ntrees),
		as.integer(concurrent),as.integer(nconcurrent),as.integer(historic),as.integer(nhistoric),
		as.double(delta),as.double(delta0),as.integer(length(delta)),
		as.integer(type),as.integer(ncat),as.double(quantiles),as.integer(length(quantiles)),
		as.double(as.matrix(trees)),as.integer(nrow(trees)),as.integer(as.matrix(oob)),as.integer(nperm),as.integer(rf),pred=double(npred))

	h=matrix(h$pred,ncol=ncol(x),byrow=FALSE)

	h
}



predict_viauxP=function(X,trees,x,oob,id,time,yindx,concurrent,nconcurrent,historic,nhistoric,delta,delta0,type,quantiles,nperm,rf=1,ncat=-10)
{

#	h_vi=predict_viaux(combt,x,oob,id,time,yindx,
#			concurrent,nconcurrent,historic,nhistoric,
#			delta,delta0,type,quantiles,nperm,ncat=object$ncat)



	NCOL_TREES=16
	ntrees=sum(trees[,2]==(-99))


	npred=nrow(x)*ncol(x)
	if(ncat>0)
		npred=npred*ncat

	h=.C("ctree_varimp",as.double(as.matrix(x)),as.integer(id),as.double(time),as.integer(nrow(x)),
		as.integer(ncol(x)),as.integer(yindx-1),as.integer(ntrees),
		as.integer(concurrent),as.integer(nconcurrent),as.integer(historic),as.integer(nhistoric),
		as.double(delta),as.double(delta0),as.integer(length(delta)),
		as.integer(type),as.integer(ncat),as.double(quantiles),as.integer(length(quantiles)),
		as.double(as.matrix(trees)),as.integer(nrow(trees)),as.integer(as.matrix(oob)),as.integer(nperm),as.integer(rf),pred=double(npred))

	h=matrix(h$pred,ncol=ncol(x),byrow=FALSE)

	h
}








varimp_hrf=function(object,nperm=20,parallel=TRUE)
{

	x=object$x
	id=object$id
	time=object$time
	yindx=object$yindx
	control=object$control
	#nti=ceiling(ntrees/control$ncores)

	x=format_data(x=x,flist=object$flist)

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


	## --- ... 
	h=rep(0,nrow(x))
	combt<-NULL
	hh=paste("combt=rbind(",paste("object$h[[",1:length(object$h),"]]$trees",collapse=","),")",sep="")
	eval(parse(text=hh))


	oob<-NULL
	hh=paste("oob=cbind(",paste("object$h[[",1:length(object$h),"]]$oob",collapse=","),")",sep="")
	eval(parse(text=hh))


	if(parallel){
		if(is.null(object$control$ncores))
		{
			object$control$ncores=detectCores()
		}
		cluster=makeCluster(object$control$ncores)
		on.exit(stopCluster(cluster))
		npi=ceiling(nperm/object$control$ncores)

		h=parLapply(cl=cluster, X=1:object$control$ncores,predict_viauxP,combt,x,oob,id,time,yindx,
			concurrent,nconcurrent,historic,nhistoric,
			delta,delta0,type,quantiles,npi,ncat=object$ncat)

		h_vi=array(0,dim=dim(h[[1]]))
		for(k in 1:length(h))
			h_vi=h_vi+h[[k]]
	
		h_vi=h_vi/length(h)
	
	}else{
		h_vi=predict_viaux(combt,x,oob,id,time,yindx,
			concurrent,nconcurrent,historic,nhistoric,
			delta,delta0,type,quantiles,nperm,ncat=object$ncat)

	}

	## predict_viaux=function(trees,x,oob,,id,time,yindx,concurrent,nconcurrent,historic,nhistoric,delta,delta0,type,quantiles,nperm)
	pa=predict_hrf(object,x=object$x,id=object$id,time=object$time,all.trees=T)

	if(object$ncat<0){
	## oob pred
	oobm<-NULL
	hh=paste("oobm=cbind(",paste("object$h[[",1:length(object$h),"]]$oob",collapse=","),")",sep="")
	eval(parse(text=hh))

	# oobm=1-oobm	
	pao=rowSums(oobm*pa)
	noob=rowSums(oobm)
	pred_oob=pao/noob

	zscore=function(z){
		error_perm=mean((x[[object$yindx]]-z)^2)
		error_oob=mean((x[[object$yindx]]-pred_oob)^2)
		d=abs(x[[object$yindx]]-z)-abs(x[[object$yindx]]-pred_oob)
		da=aggregate(d,by=list(id=object$id),mean)[,2]
		zs=mean(da)/(sd(da)/sqrt(length(da)))
		cret=c(error_oob,error_perm,(error_perm-error_oob)/error_oob,zs)
		return(cret)
	}
	h=matrix(unlist(lapply(as.data.frame(h_vi),zscore)),ncol=4,byrow=TRUE)
		dat=data.frame(vn=names(x),error_perm=round(h[,2],4),error_oob=round(h[,1],4),rchange=round(h[,3],3),zscore=round(h[,4],3))
	

	names(dat)=c("Predictor","Marginalized error","Model error","Relative change","Z-value") 
	dat=dat[order(dat[,2],decreasing=TRUE),]
	}else{

		## --- variable importance for classification (based on gini-error) 
		oobm<-NULL
		hh=paste("oobm=cbind(",paste("object$h[[",1:length(object$h),"]]$oob",collapse=","),")",sep="")
		eval(parse(text=hh))
		oobm=oobm[sort(rep(1:nrow(oobm),object$ncat)),]

		# oobm=1-oobm	
		pao=rowSums(oobm*pa)
		noob=rowSums(oobm)
		pred_oob=pao/noob
		pred_oob=matrix(pred_oob,ncol=2,byrow=T)
		pred_oob=pred_oob/rowSums(pred_oob)

		zscore=function(z){
			## GINI error 
			y=object$x[[object$yindx]]
			zm=matrix(z,ncol=object$ncat,byrow=TRUE)
			zm=zm/rowSums(zm)
			z_perm=gini_error(y,zm,object=object)
			error_perm=mean(z_perm) ##(mean((object$x[[object$yindx]]-z)^2)
			z_oob=gini_error(y,pred_oob,object=object)
			error_oob=mean(z_oob)
			d=z_perm-z_oob
			da=aggregate(d,by=list(id=object$id),mean)[,2]
			zs=mean(da)/(sd(da)/sqrt(length(da)))
			cret=c(error_oob,error_perm,(error_perm-error_oob)/error_oob,zs)
			return(cret)
		}
		h=matrix(unlist(lapply(as.data.frame(h_vi),zscore)),ncol=4,byrow=TRUE)
		dat=data.frame(vn=names(x),error_perm=round(h[,2],4),error_oob=round(h[,1],4),rchange=round(h[,3],3),zscore=round(h[,4],3))

		names(dat)=c("Predictor","Marginalized error","Model error","Relative change","Z-value") 
		dat=dat[order(dat[,2],decreasing=TRUE),]

	}
	dat

}










gini_error=function(y,p,object)
{
	m=0	
	for(k in 1:object$ncat)
	{
		pp=p[,k]
		yk=as.numeric(y==(k-1))
		m=m+((yk-pp)^2)
	}
	m=m/object$ncat
	m
}





method_convert=function(method)
{
	if(!is.element(method,c("freq","freqw","frac","fracw","mean0","meanw0")))
	{
		cat(paste(" method: '",method,"' not matched, setting to 'freq'. \n",sep=""))
		type=2
	}

	if(method=="freq")
		type=3
	if(method=="freqw")
		type=6
	if(method=="frac")
		type=2
	if(method=="fracw")
		type=5
	if(method=="mean0")
		type=1
	if(method=="meanw0")
		type=4

	return(type)

}





plot.means <- function(means = NULL,
                           SEs = NULL,
                           categories = NULL,
                           x.axis.label){
  
 ## Modified from  https://rpubs.com/brouwern/plotmeans2
 n.means=length(means)  
   
  if(is.null(SEs)){
	names(means)=categories
	barplot(means,cex.names=1.5,main="Marginal effect",cex.main=2,xlab=x.axis.label,cex.lab=2)
 }else{
  # calculate values for plotting limits 
  if(is.null(SEs)==F){
  y.max <- max(means+2*SEs) 
  y.min <- min(means-2*SEs) 
  }else{
	m.max=max(means)
	m.min=min(means)
	r=(m.max-m.min)/2 
	
  y.max <- max(means+r) 
  y.min <- min(means-r) 
  

  }

  # calculate values for plotting limits 

  
  #determine where to plot points along x-axis
  x.values <- 1:n.means
  
  #set x axis min/max
  x.axis.min <- min(x.values)-0.25
  x.axis.max <- max(x.values)+0.25
  
  x.limits <- c(x.axis.min,x.axis.max)
  
  #Plot means
  plot(means ~ x.values,
       xlim = x.limits,
       ylim = c(y.min,y.max),
       xaxt = "n",
       xlab = x.axis.label,
       ylab = "",
       cex = 1.25,
       pch = 16,main="Marginal effect",cex.main=2,cex.lab=1.5)
  
  #Add x labels
  axis(side = 1, 
       at = x.values,
       labels = categories,cex.axis=1.5
      )
  
   
  if(!is.null(SEs)){
  lwd. <- 2
  arrows(y0 = means,
         x0 = x.values,
         y1 = means+2*SEs,
         x1 = x.values,
         length = 0,
         lwd = lwd.)
  
  #Plot lower error bar
  arrows(y0 = means,
         x0 = x.values,
         y1 = means-2*SEs,
         x1 = x.values,
         length = 0,
         lwd = lwd.) 
  }
  
 }  

  
}


