
#
# Basic 
#


proc_control=function(x,yindx,id,time,control)
{

	EPSILON=.01
	## Process control list 
	if(is.null(control$nodesize))
	{
		nodesize=5
	}else{
		nodesize=control$nodesize
	}

	nsplit=ceiling(nrow(x)/nodesize)
	if(!is.null(control$nsplit))
		nsplit=control$nsplit
	
	if(is.null(control$vh)|(!control$historical))
	{
		if((is.null(id)|is.null(time))|(!control$historical))
		{
			nhistoric=0
			historic=1
		}else{
			## id=x$id;time=x$day
			ff=function(x){mean(abs(x-mean(x)))}
			aa=apply(abs(aggregate(x,by=list(id=id),ff))[,-1],2,sum)
			historic=which(aa!=0)-1
			nhistoric=length(historic)

		}
	}else{
		if(is.null(id)|is.null(time))
		{
			nhistoric=0
			historic=1
			
		}else{
		historic=control$vh
		nhistoric=length(historic)
		if(class(control$vh)=="character")
		{
			historic=(which(is.element(names(x),control$vh))-1)
			nhistoric=length(historic)
		}else{
			historic=control$vh-1

		}
		
		}
	}

	
	if(is.null(control$vc))
	{
		vc=c(0:(ncol(x)-1))
		concurrent=vc[-yindx]

	}else{
		if(class(control$vc)=="character")
		{
			concurrent=(which(is.element(names(x),control$vc))-1)
			nconcurrent=length(concurrent)
		}else{
			concurrent=control$vc-1
		}

	}
	if(sum(concurrent==(yindx-1))>0)
		concurrent=concurrent[-which(concurrent==(yindx-1))]

	nconcurrent=length(concurrent)

	if(is.null(id))
		id=c(1:nrow(x))

	if(is.null(time))
		time=rep(1,nrow(x))


	delta=0
	delta0=0
	dtry=1
	if(nhistoric>0){
		if(is.null(control$delta))
		{
		d=proc_delta(id=id,time=time,ndelta=control$ndelta)
		delta=d$delta
		delta0=d$delta0
		ndelta=length(delta)
		}else{
			delta=control$delta
			if(is.null(control$delta0))
			{
				delta0=rep(min(delta)/2,length(delta))
			}

			ndelta=length(delta)

		}
		delta=delta+EPSILON

		if(is.null(control$dtry))
			control$dtry=5 #min(c(5,ndelta))			

		dtry=min(c(control$dtry,ndelta))

	}

	
	type=1
	if(!is.null(control$type))
		type=control$type

	qtry=1
	nquantiles=0
	quantiles=0
	if(type==2|type==5)
	{
		quantiles=seq(.1,1,by=.1)
		if(!is.null(control$quantiles))
			quantiles=control$quantiles

		nquantiles=length(quantiles)
		qtry=nquantiles
		if((!is.null(control$qtry)))
			qtry=control$qtry

		if(qtry>nquantiles)
			qtry=nquantiles
	}

	if(is.null(control$sample_fraction)){
		sample_fraction=.5
	}else{
		sample_fraction=control$sample_fraction
		if(sample_fraction>1|sample_fraction<=0)
			sample_fraction=.5
	}

	update_list=list(id=id,time=time,concurrent=concurrent,nconcurrent=nconcurrent,historic=historic,nhistoric=nhistoric,
		delta=delta,delta0=delta0,dtry=dtry,nodesize=nodesize,nsplit=nsplit,type=type,quantiles=quantiles,qtry=qtry,sample_fraction=sample_fraction)


	for(n in names(update_list))
		control[[n]]=update_list[[n]]

	control
}



proc_delta=function(id,time,ndelta)
{
	if(is.null(ndelta))
		ndelta=20

	## ---- 
	df=unique(data.frame(id=id,time=time))
	df=as.data.frame(na.omit(df))
	time=df$time
	id=df$id
	EPSILON=.01
	ff=function(x){
		x=sort(x);vv=NULL;
		if(length(x)>1) {
			for(k in 2:length(x))vv=c(vv,x[k]-x[(k-1):1]);}
		return(unique(vv))}
	aa=unlist(aggregate(time,by=list(id=id),ff)[[2]])
	aa=aa[which(!is.na(aa))]
	
	delta=unique(as.numeric(quantile(aa,prob=seq(0,1,length=ndelta)))) 
	delta0=rep(min(delta[which(delta>0)],na.rm=TRUE),length(delta))

	list(delta=delta,delta0=delta0)
}





get_ctree=function(object)
{
	trees=object$trees
	colnames(trees)=c("tree","pnode","node","pred","svar","tau","spoint","delta","rL","rR","derror","size","delta0")
	trees=trees[,c(1:7,13,8:12)]
	trees=as.data.frame(trees)
	vlist=list()
	vlist[as.character(c(1:ncol(object$x))-1)]=names(object$x)
	vlist[["-1"]]="nobs"
	vlist[["-99"]]="tnode"
	trees$svar=unlist(vlist[as.character(trees$svar)])

	trees
}







cfaux=function(X=1,x,yindx,id=NULL,time=NULL,ntrees,mtry,nsplit=100,nodesize,concurrent,nconcurrent,historic,nhistoric,delta,delta0,type,quantiles,dtry,qtry,
	rf=1,lambda=1,fam=1,boot=FALSE,cvfold=0,ncat=-10,sample_fraction=1)
{

	row_train=NULL
	if(cvfold>0)
	{
		## sample cv-training set
		uid=unique(id)
		sui=sample(uid,round(length(uid)*(1-1/cvfold)))
		row_train=which(is.element(id,sui))
		x=x[row_train,]
		id=id[row_train]
		time=time[row_train]
	}


	boot_id=NULL
	if(boot)
	{
		## bootstrap sample of subjects 
		h=boot_aux(x=x,id=id,time=time)
		x=h$x
		id=h$id
		time=h$time
		boot_id=h$boot_id
	}

	NCOL_TREE=16
	multiplier=1
	if(ncat>0){
		NCOL_TREE=NCOL_TREE+ncat
		multiplier=ncat
	}
	

	h=.C("ctree",as.double(as.matrix(x)),as.integer(id),as.double(time),as.integer(nrow(x)),
		as.integer(ncol(x)),as.integer(yindx-1),as.integer(nsplit),as.integer(ntrees),as.integer(mtry),as.integer(nodesize),
		as.integer(concurrent),as.integer(nconcurrent),as.integer(historic),as.integer(nhistoric),
		as.double(delta),as.double(delta0),as.integer(length(delta)),
		as.integer(type),as.double(quantiles),as.integer(length(quantiles)),as.integer(dtry),
		as.integer(qtry),as.integer(rf),as.double(lambda),as.integer(fam),as.integer(ncat),as.double(sample_fraction),
		te=double(ntrees),pred=double(nrow(x)*multiplier),trees=double(ntrees*((nsplit*2+1)*NCOL_TREE)),oob_tree=integer(ntrees*nrow(x)),
		pred_tree=double(ntrees*nrow(x)),trows=integer(1))

	trees=matrix(h$trees,ncol=NCOL_TREE,byrow=FALSE)[1:h$trows,]

	if(ncat<0){
		colnames(trees)[1:16]=c("tree","pnode","node","pred","svar","tau","spoint","delta","rL","rR","derror","size","delta0","hindx","dcol","dcol2")	
	}else{
		cn=c("tree","pnode","node","pred","svar","tau","spoint","delta","rL","rR","derror","size","delta0","hindx","dcol","dcol2")	
		cp=paste("num_",0:(ncat-1),sep="")
		colnames(trees)=c(cn,cp)
	}
	
	list(trees=trees,oob=matrix(h$oob_tree,byrow=F,ncol=ntrees),pred=matrix(h$pred_tree,byrow=F,ncol=ntrees),boot_id=boot_id,pr=h$pred,error_train=h$te,row_train=row_train)
}


boot_aux=function(x,id,time)
{

		
		start=1
		end=length(id)
		jj=which(diff(id)!=0) 
		start=c(start,jj)
		end=c(jj-1,end)
		uid=unique(id)
		iList=list()
		nList=list()

		for(k in 1:length(uid))
		{
			iList[[as.character(uid[k])]]=c(start[k]:end[k])
			nList[[as.character(uid[k])]]=end[k]-start[k]+1
		}

		
		boot_id=sample(uid,length(uid),replace=T)

		new_id=NULL
		for(j in 1:length(boot_id))
			new_id=c(new_id,rep(j,nList[[as.character(boot_id[j])]]))

		ii=as.numeric(unlist(iList[as.character(boot_id)]))

		x=x[ii,]
		id=new_id
		time=time[ii]

	list(x=x,id=id,time=time,boot_id=boot_id)
}












predict_cfn=function(object,x,id,time,all.trees=FALSE,se=FALSE,ntrees=NULL)
{

	yindx=object$yindx
	control=object$control


	#pc=proc_control(x=x,yindx=yindx,id=id,time=time,control=control)
	pc=object$control	

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
	ncat=object$ncat

	h=rep(0,nrow(x))
	boot=FALSE

	if(object$rf==1){
	if(se==TRUE&length(object$boot)==0)
	{
		cat(" Need run 'hrf' with 'se=TRUE' to get standard errors. \n")
		
	}

	if(se==TRUE&length(object$boot)>0){
		boot=TRUE
		all.trees=FALSE
	}
	}

	if(object$rf==1){
		hh=paste("combt=rbind(",paste("object$h[[",1:length(object$h),"]]$trees",collapse=","),")",sep="")
		
		eval(parse(text=hh))
	}else{
		combt=object$trees
		if(!is.null(ntrees))
			combt=combt[which(combt[,1]<ntrees),]
	}

	combt=cbind(combt,rep(1,nrow(combt)))


	h=predict_cfaux(X=1,combt,x,id,time,yindx,
			concurrent,nconcurrent,historic,nhistoric,
			delta,delta0,type,quantiles,ncat=ncat,all.trees)


	nt=sum(combt[,2]==-99)
	
	if((!all.trees)&(object$rf==1))
		h=h/nt
	if(ncat>0&(!all.trees))
	{
		s=apply(h,1,sum)
		h=h/s
	}
	
	pred=h

	if(boot)
	{
		hh=paste("combt=rbind(",paste("object$boot[[",1:length(object$boot),"]]$trees",collapse=","),")",sep="")
		eval(parse(text=hh))
		combt=cbind(combt,rep(1,nrow(combt)))
		all.trees=TRUE
		M=predict_cfaux(X=1,combt,x,id,time,yindx,
				concurrent,nconcurrent,historic,nhistoric,
				delta,delta0,type,quantiles,all.trees=all.trees)

		B=length(object$boot);
		R=length(unique(object$boot[[1]]$trees[,1]))

		start=seq(1,((B-1)*R+1),by=R)
		end=start-1+R

		mean_v=mapply(function(x,y){ rowMeans(M[,c((x):(y))]) },x=start,y=end)
		var_v=mapply(function(x,y){ rowSums(M[,c((x):(y))]^2) },x=start,y=end)
		var_v=rowSums((var_v-R*mean_v^2)/(R-1))


		bias=var_v/B*(1/R) 
		var_biased=apply(mean_v,1,var)
		var_hat=var_biased-bias
		mean_var=mean(var_hat,na.rm=T)
		var_hat[var_hat<0]=mean_var
		se=sqrt(var_hat)


		h=cbind(pred,se) 

	}

	h
}










predict_cfaux=function(X,trees,x,id,time,yindx,concurrent,nconcurrent,historic,nhistoric,delta,delta0,type,quantiles,ncat=-10,all.trees=FALSE)
{

	NCOL_TREES=16
	if(ncat>0)
		NCOL_TREES=NCOL_TREES+ncat

	
	trees=trees[which(trees[,ncol(trees)]==X),]
	trees=trees[,-(ncol(trees))]
	ntrees=sum(trees[,2]==(-99))
	
	mt=1
	if(ncat>0)
		mt=ncat	

	npred=nrow(x)*mt
	allt=0
	if(all.trees){
		npred=nrow(x)*ntrees*mt
		allt=as.numeric(all.trees)
	}

	
	h=.C("ctree_predict",as.double(as.matrix(x)),as.integer(id),as.double(time),as.integer(nrow(x)),
		as.integer(ncol(x)),as.integer(yindx-1),as.integer(ntrees),
		as.integer(concurrent),as.integer(nconcurrent),as.integer(historic),as.integer(nhistoric),
		as.double(delta),as.double(delta0),as.integer(length(delta)),
		as.integer(type),as.integer(ncat),as.double(quantiles),as.integer(length(quantiles)),
		as.double(as.matrix(trees)),as.integer(nrow(trees)),as.integer(allt),pred=double(npred))

	
	if(ncat<0){
	if(all.trees)
	{
		retl=matrix(h$pred,ncol=ntrees,byrow=FALSE)
	}else{
		retl=h$pred
	}
	}else{
	if(all.trees)
	{
		retl=matrix(h$pred,ncol=ntrees,byrow=FALSE)
		
	}else{
		retl=matrix(h$pred,ncol=ncat,byrow=TRUE)
	}


	}
	retl
}

