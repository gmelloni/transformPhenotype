	#### LOAD FUNCTIONS###
mymapvalues <- function (x, from, to, warn_missing = FALSE) 
{
    if (length(from) != length(to)) {
        stop("`from` and `to` vectors are not the same length.")
    }
    if (!is.atomic(x)) {
        stop("`x` must be an atomic vector.")
    }
    if (is.factor(x)) {
        levels(x) <- mapvalues(levels(x), from, to, warn_missing)
        return(x)
    }
    mapidx <- match(x, from)
    mapidxNA <- is.na(mapidx)
    from_found <- sort(unique(mapidx))
    if (warn_missing && length(from_found) != length(from)) {
        message("The following `from` values were not present in `x`: ", 
            paste(from[!(1:length(from) %in% from_found)], collapse = ", "))
    }
    x[!mapidxNA] <- to[mapidx[!mapidxNA]]
    x
}

# Given a string like filt="<10,>20"
# And a dataframe with ID and trait
# It returns a vector of sample IDs that are outliers with respect to filt
applyFilters <- function(filt,dfo){
		filtList<-unlist(strsplit(filt,","))
		outliers <- character()
		if(length(filtList)>1){
			for(x in 1:length(filtList)){
				cf <- filtList[x]
				form <- paste("dfo$ID[dfo$trait",cf,"]",sep="")
				outliers <- c(outliers,as.character(eval(parse(text=form))))
			}
		}else{
			form <- paste("dfo$ID[dfo$trait",filtList,"]",sep="")
			outliers <- c(outliers,as.character(eval(parse(text=form))))
		}
		return(outliers)
	}

# This function takes the covariates, the original trait and the normalized trait
# It applies a linear model of the type normx ~ sex + age + ...
applySignifCovariates<-function(covs,x,normx){
	covsUsed<-tolower(unlist(strsplit(covs,",")))
	age2Flag<-match("age2",covsUsed)
	cl<-list()
	if(!is.na(age2Flag)){
		age2 <- x$age^2
		l<-list(name="age2",pval=NA)
		cl<-c(cl,list(l))
		covsUsed<-covsUsed[-match("age2",covsUsed)]
	}
	if(length(covsUsed)>0){
		for(xx in 1:length(covsUsed)){
			l<-list(name=covsUsed[xx],pval=NA)
			assign(covsUsed[xx],as.numeric(x[,grep(paste("^",covsUsed[xx],"$",sep=""),tolower(names(x)))]))
			cl<-c(cl,list(l))
		}
	}
	# Format regression string
	for(yy in 1:length(cl)){
		if(yy==1){
			f<-paste("normx~",cl[[yy]]$name,sep="")
			
		}else{
			f<-paste(f,"+",cl[[yy]]$name,sep="")
		}
	}
	covStr<- as.formula(f)
	model <- summary(lm(covStr))
	res<-as.numeric(model$residuals)
	retList <- list("residuals" = res, "covStr" = covStr)
	return(retList)
}

# Return the summary of the normilized data regressed over covariates
checkCovariates2<-function(covs,x,normx){
	returnStr <- character()
	covListInternal<-tolower(unlist(strsplit(covs$covariates,",")))
	cl<-list()
	for(xx in covListInternal){
			l<-list(name=xx,pval=NA)
			# assign(covListInternal[xx],x[,grep(paste("^",covListInternal[xx],"$",sep=""),tolower(names(x)))])
			assign(xx , x[ , xx])
			cl<-c(cl,list(l))
	}
	# Format regression string
	for(yy in 1:length(cl)){
		if(yy==1){
			f<-paste("normx~",cl[[yy]]$name,sep="")	
		}else{
			f<-paste(f,"+",cl[[yy]]$name,sep="")
		}
	}
	covStr<- as.formula(f)
	tst <- summary(lm(covStr))
	return(tst)
}

# Given covariates and normalized data, return a vector of covariates
# significant in the linear model
checkCovariates<-function(covs,x,normx,sxf){
	returnStr <- character()
	covListInternal<-tolower(unlist(strsplit(covs$covariates,",")))
	cl<-list()
	
	# if(!is.na(covs$age2Flag)){
	# 	age2 <- x$age^2
	# 	l<-list(name="age2",pval=NA)
	# 	cl<-c(cl,list(l))
	# 	# covListInternal<-covListInternal[-match("age2",covListInternal)]
	# }
	for(xx in covListInternal){
			l<-list(name=xx,pval=NA)
			# assign(covListInternal[xx],x[,grep(paste("^",covListInternal[xx],"$",sep=""),tolower(names(x)))])
			assign(xx , x[ , xx])
			cl<-c(cl,list(l))
	}
	# Format regression string
	for(yy in 1:length(cl)){
		if(yy==1){
			f<-paste("normx~",cl[[yy]]$name,sep="")	
		}else{
			f<-paste(f,"+",cl[[yy]]$name,sep="")
		}
	}
	covStr<- as.formula(f)
	tst <- summary(lm(covStr))
	# print(str(tst))
	# lc<-length(tst$coefficients)
	# print(lc)
	# print(tst$coefficients)
	# numCovs<-(lc/4)-1
	for(zz in 1:length(cl)){
		# pv<-tst$coefficients[(lc-numCovs)+zz]
		pv<-tst$coefficients[cl[[zz]]$name , "Pr(>|t|)"]
		if(pv<0.05){
			returnStr<-c(returnStr,cl[[zz]]$name)
		}
	}
	if(length(returnStr)==0){
		returnStr <- NA
	}
	return(returnStr)
}

checkResiduals<-function(x,dataObject,sex){
	ztestres<-(x-mean(x, na.rm=T)) / sd(x, na.rm=T)
	pv<-as.numeric(unlist(shapiro.test(ztestres))[2])
	# if(is.na(initialResidualsPvalGL)){
	# 	initialResidualsPvalGL <<- paste(sex,format(pv,digits=3,sci=TRUE))
	# }else{
	# 	initialResidualsPvalGL <<- c(initialResidualsPvalGL,paste(",",sex,format(pv,digits=3,sci=TRUE)))
	# }
	# if(pv < 0.05){
	# 	resList<-promptUserResiduals(x,sex,pv,ztestres)
	# 	# Ask user whether they want to renormalize the residuals
	# 	x=resList$x
	# 	ztestres=resList$zx
	# }
	ready_residuals<-cbind(ID=as.character(dataObject[,1]),res=x,zres=ztestres)
	retList <- list("outputData" = ready_residuals, "normResiduals" = x)
	return(retList)
}

covariateAvailabilityCheck<-function(traitNames,cvs){
	age2Flag<-NA
	okCovs<-NA
	cl<-unlist(strsplit(cvs,","))
	matched<-which(tolower(cl) %in% tolower(traitNames))
	if(length(matched)>0)okCovs<-cl[matched]
	a2<-match("age2",cl)
	if(!is.na(a2)){
		okCovs<- c(okCovs,"age2")
		matched<-c(matched,a2)
		age2Flag<-TRUE
	}
	
	print("These covariates are not valid:")
	if(length(matched)>0){
		print(cl[-matched])
	}else{
		print(cl)
	}
	retList<-list("covariates"=okCovs,"age2Flag"=age2Flag)
	return(retList)
}
	
createDF <- function(traitName,cvl,rd){
	xDF <- data.frame(ID=rd$sample_id
					, sex=rd$sex
					, trait=rd[ , as.character(traitName)])
	if(!is.na(cvl$covariates[1])){
		for(i in cvl$covariates){
			if(i!="sex" & i %in% colnames(rd))
				xDF[ , i] <- rd[ , i]
			else if(i=="age2")
				xDF[ , "age2"] <- rd[ , "age"]^2
		}
	}
	# remove NAs from trait column on
	xDF <- xDF[ !apply(xDF[ , 3:ncol(xDF) , drop=FALSE] , 1 , function(x) any(is.na(x))) , ]
	# Make sure the covariate data is clean i.e. no values < or > 5sds
	if(dim(xDF)[2]>3){
		for(yy in 4:dim(xDF)[2]){
			lwr<-mean(xDF[,yy], na.rm=T)-5*sd(xDF[,yy], na.rm=T)
			upr<-mean(xDF[,yy], na.rm=T)+5*sd(xDF[,yy], na.rm=T)
			outliers<-numeric()
			outliers<-which(xDF[,yy]<lwr)
			outliers<-c(outliers,which(xDF[,yy]>upr))
			if(length(outliers)>0)xDF<-xDF[-outliers,]
		}
	}
	return(xDF)
}


# given SD number , direction and traitObject dataframe
# return a vector of outliers sample ids
findSdOutliers <- function(sds,sdd,dfo){
		negative<-mean(dfo$trait, na.rm=T)-sds*sd(dfo$trait, na.rm=T)
		positive<-mean(dfo$trait, na.rm=T)+sds*sd(dfo$trait, na.rm=T)
		if(sdd=="both"){
			outliersIds <- as.character(dfo$ID)[dfo$trait<negative]
			outliersIds <- c(outliersIds , as.character(dfo$ID)[dfo$trait>positive])
		}
		if(sdd=="lt"){
			outliersIds <- as.character(dfo$ID)[dfo$trait<negative]
		}
		if(sdd=="gt"){
			outliersIds <- as.character(dfo$ID)[dfo$trait>positive]
		}
		return(outliersIds)
}
keyFieldsCheck <- function(traitMatrix){
	keyFields <- c("age","sex","sample_id")

	kfIndex <- which(keyFields %in% names(rawData))
	print("These fields have not been found: ")
	print(keyFields[-kfIndex])
	if(length(keyFields[-kfIndex])>0){
	  return(FALSE)
	 }else{
	  return(TRUE)
	 }
}

# example input traitObject$trait[males],traitLabel,traitLabelFull,"Males","Raw"
plotRawData<-function(x,tl,tlf,sex,type){
	if(length(x)>5000){
		x <- x[1:5000]
	}
	# setwd(file.path(workDir, paste(type,"Data-Plots",sep="")))
	# out<-paste(tl,"_",sex,"Plots",type,".png",sep="")
	# png(out,height=600,width=600)
	par(mfrow=c(2,2))
	plot(1:length(x),x,xlab="Index",ylab=tlf)
	boxplot(x,main="Boxplot")
	mymin=round(min(x),2)
	mymax=round(max(x),2)
	mymean=round(mean(x),2)
	mysd=round(sqrt(var(x)),2)
	mymain=sprintf("min=%s;max=%s;mean=%s;sd=%s",mymin,mymax,mymean,mysd)
	hist(x,main=mymain,xlab=tlf,prob=TRUE)
	m <- mean(x, na.rm=TRUE)
	std<-sd(x,na.rm=TRUE)	
	curve(dnorm(x,mean=m,sd=std),add=T, col="blue", lwd=1)
	#pv<-as.numeric(unlist(shapiro.test(x))[2])
	pv<-as.numeric(unlist(ad.test(x))[2])
	mymain=sprintf("Normal Q-Q plot (Anderson-Darling pval=%s)",format(pv,digits=3,sci = T))
	qqnorm(x,main=mymain)
	qqline(x)
	par(mfrow=c(1,1))
	invisible()
	dev.off()
	setwd(workDir)
	
	return(pv)
}

plotResidualData<-function(tl,tlf,cvs,res,normx){

		# out<-paste(tl,"_",project,"_",sex,"plotsRes.png",sep="")
		# png(out,height=600,width=600)
		x<-(res-mean(res, na.rm=T)) / sd(res, na.rm=T)
		#x<-res
		par(mfrow=c(2,2))
		plot(1:length(x),x,xlab="Index",ylab=tl)
		mymin=round(min(x),2)
		mymax=round(max(x),2)
		mymean=round(mean(x),2)
		mysd=round(sqrt(var(x)),2)
		mymain=sprintf("min=%s;max=%s;mean=%s;sd=%s",mymin,mymax,mymean,mysd)	
		hist(x,main=mymain,xlab=tl,prob=TRUE)
		m <- mean(x, na.rm=TRUE)
		std<-sd(x,na.rm=TRUE)
		curve(dnorm(x,mean=m,sd=std),add=T, col="forestgreen", lwd=2)
		# pv<-as.numeric(unlist(shapiro.test(x))[2])
		pv<-as.numeric(unlist(ad.test(x))[2])
		mymain=sprintf("Normal Q-Q plot (Anderson-Darling pval=%s)",format(pv,digits=3,sci = TRUE))
		qqnorm(x,main=mymain)
		qqline(x , col="red" , lwd=2)
		if(!is.null(cvs)){
			test<-lm(cvs)
			plot(test,which=1)
		} else {
			plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      		text(x = 0.5, y = 0.5, "No Significant Covariate",cex = 1.6, col = "black")
		}
		# par(mfrow=c(1,1))
		# invisible()
		# dev.off()
}

plotResidualDataBySex <- function(tl,tlf,myList){
	par(mfrow=c(4,2))
	for(i in names(myList)){
		color <- if(i=="Males") "navy" else if(i=="Females") "red"
		cvs <- myList[[i]]$covString
		res <- myList[[i]]$normResiduals
		normx <- myList[[i]]$normData
		x <- (res-mean(res, na.rm=T)) / sd(res, na.rm=T)		
	# Plot standardized dots in their original order
		plot(1:length(x),x,xlab="Index",ylab=tl,main=i,col=color)
	# Plot distribution of standardized	points
		mymin=round(min(x),2)
		mymax=round(max(x),2)
		mymean=round(mean(x),2)
		mysd=round(sqrt(var(x)),2)
		mymain=sprintf("min=%s;max=%s;mean=%s;sd=%s",mymin,mymax,mymean,mysd)
		mymain <- paste(i , mymain , sep="\n")
		hist(x,main=mymain,xlab=tl,prob=TRUE)
		m <- mean(x, na.rm=TRUE)
		std<-sd(x,na.rm=TRUE)
		curve(dnorm(x,mean=m,sd=std),add=T, col="forestgreen", lwd=2)
	# Plot standardized quantiles compared with normal quantiles
		pv<-as.numeric(unlist(ad.test(x))[2])
		mymain=sprintf("Normal Q-Q plot (Anderson-Darling pval=%s)",format(pv,digits=3,sci = TRUE))
		mymain <- paste(i , mymain , sep="\n")
		qqnorm(x,main=mymain)
		qqline(x , col=color , lwd=2)
	# In case of significant covariates, plot fitted residuals
		if(!is.null(cvs)){
			test<-lm(cvs)
			plot(test,which=1)
		} else {
			plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      		text(x = 0.5, y = 0.5, "No Significant Covariate",cex = 1.6, col = "black")
		}
	}
}

sexStratification <- function(X,m,f,label,labelFull){
	denmales<-density(as.numeric(X$trait[m]))
	denfemales<-density(as.numeric(X$trait[f]))
	# range(denmales$x)
	# range(denfemales$x)
	# range(denmales$y)
	# range(denfemales$y)
	# Check whether males and females have a different distribution
	# ks.test(X$trait[m], X$trait[f],alternative="two.sided")
	# wilcox.test(X$trait[m], X$trait[f],alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, conf.int = TRUE, conf.level = 0.95)
	# After analysing just used two.sided tests (i.e. testing whether the distributions are equal to each other)
	xtrial <- c(range(denmales$x),range(denfemales$x))
	x <- seq(min(xtrial),max(xtrial),length=length(denmales$x))
	ytrial <- c(range(denmales$y),range(denfemales$y))
	y <- seq(min(ytrial),max(ytrial),length=length(denmales$y))
	z<-wilcox.test(X$trait[m], X$trait[f],alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, conf.int = TRUE, conf.level = 0.95)
	pv<-signif(as.numeric(unlist(z)[[2]]),3)
	return(pv)
}
sexStratificationPlot <- function(X,m,f,label,labelFull,project){
	denmales<-density(as.numeric(X$trait[m]))
	denfemales<-density(as.numeric(X$trait[f]))
	# range(denmales$x)
	# range(denfemales$x)
	# range(denmales$y)
	# range(denfemales$y)
	# Check whether males and females have a different distribution
	# ks.test(X$trait[m], X$trait[f],alternative="two.sided")
	# wilcox.test(X$trait[m], X$trait[f],alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, conf.int = TRUE, conf.level = 0.95)
	# After analysing just used two.sided tests (i.e. testing whether the distributions are equal to each other)
	xtrial <- c(range(denmales$x),range(denfemales$x))
	x <- seq(min(xtrial),max(xtrial),length=length(denmales$x))
	ytrial <- c(range(denmales$y),range(denfemales$y))
	y <- seq(min(ytrial),max(ytrial),length=length(denmales$y))
	z<-wilcox.test(X$trait[m], X$trait[f]
		,alternative = "two.sided"
		, mu = 0
		, paired = FALSE
		, exact = NULL
		, correct = TRUE
		, conf.int = TRUE
		, conf.level = 0.95)
	pv<-signif(as.numeric(unlist(z)[[2]]),3)
	# w<-ks.test(X$trait[m], X$trait[f],alternative="two.sided")
	# Use a pval cut off of 0.05 to make the decision as whether to process males and females separately
	# pv2<-signif(as.numeric(unlist(w)[[2]]),3)
	# setwd(file.path(workDir, "SexComparison-Plots"))
	# file<-paste(label,"_",project,"_densityMalesFemales","raw.png",sep="")
	# png(file,height=600,width=600)
	# mymain=paste("UK10K Replication - ",project," ",label," Wilcoxon pval=",pv,sep="")
	mymain <- paste(project,":","Wilcoxon pval=",pv)
	plot(x,y,type="n",ylab="Density",xlab=labelFull,main=mymain)
	lines(denmales$x,denmales$y,type="l",col="navy",lwd=4)
	lines(denfemales$x,denfemales$y,type="l",col="red",lwd=4)
	# CHANGE THE FIRST TWO NUMBERS - MAKE THE AXES APPROPRIATE FOR THE DATA
	legend("topleft", c("males","females"), fill=c("navy","red"),bty="n")
	# dev.off()
	# setwd(workDir)
	# return(pv)
}

traitAvailabilityCheck<-function(traits,protocol){
	# First check the trait names to ensure the data matches
	traitsList<-protocol[-1,1]
	traitNames<-vector()
	for(a in levels(traitsList)){
	  traitNames <- c(traitNames,gsub("^\\s+|\\s+$", "", a)) 
	}
	matched<-which(tolower(traitNames) %in% tolower(traits))
	print("The following traits will be processed: ")
	print(traitNames[matched])
	print("These are not available: ")
	print(traitNames[-matched])
	processTraits <- traitNames[matched]
	toDo <-which(protocol[,1] %in% processTraits)
	# FIX, small error, This function calls protocolFile that is protocol
	# protocolFile is in global env
	traitDetails<-protocol[toDo,]
	# traitDetails<-protocolFile[toDo,]
	return(traitDetails)
}
	## END OF FUNCTIONS

#---------------------#
# DISCARDED FUNCTIONS #
#---------------------#

# boxcox <- function(y, lambda){
# 	if(lambda==0){return(log(y))}
# 	else{return((y^lambda-1)/lambda)}
# }

# normaldistance <- function(y){
# 	q=qqnorm((y-mean(y))/sd(y), plot.it=F);
# 	return(mean(abs(q$x-q$y)/sqrt(2)));
# }

# bogusoptim <- function(lambda, y){
# 	return(normaldistance(boxcox(y, lambda)));
# }

# boxit <- function(y, tolerance=1){
# 	return(optimise(bogusoptim, c(-5,5), tol=tolerance, y));
# }


# # normalizeTraitData <- function(x,tm,sex){
# normalizeTraitData <- function(x,tm){
# 	if(tm == "box_cox"){
# 		power_vals <- c(-5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5)
# 		power_opt <- c()
# 		for (i in 1:length(power_vals)){
# 			power_opt <- c(power_opt, as.numeric(bogusoptim(power_vals[i], x)))
# 		}
# 		goodness <- min(power_opt)
# 		power <- power_vals[which.min(power_opt)]
# 		norm_data <- boxcox(x, power)
# 	}else if(tm == "inverse_normal"){
# 		norm_data <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
# 	}else if(tm == "log"){
# 		norm_data <- log(x)

# 	}else if(tm == "square"){
# 		norm_data <- x^2
# 	}else if(tm == "untransformed"){
# 		norm_data<-x
# 	}else{
# 		norm_data<-x
# 	}
# 	return(norm_data)
# }

# findSdOutliers <- function(sds,sdd,dfo){
# 		line1<-mean(dfo$trait, na.rm=T)-sds*sd(dfo$trait, na.rm=T)
# 		line2<-mean(dfo$trait, na.rm=T)+sds*sd(dfo$trait, na.rm=T)

# 		if(sdd == "both"){
# 			outliersIds <- as.character(dfo$ID[which(dfo$trait<line1)])
# 			outliersIds <- c(outliersIds, as.character(dfo$ID[which(dfo$trait>line2)]))
# 				sdGL <<-paste(sds,"SDs",sep="")
# 		}else if(sdd == "gt"){
# 			outliersIds <- as.character(dfo$ID[which(dfo$trait>line2)])
# 				sdGL <<-paste(">",sds,"SDs",sep="")

# 		}else if(sdd == "lt"){
# 			outliersIds <- as.character(dfo$ID[which(dfo$trait<line1)])
# 				sdGL <<-paste("<",sds,"SDs",sep="")
# 		}else{
# 			if(is.na(commentsGL)){
# 				commentsGL <<- "An incorrect value has been entered in the sd_dir column of the protocol file. Should be both, gt or ln."
# 			}else{
# 				commentsGL <<-c(commentsGL,paste(".An incorrect value has been entered in the sd_dir column of the protocol file. Should be both, gt or ln."))
# 			}
# 		}
# 		return(outliersIds)
# 	}


# promptUserCovariates<-function(cvr,cvs,type){
# 	print(paste("For",traitLabel,"in",type,"not all covariates are significant. These were tested:"))
# 	print(cvr)
# 	print("These are significant:")
# 	print(cvs)
# 	print(paste("Enter 1 to use just the significant covariates or 2 to override with",paste(cvr,collapse=" ")))
# 	num <- readLines(con="stdin", 1)
# 	if(as.numeric(num) == 2){
# 		cvs <- cvr
# 	}else if (as.numeric(num) != 1){
# 		print("Invalid option entered - using just significant covariates")
# 	}
# 	return(cvs)
# }

# promptUserNorm<-function(meth,x){
# 	if(is.na(commentsGL)){
# 		commentsGL <<- paste(meth," didn't work")
# 	}else{
# 		commentsGL <<-c(commentsGL,paste(meth," didn't work"))
# 	}
# 	print(paste(meth," was used but the data isn't normally distributed."))
# 	print(paste("Enter 1 to inverse normalize the data or 2 to leave the normalization method as", meth))
# 	num <- readLines(con="stdin", 1)
# 	if(as.numeric(num) == 1){
# 		x <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
# 		meth<- "inverse_normal"
# 	}else if (as.numeric(num) != 2){
# 		print("Invalid option entered - leaving data as is.")
# 	}
# 	retList <- list("nd" = x, "tm" = meth)
# 	commentsGL <<-c(commentsGL,paste(" used",meth))
# 	return(retList)
# }

# promptUserResiduals<-function(x,sex,p,z){
# 	print(paste("The standardized residuals for ",sex," are not normally distributed. Shapiro p-value=",p))
# 	print("Enter 1 to inverse normalize the residuals or 2 to leave the data as is")
# 	num <- readLines(con="stdin", 1)
# 	if(as.numeric(num) == 1){
# 		x <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
# 		z<-(x-mean(x, na.rm=T)) / sd(x, na.rm=T)
		
# 		if(sex == ""){
# 			sex <- "Y"
# 		}
# 		if(is.na(renormalizeResidualsGL)){
# 			renormalizeResidualsGL <<- sex
# 		}else{
# 			renormalizeResidualsGL <<- c(renormalizeResidualsGL,sex)
# 		}
# 	}else if (as.numeric(num) != 2){
# 		print("Invalid option entered - leaving data as is.")
# 	}
# 	retList <- list("x" = x, "zx" = z)
# 	return(retList)
# }