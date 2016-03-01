# Normalization Functions
# This R source script can be modified by adding new functions

#----------------------------------#
# RULES FOR CREATING NEW FUNCTIONS #
#----------------------------------#
# The function must be added as a new element of this list
# It can accept multiple arguments, but the first one must be 'x', a numerical vector
# The return must be a list of at least 1 named element called 'norm_data' list(norm_data=xnormalized)
# 'norm_data' must be the same length of 'x' and ordered like 'x'
# Shiny doesn't now what you want to do with the other elements of your return
# they will just be output below Data Plot (see Box Cox as an example of this case)
	# e.g. mynewfunction = mynewfunction <- function(x , ...) {
	# 			# do awesome things 
	# 			return(list(norm_data=mynormalizedata , optionalelements=otheroutput))
	# }
# If possible, use the syntax ... after function(x) to allow new functions to have an unlimited number of arguments

normalizationFunctionsList <- list(
	"Untransformed" = function(x , ...) {
		norm_data <- x
		return(list(norm_data=norm_data))
	}
	,"Log base e" = function(x , ...) {
		norm_data <- log(x)
		return(list(norm_data = norm_data))
		}
	,"Log base 10" = function(x , ...) {
		norm_data <- log10(x)
		return(list(norm_data = norm_data))
		}
	,"Inverse Normal" = function(x , ...) {
		norm_data <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
		return(list(norm_data=norm_data))
		}
	,"Square" = function(x , ...) {
		norm_data <- x^2
		return(list(norm_data=norm_data))
		}
	,"Box Cox" = function(x , ...) {
		.normaldistance <- function(y){
			q=qqnorm((y-mean(y))/sd(y), plot.it=F)
			return(mean(abs(q$x-q$y)/sqrt(2)))			
		}
		.bogusoptim <- function(lambda, y){
			return(.normaldistance(.boxcox(y, lambda)))				
		}
		.boxcox <- function(y, lambda){
			if(lambda==0){return(log(y))}
			else{return((y^lambda-1)/lambda)}
		}
		# boxit <- function(y, tolerance=1){
		# 	return(optimise(.bogusoptim, c(-5,5), tol=tolerance, y));				
		# }
		power_vals <- c(-5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5)
		power_opt <- c()
		for (i in 1:length(power_vals)){
			power_opt <- c(power_opt, as.numeric(.bogusoptim(power_vals[i], x)))
		}
		goodness <- min(power_opt)
		power <- power_vals[which.min(power_opt)]
		norm_data <- .boxcox(x, power)
		return(list(norm_data=norm_data , opt_power=power))
	}
	,"Box Cox Guerrero Optimization" = function(x , ...){
		optimLambda <- BoxCox.lambda(x , method="guerrero")
		norm_data = BoxCox(x , optimLambda)
		return(list(norm_data=norm_data , optimLambda=optimLambda))
	}
	,"Box Cox LogLike Optimization" = function(x , ...){
		optimLambda <- BoxCox.lambda(x , method="loglik")
		norm_data = BoxCox(x , optimLambda)
		return(list(norm_data=norm_data , optimLambda=optimLambda))
	}
)

normalizeTraitData <- function(trait 
							, tm=NULL 
							, funcList=normalizationFunctionsList
							, ...){
	# browser()
	if(is.null(tm))
		tm <- "Untransformed"
	applyFUN <- funcList[[tm]]
	out <- applyFUN(x=trait , ...)
	return(out)
}

# To make the angela's script able to run, we have to change the name of the transformations
equivalence_table <- c("Untransformed"
,"Log base e"
,"Log base 10"
,"Inverse Normal"
,"Square"
,"Box Cox"
,"Box Cox Guerrero Optimization"
,"Box Cox LogLike Optimization")
angela_transformations <- c("untransformed" , "log" , "log" , "inverse_normal" , "square" , "box_cox" , "box_cox" , "box_cox")
names(angela_transformations) <- equivalence_table






#------------------#
# ORIGINAL VERSION #
#------------------#

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

# Another version of the function using switch.
# It works, but is a little more complicated to expand

# normalizeTraitData <- function(x , tm) {
# 	return({
# 		switch(EXPR=tm
# 		,untransformed = x
# 		,log = log(x)
# 		,log10 = log10(x)
# 		,box_cox = {
# 			.normaldistance <- function(y){
# 				q=qqnorm((y-mean(y))/sd(y), plot.it=F);
# 				return(mean(abs(q$x-q$y)/sqrt(2)));				
# 			}
# 			.bogusoptim <- function(lambda, y){
# 				return(.normaldistance(.boxcox(y, lambda)));				
# 			}
# 			.boxcox <- function(y, lambda){
# 				if(lambda==0){return(log(y))}
# 				else{return((y^lambda-1)/lambda)}
# 			}
# 			# boxit <- function(y, tolerance=1){
# 			# 	return(optimise(.bogusoptim, c(-5,5), tol=tolerance, y));				
# 			# }
# 			power_vals <- c(-5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5)
# 			power_opt <- c()
# 			for (i in 1:length(power_vals)){
# 				power_opt <- c(power_opt, as.numeric(.bogusoptim(power_vals[i], x)))
# 			}
# 			goodness <- min(power_opt)
# 			power <- power_vals[which.min(power_opt)]
# 			norm_data <- .boxcox(x, power)
# 			return(norm_data)
# 			}
# 		,inverse_normal = qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
# 		,square=x^2
# 	)})
# }