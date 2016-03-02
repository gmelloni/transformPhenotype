
# Implementation of a special for a "not in"
'%notin%' <- function(x, table) 
				match(x, table, nomatch = 0L) == 0L

# Implementation of two specials for comparison of vectors with NA
# These functions differ from "==" and "!=" in the sense that they can also compare NA
# NA==NA is TRUE and "foo"==NA is FALSE (generally they both return NA)
"%eq%" <- function(x , y) {
	out <- x==y
	nacomp <- is.na(x)==is.na(y)
	out[is.na(out)] <- nacomp[is.na(out)]
	return(out)
}
"%noteq%" <- Negate("%eq%")

# Implementation of a '+' for strings like in python
"%+%" <- function(a, b) paste0(a, b)
"%++%" <- function(a, b) paste(a, b)

#change all factor columns to character
.changeFactor <- function(df) {
  for( i in colnames(df)){
    if(is.factor(df[,i]))
      df[ , i] <- as.character(df[,i])
  }
  return(df)
}

# divide a vector in chunks of the specified length
.subsetter <- function(x , len) {
  subs <- seq(1,length(x),len)
  xBlocks <- lapply(2:length(subs) , function(i) {
                    start <- subs[i-1]
                    end <- subs[i]-1
                    return(x[start:end])
                    })
  xBlocks[[length(xBlocks)+1]] <- x[subs[length(subs)]:length(x)]
  return(xBlocks)
}

# fast trimmer of trailing spaces
#Trims both before and after the submitted object
.myTrimmer <- function (object, ...) 
{
   s <- sub("^[\t\n\f\r ]*", "", as.character(object))
   s <- sub("[\t\n\f\r ]*$", "", s)
   s
}

# Merge mantaining the order of the first dataframe
.mergeOrder <- function(df1 , df2 , ...){
  if("INDEXVARTOBEDESTROYED" %in% colnames(df1)| "INDEXVARTOBEDESTROYED" %in% colnames(df2)){
    stop("There should not be any column of df1 with the name INDEXVARTOBEDESTROYED")
  }
  df1$INDEXVARTOBEDESTROYED <- 1:nrow(df1)
  out <- merge(x=df1 , y=df2 , ...)
  out <- out[ order(out$INDEXVARTOBEDESTROYED , na.last=TRUE) , ]
  out$INDEXVARTOBEDESTROYED <- NULL
  return(out)
}

# find all duplicated elements in a vector
.superdup <- function(x) {
  duplicated(x) | duplicated(x , fromLast=TRUE)
}

# Substitute all NA of a vector with the specified value, default is ""
.noNA <- function(x , subs="") {
  ifelse(is.na(x) , subs , x)
}

# set a variable to numeric if possible, character otherwise
# .NumCheck <- function(x)
# {
#     x <- as.character(x)
#    Result <- suppressWarnings(as.numeric(x))
#    if ( all(is.na(Result) ) ) {
#         return(x)
#    } else {
#         return(Result)
#    }
# }

# You can load library even if they are not installed. 
# In case is not present, this function will check on CRAN first and the Bioconductor
# USAGE = .loadLibrary("stringr")

# .loadLibrary <- function(x) {
#   if(suppressWarnings(require(x , character.only = TRUE , quietly=TRUE))){
#       print( paste(x , "is loaded correctly") )
#   } else {
      
#       tryCatch( {
#           print( paste("Trying to install from CRAN:" , x) )
#           install.packages(x) 
#           }
#           , error=function() {
#             print("The library is not on CRAN")
#             print("Trying to install from Bioconductor")
#             source("http://bioconductor.org/biocLite.R")
#           biocLite(x , suppressUpdates=TRUE , suppressAutoUpdate=TRUE)
#           }
#           , warning=function(w) {
#             if( grepl('not available' , w) ) {
#               print("The library is not on CRAN")
#               print("Trying to install from Bioconductor")
#               source("http://bioconductor.org/biocLite.R")
#             biocLite(x , suppressUpdates=TRUE , suppressAutoUpdate=TRUE)
#             } else {
#             print(w)
#             }}
#           )
#       if(require(x , character.only = TRUE , quietly=TRUE)){
#           print( paste(x , "installed and loaded") )
#       } else {
#           stop( paste("Could not install" , x) )
#       }
#   }
# }

# MFROW autosetter
.mfrow <- function(nplots) {
  if(nplots <= 3) {
    return(c(1, nplots))
  }
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }
  sqrt_nplots <- sqrt(nplots)
  if(is.wholenumber(sqrt_nplots)){
    # par(mfrow=c(sqrt_nplots , sqrt_nplots))
    return(c(sqrt_nplots , sqrt_nplots))
  } else {
    if(round((sqrt_nplots))^2 > nplots){
      # par(mfrow=c(round(sqrt_nplots) , round(sqrt_nplots)))
      return(c(round(sqrt_nplots) , round(sqrt_nplots)))
    }else{
      # par(mfrow=c(round(sqrt_nplots) , round(sqrt_nplots) + 1))
      return(c(round(sqrt_nplots) , round(sqrt_nplots) + 1))
    }
  }
}

.mymapvalues <- function (x, from, to, warn_missing = FALSE) 
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