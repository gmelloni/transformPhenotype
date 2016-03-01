
library(shiny)
library(magrittr)
library(nortest)
library(forecast)
library(DT)
options(DT.options = list(pageLength = 20))

# Internal functions of the app
source(file.path("data" , "custom_functions" , "transformFunctions.R"))

shinyServer(function(input, output,session) {
#Put The Logo
output$myImage <- renderImage({
    width  <- session$clientData$output_myImage_width
    height <- as.character(as.numeric(session$clientData$output_myImage_height)/2)
    # For high-res displays, this will be greater than 1
    pixelratio <- session$clientData$pixelratio
    # filename <- normalizePath(file.path('data', 'Lognormal_distribution.png'))
    filename <- normalizePath(file.path('data', 'logo.png'))
    list(src = filename, alt="Coolest Logo Ever!", width=width , height=height)
  }, deleteFile = FALSE)

# Reactive object with the full dataset
# This object is activated when you upload some data
rawData <- reactive({
  # If you click on Load example dataset, it will upload an example
  if(!is.null(input$myfile)){
    inFile <- input$myfile
  } else if(as.logical(input$example)){
    inFile <- list(name="exampleDataset.txt" 
      , datapath=file.path("data" , "exampleDataset.txt"))
  } else {
    return(NULL)
  }
  out <- read.table(inFile$datapath
            , header=TRUE
            , sep="\t"
            , fill=TRUE
            , strip.white=TRUE
            , as.is=TRUE
            )
  colnames(out) <- tolower(colnames(out))
  if("sample_id" %in% colnames(out)){
    out$sample_id <- as.character(out$sample_id)
    out <- out[ , c("sample_id" , colnames(out)[colnames(out)!="sample_id"])]
  } else if( length(unique(out[ , 1]))==nrow(out) ){
    out[ , 1] <- as.character(out[ , 1])
    colnames(out)[1] <- "sample_id"
  } else {
    stop("FORMAT OF THE FILE NOT SUPPORTED. NEED A SAMPLE_ID COLUMN WITH UNIQUE IDENTIFIERS")
  }
  return(out)
})
# reactive column of the chosen trait
mytrait <- reactive({
  rawData()[[input$trait]]
})
# Show the uploaded table
output$contents <- DT::renderDataTable({
  rawData()
} ,rownames=FALSE
  ,extensions=list(FixedColumns=list(leftColumns = 1))
  ,options = list(searchHighlight = TRUE
    					,scrollX = TRUE
    					,scrollCollapse = FALSE
              ,pageLength = 10)
  ,filter = 'bottom'
)
#Change Trait selector according to the colnames of user table
observe({
  updateSelectInput(session , "trait" , "Choose Your Trait:"
      , choices= if(is.null(rawData())){
             				""
           			  }else{
           				 notTraits2 <- c(notTraits , colnames(rawData())[sapply(rawData() , class)=="character"])
                   colnames(rawData())[!colnames(rawData()) %in% notTraits2]
           			  }
      , selected=if(is.null(rawData())){
                    ""
                  }else{
                   notTraits2 <- c(notTraits , colnames(rawData())[sapply(rawData() , class)=="character"])
                   colnames(rawData())[!colnames(rawData()) %in% notTraits2][1]
                  }
  )
})
# Modify filter selector, based on min and max value of the trait chosen
observe({
  if(input$trait!=""){
    filtermin <- min(mytrait() , na.rm=TRUE)
    filtermax <- max(mytrait() , na.rm=TRUE)
    updateSliderInput(session , "filters" , "Additional filter on raw values:"
                ,min = filtermin
                ,max = filtermax
                ,value = c(filtermin , filtermax)
  )
  }
})
# Modify Covariates selector adding all variables except the one under analysis
# If age is among the variable, age2 is added too
observe({
  if(!is.null(rawData())){
    newchoices <- colnames(rawData()) %>% .[.!="sample_id"]
    newchoices <- c(NA , newchoices)
    if("age" %in% newchoices)
      newchoices <- c(newchoices , "age2") %>% sort
    newchoices <- newchoices[ newchoices!=input$trait ]
    updateSelectInput(session , "covariates_tested" , "Choose one or more covariate:"
        , choices= newchoices
        , selected=if("sex" %in% newchoices) "sex" else NA
    )
  }
})

# Generate protocol file from user choices like filters, transformation type etc.
protocolFile <- reactive({
  # filter must be treated separately because of the original format of protocol file
  # The slider return a vector with a sort range c(10,20) that becomes <10,>20
if(input$trait!=""){
  myfilter <- input$filters %>%
              paste(c("<" , ">") , . , sep="") %>%
              paste(. , collapse=",")
  covariates <- if(any(input$covariates_tested %in% c("NA","") )) {
                  NA
                } else {
                  paste(input$covariates_tested , collapse=",")
                }
  data.frame(trait=as.factor(input$trait)
            ,units=input$units
            ,transformation_method=input$transformation_method
            ,covariates_tested=covariates
            ,sd_num=if(input$sd_num=="NA") NA else as.numeric(input$sd_num)
            ,sd_dir=input$sd_dir
            ,filters=myfilter
            ,stringsAsFactors=FALSE
            )
} else {
  return(NULL)
}
})
# Reduce dataset according to protocol specification
  # remove missing sex specification
  # remove absolute outliers
  # remove tails according to SD specifications
    # note: the SD filter is influenced by stratification by sex
traitObject <- reactive({
  if(is.null(rawData())){
      return(NULL)
  }
  # Variable from protocol
  altFilter <- protocolFile()$filters
  currTrait <- protocolFile()$trait
  numSDs <- protocolFile()$sd_num
  sdDir <- protocolFile()$sd_dir
  covariates <- as.character(protocolFile()$covariates_tested)
  covariatesSplit <- unlist(strsplit(covariates , ","))
  covList<- list(covariates=covariatesSplit 
            , age2Flag=if("age2" %in% covariatesSplit) TRUE else NA)
  # Remove missing sex
  missingSex <- which(is.na(rawData()$sex))
  dataset <- rawData()
  if(length(missingSex>0)){
    dataset <- dataset[-missingSex,]
  }
  traitObject <- createDF(currTrait,covList,dataset)
  # Apply hard filter for oulier
  if(!is.na(altFilter)){
    outliersFilter <- applyFilters(altFilter,traitObject)
    excl <- which(traitObject$ID %in% outliersFilter)
    if(length(excl)>0)traitObject <- traitObject[-excl,]
  }
  females<-which(traitObject$sex==2)
  males<-which(traitObject$sex==1)
  if(input$sexStratFlag=="Yes"){
    if(!is.na(numSDs)){
      maleSdOutliers <-findSdOutliers(numSDs,sdDir,traitObject[males,])
      femaleSdOutliers <-findSdOutliers(numSDs,sdDir,traitObject[females,])
      sdOutliers <- c(maleSdOutliers,femaleSdOutliers)
      excl<-which(traitObject$ID %in% sdOutliers)
      if(length(excl)>0){
        traitObjforRawPlot <- traitObject[-excl,]
      } else {
        traitObjforRawPlot <- traitObject  
      }
    } else {
      traitObjforRawPlot <- traitObject
    }
  } else {
    if(!is.na(numSDs)){
      sdOutliers <- findSdOutliers(numSDs,sdDir,traitObject)
      excl<-which(traitObject$ID %in% sdOutliers)
      if(length(excl)>0){
        traitObjforRawPlot <- traitObject[-excl,]
      } else {
        traitObjforRawPlot <- traitObject  
      }
    } else {
      traitObjforRawPlot <- traitObject
    }
  }
  return(traitObjforRawPlot)
})

# Display protocol file
output$protocolFile <- renderTable({
  if(is.null(protocolFile()))
    return(NULL)
  protocolFile()
  } , row.names=FALSE)

# Follow all attempts of the user
values <- reactiveValues(df_data=NULL)
observeEvent(input$storeattempt , {
  temp <- rbind(values$df_data , protocolFile())
  temp$transformation_method <- mymapvalues(temp$transformation_method 
                                          , from=names(angela_transformations) 
                                          , to=unname(angela_transformations))
  values$df_data <- unique(temp)
  })
output$cumulativeProtocolFile <- DT::renderDataTable({
  if(is.null(protocolFile()))
    return(NULL)
  else
    return(values$df_data)
  })
observeEvent(input$deleteattempt , {
  if(!is.null(input$cumulativeProtocolFile_rows_selected)){
    temp <- values$df_data[ -input$cumulativeProtocolFile_rows_selected , ]
    values$df_data <- temp
  }
})



# Plot density curves of trait between males and females
output$sexplot <- renderPlot({
  if(is.null(rawData())){
    return({
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, "No data yet"
        ,cex = 3, col = "black")
    })
  }
  # Extract covariate field: if it is NA, abort the plot
  covariates <- as.character(protocolFile()$covariates_tested)
  if(is.na(covariates[1])){
    return({
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, "Select sex in your covariates\nif you want to see this plot"
        ,cex = 3, col = "black")
    })
  }
  # Now that we now that is not NA, let's split it and check if sex is among covariates
  covariatesSplit <- unlist(strsplit(covariates , ","))
  if(!"sex" %in% covariatesSplit){
    return({
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, "Select sex in your covariates\nif you want to see this plot"
        ,cex = 3, col = "black")
    }) 
  }
  currTrait <- protocolFile()$trait
  traitUnit <- protocolFile()$units
  # transformMethod <- protocolFile()$transformation_method
  # if(length(bmiCheck)>0) currTrait <- paste(currTrait,"BMIadj",sep="")
  traitLabelFull <- paste(currTrait," (",traitUnit,")",sep="")
  # Start checking between sexes    
  females<-which(traitObject()$sex==2)
  males<-which(traitObject()$sex==1)

  # 6. Run Wilcoxon test to see if sexes differ and print out density plot
  # sexPval<-sexStratification(traitObject,males,females,traitLabel,traitLabelFull)
  project <- as.character(currTrait)
  sexStratificationPlot(X=traitObject(),m=males,f=females
                      ,label=currTrait,labelFull=traitLabelFull
                      ,project=project)
})

# Plot data according to normalization criterion
output$transformedplot <- renderPlot({
  if(is.null(rawData())){
    return({
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1)
          , ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, "No data yet"
        ,cex = 3, col = "black")
    })
  }
  # Protocol Variables
  currTrait <- protocolFile()$trait
  traitUnit <- protocolFile()$units
  traitLabelFull <- paste(currTrait," (",traitUnit,")",sep="")
  transformMethod <- protocolFile()$transformation_method
  # If radiobutton stratified by sex is on Yes, the plot is run twice for m and f
    return({
      if(input$sexStratFlag=="Yes"){
        females<-which(traitObject()$sex==2)
        males<-which(traitObject()$sex==1)
        loop <- c("Males" , "Females")
        par(mfrow=c(4,2))
      } else {
        loop <- "No Stratification"
        par(mfrow=c(2,2))
      }
      for(i in loop){
        subs <- if(i=="Males") {
                  males 
                } else if(i=="Females"){
                  females
                } else {
                  1:length(traitObject()$trait)
                }
        x <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=transformMethod)$norm_data
        if(length(x)>5000){
          x <- x[1:5000]
        }
        plot( x=1:length(x)
            , y=x
            , xlab="Index"
            , ylab=traitLabelFull
            , main=paste(i , "ScatterPlot" , paste("N =" , length(x)) , sep="\n")
            , col=if(i=="Males") "navy" else if(i=="Females") "red" else "black")
        boxplot( x=x
                , main=paste(i , "Boxplot" , sep="\n")
                , col="cyan")
        mymin=round(min(x),2)
        mymax=round(max(x),2)
        mymean=round(mean(x),2)
        mysd=round(sqrt(var(x)),2)
        mymain=sprintf("min=%s;max=%s;mean=%s;sd=%s",mymin,mymax,mymean,mysd)
        hist( x=x
            , main=paste(i , mymain , sep="\n")
            , xlab=traitLabelFull
            , prob=TRUE)
        m <- mean(x, na.rm=TRUE)
        std <- sd(x,na.rm=TRUE) 
        curve(dnorm(x,mean=m,sd=std)
            , add=T
            , col="forestgreen"
            , lwd=4)
        #pv<-as.numeric(unlist(shapiro.test(x))[2])
        pv<-as.numeric(unlist(ad.test(x))[2])
        mymain=sprintf("Normal Q-Q plot (Anderson-Darling pval=%s)"
          ,format(pv,digits=3,sci = T))
        qqnorm(x
              ,main=paste(i , mymain , sep="\n"))
        qqline(x 
              ,col=if(i=="Males") "navy" else if(i=="Females") "red" else "black"
              , lwd=4)
      }
    })
})

# The output of every normalization function is a list
# This list contains norm_data as the first element
# The other elements are optional
# In case they are present they will be displayed as a html box at the bottom
# An example is box_cox that report the optimal cut off
output$normalizationSideEffect <- renderPrint({
  if(is.null(traitObject())){
    return(NULL)
  }
  transformMethod <- protocolFile()$transformation_method
  outputList <- list(paste(transformMethod , "further Normalization Info:"))
  if(input$sexStratFlag=="Yes"){
    females<-which(traitObject()$sex==2)
    males<-which(traitObject()$sex==1)
    loop <- c("Males" , "Females")
  } else {
    loop <- "No Stratification"
  }
  for(i in loop){
    subs <- if(i=="Males") {
                males 
              } else if(i=="Females"){
                females
              } else {
                1:length(traitObject()$trait)
              }
    x <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=transformMethod)
    if(length(x)>1){
      outputList <- c(outputList , list(i) , list(x[-1]))
    } else {
      outputList <- c(outputList , list(i) , list("No other Info"))
    }
  }
  return(outputList)
})

# Apply all normalization and calculate normality test
# The output is a table with the pvalue for every transformation
# Yellow color is a NON significant pvalue (that's what we are interested in)
output$normalTable <- DT::renderDataTable({
  if(is.null(rawData())){
    return(NULL)
  }
  if(input$sexStratFlag=="Yes"){
    Transformation <- names(normalizationFunctionsList)
    females<-which(traitObject()$sex==2)
    males<-which(traitObject()$sex==1)
    loop <- c("Males" , "Females")
    normalTable <- lapply(loop , function(sex) {
      lapply(Transformation , function(trans) {
          subs <- if(sex=="Males") males else if(sex=="Females") females else stop("Error in loop")
          x <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=trans)$norm_data
          ad <- tryCatch(ad.test(x)$p.value , error=function(e) e)
          sw <- tryCatch(shapiro.test(x)$p.value, error=function(e) e)
          # ks is annoying with warnings in case of ties, so I suppress them
          ks <- suppressWarnings({
                  tryCatch(ks.test(x,rnorm(50))$p.value, error=function(e) e)
                })
          return(c(sex , trans , ad , sw , ks))
        })
      })
    normalTable <- lapply(normalTable , function(x) do.call("rbind" , x))
    normalTable <- do.call("rbind" , normalTable)
    colnames(normalTable) <- c("Sex" 
                              , "Transformation" 
                              , "Anderson-Darling Test" 
                              , "Shapiro-Wilks Test" 
                              , "Kolmogorov-Smirnov Test")
  } else if(input$sexStratFlag=="No") {
    Transformation <- names(normalizationFunctionsList)
    normalTable <- lapply(Transformation , function(trans) {
          x <- normalizeTraitData(trait=traitObject()$trait , tm=trans)$norm_data
          ad <- tryCatch(ad.test(x)$p.value , error=function(e) e)
          sw <- tryCatch(shapiro.test(x)$p.value, error=function(e) e)
          # ks is annoying with warnings in case of ties, so I suppress them
          ks <- suppressWarnings({
                  tryCatch(ks.test(x,rnorm(50))$p.value, error=function(e) e)
                })
          return(c(trans , ad , sw , ks))
        })
    normalTable <- do.call("rbind" , normalTable)
    colnames(normalTable) <- c("Transformation" 
                              , "Anderson-Darling Test" 
                              , "Shapiro-Wilks Test" 
                              , "Kolmogorov-Smirnov Test")
  }
  return({
      DT::datatable(normalTable) %>% 
            formatStyle("Anderson-Darling Test" 
                , backgroundColor = styleInterval(0.05, c('lightgray', 'yellow'))) %>%
            formatStyle("Shapiro-Wilks Test" 
                , backgroundColor = styleInterval(0.05, c('lightgray', 'yellow'))) %>%
            formatStyle("Kolmogorov-Smirnov Test" 
                , backgroundColor = styleInterval(0.05, c('lightgray', 'yellow')))
            })
} 
)

# JUST FOR CONTROL, GOTTA REMOVE IT IN RELEASE VERSION
# output$traitObject2 <- renderTable({
#   head(traitObject())
#   })

# Linear Model on covariates results
output$linearCovariates <- renderPrint({
  if(is.null(rawData())){
    return(NULL)
  }
  covariates <- protocolFile()$covariates_tested
  covariates <- if(is.na(covariates) | covariates=="") NA else covariates
  if(is.na(covariates)){
    return(NULL)
  }
  currTrait <- protocolFile()$trait
  traitUnit <- protocolFile()$units
  traitLabelFull <- paste(currTrait," (",traitUnit,")",sep="")
  transformMethod <- protocolFile()$transformation_method
  cvr <- unlist(strsplit(covariates , ","))
  covList<- list(covariates=cvr
              , age2Flag=if("age2" %in% cvr) TRUE else NA)
  if(input$sexStratFlag=="Yes"){
    females<-which(traitObject()$sex==2)
    males<-which(traitObject()$sex==1)
    loop <- c("Males" , "Females")
    normData <- lapply(list("Males"=males , "Females"=females) , function(subs) {
        normData <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=transformMethod)$norm_data
        signifCovs<-checkCovariates2(covs=covList
                                ,x=traitObject()[subs , ]
                                ,normx=normData)
      })
    return(invisible(sapply(names(normData) , function(x) {
        return(cat(toupper(x) , show(normData[[x]]) , sep="\n"))
        }))
    )
  } else {
    # subs <- 1:nrow(traitObject())
    normData <- normalizeTraitData(trait=traitObject()$trait , tm=transformMethod)$norm_data
    signifCovs<-checkCovariates2(covs=covList
                                ,x=traitObject()
                                ,normx=normData)
    return(show(signifCovs))
  }
})

# This object is the residual table
# It is updated every time the residualPlot changes
residual <- reactiveValues(df_res=NULL)

# Normalize data and check for residual of the covariates
output$residualsPlot <- renderPlot({
  if(is.null(rawData())){
    return({
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, "No data yet"
        ,cex = 3, col = "black")
    })
  }
  covariates <- protocolFile()$covariates_tested
  covariates <- if(is.na(covariates) | covariates=="") NA else covariates
  if(is.na(covariates)){
    return({
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, "No covariates selected"
        ,cex = 3, col = "black")
    })
  }
  currTrait <- protocolFile()$trait
  traitUnit <- protocolFile()$units
  traitLabelFull <- paste(currTrait," (",traitUnit,")",sep="")
  transformMethod <- protocolFile()$transformation_method
  if(input$sexStratFlag=="Yes"){
    cvr <- unlist(strsplit(covariates , ",")) %>% .[.!="sex"]
    if(length(cvr)==0){
      return({
        par(mar = c(0,0,0,0))
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, "If you stratify by sex, you can't have sex...\n\n... as the only covariate :)"
          ,cex = 3, col = "black")
      })
    }
    covList<- list(covariates=cvr
             , age2Flag=if("age2" %in% cvr) TRUE else NA)
    females<-which(traitObject()$sex==2)
    males<-which(traitObject()$sex==1)
    loop <- list("Males"=males , "Females"=females)
    # loop <- c("Males" , "Females")
    normDataListForPlot <- lapply(loop , function(subs) {
      normData <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=transformMethod)$norm_data
      signifCovs<-checkCovariates(covs=covList
                                  ,x=traitObject()[subs,]
                                  ,normx=normData
                                  ,sxf=sexStratFlag)
      if(!is.na(signifCovs[1])){
        resList<-applySignifCovariates(signifCovs,traitObject()[subs,],normData)
        normResiduals<-resList$residuals
        covString <- resList$covStr
        retData <- checkResiduals(normResiduals,traitObject()[subs,],"")
        outputData <- retData$outputData
        normResiduals <- retData$normResiduals
      } else {
        normResiduals <- normData
        covString <- NULL
        retData <- checkResiduals(normResiduals,traitObject()[subs,],"")
        outputData <- retData$outputData
        normResiduals <- retData$normResiduals
      }
      return(list(covString=covString , normResiduals=normResiduals , normData=normData , outputData=outputData))
      })
    # Add a side effect for this function, by reporting the actual residual table
    outputData <- do.call("rbind" , lapply(normDataListForPlot , '[[' , "outputData") )
    observe({
      residual$df_res <- outputData
    })
    return(plotResidualDataBySex(
                    tl=currTrait
                    ,tlf=traitLabelFull
                    ,myList=normDataListForPlot
                    ))
  } else {
    subs <- 1:nrow(traitObject())
    normData <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=transformMethod)$norm_data
    cvr <- unlist(strsplit(covariates , ","))
    covList<- list(covariates=cvr
              , age2Flag=if("age2" %in% cvr) TRUE else NA)
    signifCovs<-checkCovariates(covs=covList
                                ,x=traitObject()
                                ,normx=normData
                                ,sxf=sexStratFlag)
    if(!is.na(signifCovs[1])){
      resList<-applySignifCovariates(signifCovs,traitObject(),normData)
      normResiduals<-resList$residuals
      covString <- resList$covStr
      retData <- checkResiduals(normResiduals,traitObject(),"")
      outputData <- retData$outputData
      normResiduals <- retData$normResiduals
    } else {
      normResiduals <- normData
      covString <- NULL
      retData <- checkResiduals(normResiduals,traitObject(),"")
      outputData <- retData$outputData
      normResiduals <- retData$normResiduals
    }
    # Side effect, update residual final table
    observe({
      residual$df_res <- outputData
    })
    return(plotResidualData(tl=currTrait
                    ,tlf=traitLabelFull
                    # ,"No Stratification"
                    ,cvs=covString
                    ,res=normResiduals
                    ,normx=normData))
  }
})

# This observer will wait for the residual table to be present
# When this event is observed, the download button is enabled
observeEvent(!is.null(residual$df_res) , {
    # notify the browser that the residual table is ready to be downloaded
    session$sendCustomMessage("download_ready_res" , list(mex=""))
})

# Download residuals table
output$downloadResidual <- downloadHandler(
    filename=function() {
      timeTag <- Sys.time() %>% 
            sub(" GMT$" , "" , .) %>% 
            gsub(":" , "_" , .) %>%
            gsub("-" , "" , .) %>%
            sub(" " , "." , .)
      paste(protocolFile()$trait , "residuals" , timeTag , "stand_residuals.txt" , sep=".")
    }
    , content=function(file) {
        write.table(residual$df_res
                , file=file
                , sep = "\t"
                , col.names = TRUE
                , row.names = FALSE
                , quote=FALSE)
})

# This observer will wait for the store button to be pressed at least once
# When this event is observed, the download button is enabled
observeEvent(input$storeattempt , {
    # notify the browser that the the protocol is ready to be downloaded
    session$sendCustomMessage("download_ready" , list(mex=""))
})

output$downloadFinalProtocol <- downloadHandler(
    filename=function() {
      timeTag <- Sys.time() %>% 
      			sub(" GMT$" , "" , .) %>% 
      			gsub(":" , "_" , .) %>%
      			gsub("-" , "" , .) %>%
      			sub(" " , "." , .)
      paste("ProtocolFile" , timeTag , "txt" , sep=".")
    }
    , content=function(file) {
        write.table(values$df_data
                , file=file
                , sep = "\t"
                , col.names = TRUE
                , row.names = FALSE
                , quote=FALSE)
})

#Close the R process if the browser is closed
session$onSessionEnded(function() {
      stopApp()
})

})



#---------------#
# CODE CEMETERY #
#---------------#

##### OLD VERSION DEPRECATED #####
# NOTE: this version is a two step filter like the original code
# In the original code SD filter is used after sex stratification. Now it is all in 1 step
# Filtered Data Reactive Object
  # This object take the result from traiObject and apply SD based filter
# traitObject2 <- reactive({
#   numSDs <- protocolFile()$sd_num
#   sdDir <- protocolFile()$sd_dir
#   currTrait <- protocolFile()$trait
#   traitUnit <- protocolFile()$units
#   traitLabelFull <- paste(currTrait," (",traitUnit,")",sep="")
#   transformMethod <- protocolFile()$transformation_method
#   females<-which(traitObject()$sex==2)
#   males<-which(traitObject()$sex==1)
#   if(input$sexStratFlag=="Yes"){
#     if(!is.na(numSDs)){
#       maleSdOutliers <-findSdOutliers(numSDs,sdDir,traitObject()[males,])
#       femaleSdOutliers <-findSdOutliers(numSDs,sdDir,traitObject()[females,])
#       sdOutliers <- c(maleSdOutliers,femaleSdOutliers)
#       excl<-which(traitObject()$ID %in% sdOutliers)
#       if(length(excl)>0){
#         traitObjforRawPlot <- traitObject()[-excl,]
#       } else {
#         traitObjforRawPlot <- traitObject()  
#       }
#     } else {
#       traitObjforRawPlot <- traitObject()
#     }
#   } else {
#     if(!is.na(numSDs)){
#       sdOutliers <- findSdOutliers(numSDs,sdDir,traitObject())
#       excl<-which(traitObject()$ID %in% sdOutliers)
#       if(length(excl)>0){
#         traitObjforRawPlot <- traitObject()[-excl,]
#       } else {
#         traitObjforRawPlot <- traitObject()  
#       }
#     } else {
#       traitObjforRawPlot <- traitObject()
#     }
#   }
# })

# rm(list = ls())
# library(shiny)

# Logged = FALSE;
# my_username <- "test"
# my_password <- "test"

# ui1 <- function(){
#   tagList(
#     div(id = "login",
#         wellPanel(textInput("userName", "Username"),
#                   passwordInput("passwd", "Password"),
#                   br(),actionButton("Login", "Log in"))),
#     tags$style(type="text/css", "#login {font-size:10px;   text-align: left;position:absolute;top: 40%;left: 50%;margin-top: -100px;margin-left: -150px;}")
#   )}

# ui2 <- function(){tagList(tabPanel("Test"))}

# ui = (htmlOutput("page"))
# server = (function(input, output,session) {
#   USER <<- reactiveValues(Logged = Logged)
#   observe({ 
#     if (USER$Logged == FALSE) {
#       if (!is.null(input$Login)) {
#         if (input$Login > 0) {
#           Username <- isolate(input$userName)
#           Password <- isolate(input$passwd)
#           Id.username <- which(my_username == Username)
#           Id.password <- which(my_password == Password)
#           if (length(Id.username) > 0 & length(Id.password) > 0) {
#             if (Id.username == Id.password) {
#               USER$Logged <<- TRUE
#             } 
#           }
#         } 
#       }
#     }    
#   })
#   observe({
#     if (USER$Logged == FALSE) {

#       output$page <- renderUI({
#         div(class="outer",do.call(bootstrapPage,c("",ui1())))
#       })
#     }
#     if (USER$Logged == TRUE) 
#     {
#       output$page <- renderUI({
#         div(class="outer",do.call(navbarPage,c(inverse=TRUE,title = "Contratulations you got in!",ui2())))
#       })
#       print(ui)
#     }
#   })
# })
# runApp(list(ui = ui, server = server))
