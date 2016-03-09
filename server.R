
library(shiny)
library(magrittr)
library(nortest)
library(forecast)
library(DT)
options(DT.options = list(pageLength = 20))

# Internal functions of the app
source(file.path("data" , "custom_functions" , "transformFunctions.R"))
# Serie of cool utility functions used from time to time
source(file.path("data" , "custom_functions" , "utilitiesFunctions.R"))

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

#---------------------#
#     LOAD DATASET    #
#---------------------#

# Create a reactive value for error
dferror <- reactiveValues(myerror=c("Welcome to the Transform Phenotype App!"
                                    ,"If you don't have your own dataset or you are not sure about the format"
                                    ,"Click on Load Example Dataset ;)"))

# Reactive object with the full dataset
# This object is activated when you upload some data
rawData <- reactive({
  # If you click on Load example dataset, it will upload an example
  if(!is.null(input$myfile)){
    inFile <- input$myfile
  } else if(as.logical(input$example)){
    inFile <- data.frame(name="exampleDataset.txt" 
      , datapath=file.path("data" , "exampleDataset.txt")
      ,stringsAsFactors=FALSE)
  } else {
    return(NULL)
  }
  if(nrow(inFile)==1){
    out <- read.table(inFile$datapath
              , header=TRUE
              , sep="\t"
              , fill=TRUE
              , strip.white=TRUE
              , as.is=TRUE
              )
  } else {
    checkCols <- sapply(inFile$datapath , function(x) readLines(x , n=1))
    if(any(checkCols!=checkCols[1])){
      dferror$myerror <- "You select multiple files with different column names"
      return(NULL)
    }
    fileNames <- if(any(duplicated(inFile$name))){
                  paste("Cohort" , 1:nrow(inFile) , sep=".")
                 } else {
                  inFile$name
                }
    out <- lapply(1:length(inFile$datapath) , function(x) {
            outInternal <- read.table(inFile$datapath[[x]]
              , header=TRUE
              , sep="\t"
              , fill=TRUE
              , strip.white=TRUE
              , as.is=TRUE
              )
            if(nrow(outInternal)==0)
              return(NULL)
            if("InputFileNum" %notin% colnames(outInternal)){
              outInternal$InputFileNum <- fileNames[x]
            } else {
              return(NULL)
            }
            return(outInternal)
            }) %>% do.call("rbind" , .)
    if(is.null(out)){
      dferror$myerror <- "Either all files are empty or InputFileNum is among the column names"
      return(NULL)
    }
  }
  colnames(out) <- tolower(colnames(out))
  rawRows <- nrow(out)
  out <- unique(out)
  if(nrow(out)!=rawRows)
    warning("FOUND DUPLICATED ROWS. REMOVED")
  #### Check for Sanger Phenotype Database format
    # There is currently a typo in the database (phentoype instead of phenotype)
    # We accept both forms
  if( all(c("genotype" , "phenotype") %in% colnames(out)) | all(c("genotype" , "phentoype") %in% colnames(out)) ){
    out$sample_id <- out$genotype
    out$phenotype <- NULL
    out$phentoype <- NULL
    out$genotype <- NULL
  }
  # If sample_id is among the colnames it is going to be used as sample identifier
  # If sample_id doesn't exist, shiny tries with the first column, otherwise dies
  if("sample_id" %in% colnames(out)){
    out$sample_id <- as.character(out$sample_id)
    if( length(unique(out$sample_id))==nrow(out) ){
      out <- out[ , c("sample_id" , colnames(out)[colnames(out)!="sample_id"])]
    } else {
      dferror$myerror <- "sample_id COLUMN HAS DUPLICATED VALUES"
      return(NULL)
    }
  } else if( length(unique(out[ , 1]))==nrow(out) ){
    out[ , 1] <- as.character(out[ , 1])
    colnames(out)[1] <- "sample_id"
  } else {
    dferror$myerror <- "FORMAT OF THE FILE NOT SUPPORTED. NEED A SAMPLE_ID COLUMN WITH UNIQUE IDENTIFIERS"
    return(NULL)
  }
  ##### Check sex
    # At least sex or gender must be present
  if( all(c("sex" , "gender") %notin% colnames(out)) ){
    dferror$myerror <- "SEX OR GENDER REQUIRED"
    return(NULL)
  }
  if( all(c("sex" , "gender") %in% colnames(out)) ){
    if( all(out$sex %eq% out$gender) ){
      out$gender <- NULL
      warning("FOUND BOTH SEX AND GENDER COLUMNS. KEPT SEX")
    } else {
      dferror$myerror <- "FOUND BOTH SEX AND GENDER COLUMNS AND THEY ARE DIFFERENT. GET RID OF ONE OF THEM"
      return(NULL)
    }
  } else if("gender" %in% colnames(out)) {
    warning("GENDER COLUMN CHANGED TO SEX")
    out$sex <- out$gender
    out$gender <- NULL
  }
    # We can accept sex coded as 0/1 or M/F or m/f
    # Our standard will be 1/2 (plink style)
  out$sex[out$sex==""] <- NA
  if( all(na.omit(out$sex) %in% c(0,1)) ){
    out$sex <- out$sex + 1
  } else if( all( na.omit(tolower(out$sex)) %in% c("m","f")) ){
    out$sex <- .mymapvalues(tolower(out$sex) , from=c("m" , "f") , to=c(1,2)) %>%
                as.integer
  } else if( any( na.omit(out$sex) %notin% c(1,2) ) ){
    dferror$myerror <- "UNRECOGNIZED OR MIXED SEX CODIFICATION. REVERT TO 1==MALE , 2==FEMALE"
    return(NULL)
  }
  out$sex <- suppressWarnings( as.numeric(out$sex) )
  if("age" %in% colnames(out))
    out$age <- suppressWarnings( as.numeric(out$age) )
  numTraits <- sapply(out[ , colnames(out) %notin% c("age" , "sex")] , class) %in% c("numeric" , "integer") %>%
                which %>%
                length
  # numStrats <- (length(which(sapply(out , class)=="character")) - 1)
  rule1 <- sapply(out , class) %in% c("character" , "integer")
  rule2 <- sapply(out , function(x) any(duplicated(x)))
  numStrats <- length(which(rule1 & rule2))
  dferror$myerror <- c(
    "DATA LOADED CORRECTLY:"
    ,"Numeric Traits (excluding sex and age):" %++% numTraits
    ,"Sex detected:" %++% c("No" , "Yes")[as.numeric(any(colnames(out)=="sex"))+1]
    ,"Age detected:" %++% c("No" , "Yes")[as.numeric(any(colnames(out)=="age"))+1]
    ,"BMI detected:" %++% c("No" , "Yes")[as.numeric(any(colnames(out)=="bmi"))+1]
    ,"Stratification Variables:" %++% numStrats
                          )
  return(out)
})

# What was wrong with my data?
# Error message about format is stored together with rawData() reactive object
output$dferror <- renderPrint(cat(dferror$myerror , sep="\n"))

# Reactive column of the chosen trait
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

#---------------------------------#
#       SELECTOR MODIFIERS        #
#---------------------------------#

#Change Trait selector according to the colnames of user table
observe({
  updateSelectInput(session , "trait" , "Choose Your Trait:"
      , choices= if(is.null(rawData())){
             				""
           			  }else{
           				 notTraits2 <- c(notTraits 
                                , colnames(rawData())[sapply(rawData() , class)=="character"])
                   colnames(rawData())[!colnames(rawData()) %in% notTraits2]
           			  }
      , selected=if(is.null(rawData())){
                    ""
                  }else{
                   notTraits2 <- c(notTraits 
                                , colnames(rawData())[sapply(rawData() , class)=="character"])
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
    # covariates must be numerical or integer
    newchoices <- colnames(rawData())[sapply(rawData() , class)!="character"]
    # newchoices <- c(NA , newchoices)
    if("age" %in% newchoices)
      newchoices <- c(newchoices , "age2") %>% sort
    newchoices <- newchoices[ newchoices!=input$trait ]
    updateSelectInput(session , "covariates_tested" , "Choose one or more covariate:"
        , choices= newchoices
        , selected=if("sex" %in% newchoices) "sex" else NA
    )
  }
})
# Adjust the stratifier
# Stratifiers are all the character columns, except for sample_id column
observe({
  if(!is.null(rawData())){
    # Rules to be a stratifier:
      # being a vector of type character or integer
      # contain at least a duplicated value (otherwise the splitting makes no sense)
      # sex cannot be a stratifier (it is treated separately)
    rule1 <- sapply(rawData() , class) %in% c("character" , "integer")
    rule2 <- sapply(rawData() , function(x) any(duplicated(x)))
    newstrats <- colnames(rawData())[rule1 & rule2]
    stratifiers <- newstrats
    stratifiers <- stratifiers[stratifiers %notin% c("sample_id" , "sex")]
    updateSelectInput(session , "stratifier" , "Choose one stratification variable:"
        , choices=stratifiers
        , selected=NULL
    )
  }
})

#---------------------#
# PROTOCOL GENERATION #
#---------------------#

# Generate protocol file from user choices like filters, transformation type etc.
protocolFile <- reactive({
  # filter must be treated separately because of the original format of protocol file
  # The slider return a vector with a sort range c(10,20) that becomes <10,>20
if(input$trait!=""){
  myfilter <- input$filters %>%
              paste(c("<" , ">") , . , sep="") %>%
              paste(. , collapse=",")
  covariates <- if(is.null(input$covariates_tested)) {
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

#------------------------------------------#
#   CREATING TRAITOBJECT USING PROTOCOL    #
#------------------------------------------#

# Reduce dataset according to protocol specification
  # remove missing sex specification
  # remove absolute outliers
  # remove tails according to SD specifications
    # note: the SD filter is influenced by stratification by sex
# If input$stratifier is not set or is NA, the result is a dataframe
# If input$stratifier is set the result is a list of dataframe
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
  if(is.null(input$stratifier)){
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
  } else {
    datasetsplit <- split(dataset , as.list(dataset[ , input$stratifier , drop=FALSE]))
    traitObjectSplit <- lapply(datasetsplit , function(dataset){
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
    names(traitObjectSplit) <- names(datasetsplit)
    return(traitObjectSplit)
  }
})

#-------------------------------#
#   PROTOCOL FILE ADJUSTMENT    #
#-------------------------------#

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
  temp$transformation_method <- .mymapvalues(temp$transformation_method 
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

#--------------------------#
#   SEX DIFFERENCE PLOT    #
#--------------------------#

# This closure is a generator of empty plots with a message
emptyPlotter <-function(message){
  renderPlot({
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, message
        ,cex = 3, col = "black")
  })
}

# Plot density curves of trait between males and females
# output$sexplot <- renderPlot({
output$sexplot <- renderUI({
  if(is.null(rawData())){
    output$sexplotempty <- emptyPlotter("No Data Yet")
    return(
    plotOutput("sexplotempty" , height="800px")
    )
  }
  # Global parameters, indipendent from traitObject
  currTrait <- protocolFile()$trait
  traitUnit <- protocolFile()$units
  traitLabelFull <- paste(currTrait," (",traitUnit,")",sep="")
  # Start checking between sexes    
  project <- as.character(currTrait)
  # Extract covariate field: if it is NA, abort the plot
  covariates <- as.character(protocolFile()$covariates_tested)
  covariatesSplit <- unlist(strsplit(covariates , ","))
  if(!"sex" %in% covariatesSplit){
    output$sexplotempty <- emptyPlotter("Select sex in your covariates\nif you want to see this plot")
    return(
      plotOutput("sexplotempty" , height="800px")
    )
  }
  if(is.null(input$stratifier)){
    # Check if sex is among covariates
    # Run Wilcoxon test to see if sexes differ and print out density plot
    females<-which(traitObject()$sex==2)
    males<-which(traitObject()$sex==1)
    output$monosexplot <- renderPlot(sexStratificationPlot(X=traitObject(),m=males,f=females
                                    ,label=currTrait,labelFull=traitLabelFull
                                    ,project=project))
    return(plotOutput("monosexplot" , height="800px"))
  } else {
    tabs <- lapply(names(traitObject()) , function(strat){
              females<-which(traitObject()[[strat]]$sex==2)
              males<-which(traitObject()[[strat]]$sex==1)
              output[[paste0(which(strat==names(traitObject())) , "sexstrat_")]] <- renderPlot({
                    sexStratificationPlot(X=traitObject()[[strat]],m=males,f=females
                                      ,label=currTrait,labelFull=traitLabelFull
                                      ,project=project)
                  })
              return(tabPanel(strat , plotOutput(paste0(which(strat==names(traitObject())) , "sexstrat_") , height="800px")))
            })
    return(do.call(tabsetPanel , tabs))
  }
})

# mytransform <- reactiveValues(transformer = NULL)
##### render of the tab Data Plot tab
output$transformer <- renderUI({
  if(is.null(rawData())){
    output$transformerempty <- emptyPlotter("No Data Yet")
    return(
    fludifPage(plotOutput("transformerempty" , height="800px"))
    )
  }
  # Protocol Variables
  currTrait <- protocolFile()$trait
  traitUnit <- protocolFile()$units
  traitLabelFull <- paste(currTrait," (",traitUnit,")",sep="")
  transformMethod <- protocolFile()$transformation_method
  if(is.null(input$stratifier)){
    # Initialize transformation side effect output
    outputList <- list(paste(transformMethod , "further Normalization Info:"))
    if(input$sexStratFlag=="Yes"){
      females<-which(traitObject()$sex==2)
      males<-which(traitObject()$sex==1)
      loop <- c("Males" , "Females")
    } else {
      loop <- "No Stratification"
    }
    forPlotandTable <- lapply(loop , function(i){
         subs <- if(i=="Males") {
                    males 
                  } else if(i=="Females"){
                    females
                  } else {
                    1:length(traitObject()$trait)
                  }
          x <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=transformMethod)
          if(length(x)>1){
            outputList <<- c(outputList , list(i) , list(x[-1]))
          } else {
            outputList <<- c(outputList , list(i) , list("No other Info"))
          }
          x <- x$norm_data
          return(x)
    })
    names(forPlotandTable) <- loop
    # Store all the output of the transformation, I am sure is going to be useful
    # mytransform$transformer <- forPlotandTable
    # If radiobutton stratified by sex is on Yes, the plot is run twice for m and f
    output$transformerplot <- renderPlot({
      return({
        if(input$sexStratFlag=="Yes"){
          # loop <- c("Males" , "Females")
          par(mfrow=c(4,2))
        } else {
          # loop <- "No Stratification"
          par(mfrow=c(2,2))
        }
        for(i in loop){
          normalizationPlot( x = forPlotandTable[[i]] , i = i , traitLabelFull = traitLabelFull)
        }
        })
    })
    # Display all the residual output from transformation (see box cox example)
    output$normalizationSideEffect <- renderPrint({
      if(is.null(traitObject())){
        return(NULL)
      }
      return(outputList)
    })
    return({
      fluidPage(
        plotOutput("transformerplot",height = "800px")
        ,tags$hr()
        ,verbatimTextOutput("normalizationSideEffect")
      )
    })
  } else {
    # Initialize transformation side effect output
    tabs <- lapply(names(traitObject()) , function(strat){
      outputList <- list(paste(transformMethod , "further Normalization Info:"))
      if(input$sexStratFlag=="Yes"){
        females<-which(traitObject()[[strat]]$sex==2)
        males<-which(traitObject()[[strat]]$sex==1)
        loop <- c("Males" , "Females")
      } else {
        loop <- "No Stratification"
      }
      forPlotandTable <- lapply(loop , function(i){
           subs <- if(i=="Males") {
                      males 
                    } else if(i=="Females"){
                      females
                    } else {
                      1:length(traitObject()[[strat]]$trait)
                    }
            x <- normalizeTraitData(trait=traitObject()[[strat]]$trait[subs] , tm=transformMethod)
            if(length(x)>1){
              outputList <<- c(outputList , list(i) , list(x[-1]))
            } else {
              outputList <<- c(outputList , list(i) , list("No other Info"))
            }
            x <- x$norm_data
            return(x)
      })
      names(forPlotandTable) <- loop
      # Store all the output of the transformation, I am sure is going to be useful
      # mytransform$transformer <- forPlotandTable
      # If radiobutton stratified by sex is on Yes, the plot is run twice for m and f
      output[[paste0(which(strat==names(traitObject())) , "_transformerplot")]] <- renderPlot({
        return({
          if(input$sexStratFlag=="Yes"){
            # loop <- c("Males" , "Females")
            par(mfrow=c(4,2))
          } else {
            # loop <- "No Stratification"
            par(mfrow=c(2,2))
          }
          for(i in loop){
            normalizationPlot( x = forPlotandTable[[i]] , i = i , traitLabelFull = traitLabelFull)
          }
          })
      })
      # Display all the residual output from transformation (see box cox example)
      output[[paste0(which(strat==names(traitObject())) , "_normalizationSideEffect")]] <- renderPrint({
        if(is.null(traitObject())){
          return(NULL)
        }
        return(outputList)
      })
      return({
        tabPanel(strat , fluidPage(
          plotOutput(paste0(which(strat==names(traitObject())) , "_transformerplot"),height = "800px")
          ,tags$hr()
          ,verbatimTextOutput(paste0(which(strat==names(traitObject())) , "_normalizationSideEffect"))
        ))
      })
    })
  return(do.call(tabsetPanel , tabs))
  }
})


#-------------------------------#
#   NORMALIZATION TEST TABLE    #
#-------------------------------#

# Apply all normalization and calculate normality test
# The output is a table with the pvalue for every transformation
# Yellow color is a NON significant pvalue (that's what we are interested in)
# output$normalTable <- DT::renderDataTable({
output$normalTable <- renderUI({
  if(is.null(rawData())){
    output$normalTableEmpty <- emptyPlotter("No Data Yet")
    return(
    plotOutput("normalTableEmpty" , height="800px")
    )
  }
  if(is.null(input$stratifier)){
    if(input$sexStratFlag=="Yes"){
      Transformation <- names(normalizationFunctionsList)
      females<-which(traitObject()$sex==2)
      males<-which(traitObject()$sex==1)
      loop <- c("Males" , "Females")
      normalTable <- lapply(loop , function(sex) {
        lapply(Transformation , function(trans) {
            subs <- if(sex=="Males") males else if(sex=="Females") females else stop("Error in loop")
            x <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=trans)$norm_data
            ad <- tryCatch(ad.test(x)$p.value , error=function(e) NA)
            sw <- tryCatch(shapiro.test(x)$p.value, error=function(e) NA)
            # ks is annoying with warnings in case of ties, so I suppress them
            ks <- suppressWarnings({
                    tryCatch(ks.test(x,rnorm(50))$p.value, error=function(e) NA)
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
            ad <- tryCatch(ad.test(x)$p.value , error=function(e) NA)
            sw <- tryCatch(shapiro.test(x)$p.value, error=function(e) NA)
            # ks is annoying with warnings in case of ties, so I suppress them
            ks <- suppressWarnings({
                    tryCatch(ks.test(x,rnorm(50))$p.value, error=function(e) NA)
                  })
            return(c(trans , ad , sw , ks))
          })
      normalTable <- do.call("rbind" , normalTable)
      colnames(normalTable) <- c("Transformation" 
                                , "Anderson-Darling Test" 
                                , "Shapiro-Wilks Test" 
                                , "Kolmogorov-Smirnov Test")
    }
    output$singleNormalTable <- DT::renderDataTable({
              DT::datatable(normalTable) %>% 
              formatStyle("Anderson-Darling Test" 
                  , backgroundColor = styleInterval(0.05, c('lightgray', 'yellow'))) %>%
              formatStyle("Shapiro-Wilks Test" 
                  , backgroundColor = styleInterval(0.05, c('lightgray', 'yellow'))) %>%
              formatStyle("Kolmogorov-Smirnov Test" 
                  , backgroundColor = styleInterval(0.05, c('lightgray', 'yellow')))
              })
    return({
      fluidPage(DT::dataTableOutput("singleNormalTable")
              ,helpText("If you see a yellow box, p-value is over 0.05 and normality test is OK ;=)")
              ,helpText("Empty cells means that the test could not be evaluated")) 
              })
  } else {
    tabs <- lapply(names(traitObject()) , function(strat) {
      if(input$sexStratFlag=="Yes"){
        Transformation <- names(normalizationFunctionsList)
        females<-which(traitObject()[[strat]]$sex==2)
        males<-which(traitObject()[[strat]]$sex==1)
        loop <- c("Males" , "Females")
        normalTable <- lapply(loop , function(sex) {
          lapply(Transformation , function(trans) {
              subs <- if(sex=="Males") males else if(sex=="Females") females else stop("Error in loop")
              x <- normalizeTraitData(trait=traitObject()[[strat]]$trait[subs] , tm=trans)$norm_data
              ad <- tryCatch(ad.test(x)$p.value , error=function(e) NA)
              sw <- tryCatch(shapiro.test(x)$p.value, error=function(e) NA)
              # ks is annoying with warnings in case of ties, so I suppress them
              ks <- suppressWarnings({
                      tryCatch(ks.test(x,rnorm(50))$p.value, error=function(e) NA)
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
              x <- normalizeTraitData(trait=traitObject()[[strat]]$trait , tm=trans)$norm_data
              ad <- tryCatch(ad.test(x)$p.value , error=function(e) NA)
              sw <- tryCatch(shapiro.test(x)$p.value, error=function(e) NA)
              # ks is annoying with warnings in case of ties, so I suppress them
              ks <- suppressWarnings({
                      tryCatch(ks.test(x,rnorm(50))$p.value, error=function(e) NA)
                    })
              return(c(trans , ad , sw , ks))
            })
        normalTable <- do.call("rbind" , normalTable)
        colnames(normalTable) <- c("Transformation" 
                                  , "Anderson-Darling Test" 
                                  , "Shapiro-Wilks Test" 
                                  , "Kolmogorov-Smirnov Test")
      }
      output[[paste0(which(strat==names(traitObject())) , "normalTable_")]] <- DT::renderDataTable({
                DT::datatable(normalTable) %>% 
                formatStyle("Anderson-Darling Test" 
                    , backgroundColor = styleInterval(0.05, c('lightgray', 'yellow'))) %>%
                formatStyle("Shapiro-Wilks Test" 
                    , backgroundColor = styleInterval(0.05, c('lightgray', 'yellow'))) %>%
                formatStyle("Kolmogorov-Smirnov Test" 
                    , backgroundColor = styleInterval(0.05, c('lightgray', 'yellow')))
                })
      return({
        tabPanel(strat 
              , DT::dataTableOutput(paste0(which(strat==names(traitObject())) , "normalTable_"))
              ,helpText("If you see a yellow box, p-value is over 0.05 and normality test is OK ;=)")
              ,helpText("Empty cells means that the test could not be evaluated"))
                })
    })
    do.call(tabsetPanel , tabs)
  }
})

# Keep the residual for download
residual <- reactiveValues(df_res=NULL)

output$covariateAnalysis <- renderUI({
  if(is.null(rawData())){
    output$covariateAnalysisEmpty <- emptyPlotter("No Data Yet")
    return(
    plotOutput("covariateAnalysisEmpty" , height="800px")
    )
  }
  covariates <- protocolFile()$covariates_tested
  covariates <- if(is.null(covariates)) NA else covariates
  if(is.na(covariates)){
    output$covariateAnalysisEmpty <- emptyPlotter("No covariates selected")
    return(
    plotOutput("covariateAnalysisEmpty" , height="800px")
    )
  }
  currTrait <- protocolFile()$trait
  traitUnit <- protocolFile()$units
  traitLabelFull <- paste(currTrait," (",traitUnit,")",sep="")
  transformMethod <- protocolFile()$transformation_method
  cvr <- unlist(strsplit(covariates , ","))
  covList<- list(covariates=cvr
              , age2Flag=if("age2" %in% cvr) TRUE else NA)
  if(is.null(input$stratifier)){
    if(input$sexStratFlag=="Yes"){
      cvr <- unlist(strsplit(covariates , ",")) %>% .[.!="sex"]
      if(length(cvr)==0){
        mymex <- "If you stratify by sex, you can't have sex...\n\n... as the only covariate :)"
        output$covariateAnalysisEmpty <- emptyPlotter(mymex)
        return(
          plotOutput("covariateAnalysisEmpty" , height="800px")
        )
      }
      covList<- list(covariates=cvr
                , age2Flag=if("age2" %in% cvr) TRUE else NA)
      females<-which(traitObject()$sex==2)
      males<-which(traitObject()$sex==1)
      loop <- list("Males"=males , "Females"=females)
    } else {
      loop <- list("No Stratification"=1:nrow(traitObject()))
    }
    normDataListForPlot <- lapply(loop , function(subs) {
        normData <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=transformMethod)$norm_data
        signifCovs <- checkCovariates(covs=covList
                                    ,x=traitObject()[subs,]
                                    ,normx=normData
                                    )
        if(!is.na(signifCovs[["signCov"]][1])){
          resList<-applySignifCovariates(signifCovs[["signCov"]],traitObject()[subs,],normData)
          normResiduals<-resList$residuals
          covString <- resList$covStr
          retData <- checkResiduals(normResiduals,traitObject()[subs,])
          outputData <- retData$outputData
          normResiduals <- retData$normResiduals
        } else {
          normResiduals <- normData
          covString <- NULL
          retData <- checkResiduals(normResiduals,traitObject()[subs,])
          outputData <- retData$outputData
          normResiduals <- retData$normResiduals
        }
        return(list(signifCovs=signifCovs 
                  , covString=covString 
                  , normResiduals=normResiduals 
                  , normData=normData 
                  , outputData=outputData))
        })
          # First checkpoint, show the summary of the linear model
      output$linearCovariates <- renderPrint(
        invisible(sapply(names(normDataListForPlot) , function(x) {
                      return(cat(toupper(x) , show(normDataListForPlot[[x]][["signifCovs"]][["summary"]]) , sep="\n"))
                    }))
        )
      # Add a side effect for this function, by reporting the actual residual table
      outputData <- do.call("rbind" , lapply(normDataListForPlot , '[[' , "outputData") )
      observe({
        residual$df_res <- outputData
      })
      output$residualsPlot <- renderPlot(
        plotResidual2(
                      tl=currTrait
                      ,tlf=traitLabelFull
                      ,myList=normDataListForPlot
                      ,sexStratFlag=input$sexStratFlag
                      )
      )
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
    return(fluidPage(
            downloadButton("downloadResidual" , label="Download Final Residuals")
            ,verbatimTextOutput("linearCovariates")
            ,plotOutput("residualsPlot", height="800px") )
    )
  } else {
    tabs <- lapply(names(traitObject()) , function(strat) {
      if(input$sexStratFlag=="Yes"){
        cvr <- unlist(strsplit(covariates , ",")) %>% .[.!="sex"]
        if(length(cvr)==0){
          mymex <- "If you stratify by sex, you can't have sex...\n\n... as the only covariate :)"
          output$covariateAnalysisEmpty <- emptyPlotter(mymex)
          return(
            plotOutput("covariateAnalysisEmpty" , height="800px")
          )
        }
        covList<- list(covariates=cvr
                  , age2Flag=if("age2" %in% cvr) TRUE else NA)
        females<-which(traitObject()[[strat]]$sex==2)
        males<-which(traitObject()[[strat]]$sex==1)
        loop <- list("Males"=males , "Females"=females)
      } else {
        loop <- list("No Stratification"=1:nrow(traitObject()[[strat]]))
      }
      normDataListForPlot <- lapply(loop , function(subs) {
          normData <- normalizeTraitData(trait=traitObject()[[strat]]$trait[subs] , tm=transformMethod)$norm_data
          signifCovs <- checkCovariates(covs=covList
                                      ,x=traitObject()[[strat]][subs,]
                                      ,normx=normData
                                      )
          if(!is.na(signifCovs[["signCov"]][1])){
            resList<-applySignifCovariates(signifCovs[["signCov"]],traitObject()[[strat]][subs,],normData)
            normResiduals<-resList$residuals
            covString <- resList$covStr
            retData <- checkResiduals(normResiduals,traitObject()[[strat]][subs,])
            outputData <- retData$outputData
            normResiduals <- retData$normResiduals
          } else {
            normResiduals <- normData
            covString <- NULL
            retData <- checkResiduals(normResiduals,traitObject()[[strat]][subs,])
            outputData <- retData$outputData
            normResiduals <- retData$normResiduals
          }
          return(list(signifCovs=signifCovs 
                    , covString=covString 
                    , normResiduals=normResiduals 
                    , normData=normData 
                    , outputData=outputData))
          })
            # First checkpoint, show the summary of the linear model
        output[[paste0(which(strat==names(traitObject())) , "_linearCovariates")]] <- renderPrint(
          invisible(sapply(names(normDataListForPlot) , function(x) {
                        return(cat(toupper(x) , show(normDataListForPlot[[x]][["signifCovs"]][["summary"]]) , sep="\n"))
                      }))
          )
        # Add a side effect for this function, by reporting the actual residual table
        outputData <- do.call("rbind" , lapply(normDataListForPlot , '[[' , "outputData") )
        observe({
          residual[[paste0(which(strat==names(traitObject())) , "_df_res")]] <- outputData
        })
        output[[paste0(which(strat==names(traitObject())) , "_residualsPlot")]] <- renderPlot(
          plotResidual2(
                        tl=currTrait
                        ,tlf=traitLabelFull
                        ,myList=normDataListForPlot
                        ,sexStratFlag=input$sexStratFlag
                        )
        )
        output[[paste0(which(strat==names(traitObject())) , "_downloadResiduals")]] <- downloadHandler(
          filename=function() {
            timeTag <- Sys.time() %>% 
                  sub(" GMT$" , "" , .) %>% 
                  gsub(":" , "_" , .) %>%
                  gsub("-" , "" , .) %>%
                  sub(" " , "." , .)
            paste(protocolFile()$trait , "residuals" , strat , timeTag , "stand_residuals.txt" , sep=".")
          }
          , content=function(file) {
              write.table(residual[[paste0(which(strat==names(traitObject())) , "_df_res")]]
                      , file=file
                      , sep = "\t"
                      , col.names = TRUE
                      , row.names = FALSE
                      , quote=FALSE)
        })
      return(tabPanel(
              strat
              ,downloadButton(paste0(which(strat==names(traitObject())) , "_downloadResiduals") , label="Download Final Residuals")
              ,verbatimTextOutput(paste0(which(strat==names(traitObject())) , "_linearCovariates"))
              ,plotOutput(paste0(which(strat==names(traitObject())) , "_residualsPlot"), height="800px") )
      )
    })
    do.call(tabsetPanel , tabs)
  }
})

# This observer will wait for the residual table to be present
# When this event is observed, the download button is enabled
observeEvent(!is.null(residual$df_res) , {
    # notify the browser that the residual table is ready to be downloaded
    session$sendCustomMessage("download_ready_res" , list(mex=""))
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
  # This object take the result from traitObject and apply SD based filter
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


#-----------#
# Data Plot #
#-----------#

###################### OLD VERSION DEPRECATED START
# Plot data according to normalization criterion
# output$transformedplot <- renderPlot({
#   if(is.null(rawData())){
#     return({
#       par(mar = c(0,0,0,0))
#       plot(c(0, 1), c(0, 1)
#           , ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
#       text(x = 0.5, y = 0.5, "No data yet"
#         ,cex = 3, col = "black")
#     })
#   }
#   # Protocol Variables
#   currTrait <- protocolFile()$trait
#   traitUnit <- protocolFile()$units
#   traitLabelFull <- paste(currTrait," (",traitUnit,")",sep="")
#   transformMethod <- protocolFile()$transformation_method
#   # If radiobutton stratified by sex is on Yes, the plot is run twice for m and f
#     return({
#       if(input$sexStratFlag=="Yes"){
#         females<-which(traitObject()$sex==2)
#         males<-which(traitObject()$sex==1)
#         loop <- c("Males" , "Females")
#         par(mfrow=c(4,2))
#       } else {
#         loop <- "No Stratification"
#         par(mfrow=c(2,2))
#       }
#       for(i in loop){
#         subs <- if(i=="Males") {
#                   males 
#                 } else if(i=="Females"){
#                   females
#                 } else {
#                   1:length(traitObject()$trait)
#                 }
#         x <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=transformMethod)$norm_data
#         if(length(x)>5000){
#           x <- x[1:5000]
#         }
#         plot( x=1:length(x)
#             , y=x
#             , xlab="Index"
#             , ylab=traitLabelFull
#             , main=paste(i , "ScatterPlot" , paste("N =" , length(x)) , sep="\n")
#             , col=if(i=="Males") "navy" else if(i=="Females") "red" else "black")
#         boxplot( x=x
#                 , main=paste(i , "Boxplot" , sep="\n")
#                 , col="cyan")
#         mymin=round(min(x),2)
#         mymax=round(max(x),2)
#         mymean=round(mean(x),2)
#         mysd=round(sqrt(var(x)),2)
#         mymain=sprintf("min=%s;max=%s;mean=%s;sd=%s",mymin,mymax,mymean,mysd)
#         hist( x=x
#             , main=paste(i , mymain , sep="\n")
#             , xlab=traitLabelFull
#             , prob=TRUE)
#         m <- mean(x, na.rm=TRUE)
#         std <- sd(x,na.rm=TRUE) 
#         curve(dnorm(x,mean=m,sd=std)
#             , add=T
#             , col="forestgreen"
#             , lwd=4)
#         #pv<-as.numeric(unlist(shapiro.test(x))[2])
#         pv<-as.numeric(unlist(ad.test(x))[2])
#         mymain=sprintf("Normal Q-Q plot (Anderson-Darling pval=%s)"
#           ,format(pv,digits=3,sci = T))
#         qqnorm(x
#               ,main=paste(i , mymain , sep="\n"))
#         qqline(x 
#               ,col=if(i=="Males") "navy" else if(i=="Females") "red" else "black"
#               , lwd=4)
#       }
#     })
# })

# The output of every normalization function is a list
# This list contains norm_data as the first element
# The other elements are optional
# In case they are present they will be displayed as a html box at the bottom
# An example is box_cox that report the optimal cut off
# output$normalizationSideEffect <- renderPrint({
#   if(is.null(traitObject())){
#     return(NULL)
#   }
#   transformMethod <- protocolFile()$transformation_method
#   outputList <- list(paste(transformMethod , "further Normalization Info:"))
#   if(input$sexStratFlag=="Yes"){
#     females<-which(traitObject()$sex==2)
#     males<-which(traitObject()$sex==1)
#     loop <- c("Males" , "Females")
#   } else {
#     loop <- "No Stratification"
#   }
#   for(i in loop){
#     subs <- if(i=="Males") {
#                 males 
#               } else if(i=="Females"){
#                 females
#               } else {
#                 1:length(traitObject()$trait)
#               }
#     x <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=transformMethod)
#     if(length(x)>1){
#       outputList <- c(outputList , list(i) , list(x[-1]))
#     } else {
#       outputList <- c(outputList , list(i) , list("No other Info"))
#     }
#   }
#   return(outputList)
# })
###################### OLD VERSION DEPRECATED END

#--------------------#
# COVARIATE ANALYSIS #
#--------------------#

# Linear Model on covariates results
# output$linearCovariates <- renderPrint({
#   if(is.null(rawData())){
#     return(NULL)
#   }
#   covariates <- protocolFile()$covariates_tested
#   covariates <- if(is.na(covariates) | covariates=="") NA else covariates
#   if(is.na(covariates)){
#     return(NULL)
#   }
#   currTrait <- protocolFile()$trait
#   traitUnit <- protocolFile()$units
#   traitLabelFull <- paste(currTrait," (",traitUnit,")",sep="")
#   transformMethod <- protocolFile()$transformation_method
#   cvr <- unlist(strsplit(covariates , ","))
#   covList<- list(covariates=cvr
#               , age2Flag=if("age2" %in% cvr) TRUE else NA)
#   if(input$sexStratFlag=="Yes"){
#     females<-which(traitObject()$sex==2)
#     males<-which(traitObject()$sex==1)
#     loop <- c("Males" , "Females")
#     normData <- lapply(list("Males"=males , "Females"=females) , function(subs) {
#         normData <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=transformMethod)$norm_data
#         signifCovs<-checkCovariates2(covs=covList
#                                 ,x=traitObject()[subs , ]
#                                 ,normx=normData)
#       })
#     return(invisible(sapply(names(normData) , function(x) {
#         return(cat(toupper(x) , show(normData[[x]]) , sep="\n"))
#         }))
#     )
#   } else {
#     # subs <- 1:nrow(traitObject())
#     normData <- normalizeTraitData(trait=traitObject()$trait , tm=transformMethod)$norm_data
#     signifCovs <- checkCovariates2(covs=covList
#                                 ,x=traitObject()
#                                 ,normx=normData)
#     return(show(signifCovs))
#   }
# })

# This object is the residual table
# It is updated every time the residualPlot changes
# residual <- reactiveValues(df_res=NULL)
# Normalize data and check for residual of the covariates
# output$residualsPlot <- renderPlot({
#   if(is.null(rawData())){
#     return({
#       par(mar = c(0,0,0,0))
#       plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
#       text(x = 0.5, y = 0.5, "No data yet"
#         ,cex = 3, col = "black")
#     })
#   }
#   covariates <- protocolFile()$covariates_tested
#   covariates <- if(is.na(covariates) | covariates=="") NA else covariates
#   if(is.na(covariates)){
#     return({
#       par(mar = c(0,0,0,0))
#       plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
#       text(x = 0.5, y = 0.5, "No covariates selected"
#         ,cex = 3, col = "black")
#     })
#   }
#   currTrait <- protocolFile()$trait
#   traitUnit <- protocolFile()$units
#   traitLabelFull <- paste(currTrait," (",traitUnit,")",sep="")
#   transformMethod <- protocolFile()$transformation_method
#   if(input$sexStratFlag=="Yes"){
#     cvr <- unlist(strsplit(covariates , ",")) %>% .[.!="sex"]
#     if(length(cvr)==0){
#       return({
#         par(mar = c(0,0,0,0))
#         plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
#         text(x = 0.5, y = 0.5, "If you stratify by sex, you can't have sex...\n\n... as the only covariate :)"
#           ,cex = 3, col = "black")
#       })
#     }
#     covList<- list(covariates=cvr
#              , age2Flag=if("age2" %in% cvr) TRUE else NA)
#     females<-which(traitObject()$sex==2)
#     males<-which(traitObject()$sex==1)
#     loop <- list("Males"=males , "Females"=females)
#     # loop <- c("Males" , "Females")
#     normDataListForPlot <- lapply(loop , function(subs) {
#       normData <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=transformMethod)$norm_data
#       signifCovs<-checkCovariates(covs=covList
#                                   ,x=traitObject()[subs,]
#                                   ,normx=normData
#                                   ,sxf=sexStratFlag)
#       if(!is.na(signifCovs[1])){
#         resList<-applySignifCovariates(signifCovs,traitObject()[subs,],normData)
#         normResiduals<-resList$residuals
#         covString <- resList$covStr
#         retData <- checkResiduals(normResiduals,traitObject()[subs,],"")
#         outputData <- retData$outputData
#         normResiduals <- retData$normResiduals
#       } else {
#         normResiduals <- normData
#         covString <- NULL
#         retData <- checkResiduals(normResiduals,traitObject()[subs,],"")
#         outputData <- retData$outputData
#         normResiduals <- retData$normResiduals
#       }
#       return(list(covString=covString , normResiduals=normResiduals , normData=normData , outputData=outputData))
#       })
#     # Add a side effect for this function, by reporting the actual residual table
#     outputData <- do.call("rbind" , lapply(normDataListForPlot , '[[' , "outputData") )
#     observe({
#       residual$df_res <- outputData
#     })
#     return(plotResidualDataBySex(
#                     tl=currTrait
#                     ,tlf=traitLabelFull
#                     ,myList=normDataListForPlot
#                     ))
#   } else {
#     # subs <- 1:nrow(traitObject())
#     # normData <- normalizeTraitData(trait=traitObject()$trait[subs] , tm=transformMethod)$norm_data
#     normData <- normalizeTraitData(trait=traitObject()$trait , tm=transformMethod)$norm_data
#     cvr <- unlist(strsplit(covariates , ","))
#     covList<- list(covariates=cvr
#               , age2Flag=if("age2" %in% cvr) TRUE else NA)
#     signifCovs<-checkCovariates(covs=covList
#                                 ,x=traitObject()
#                                 ,normx=normData
#                                 ,sxf=sexStratFlag)
#     if(!is.na(signifCovs[1])){
#       resList<-applySignifCovariates(signifCovs,traitObject(),normData)
#       normResiduals<-resList$residuals
#       covString <- resList$covStr
#       retData <- checkResiduals(normResiduals,traitObject(),"")
#       outputData <- retData$outputData
#       normResiduals <- retData$normResiduals
#     } else {
#       normResiduals <- normData
#       covString <- NULL
#       retData <- checkResiduals(normResiduals,traitObject(),"")
#       outputData <- retData$outputData
#       normResiduals <- retData$normResiduals
#     }
#     # Side effect, update residual final table
#     observe({
#       residual$df_res <- outputData
#     })
#     return(plotResidualData(tl=currTrait
#                     ,tlf=traitLabelFull
#                     # ,"No Stratification"
#                     ,cvs=covString
#                     ,res=normResiduals
#                     ,normx=normData))
#   }
# })

#-------------------#
# DOWNLOAD RESIDUAL #
#-------------------#
# OLD AND DEPRECATED
# Download residuals table
# output$downloadResidual <- downloadHandler(
#     filename=function() {
#       timeTag <- Sys.time() %>% 
#             sub(" GMT$" , "" , .) %>% 
#             gsub(":" , "_" , .) %>%
#             gsub("-" , "" , .) %>%
#             sub(" " , "." , .)
#       paste(protocolFile()$trait , "residuals" , timeTag , "stand_residuals.txt" , sep=".")
#     }
#     , content=function(file) {
#         write.table(residual$df_res
#                 , file=file
#                 , sep = "\t"
#                 , col.names = TRUE
#                 , row.names = FALSE
#                 , quote=FALSE)
# })


# SOME ATTEMPTS TO A BETTER AESTHETIC
#    ,tags$head(
#   tags$style(type="text/css", "label.radio { display: inline-block; }", ".radio input[type=\"radio\"] { float: none; }"),
#   tags$style(type="text/css", "select { max-width: 200px; }"),
#   tags$style(type="text/css", "textarea { max-width: 4000px; }"),
#   tags$style(type="text/css", ".jslider { max-width: 200px; }"),
#   tags$style(type='text/css', ".well { padding: 12px; margin-bottom: 5px; max-width: 280px; }"),
#   tags$style(type='text/css', ".span4 { max-width: 10px; }"),
#   tags$style(type='text/css', ".span8 { max-width: 4000px; min-width: 4000px;}")
# )
# tags$head(HTML('<style>.span2 {min-width: 265px; max-width: 265px; }</style>'))