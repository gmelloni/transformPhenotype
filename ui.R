library(shinythemes)
library(shinydashboard)

# This tag is a button disabler (triggered by some event)
tagfordisable <- singleton(tags$head(HTML(
'
  <script type="text/javascript">
    $(document).ready(function() {
      // disable download at startup. downloadFinalProtocol is the id of the downloadButton
      $("#downloadFinalProtocol").attr("disabled", "true").attr("onclick", "return false;");
      Shiny.addCustomMessageHandler("download_ready", function(message) {
        $("#downloadFinalProtocol").removeAttr("disabled").removeAttr("onclick").html(
          "<i class=\\"fa fa-download\\"></i> Download of final protocol file is ready" + message.mex);
      });
    })
  </script>
'
)))


shinyUI(
  dashboardPage(skin="purple"
    ,header=dashboardHeader(title = "Transform Phenotype App" , titleWidth=350)
    ,sidebar=dashboardSidebar(width=200
      #Logo
      ,imageOutput("myImage")
      ,sidebarMenu(
        menuItem("Analysis" , tabName="analysispan" , icon=icon("area-chart"))
        ,menuItem("Help" , tabName="help" , icon=icon("info-sign" , lib="glyphicon"))
        ,menuItem("Contacts" , tabName="contact" , icon=icon("list"))
      )
      ,disable=FALSE
    )
    ,body=dashboardBody(
       tags$head(tags$style(HTML("
          @import url('//fonts.googleapis.com/css?family=Lobster|Cabin:400,700');
        ")))
      ,tags$head(tags$style(HTML('
          .main-header .logo {
            font-family: "Lobster", cursive;
            font-weight: bold;
          }
      ')))
      ,tabItems(
        tabItem("analysispan" , 
          fluidPage(
            tagfordisable
            ,sidebarLayout(
              sidebarPanel(id="sideID"
              ,style = "overflow-y:scroll; max-height: 1000px; position:relative;"
              ,fileInput('myfile'
                ,'Choose a tab delimited file'
                # , accept=c('text/comma-separated-values' , 'text/plain') 
                , multiple=TRUE)
              ,actionButton("example" , "Load Example Dataset")
              ,helpText("Transform Phenotype accepts a tab-delimited file with header")
              ,tags$hr()
              ,actionButton("storeattempt", "Store this transformation")
              ,actionButton("deleteattempt", "Delete this transformation")
              ,selectInput("trait" , "Choose A Trait to analyze:"
                ,choices=c("")
                ,multiple=FALSE
                ,selected=""
              )
              ,textInput("units", label = "Input Unit of measurement", value = "nmol/l")
              ,selectInput("transformation_method", 
                "Choose a function for the transformation:", 
                ,choices=names(normalizationFunctionsList)
                ,multiple=FALSE
                ,selected="untransformed")
              ,selectInput("sd_num"
                ,"Filter Outliers based on SD value:"
                ,choices=c( 1 , 2 , 3 , 4 , 5 , 6 , NA )
                ,multiple=FALSE
                ,selected=NA)
              ,selectInput("sd_num_sex"
                           ,"Filter Outliers based on SD value (Sex stratified):"
                           ,choices=c( 1 , 2 , 3 , 4 , 5 , 6 , NA )
                           ,multiple=FALSE
                           ,selected=NA)
              ,radioButtons('sd_dir', 'Set Direction of the SD filter'
                   ,c("Greater than"='gt'
                     ,"Less than"='lt'
                     ,"Both"='both')
                    ,'both')
              ,sliderInput("filters", "Additional filter on raw values:"
                ,min = 0
                ,max = 1
                ,value = c(0,1))
              ,selectInput("covariates_tested", 
                "Choose one or more covariate:", 
                ,choices=NULL
                ,multiple=TRUE
                # ,selected=NA
                )
              # ,helpText("If you select NA among covariates, it will overwrite any other choice")
              ,tags$hr()
              ,radioButtons('sexStratFlag', 'Stratify by sex?'
                   ,c("Yes"="Yes"
                     ,"No"="No"
                     )
                    ,"No")
              ,selectInput("stratifier"
                ,"Choose one stratification column:"
                ,choices=NULL
                ,multiple=TRUE
                # ,selected=NA
                )
            ,width=3)
            ,mainPanel(
              tabsetPanel(
                type = "tabs"
                ,tabPanel("Uploaded Data" , fluidPage(DT::dataTableOutput("contents")
                                          , verbatimTextOutput("dferror")
                                          ))
                , tabPanel("Gender Difference" , fluidPage(
                    uiOutput("sexplot")
                    )
                  )
                , tabPanel("Transform Plot" 
                    ,uiOutput("transformer")
                  )
                , tabPanel("Normality Summary" , fluidPage(
                    uiOutput("normalTable")
                    )
                  )
                , tabPanel("Covariate Analysis" , uiOutput("covariateAnalysis")
                  )
                , tabPanel("Protocol" , fluidPage(
                      p(h3(span("Current Trait" , style = "color:blue")))
                      ,tableOutput("protocolFile")
                      ,tagfordisable
                      ,p(h3(span("Cumulative Protocol File" , style = "color:blue")))
                      ,downloadButton("downloadFinalProtocol" , label="Download Protocol File")
                      ,helpText("Download will be available after first click to Store this transformation.")             
                      ,DT::dataTableOutput("cumulativeProtocolFile")
                      ))
              )
            ,width=9)
          ))
        )
        #### HELP SECTION IN DEVELOPMENT ####
        # ,tabItem("help" 
        #   , fluidRow(
        #                   column(4,
        #                       # h3("What's this app?"),
        #                       imageOutput("cheatsheet")
        #                       ),
        #                   column(8,
        #                       htmlOutput("myVideo")
        #                       )
        #                       )
        #         )
        ,tabItem("contact" 
          ,fluidPage(sidebarLayout(
            sidebarPanel(h2("TransformPhenotype Contacts"))
            ,mainPanel(
              p(h4("The TransformPhenotype App was created by Giorgio Melloni at the Wellcome Trust Sange Institute")),
              p(h4("The original code was written by Angela Matchan and improved by Arthur Gilly, Rachel Moore, Loz Southam and Vincent Appel")),
              p(h4("For any question or suggestion regarding the app, write to" 
                  , span("melloni.giorgio [a] gmail.com" , style = "color:blue"))),
              br(),
              p(h3(strong("Thank you for using our app :)")))
            ))))
      )
    )
  )
)


#-----------------------------#
#       CODE CEMETERY         #
#-----------------------------#
# FAILED AESTHESTICS IMPROVEMENTS
# tags$head(HTML('<style>.span2 {min-width: 265px; max-width: 265px; }</style>'))
         #    ,tags$head(
          #   tags$style(type="text/css", "label.radio { display: inline-block; }", ".radio input[type=\"radio\"] { float: none; }"),
          #   tags$style(type="text/css", "select { max-width: 200px; }"),
          #   tags$style(type="text/css", "textarea { max-width: 4000px; }"),
          #   tags$style(type="text/css", ".jslider { max-width: 200px; }"),
          #   tags$style(type='text/css', ".well { padding: 12px; margin-bottom: 5px; max-width: 280px; }"),
          #   tags$style(type='text/css', ".span4 { max-width: 10px; }"),
          #   tags$style(type='text/css', ".span8 { max-width: 4000px; min-width: 4000px;}")
            # )
            # , conditionalDisabledPanel(element="goDDButton" , tagmessage="go_with_the_analysis")

# UNUSED DISABLING TAG
# tagfordisable2 <- singleton(tags$head(HTML(
# '
#   <script type="text/javascript">
#     $(document).ready(function() {
#       // disable download at startup. downloadResidual is the id of the downloadButton
#       $("#downloadResidual").attr("disabled", "true").attr("onclick", "return false;");
#       Shiny.addCustomMessageHandler("download_ready_res", function(message) {
#         $("#downloadResidual").removeAttr("disabled").removeAttr("onclick").html(
#           "<i class=\\"fa fa-download\\"></i> Download of residual table file is ready" + message.mex);
#       });
#     })
#   </script>
# '
# )))