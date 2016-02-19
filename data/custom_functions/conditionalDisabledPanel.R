conditionalDisabledPanel <- function(element , tagmessage) {
  scriptOpen <- '<script type="text/javascript">'
  startFunc <- '$(document).ready(function() {'
  initialAttr <- paste0('$("#' , element , '").attr("disabled", "true").attr("onclick", "return false;");')
  condition <- paste0('Shiny.addCustomMessageHandler("' , tagmessage , '", function(message) {')
  finalAttr <- paste0('$("#' , element , '").removeAttr("disabled").removeAttr("onclick").html(message.mex);')
  closescript <- paste('});' , '})' , '</script>' , sep="\n")
  singleton(tags$head(HTML(
    paste(scriptOpen
          ,startFunc
          ,initialAttr
          ,condition
          ,finalAttr
          ,closescript
          ,sep="\n")
  )))
}