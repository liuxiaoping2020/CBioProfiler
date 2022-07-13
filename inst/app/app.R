library(shiny)
# library(shinyjs)
# jsResetCode <- "shinyjs.reset = function() {history.go(0)}" 
ui <- navbarPage(
    source(file.path("ui", "ui.R"),  local = TRUE)$value
)

server <- function(input, output, session) {

    source(file.path("server", "server.R"),  local = TRUE)$value
}

shinyApp(ui = ui, server = server)



 