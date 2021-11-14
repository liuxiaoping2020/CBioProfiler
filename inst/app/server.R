
source("funcs.R")
dataset<-read.csv("dataset.csv",header=T,row.names=1)


CBioExplorerServer <- function(input, output, session) {

    session$onSessionEnded(stopApp)
    server.path <- ifelse(system.file("app", package = "CBioExplorer") == "",
                          "server",
                          file.path(system.file("app", package = "CBioExplorer"),"server"))
    eset.path <- ifelse(system.file("app", package = "CBioExplorer") == "",
                          "eset",
                          file.path(system.file("app", package = "CBioExplorer"),"eset"))


    source(file.path(server.path, "datainput.R"),  local = TRUE)$value
    source(file.path(server.path, "ClinicalCorrelation.R"),  local = TRUE)$value
    source(file.path(server.path, "KM.R"),  local = TRUE)$value
    source(file.path(server.path, "CoxPH.R"),  local = TRUE)$value
    source(file.path(server.path, "SurvROC.R"),  local = TRUE)$value
    source(file.path(server.path, "genecor.R"),  local = TRUE)$value
    source(file.path(server.path, "mcorgene.R"),  local = TRUE)$value
    source(file.path(server.path, "genediff.R"),  local = TRUE)$value
    source(file.path(server.path, "unicox.R"),  local = TRUE)$value
    source(file.path(server.path, "DEG.R"),  local = TRUE)$value
    source(file.path(server.path, "WGCNA.R"),  local = TRUE)$value
    source(file.path(server.path, "nestr.R"),  local = TRUE)$value
    source(file.path(server.path, "predmo.R"),  local = TRUE)$value
    source(file.path(server.path, "valmo.R"),  local = TRUE)$value
    source(file.path(server.path, "nomo.R"),  local = TRUE)$value
    source(file.path(server.path, "Bioann.R"),  local = TRUE)$value
    source(file.path(server.path, "immune.R"),  local = TRUE)$value
    source(file.path(server.path, "stemness.R"),  local = TRUE)$value

    observe({
        if(!is.null(input$test)) stopApp()  # stop shiny
    })

}

