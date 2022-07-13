
source("funcs.R")
dataset<-read.csv("dataset.csv",header=T,row.names=1)

# PCBC_stemSig<-readRDS("PCBC_stemSig.rds")
common_genes<-readRDS("common_genes.rds")
oncopath<-readRDS("oncopath.rds")
metabolism<-readRDS("metabolismpath.rds") 

CBioExplorerServer <- function(input, output, session) {
    server.path <- ifelse(system.file("app", package = "CBioExplorer") == "", "server",
                          file.path(system.file("app", package = "CBioExplorer"),"server"))
    eset.path <- ifelse(system.file("app", package = "CBioExplorer") == "","eset",
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
    source(file.path(server.path, "estimate.R"),  local = TRUE)$value
    source(file.path(server.path, "ICB.R"),  local = TRUE)$value
    source(file.path(server.path, "IFN.R"),  local = TRUE)$value
    source(file.path(server.path, "CYA.R"),  local = TRUE)$value
    source(file.path(server.path, "oncopath.R"),  local = TRUE)$value
    source(file.path(server.path, "metabol.R"),  local = TRUE)$value
    source(file.path(server.path, "hall.R"),  local = TRUE)$value
    source(file.path(server.path, "drug.R"),  local = TRUE)$value

    
    source(file.path(server.path, "cba.R"),  local = TRUE)$value
    source(file.path(server.path, "cbaclinical.R"),  local = TRUE)$value
    source(file.path(server.path, "cbakm.R"),  local = TRUE)$value
    source(file.path(server.path, "cbaCoxPH.R"),  local = TRUE)$value
    source(file.path(server.path, "cbaSurvROC.R"),  local = TRUE)$value
    source(file.path(server.path, "cbaDEG.R"),  local = TRUE)$value
    source(file.path(server.path, "cbaimm.R"),  local = TRUE)$value
    source(file.path(server.path, "cbastem.R"),  local = TRUE)$value
    source(file.path(server.path, "cbaestimate.R"),  local = TRUE)$value
    source(file.path(server.path, "cbaICB.R"),  local = TRUE)$value
    source(file.path(server.path, "cbaIFN.R"),  local = TRUE)$value
    source(file.path(server.path, "cbaCYA.R"),  local = TRUE)$value
    source(file.path(server.path, "cbaoncopath.R"),  local = TRUE)$value
    source(file.path(server.path, "cbametabol.R"),  local = TRUE)$value
    source(file.path(server.path, "cbahall.R"),  local = TRUE)$value
    source(file.path(server.path, "cbadrug.R"),  local = TRUE)$value
    # observe({
    #     if(!is.null(input$test)) stopApp()  # stop shiny
    # })
    observeEvent(input$reset_button,
                 {
                   withProgress(message = "Restarting CBioExplorer",
                                detail = "This may take a while...",
                                value = 3,
                                {
                   session$reload()
                                })
                 }
    )
    
    # hide("loading-content", TRUE, "fade")
}
