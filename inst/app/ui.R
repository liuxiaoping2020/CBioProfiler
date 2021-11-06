suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyBS)
  library(table1)
  library(survminer)
  library(survival)
  library(DT)
  library(psych)
  library(survivalROC)
  library(tidyverse)
  library(WGCNA)
  library(ggpubr)
  library(limma)
  library(colourpicker)
  library(randomForestSRC)
  library(mlr)
  library(shinyjs)
  library(enrichplot)
  library(rms)
  library(TCGAbiolinks)
  library(ReactomePA)
  library(clusterProfiler)
  library(msigdbr)
  library(caret)
  library(Biobase)
  library(CuratedCancerPrognosisData)
})


ui.path <- ifelse(system.file("app", package = "CBioExplorer") == "",
                  "ui",
                  file.path(system.file("app", package = "CBioExplorer"),"ui"))


source("funcs.R")
dataset<-read.csv("dataset.csv",header=T,row.names=1)

header <- dashboardHeader(
    title = "CBioExplorer",
    titleWidth = 250
)

sidebar <-  dashboardSidebar(
    width = 250,
    sidebarMenu(
      id = "tabs",

      tags$div("Introduction",
               style= "font-size: 1.5em;
                 margin-top: 6px;
                 padding:0 0.25em;
                 text-align: center;
                 background: rgba(255, 255, 255, 0);
                 color: orange"),
      menuItem("Introduction", tabName = "welcome1", icon = icon("home"), selected = T),

      tags$div("Data",
               style= "font-size: 1.5em;
                 margin-top: 6px;
                 padding:0 0.25em;
                 text-align: center;
                 background: rgba(255, 255, 255, 0);
                 color: #FF5151"),
      menuSubItem("Data input", tabName = "dataset",icon=icon("database")),
      tags$div("Analysis",
               style = "margin-top: 6px;
                          font-size: 1.5em;
                          padding: 0 1.25em;
                          text-align: center;
                          background: rgba(255, 255, 255, 0);
                          color: #0AC71B"),

      menuItem("Dimensionality reduction", icon=icon("tree"),
               menuSubItem("WGCNA",tabName = "wgcna",icon=NULL),
               menuSubItem("Survival related genes",tabName = "msurv",icon=NULL),
               menuSubItem("Differentially expressed genes",tabName = "DEG",icon=NULL)
      ),
      menuItem("Benchmark experiment", icon=icon("laptop"),

               menuSubItem("Benchmark experiment",tabName = "nestr",icon=NULL)),

      menuItem("Prediction model", icon=icon("notes-medical"),
               menuSubItem("Construct model",tabName = "bmpm",icon=NULL),
               menuSubItem("Validate model",tabName = "valmo",icon=NULL),
               menuSubItem("Nomogram",tabName = "nomo",icon=NULL)
      ),
      menuItem("Clinical annotation", icon=icon("layer-group"),
               menuSubItem("Correlation with clinical features",tabName = "clinical",icon=NULL),
               menuSubItem("Kaplan-Meier curve",tabName = "KM",icon=NULL),
               menuSubItem("CoxPH model",tabName = "CoxPH",icon=NULL),
               menuSubItem("Time-dependent ROC",tabName = "SurvROC",icon=NULL),
               menuSubItem("Most correlated genes",tabName = "mcorgene",icon=NULL),
               menuSubItem("Correlation with specific gene",tabName = "gene",icon=NULL),
               menuSubItem("Gene expression in different groups",tabName = "genediff",icon=NULL),
               menuSubItem("Correlation with Immune infiltration",tabName = "immune",icon=NULL),
               menuSubItem("Correlation with stemness score",tabName = "stemness",icon=NULL)
      )

      ,
      menuItem("Biological annotation", icon=icon("dna"),
               menuSubItem("Biological annotation",tabName = "bioan",icon=NULL)
      ) ,

      tags$div("Source",
               style= "margin-top: 6px;
                         font-size: 1.5em;
                         padding: 0 1.25em;
                         text-align: center;
                         background: rgba(255, 255, 255, 0);
                        color: #909CFF"),

      menuSubItem("Dataset soruce", href = "https://liuxiaoping2020.github.io/CBioExplorerDatasource/",icon=icon("external-link")),
      menuSubItem("CBioExplorer standalone app", href = "https://liuxiaoping2020.github.io/CBioExplorer_reference/",icon=icon("external-link")),
      menuSubItem("R package reference", href = "https://liuxiaoping2020.github.io/CBioExplorer_reference/",icon=icon("external-link")),

      tags$div("Tutorial",
               style= "margin-top: 6px;
                         font-size: 1.5em;
                         padding: 0 1.25em;
                         text-align: center;
                         background: rgba(255, 255, 255, 0);
                        color: deepskyblue"),

      menuSubItem("Tutorial" , href = "https://github.com/liuxiaoping2020/CBioExplorer/blob/main/CBioExplorer%20user%20guide.pdf", icon = icon("book")),
      tags$div("Contact",
               style= "margin-top: 6px;
                           font-size: 1.5em;
                           padding: 0 1.25em;
                           text-align: center;
                           background: rgba(255, 255, 255, 0);
                          color: lightgreen"),
      tags$div("First author: Xiao-Ping Liu",
               style= "margin-top: 6px;
                           font-size: 1em;
                           padding: 0 1.25em;
                           text-align: left;
                           background: rgba(255, 255, 255, 0);
                          color: white"),

      tags$div("Email: liuxiaoping@whu.edu.cn",
               style= "margin-top: 6px;
                           font-size: 1em;
                           padding: 0 1.25em;
                           text-align: left;
                           background: rgba(255, 255, 255, 0);
                          color: white")


        )
)

body <-  dashboardBody(

  tabItems(
    source(file.path(ui.path, "intro.R"),local=T)$value,
    source(file.path(ui.path, "datainput.R"),local=T)$value,
    source(file.path(ui.path, "ClinicalCorrelation.R"),local=T)$value,
    source(file.path(ui.path, "KM.R"),local=T)$value,
    source(file.path(ui.path, "CoxPH.R"),local=T)$value,
    source(file.path(ui.path, "SurvROC.R"),local=T)$value,
    source(file.path(ui.path, "genecor.R"),local=T)$value,
    source(file.path(ui.path, "mcorgene.R"),local=T)$value,
    source(file.path(ui.path, "genediff.R"),local=T)$value,
    source(file.path(ui.path, "unicox.R"),local=T)$value,
    source(file.path(ui.path, "DEG.R"),local=T)$value,
    source(file.path(ui.path, "WGCNA.R"),local=T)$value,
    source(file.path(ui.path, "nestr.R"),local=T)$value,
    source(file.path(ui.path, "predmo.R"),local=T)$value,
    source(file.path(ui.path, "valmo.R"),local=T)$value,
    source(file.path(ui.path, "nomo.R"),local=T)$value,
    source(file.path(ui.path, "Bioann.R"),local=T)$value,
    source(file.path(ui.path, "immune.R"),local=T)$value,
    source(file.path(ui.path, "stemness.R"),local=T)$value


  )
)


shinyUI(
    bootstrapPage(

        dashboardPage(title="CBioExplorer",
                      skin = "green",
                      header,
                      sidebar,
                      body)

    )
)
