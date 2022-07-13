suppressPackageStartupMessages(
{
  library(shiny)
  library(shinydashboard)
  library(shinyBS)
  library(ggpubr)
  library(shinyjs)
  
  # library(table1)
  # library(survminer)
  # library(survival)
  # library(DT)
  # library(psych)
  # library(tibble)
  # library(survivalROC)
  # library(tidyverse)
  # library(WGCNA)
  # library(limma)
  # library(colourpicker)
  # library(randomForestSRC)
  # library(mlr)
  # library(enrichplot)
  # library(rms)
  # library(TCGAbiolinks)
  # library(ReactomePA)
  # library(clusterProfiler)
  # library(msigdbr)
  # library(caret)
  # library(Biobase)
  # library(CuratedCancerPrognosisData)
  # library(ggplotify)
  # library(ConsensusTME)
  # library(M3C)
  # library(ComplexHeatmap)
  # library(estimate)
  # library(ConsensusClusterPlus)
}
)

ui.path <- ifelse(system.file("app", package = "CBioProfiler") == "", "ui",
                  file.path(system.file("app", package = "CBioProfiler"),"ui"))

source("funcs.R")
dataset<-read.csv("dataset.csv",header=T,row.names=1)
# PCBC_stemSig<-readRDS("PCBC_stemSig.rds")

# common_genes<-readRDS("common_genes.rds")

header <- dashboardHeader(
    title = "CBioProfiler",
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
      tags$div("Biomarker Analysis",
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
               menuSubItem("Correlation with immune infiltration",tabName = "immune",icon=NULL),
               menuSubItem("Correlation with stemness score",tabName = "stemness",icon=NULL),
               menuSubItem("Correlation with ESTIMATE score",tabName = "estimate",icon=NULL),
               menuSubItem("Correlation with immune checkpoint",tabName = "ICB",icon=NULL),
               menuSubItem("Correlation with IFN-gamma score",tabName = "IFN",icon=NULL),
               menuSubItem("Correlation with cytolytic activity",tabName = "CYA",icon=NULL),
               menuSubItem("Correlation with cancer pathway",tabName = "oncopath",icon=NULL),
               menuSubItem("Correlation with metabolism pathway",tabName = "metapath",icon=NULL),
               menuSubItem("Correlation with hallmark signature",tabName = "hallpath",icon=NULL),
               menuSubItem("Correlation with drug response",tabName = "drug",icon=NULL)
      ),
      menuItem("Biological annotation", icon=icon("dna"),
               menuSubItem("Biological annotation",tabName = "bioan",icon=NULL)
      ),
      # menuItem("Meta analysis", icon=icon("align-center"),
      #          menuSubItem("Meta-analysis of biomarker",tabName = "meta",icon=NULL)
      # ),
      
      
      tags$div("Subtype Analysis",
               style = "margin-top: 6px;
                          font-size: 1.5em;
                          padding: 0 1.25em;
                          text-align: center;
                          background: rgba(255, 255, 255, 0);
                          color: #909CFF"),
      menuItem("Subtype identification", icon=icon("sitemap"),
               menuSubItem("Subtype identification",tabName = "cba",icon=NULL)),
      menuItem("Subtype characterization", icon=icon("fingerprint"),
               menuSubItem("Correlation with clinical features",tabName = "cbaclinical",icon=NULL),
               menuSubItem("Kaplan-Meier curve",tabName = "cbaKM",icon=NULL),
               menuSubItem("CoxPH model",tabName = "cbaCoxPH",icon=NULL),
               menuSubItem("Time-dependent ROC",tabName = "cbaSurvROC1",icon=NULL),
               # menuSubItem("Most correlated genes",tabName = "mcorgene",icon=NULL),
               # menuSubItem("Correlation with specific gene",tabName = "gene",icon=NULL),
               menuSubItem("Differentially expressed genes",tabName = "cbaDEG",icon=NULL),
               menuSubItem("Immune infiltration among subtypes",tabName = "cbaimmune",icon=NULL),
               menuSubItem("Stemness score among subtypes",tabName = "cbastemness",icon=NULL),
               menuSubItem("ESTIMATE score among subtypes",tabName = "cbaestimate",icon=NULL),
               menuSubItem("Immune checkpoints among subtypes",tabName = "cbaICB",icon=NULL),
               menuSubItem("Interferon-gamma among subtypes",tabName = "cbaIFN",icon=NULL),
               menuSubItem("Cytolytic activity among subtypes",tabName = "cbaCYA",icon=NULL),
               menuSubItem("Cancer pathway score among subtypes",tabName = "cbaoncopath",icon=NULL),
               menuSubItem("Metabolisim score among subtypes",tabName = "cbametapath",icon=NULL),
               menuSubItem("Hallmark signature among subtypes",tabName = "cbahall",icon=NULL),
               menuSubItem("Drug response among subtypes",tabName = "cbadrug",icon=NULL)
               
      ),

      tags$div("Source",
               style= "margin-top: 6px;
                         font-size: 1.5em;
                         padding: 0 1.25em;
                         text-align: center;
                         background: rgba(255, 255, 255, 0);
                        color: yellow"),

      menuSubItem("Dataset source", href = "https://liuxiaoping2020.github.io/CBioProfilerDatasource/",icon=icon("external-link")),
      menuSubItem("CBioProfiler standalone app", href = "https://gitee.com/liuxiaoping2020/CBioProfiler",icon=icon("external-link")),
      menuSubItem("R package reference", href = "https://liuxiaoping2020.github.io/CBioProfiler_reference/",icon=icon("external-link")),

      tags$div("Tutorial",
               style= "margin-top: 6px;
                         font-size: 1.5em;
                         padding: 0 1.25em;
                         text-align: center;
                         background: rgba(255, 255, 255, 0);
                        color: deepskyblue"),

      menuSubItem("Tutorial" , href = "https://github.com/liuxiaoping2020/CBioProfiler_tutorial/blob/main/CBioProfiler_tutorial.pdf", icon = icon("book")),
      tags$div("Contact",
               style= "margin-top: 6px;
                           font-size: 1.5em;
                           padding: 0 1.25em;
                           text-align: center;
                           background: rgba(255, 255, 255, 0);
                          color: lightgreen"),
      tags$div("Xing-Huan Wang",
               style= "margin-top: 6px;
                           font-size: 1em;
                           padding: 0 1.25em;
                           text-align: left;
                           background: rgba(255, 255, 255, 0);
                          color: white"),

      tags$div("Email: wangxinghuan@whu.edu.cn",
               style= "margin-top: 6px;
                           font-size: 1em;
                           padding: 0 1.25em;
                           text-align: left;
                           background: rgba(255, 255, 255, 0);
                          color: white"),
      br(),
      tags$div("Xiao-Ping Liu",
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
                          color: white"),
      br(),
 tags$div("Sheng Li",
               style= "margin-top: 6px;
                           font-size: 1em;
                           padding: 0 1.25em;
                           text-align: left;
                           background: rgba(255, 255, 255, 0);
                          color: white"),

      tags$div("Email: lisheng-znyy@whu.edu.cn",
               style= "margin-top: 6px;
                           font-size: 1em;
                           padding: 0 1.25em;
                           text-align: left;
                           background: rgba(255, 255, 255, 0);
                          color: white"),


# useShinyjs(),                                          
# actionButton("reset_button",
#              "Restart CBioProfiler",
#              type="success",
#              style = "background-color: forestgreen;
#                                             color: floralwhite;
#                                             margin-left: auto;
#                                             margin-right: auto;
#                                             width: 80%",
#              icon = icon("refresh"))
br(),

actionButton("reset_button",
             "Restart CBioProfiler",
              type="success",
              style = "background-color: #000080;
                      color: floralwhite;
                      margin-left: auto;
                      margin-right: auto;
                      width: 80%",
                      icon = icon("refresh"))
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
    source(file.path(ui.path, "stemness.R"),local=T)$value,
    source(file.path(ui.path, "estimate.R"),local=T)$value,
    source(file.path(ui.path, "ICB.R"),local=T)$value,
    source(file.path(ui.path, "IFN.R"),local=T)$value,
    source(file.path(ui.path, "CYA.R"),local=T)$value,
    source(file.path(ui.path, "oncopath.R"),local=T)$value,
    source(file.path(ui.path, "metabol.R"),local=T)$value,
    source(file.path(ui.path, "hall.R"),local=T)$value,
    source(file.path(ui.path, "drug.R"),local=T)$value,


    source(file.path(ui.path, "cba.R"),local=T)$value,
    source(file.path(ui.path, "cbaclinical.R"),local=T)$value,
    source(file.path(ui.path, "cbakm.R"),local=T)$value,
    source(file.path(ui.path, "cbaCoxPH.R"),local=T)$value,
    source(file.path(ui.path, "cbaSurvROC2.R"),local=T)$value,
    source(file.path(ui.path, "cbaDEG.R"),local=T)$value,
    source(file.path(ui.path, "cbaimm.R"),local=T)$value,
    source(file.path(ui.path, "cbastem.R"),local=T)$value,
    source(file.path(ui.path, "cbaestimate.R"),local=T)$value,
    source(file.path(ui.path, "cbaICB.R"),local=T)$value,
    source(file.path(ui.path, "cbaIFN.R"),local=T)$value,
    source(file.path(ui.path, "cbaCYA.R"),local=T)$value,
    source(file.path(ui.path, "cbaoncopath.R"),local=T)$value,
    source(file.path(ui.path, "cbametabol.R"),local=T)$value,
    source(file.path(ui.path, "cbahall.R"),local=T)$value,
    source(file.path(ui.path, "cbadrug.R"),local=T)$value
  )
)

shinyUI(
    bootstrapPage(
      # useShinyjs(),
      # div(id = "loading-content",
      #     img(src = "zcool.gif",align = "center"),
      #     # img(src = "progress.gif"),
      #     
      #     tags$br(),
      #     "Loading CBioProfiler..."
      # ),
        dashboardPage(title="CBioProfiler",
                      skin = "green",
                      header,
                      sidebar,
                      body)

    )
)






