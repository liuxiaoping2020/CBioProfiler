tabItem (tabName = "SurvROC",
  fluidRow(column(
    9,
    # box(
    #   title = "Time-dependent ROC analysis",
    #   width = NULL,
    #   height = NULL,
    #   solidHeader = T,
    #   collapsible = F,
    #   status = "success",
    bsAlert("survrocmess"),
    bsCollapse(id = "collapseSurvROC", open = "Descriptions and parameters",
               bsCollapsePanel("Descriptions and parameters",  includeHTML("SurvROC1.html"), style = "default"),
               bsCollapsePanel("Time-dependent ROC analysis",style = "default",
     

      uiOutput("SurvROCplot", align="center"),
      fluidRow(column(4,align="left",
                      useShinyjs(),
                      hidden(div(id = "survROC_wrapper",
                                 splitLayout(
                                   numericInput("survROCwidth.dl","Figure width",value = 10),
                                   numericInput("survROCheight.dl","Figure height",value = 10)),
                                   downloadButton('downloadsurvROC', 'Download figure', class = "butt2")
                      ))
      ))
    ))
  ),
  column(
    3,
    box(
      title = "Survival ROC plot",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      selectizeInput(
        "SurvROCgene",
        label = "Official gene symbol",
        choices = NULL,
        multiple = F,
        selected= "TP53"
      ),
      bsTooltip("SurvROCgene", "Input a gene with official gene symbol", "left"),
      selectizeInput(
        "SurvROCtime",
        label = "Survival time",
        choices = NULL,
        multiple = F,
        selected= "OS.time"
      ),
      bsTooltip(
        "SurvROCtime",
        "Select survival time column. Example: OS.time, RFS.time, PFS.time",
        "left"
      ),

      selectizeInput(
        "SurvROCstatus",
        label = "Survival status",
        choices = NULL,
        multiple = F,
        selected= "OS"
      ),
      bsTooltip(
        "SurvROCstatus",
        "Select survival status column. Example: OS, RFS, PFS",
        "left"
      ),
      selectizeInput(
        "predictyear",
        label = "Prediction years",
        choices = NULL,
        multiple = T
      ),
      bsTooltip(
        "predictyear",
        "Define the time points in years you want to predict based on time dependent ROC analysis. the longest time point should not the max of survival (relapse) duration",
        "left"
      ),
      selectInput('SurvROCmethod', 'Method for survival ROC', choices =c("KM","NNE"), selected = "NNE"),
      checkboxInput("cutoff", "Show optimal cutoff ?", value = F, width = NULL),

      box(
        title = "Size control",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = T,
        sliderInput("SurvROCwidth","Plot Width (%)",min = 0,max = 100,value = 50),
        sliderInput("SurvROCheight",
          "Plot Height (px)",  min = 0,max = 1000,value = 410)
      ),
      useShinyjs(),
      actionButton(
        "SurvROCbt",
        "Submit",
        style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
        icon = icon("picture-o")

      )
     )
   )
  ),
  fluidRow(
    div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
        actionButton(inputId = 'page_before_SurvROC',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>CoxPH model</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Most correlated genes</i>'),
        actionButton(inputId = 'page_after_SurvROC',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )

 )
