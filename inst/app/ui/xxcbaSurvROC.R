tabItem (tabName = "cbaSurvROC1", 
         # fluidRow(column( 
         #   9,
         #   # bsAlert("cbasurvrocmess"),
         #   # bsCollapse(id = "cbacollapseSurvROC", open = "Descriptions and parameters",
         #   #            bsCollapsePanel("Descriptions and parameters",  includeHTML("cbaSurvROC.html"), style = "default"),
         #   #            bsCollapsePanel("Time-dependent ROC analysis",style = "default",
         #   #                            
         #   #                            
         #   #                            uiOutput("cbaSurvROCplot", align="center"),
         #   #                            fluidRow(column(4,align="left",
         #   #                                            useShinyjs(),
         #   #                                            hidden(div(id = "cbasurvROC_wrapper",
         #   #                                                       splitLayout(
         #   #                                                         numericInput("cbasurvROCwidth.dl","Figure width",value = 10),
         #   #                                                         numericInput("cbasurvROCheight.dl","Figure height",value = 10)),
         #   #                                                       downloadButton('cbadownloadsurvROC', 'Download figure', class = "butt2")
         #   #                                            ))
         #   #                            ))
         #   #            ))
         # ),
         # column(
         #   3,
         #   # box(
         #   #   title = "Survival ROC plot",
         #   #   width = NULL,
         #   #   status = "danger",
         #   #   solidHeader = T,
         #   #   collapsible = T,
         #   #   selectizeInput(
         #   #     "cbasurvROCcor",
         #   #     label = "Select the cohort",
         #   #     choices = c("Training set","Validation set"),
         #   #     selected="Training set",
         #   #     multiple = F
         #   #   ),
         #   #   bsTooltip("cbasurvROCcor", "Select the cohort (training set or validation set) to characterize the subtype","left"),
         #   #   
         #   #   selectizeInput(
         #   #     "cbaSurvROCtime",
         #   #     label = "Survival time",
         #   #     choices = NULL,
         #   #     multiple = F,
         #   #     selected= "OS.time"
         #   #   ),
         #   #   bsTooltip("cbaSurvROCtime","Select survival time column. Example: OS.time, RFS.time, PFS.time", "left"),
         #   #   
         #   #   selectizeInput(
         #   #     "cbaSurvROCstatus",
         #   #     label = "Survival status",
         #   #     choices = NULL,
         #   #     multiple = F,
         #   #     selected= "OS"
         #   #   ),
         #   #   bsTooltip(
         #   #     "cbaSurvROCstatus",
         #   #     "Select survival status column. Example: OS, RFS, PFS",
         #   #     "left"
         #   #   ),
         #   #   selectizeInput(
         #   #     "cbapredictyear",
         #   #     label = "Prediction years",
         #   #     choices = NULL,
         #   #     multiple = T
         #   #   ),
         #   #   bsTooltip(
         #   #     "cbapredictyear",
         #   #     "Define the time points in years you want to predict based on time dependent ROC analysis. the longest time point should not the max of survival (relapse) duration",
         #   #     "left"
         #   #   ),
         #   #   selectInput('cbaSurvROCmethod', 'Method for survival ROC', choices =c("KM","NNE"), selected = "NNE"),
         #   #   checkboxInput("cbacutoff", "Show optimal cutoff ?", value = F, width = NULL),
         #   #   
         #   #   box(
         #   #     title = "Size control",
         #   #     width = NULL,
         #   #     solidHeader = TRUE,
         #   #     collapsible = TRUE,
         #   #     collapsed = T,
         #   #     sliderInput("cbaSurvROCwidth","Plot Width (%)",min = 0,max = 100,value = 50),
         #   #     sliderInput("cbaSurvROCheight",
         #   #                 "Plot Height (px)",  min = 0,max = 1000,value = 410)
         #   #   ),
         #   #   useShinyjs(),
         #   #   actionButton(
         #   #     "cbaSurvROCbt",
         #   #     "Submit",
         #   #     style = "background-color: #000080;
         #   #                                  color: #FFFFFF;
         #   #                                  margin-left: auto;
         #   #                                  margin-right: auto;
         #   #                                  width: 100%",
         #   #     icon = icon("picture-o")
         #   #     
         #   #   )
         #   # )
         #   
         # )
         # ),
         # fluidRow(
         #   div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
         #       actionButton(inputId = 'page_before_cbaSurvROC',label = '',icon = icon('arrow-left'),
         #                    style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
         #       HTML('<i>CoxPH model</i>')
         #   ),
         #   div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
         #       HTML('<i>Differentially expressed genes</i>'),
         #       actionButton(inputId = 'page_after_cbaSurvROC',label = '',icon = icon('arrow-right'),
         #                    style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
         #   )
         # )
         
)
