tabItem(tabName = "cbaCoxPH",
        fluidRow(
          column(9,
                 # box(title="Cox proportional hazard regression model",width = NULL,
                 #       solidHeader=T,collapsible=T,status="success",#align="center",
                 bsAlert("cbaCoxPHmess"),
                 bsCollapse(id = "cbacollapseCoxPH", open = "Descriptions and parameters",
                            bsCollapsePanel("Descriptions and parameters",  includeHTML("cbaCoxPH.html"), style = "default"),
                            bsCollapsePanel("Cox proportional hazard regression model",style = "default",
                                            
                                            tabBox(width = NULL,
                                                   tabPanel("Forestplot", align="center",
                                                            
                                                            uiOutput('cbaCoxforestplot'),
                                                            useShinyjs(),
                                                            fluidRow(column(4,align="left",
                                                                            hidden(div(id = "cbamybox_wrapper",
                                                                                       splitLayout(
                                                                                         numericInput("cbaFPwidth","Figure width",value = 10),
                                                                                         numericInput("cbaFPheight","Figure height",value = 10)),
                                                                                       downloadButton('cbasaveforest', 'Download figure', class = "butt2")
                                                                            ))
                                                            ))
                                                   ),
                                                   
                                                   tabPanel("Table", DT::dataTableOutput('cbaCoxtable'),
                                                            useShinyjs(),
                                                            fluidRow(column(4,
                                                                            hidden(div(id = "cbamybox_wrapper.table",
                                                                                       downloadButton('cbasaveforesttable', 'Download table', class = "butt2")
                                                                            ))
                                                            ))
                                                   ),
                                                   tabPanel("Summary", verbatimTextOutput("cbaCoxsummary"),align="left",
                                                            
                                                            
                                                            tags$head(
                                                              tags$style(
                                                                "#Coxsummary{ font-size:12px; font-style:Arial;width: 1000px; max-width: 215%;background: ghostwhite;}"
                                                              )
                                                            )
                                                   )
                                            )
                            )
                 )
          ),
          column(3,box(
            title = "Cox proportional hazards regression model",
            width = NULL,
            status = "danger",
            solidHeader = T,
            collapsible = T,
            
            selectizeInput(
              "cbaCoxphcor",
              label = "Select the cohort",
              choices = c("Training set","Validation set"),
              selected="Training set",
              multiple = F
            ),
            bsTooltip("cbaCoxphcor", "Select the cohort (training set or validation set) to characterize the subtype","left"),
            selectizeInput(
              "cbacoxclinvar",
              label = "Select clinical features",
              choices = NULL,
              multiple = T
            ),
            bsTooltip("cbacoxclinvar", "Select clinical variables to included to CoxPH model","left"),
            # selectizeInput(
            #   "coxgene",
            #   label = "Official gene symbol",
            #   choices = NULL,
            #   multiple = T
            # ),
            # 
            # bsTooltip("coxgene", "Input one or more genes with official gene symbol that you want to included in the CoxPH model.","left"),
            
            selectizeInput("cbaCoxPHtime",label = "Survival time",choices = NULL,multiple = F,selected = "OS.time"),
            bsTooltip("cbaCoxPHtime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),
            
            selectizeInput("cbaCoxPHstatus",label = "Survival status",choices = NULL, multiple = F,selected = "OS"),
            bsTooltip("cbaCoxPHstatus", "Select survival status column. Example: OS, RFS, PFS","left"),
            
            box(title = "Size control",width = NULL,
                solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                sliderInput("cbaCoxwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
                sliderInput("cbaCoxheight", "Plot Height (px)", min = 0, max = 1000, value = 430),
                sliderInput("cbacoxfiglg", "Figure legend position", value = c(0.5, 0.9), min = 0, max = 1),
                bsTooltip("cbacoxfiglg", "Set the coordinates of the legend box. Their values should be between 0 and 1. c(0,0) corresponds to the 'bottom left' and c(1,1) corresponds to the 'top right' position","left"),
                textInput("cbacoxPHvarname", "Variable names", placeholder="Age,Gender,Stage,Grade..."),
                bsTooltip("cbacoxPHvarname", "Define the variable names you interested,the number of variable names should be identical with the number of variables you included in the CoxPH model and use '|' to separate multiple variable names","left"),
                numericInput("cbamaxtick", "Maximum of xticks", 5, min = 1, max = 20,step=1),
                bsTooltip("cbamaxtick", "Define the maximum of xticks and clip to adjust the size of the forestplot","left")
                
            ),
            
            useShinyjs(),
            actionButton("cbaCoxPHbt",
                         "Perform CoxPH model",
                         style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                         icon = icon("picture-o"))
          )
          
          )
        ),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_cbaCoxPH',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>Kaplan-Meier curve</i>')
          ),
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>Time-dependent ROC</i>'),
              actionButton(inputId = 'page_after_cbaCoxPH',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )
)
