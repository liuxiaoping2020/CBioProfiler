tabItem(tabName = "CoxPH",
        fluidRow(
          column(9,
                 # box(title="Cox proportional hazard regression model",width = NULL,
                 #       solidHeader=T,collapsible=T,status="success",#align="center",
                     bsAlert("CoxPHmess"),
                 bsCollapse(id = "collapseCoxPH", open = "Descriptions and parameters",
                            bsCollapsePanel("Descriptions and parameters",  includeHTML("CoxPH.html"), style = "default"),
                            bsCollapsePanel("Cox proportional hazard regression model",style = "default",
                                              
            tabBox(width = NULL,
            tabPanel("Forestplot", align="center",

                     uiOutput('Coxforestplot'),
                     useShinyjs(),
                        fluidRow(column(4,align="left",
                      hidden(div(id = "mybox_wrapper",
                       splitLayout(
                       numericInput("FPwidth","Figure width",value = 10),
                       numericInput("FPheight","Figure height",value = 10)),
                       downloadButton('saveforest', 'Download figure', class = "butt2")
                        ))
                       ))
                     ),

            tabPanel("Table", dataTableOutput('Coxtable'),
                     useShinyjs(),
                     fluidRow(column(4,
                                     hidden(div(id = "mybox_wrapper.table",
                                                downloadButton('saveforesttable', 'Download table', class = "butt2")
                                     ))
                     ))
                     ),
            tabPanel("Summary", verbatimTextOutput("Coxsummary"),align="left",


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
                  "coxclinvar",
                  label = "Select clinical features",
                  choices = NULL,
                  multiple = T
                ),
            bsTooltip("coxclinvar", "Select clinical variables to included to CoxPH model","left"),
            selectizeInput(
              "coxgene",
              label = "Official gene symbol",
              choices = NULL,
              multiple = T
            ),

            bsTooltip("coxgene", "Input one or more genes with official gene symbol that you want to included in the CoxPH model.","left"),

            selectizeInput("CoxPHtime",label = "Survival time",choices = NULL,multiple = F,selected = "OS.time"),
            bsTooltip("CoxPHtime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),

            selectizeInput("CoxPHstatus",label = "Survival status",choices = NULL, multiple = F,selected = "OS"),
            bsTooltip("CoxPHstatus", "Select survival status column. Example: OS, RFS, PFS","left"),

            box(title = "Size control",width = NULL,
                solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                sliderInput("Coxwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
                sliderInput("Coxheight", "Plot Height (px)", min = 0, max = 1000, value = 430),
                sliderInput("coxfiglg", "Figure legend position", value = c(0.5, 0.9), min = 0, max = 1),
                bsTooltip("coxfiglg", "Set the coordinates of the legend box. Their values should be between 0 and 1. c(0,0) corresponds to the 'bottom left' and c(1,1) corresponds to the 'top right' position","left"),
                textInput("coxPHvarname", "Variable names", placeholder="Age,Gender,Stage,Grade..."),
                bsTooltip("coxPHvarname", "Define the variable names you interested,the number of variable names should be identical with the number of variables you included in the CoxPH model and use '|' to separate multiple variable names","left"),
                numericInput("maxtick", "Maximum of xticks", 5, min = 1, max = 20,step=1),
                bsTooltip("maxtick", "Define the maximum of xticks and clip to adjust the size of the forestplot","left")

                ),

             useShinyjs(),
              actionButton("CoxPHbt",
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
             actionButton(inputId = 'page_before_CoxPH',label = '',icon = icon('arrow-left'),
                          style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
             HTML('<i>Kaplan-Meier curve</i>')
         ),
         div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
             HTML('<i>Time-dependent ROC</i>'),
             actionButton(inputId = 'page_after_CoxPH',label = '',icon = icon('arrow-right'),
                          style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
         )
       )
)
