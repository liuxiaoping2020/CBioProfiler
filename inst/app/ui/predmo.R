tabItem(tabName = "bmpm",
        fluidRow(column(9,
                        box(
                          title = "Prediction model",
                          width = NULL,
                          solidHeader = T,
                          collapsible = F,
                          status = "success",
                          bsAlert("bmmess"),
                          align="center",
                          tabBox(width = NULL,
                            tabPanel("Survival ROC curve", uiOutput('bmsurvROC'),
                                     useShinyjs(),
                                     fluidRow(column(4,align="left",
                                                     hidden(div(id = "bmsurvROC_wrapper",
                                                                splitLayout(
                                                                    numericInput("bmsurvROCwidthdl","Figure width",value = 10),
                                                                    numericInput("bmsurvROCheightdl","Figure height",value = 10)),
                                                                downloadButton('savebmsurvROC', 'Download figure', class = "butt2")
                                                                ))
                                     ))
                                     ),
                            tabPanel("KM plot", uiOutput('bmKMplot'),
                                     fluidRow(column(4,align="left",
                                                     hidden(div(id = "bmKMplot_wrapper",
                                                                splitLayout(
                                                                    numericInput("bmKMplotwidthdl","Figure width",value = 10),
                                                                    numericInput("bmKMplotheightdl","Figure height",value = 10)),
                                                                downloadButton('savebmKMplot', 'Download figure', class = "butt2")
                                                     ))
                                     ))
                                     ),
                            tabPanel("Forestplot",uiOutput('bmCoxforest'),
                                     fluidRow(column(4,align="left",
                                                     hidden(div(id = "bmCoxforest_wrapper",
                                                                splitLayout(
                                                                    numericInput("bmCoxforestwidthdl","Figure width",value = 10),
                                                                    numericInput("bmCoxforestheightdl","Figure height",value = 10)),
                                                                downloadButton('savebmCoxforest', 'Download figure', class = "butt2")
                                                     ))
                                     ))
                                     ),
                            tabPanel("CoxPH table", dataTableOutput('bmCoxtable'),
                                     fluidRow(column(4,align="left",
                                                     hidden(div(id = "bmCoxtable_wrapper",
                                                                downloadButton('savebmCoxtable', 'Download table', class = "butt2")
                                                     ))
                                     ))
                                     )#,


                            )
                        )
        ),
        column(3,box(
          title = "Prediction model",
          width = NULL,
          status = "danger",
          solidHeader = T,
          collapsible = T,
          selectizeInput(
            "bmtime",
            label = "Survival time",
            choices = NULL,
            multiple = F,
            selected = "OS.time"
          ),
          bsTooltip("bmtime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),

          selectizeInput(
            "bmstatus",
            label = "Survival status",
            choices = NULL,
            multiple = F,
            selected = "OS"
          ),
          bsTooltip("bmstatus", "Select survival status column. Example: OS, RFS, PFS","left"),

          box(title="Survival ROC setting",solidHeader = F,collapsible = T,width = NULL,collapsed = F,
              selectizeInput(
                "bmpredictyear",
                label = "Prediction years",
                choices = NULL,
                multiple = T,
                selected="1-year"
              ),
              bsTooltip(
                "bmpredictyear",
                "Define the time points in years you want to predict based on time dependent ROC analysis. The longest time point should not exceed the max of survival (relapse) duration",
                "left"
              ),
              selectInput('bmSurvROCmethod', 'Method for survival ROC', choices =c("KM","NNE"), selected = "NNE"),
              checkboxInput("bmcutoff", "Show optimal cutoff ?", value = F, width = NULL),
              box(title = "Size control",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
              sliderInput("bmSurvROCwidth","Plot Width (%)",min = 0,max = 100, value = 50),
              sliderInput("bmSurvROCheight","Plot Height (px)",min = 0,max = 1000,value = 600))
              ),

          box(title = "KM plot setting",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
              selectizeInput(
                  "bmkmgroupby",
                  label = "Group by",
                  choices = c("Percentage","Value"),
                  multiple = F,
                  selected = "Percentage"
              ),
              conditionalPanel(
                  condition = "input.bmkmgroupby == 'Percentage'",
                  sliderInput("riskcut", "Cutoff percentage", min = 0, max = 1, value = 0.5),
                  bsTooltip("riskcut", "Input a cutoff (range:0-1) to categorize the samples into low and high risk group. Example if 0.25: 0%-25% = Low, 25%-100% high","left"),

              ),
              conditionalPanel(
                  condition = "input.bmkmgroupby == 'Value'",
                  numericInput('bmkmgpvalue', "Cutoff value", value=1,step = 1),
                  bsTooltip("bmkmgpvalue", "Input a specific cutoff value to categorize the samples into low and high risk group.","left")
              ),

          box(title="Customized setting",width=NULL,solidHeader=F,collapsible=T,status=NULL,collapsed =T,
          textInput("bmsurvxlab","Label of X axis",placeholder="Overall survival (months)",value="Overall survival (months)"),
          bsTooltip("bmsurvxlab", "Define the label of X axis","left"),
          checkboxInput("bmsurvP", "Show p-value?", value = TRUE, width = NULL),
          checkboxInput("bmsurvRT", "Show risk table?", value = TRUE, width = NULL),
          checkboxInput("bmsurvCI", "Show confidence interval?", value = F, width = NULL),
          colourpicker::colourInput("bmkmcolor1", "Color of group 1", value = "#00AFBB", showColour = "background",closeOnClick = TRUE),
          colourpicker::colourInput("bmkmcolor2", "Color of group 2", value = "#FC4E07", showColour = "background",closeOnClick = TRUE)),
          box(title="Size control",width=NULL,solidHeader=F,collapsible=T,status=NULL,collapsed =T,
          sliderInput("bmsurvwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
          sliderInput("bmsurvheight", "Plot Height (px)", min = 0, max = 1200, value = 650))
          ),
          box(title = "CoxPH model setting",solidHeader = F,collapsible = T,width = NULL,collapsed = F,
            selectizeInput("bmcoxclinvar",label = "Clinical covariates",choices = NULL,multiple = T),
            bsTooltip("bmcoxclinvar", "Select clinical variables to included to CoxPH model","left"),

            box(title = "Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                sliderInput("bmCoxwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
                sliderInput("bmCoxheight", "Plot Height (px)", min = 0, max = 2500, value = 550),
                textInput("bmcoxfeat", "Variable names", placeholder="Age|Gender|Stage|Grade"),
                bsTooltip("bmcoxfeat", "Define the variable names you interested,the number of variable names should be identical with the number of variables you included in the CoxPH model and use '|' to separate multiple variable names","left"),
                numericInput("bmmaxtick", "Maximum of xticks", 5, min = 1, max = 20,step=1),
                bsTooltip("bmmaxtick", "Define the maximum of xticks and clip to adjust the size of the forestplot","left"),
                sliderInput("bmlegend.pos", "Figure legend position", value = c(0.5, 0.9), min = 0, max = 1),
                bsTooltip("bmlegend.pos", "Set the coordinates of the legend box. Their values should be between 0 and 1. c(0,0) corresponds to the 'bottom left' and c(1,1) corresponds to the 'top right' position","left")

            )

          ),
          actionButton(
            "pmbt",
            "Submit",
            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
            icon = icon("picture-o")
          )

        ))
    ),

    fluidRow(
        div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
            actionButton(inputId = 'page_before_bmpm',label = '',icon = icon('arrow-left'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
            HTML('<i>Benchmark experiment</i>')
        ),
        div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
            HTML('<i>Prediction model: Validate model</i>'),
            actionButton(inputId = 'page_after_bmpm',label = '',icon = icon('arrow-right'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;")

        )
    )






    )




