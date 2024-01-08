tabItem(
  tabName = "valmo",
  
  fluidRow(column(9,
                  # box(
                  #   title = "Validation model",
                  #   width = NULL,
                  #   solidHeader = T,
                  #   collapsible = T,
                  #   status = "success",
                  #    bsAlert("valmess"),
                  #   align="center",
                  bsAlert("valmess"),
                  bsCollapse(id = "collapsevalmo", open = "Descriptions and parameters",
                             bsCollapsePanel("Descriptions and parameters",  includeHTML("valmo.html"), style = "default"),
                             bsCollapsePanel("Validation model",style = "default",
             tabBox(width = NULL,
               tabPanel("Survival ROC curve", align="center",uiOutput('valsurvROC'),
                        fluidRow(column(4,align="left",
                                        hidden(div(id = "valsurvROC_wrapper",
                                                   splitLayout(
                                                     numericInput("valsurvROCwidthdl","Figure width",value = 5),
                                                     numericInput("valsurvROCheightdl","Figure height",value = 5)),
                                                   downloadButton('savevalsurvROC', 'Download figure', class = "butt2")
                                        ))
                        ))
                        ),
               tabPanel("KM plot", align="center",uiOutput('valKMplot'),
                        fluidRow(column(4,align="left",
                                        hidden(div(id = "valKMplot_wrapper",
                                                   splitLayout(
                                                     numericInput("valKMplotwidthdl","Figure width",value = 5),
                                                     numericInput("valKMplotheightdl","Figure height",value = 5)),
                                                   downloadButton('savevalvalKMplot', 'Download figure', class = "butt2")
                                        ))
                        ))
                        ),
               tabPanel("Forestplot", align="center",uiOutput('valCoxforest'),
                        fluidRow(column(4,align="left",
                                        hidden(div(id = "valCoxforest_wrapper",
                                                   splitLayout(
                                                     numericInput("valCoxforestwidthdl","Figure width",value = 10),
                                                     numericInput("valCoxforestheightdl","Figure height",value = 5)),
                                                   downloadButton('savevalvalCoxforest', 'Download figure', class = "butt2")
                                        ))
                        ))
                        ),
               tabPanel("CoxPH table", align="center",DT::dataTableOutput('valCoxtable'),
                        fluidRow(column(4,align="left",
                                        hidden(div(id = "valCoxtable_wrapper",

                                                   downloadButton('savevalCoxtable', 'Download table', class = "butt2")
                                        ))
                        ))
             )
           )
    ))),
    column(3,
           box(
             title = "Valiation model",
             width = NULL,
             status = "danger",
             solidHeader = T,
             collapsible = T,
             selectizeInput(
               "Valtime",
               label = "Survival time",
               choices = NULL,
               multiple = F,
               selected = "OS.time"
             ),
             bsTooltip("Valtime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),

             selectizeInput(
               "Valstatus",
               label = "Survival status",
               choices = NULL,
               multiple = F,
               selected = "OS"
             ),
             bsTooltip("Valstatus", "Select survival status column. Example: OS, RFS, PFS","left"),

             box(title="Survival ROC setting",solidHeader = F,collapsible = T,width = NULL,collapsed = F,
                 selectizeInput(
                   "Valpredictyear",
                   label = "Prediction years",
                   choices = NULL,
                   multiple = T,
                   selected="1-year"
                 ),
                 bsTooltip(
                   "Valpredictyear",
                   "Define the time points in years you want to predict based on time dependent ROC analysis. the longest time point should not the max of survival (relapse) duration",
                   "left"
                 ),
                 selectInput('ValSurvROCmethod', 'Method for survival ROC', choices =c("KM","NNE"), selected = "NNE"),
                 checkboxInput("Valcutoff", "Show optimal cutoff ?", value = F, width = NULL),
                 box(title="Size control",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                 sliderInput("ValSurvROCwidth","Plot Width (%)",min = 0,max = 100, value = 45),
                 sliderInput("ValSurvROCheight","Plot Height (px)",min = 0,max = 1000,value = 400))
             ),

             box(title = "KM plot setting",solidHeader = F,collapsible = T,width = NULL,collapsed = T,

                 selectizeInput(
                   "Valkmgroupby",
                   label = "Group by",
                   choices = c("Percentage","Value"),
                   multiple = F,
                   selected = "Percentage"
                 ),
                 conditionalPanel(
                   condition = "input.Valkmgroupby == 'Percentage'",
                   sliderInput("Valriskcut", "Cutoff", min = 0, max = 1, value = 0.5),
                   bsTooltip("Valriskcut", "Input a cutoff (range:0-1) to categorize the samples into low and high risk groups.","left"),
                   ),
                 conditionalPanel(
                   condition = "input.Valkmgroupby == 'Value'",
                   numericInput('Valkmgpvalue', "Cutoff value", value=1,step = 1),
                   bsTooltip("Valkmgpvalue", "Input a specific cutoff value to categorize the samples into low and high risk group.","left")
                 ),

                 box(title = "Customized setting",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                 textInput("Valsurvxlab","Label of X axis",placeholder="Overall survival (months)",value="Overall survival (months)"),
                 bsTooltip("Valsurvxlab", "Define the label of X axis","left"),
                 checkboxInput("ValsurvP", "Show p-value?", value = TRUE, width = NULL),
                 checkboxInput("ValsurvRT", "Show risk table?", value = TRUE, width = NULL),
                 checkboxInput("ValsurvCI", "Show confidence interval?", value = F, width = NULL),
                 colourpicker::colourInput("Valkmcolor1", "Color of group 1", value = "#00AFBB", showColour = "background",closeOnClick = TRUE),
                 colourpicker::colourInput("Valkmcolor2", "Color of group 2", value = "#FC4E07", showColour = "background",closeOnClick = TRUE)),

                 box(title = "Size control",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                 sliderInput("Valsurvwidth", "Plot Width (%)", min = 0, max = 100, value = 45),
                 sliderInput("Valsurvheight", "Plot Height (px)", min = 0, max = 1200, value = 400))
             ),
             box(title = "CoxPH model",solidHeader = F,collapsible = T,width = NULL,collapsed = F,
                 selectizeInput("Valcoxclinvar",label = "Clinical Covariates" ,choices = NULL,multiple = T),
                 bsTooltip("Valcoxclinvar", "Select clinical variables to included to CoxPH model","left"),

                 box(title = "Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                     sliderInput("ValCoxwidth", "Plot Width (%)", min = 0, max = 100, value = 90),
                     sliderInput("ValCoxheight", "Plot Height (px)", min = 0, max = 1000, value = 200),
                     textInput("Valcoxfeat", "Variable names", placeholder="Age|Gender|Stage|Grade"),
                     bsTooltip("Valcoxfeat", "Define the variable names you interested,the number of variable names should be identical with the number of variables you included in the CoxPH model and use '|' to separate multiple variable names","left"),
                     numericInput("Valmaxtick", "maximum of xticks", 5, min = 1, max = 20,step=1),
                     bsTooltip("Valmaxtick", "Define the maximum of xticks and clip to adjust the size of the forestplot","left"),
                     sliderInput("valegend.pos", "Figure legend position", value = c(0.5, 0.9), min = 0, max = 1),
                     bsTooltip("valegend.pos", "Set the coordinates of the legend box. Their values should be between 0 and 1. c(0,0) corresponds to the 'bottom left' and c(1,1) corresponds to the 'top right' position","left")
                     )
             ),
             actionButton(
               "Valbt",
               "Submit",
               style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
               icon = icon("upload")
             )
           )
         )
       # )
     # )
  ),
  fluidRow(
    div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
        actionButton(inputId = 'page_before_valmo',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Prediction model: Construct model</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Prediction model: Nomogram</i>'),
        actionButton(inputId = 'page_after_valmo',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")

    )
   )
  )
