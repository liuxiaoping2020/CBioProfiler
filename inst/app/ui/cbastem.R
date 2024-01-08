tabItem(
  tabName = "cbastemness",
  fluidRow(
    column(9,
           bsAlert("cbastemnessmess"),
           bsCollapse(id = "cbacollapsestemness", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("cbastem.html"), style = "default"),
                      bsCollapsePanel("Correlation with stemness score",style = "default",
                                      # tabBox(width = NULL,
                                             # tabPanel("Correlation plot", align="center",
                                                      uiOutput('cbastemnessplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "cbastemnessplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("cbastempltwidth","Figure width",value = 5),
                                                                                   numericInput("cbastempltheight","Figure height",value = 5)),
                                                                                 downloadButton('cbasavestemness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             # )
                                             # ,
                                             # 
                                             # tabPanel("Stemness table", dataTableOutput('cbastemtable'),align="center",
                                             #          useShinyjs(),
                                             #          fluidRow(column(4,
                                             #                          hidden(div(id = "cbastemtable_wrapper",
                                             #                                     downloadButton('cbasavestemtable', 'Download table', class = "butt2")
                                             #                          ))
                                             #          ))
                                             # )
                                             
                                      # )
                      ))
    ),
    column(3,box(
      title = "Correlation with Stemness score",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      selectizeInput(
        "cbastemcor",
        label = "Select the cohort",
        choices = c("Training set","Validation set"),
        selected="Training set",
        multiple = F
      ),
      bsTooltip("cbastemcor", "Select the cohort (training set or validation set) to characterize the subtype","left"),
      
      selectizeInput(
        "cbastemmeth",
        label = "Group comparison method",
        choices = c("Nonparametric test","Parametric test"),
        multiple = F,
        selected="Nonparametric test"
      ),
      bsTooltip("cbastemmeth", "Whether to use a parametric or nonparametric test for comparisons between groups","left"),
      
      checkboxInput("cbastempair", "Perform paire-wised comparision ?", value = F, width = NULL),
      selectizeInput(
        "cbastemplot",
        label = "Plot type",
        choices = c("Bar plot", "Box plot", "Violin plot"),
        multiple = F,
        selected="Box plot"
      ),
      
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("cbastemwidth", "Plot Width (%)", min = 0, max = 100, value = 45),
          sliderInput("cbastemheight", "Plot Height (px)", min = 0, max = 1000, value = 350)
      ),
      actionButton("cbastembt",
                   "Submit",
                   style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                   icon = icon("upload"))
    )
   )
  ),
  fluidRow(
    div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
        actionButton(inputId = 'page_before_cbastemness',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Immune infiltration among subtypes</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>ESTIMATE score among subtypes</i>'),
        actionButton(inputId = 'page_after_cbastemness',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)
