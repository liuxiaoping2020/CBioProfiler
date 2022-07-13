tabItem(
  tabName = "cbaICB",
  fluidRow(
    column(9,
           bsAlert("cbaICBmess"),
           bsCollapse(id = "cbacollapseICB", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("cbaICB.html"), style = "default"),
                      bsCollapsePanel("Comparision of immune checkpoint distribtion",style = "default",
                                      tabBox(width = NULL,
                                             tabPanel("Immune checkpoint plot", align="center",
                                                      uiOutput('cbaICBplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "cbaICBplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("cbaICBpltwidth","Figure width",value = 10),
                                                                                   numericInput("cbaICBpltheight","Figure height",value = 10)),
                                                                                 downloadButton('cbasaveICBness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             ),
                                             tabPanel("Immune checkpoint table", DT::dataTableOutput('cbaICBtable'),align="center",
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "cbaICBtable_wrapper",
                                                                                 downloadButton('cbasaveICBtable', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             )
                                             
                                      )
                      ))
    ),
    column(3,box(
      title = "Correlation with ICB score",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      selectizeInput(
        "cbaICBcor",
        label = "Select the cohort",
        choices = c("Training set","Validation set"),
        selected="Training set",
        multiple = F
      ),
      bsTooltip("cbaICBcor", "Select the cohort (training set or validation set) to characterize the subtype","left"),
      
      selectizeInput(
        "cbaICBmeth",
        label = "Group comparison method",
        choices = c("Nonparametric test","Parametric test"),
        multiple = F,
        selected="Nonparametric test"
      ),
      bsTooltip("cbaICBmeth", "Whether to use a parametric or nonparametric test for comparisons between groups","left"),
      
      checkboxInput("cbaICBpair", "Perform paire-wised comparision ?", value = F, width = NULL),
      selectizeInput(
        "cbaICBplot",
        label = "Plot type",
        choices = c("Bar plot", "Box plot", "Violin plot"),
        multiple = F,
        selected="Box plot"
      ),
      
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("cbaICBwidth", "Plot Width (%)", min = 0, max = 100, value = 80),
          sliderInput("cbaICBheight", "Plot Height (px)", min = 0, max = 1000, value = 800)
      ),
      actionButton("cbaICBbt",
                   "Submit",
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
        actionButton(inputId = 'page_before_cbaICB',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>ESTIMATE score among subtypes</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Interferon-gamma among subtypes</i>'),
        actionButton(inputId = 'page_after_cbaICB',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)