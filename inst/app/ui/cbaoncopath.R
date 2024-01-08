tabItem(
  tabName = "cbaoncopath",
  fluidRow(
    column(9,
           bsAlert("cbaoncopathmess"),
           bsCollapse(id = "cbacollapseoncopath", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("cbaoncopath.html"), style = "default"),
                      bsCollapsePanel("Comparision of cancer pathway",style = "default",
                                      tabBox(width = NULL,
                                             tabPanel("Cancer pathway plot", align="center",
                                                      uiOutput('cbaoncopathplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "cbaoncopathplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("cbaoncopathpltwidth","Figure width",value = 10),
                                                                                   numericInput("cbaoncopathpltheight","Figure height",value = 10)),
                                                                                 downloadButton('cbasaveoncopathness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             ),
                                             tabPanel("Cancer pathway table", DT::dataTableOutput('cbaoncopathtable'),align="center",
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "cbaoncopathtable_wrapper",
                                                                                 downloadButton('cbasaveoncopathtable', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             )
                                             
                                      )
                      ))
    ),
    column(3,box(
      title = "Comparision of cancer pathway",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      selectizeInput(
        "cbaoncopathcor",
        label = "Select the cohort",
        choices = c("Training set","Validation set"),
        selected="Training set",
        multiple = F
      ),
      bsTooltip("cbaoncopathcor", "Select the cohort (training set or validation set) to characterize the subtype","left"),
      
      selectizeInput(
        "cbaoncopathmeth",
        label = "Group comparison method",
        choices = c("Nonparametric test","Parametric test"),
        multiple = F,
        selected="Nonparametric test"
      ),
      bsTooltip("cbaoncopathmeth", "Whether to use a parametric or nonparametric test for comparisons between groups","left"),
      
      checkboxInput("cbaoncopathpair", "Perform paire-wised comparision ?", value = F, width = NULL),
      selectizeInput(
        "cbaoncopathplot",
        label = "Plot type",
        choices = c("Bar plot", "Box plot", "Violin plot"),
        multiple = F,
        selected="Box plot"
      ),
      
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("cbaoncopathwidth", "Plot Width (%)", min = 0, max = 100, value = 80),
          sliderInput("cbaoncopathheight", "Plot Height (px)", min = 0, max = 1000, value = 800)
      ),
      actionButton("cbaoncopathbt",
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
        actionButton(inputId = 'page_before_cbaoncopath',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Cytolytic activity among subtypes</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Metabolisim score among subtypes</i>'),
        actionButton(inputId = 'page_after_cbaoncopath',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)