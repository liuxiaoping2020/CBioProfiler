tabItem(
  tabName = "cbametapath",
  fluidRow(
    column(9,
           bsAlert("cbametapathmess"),
           bsCollapse(id = "cbacollapsemetapath", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("cbametapath.html"), style = "default"),
                      bsCollapsePanel("Comparision of metabolisim pathway",style = "default",
                                      tabBox(width = NULL,
                                             tabPanel("Metabolism pathway plot", align="center",
                                                      uiOutput('cbametapathplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "cbametapathplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("cbametapathpltwidth","Figure width",value = 7),
                                                                                   numericInput("cbametapathpltheight","Figure height",value = 10)),
                                                                                 downloadButton('cbasavemetapathness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             ),
                                             tabPanel("Metabolisim pathway table", DT::dataTableOutput('cbametapathtable'),align="center",
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "cbametapathtable_wrapper",
                                                                                 downloadButton('cbasavemetapathtable', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             )
                                             
                                      )
                      ))
    ),
    column(3,box(
      title = "Comparision of metabolisim pathway",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      selectizeInput(
        "cbametapathcor",
        label = "Select the cohort",
        choices = c("Training set","Validation set"),
        selected="Training set",
        multiple = F
      ),
      bsTooltip("cbametapathcor", "Select the cohort (training set or validation set) to characterize the subtype","left"),
      
      selectizeInput(
        "cbametapathmeth",
        label = "Group comparison method",
        choices = c("Nonparametric test","Parametric test"),
        multiple = F,
        selected="Nonparametric test"
      ),
      bsTooltip("cbametapathmeth", "Whether to use a parametric or nonparametric test for comparisons between groups","left"),
      
      checkboxInput("cbametapathpair", "Perform paire-wised comparision ?", value = F, width = NULL),
      selectizeInput(
        "cbametapathplot",
        label = "Plot type",
        choices = c("Bar plot", "Box plot", "Violin plot"),
        multiple = F,
        selected="Box plot"
      ),
      
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("cbametapathwidth", "Plot Width (%)", min = 0, max = 100, value = 60),
          sliderInput("cbametapathheight", "Plot Height (px)", min = 0, max = 1000, value = 800)
      ),
      actionButton("cbametapathbt",
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
        actionButton(inputId = 'page_before_cbametapath',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Cancer pathway score among subtypes</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Hallmark signature among subtypes</i>'),
        actionButton(inputId = 'page_after_cbametapath',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)