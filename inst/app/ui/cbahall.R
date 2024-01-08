tabItem(
  tabName = "cbahall",
  fluidRow(
    column(9,
           bsAlert("cbahallpathmess"),
           bsCollapse(id = "cbacollapsehallpath", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("cbahall.html"), style = "default"),
                      bsCollapsePanel("Comparision of hallmark signature",style = "default",
                                      tabBox(width = NULL,
                                             tabPanel("Hallmark signature plot", align="center",
                                                      uiOutput('cbahallpathplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "cbahallpathplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("cbahallpathpltwidth","Figure width",value = 15),
                                                                                   numericInput("cbahallpathpltheight","Figure height",value = 30)),
                                                                                 downloadButton('cbasavehallpathness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             ),
                                             tabPanel("Hallmark signature table", DT::dataTableOutput('cbahallpathtable'),align="center",
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "cbahallpathtable_wrapper",
                                                                                 downloadButton('cbasavehallpathtable', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             )
                                             
                                      )
                      ))
    ),
    column(3,box(
      title = "Comparision of hallmark signature",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      selectizeInput(
        "cbahallpathcor",
        label = "Select the cohort",
        choices = c("Training set","Validation set"),
        selected="Training set",
        multiple = F
      ),
      bsTooltip("cbahallpathcor", "Select the cohort (training set or validation set) to characterize the subtype","left"),
      
      selectizeInput(
        "cbahallpathmeth",
        label = "Group comparison method",
        choices = c("Nonparametric test","Parametric test"),
        multiple = F,
        selected="Nonparametric test"
      ),
      bsTooltip("cbahallpathmeth", "Whether to use a parametric or nonparametric test for comparisons between groups","left"),
      
      checkboxInput("cbahallpathpair", "Perform paire-wised comparision ?", value = F, width = NULL),
      selectizeInput(
        "cbahallpathplot",
        label = "Plot type",
        choices = c("Bar plot", "Box plot", "Violin plot"),
        multiple = F,
        selected="Box plot"
      ),
      
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("cbahallpathwidth", "Plot Width (%)", min = 0, max = 100, value = 100),
          sliderInput("cbahallpathheight", "Plot Height (px)", min = 0, max = 1800, value = 1800)
      ),
      actionButton("cbahallpathbt",
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
        actionButton(inputId = 'page_before_cbahallpath',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Metabolisim score among subtypes</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Drug response among subtypes</i>'),
        actionButton(inputId = 'page_after_cbahallpath',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)