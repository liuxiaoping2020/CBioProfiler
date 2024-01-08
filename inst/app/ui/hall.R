tabItem(
  tabName = "hallpath",
  fluidRow(
    column(9,
           bsAlert("hallpathmess"),
           bsCollapse(id = "collapsehallpath", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("hallmark.html"), style = "default"),
                      bsCollapsePanel("Correlation with hallmark signature",style = "default",
                                      tabBox(width = NULL,
                                             tabPanel("Hallmark signature plot", align="center",
                                                      uiOutput('hallpathplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "hallpathplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("hallpathpltwidth","Figure width",value = 15),
                                                                                   numericInput("hallpathpltheight","Figure height",value = 30)),
                                                                                 downloadButton('savehallpathness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             ),
                                             tabPanel("Hallmark signature table", DT::dataTableOutput('hallpathtable'),align="center",
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "hallpathtable_wrapper",
                                                                                 downloadButton('savehallpathtable', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             )
                                             
                                      )
                      ))
    ),
    column(3,box(
      title = "Correlation with hallmark signature",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      
      selectizeInput(
        "hallpathgene",
        label = "Official gene symbol",
        choices = NULL,
        multiple = F
      ),
      bsTooltip("hallpathgene", "Input one gene (biomarker candidate) with official gene symbol that you want to analyze the correlation between the gene and the hallmark signature","left"),
      selectizeInput(
        "hallpathmethod",
        label = "Correlation method",
        choices = c("pearson", "spearman"),
        multiple = F,
        selected="spearman"
      ),
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("hallpathwidth", "Plot Width (%)", min = 0, max = 100, value = 100),
          sliderInput("hallpathheight", "Plot Height (px)", min = 0, max = 1800, value = 1800)
      ),
      actionButton("hallpathbt",
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
        actionButton(inputId = 'page_before_hallpath',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Correlation with metabolism pathway</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Correlation with drug response</i>'),
        actionButton(inputId = 'page_after_hallpath',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)