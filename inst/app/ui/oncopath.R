tabItem(
  tabName = "oncopath",
  fluidRow(
    column(9,
           bsAlert("oncopathmess"),
           bsCollapse(id = "collapseoncopath", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("oncopath.html"), style = "default"),
                      bsCollapsePanel("Correlation with cancer pathway",style = "default",
                                      tabBox(width = NULL,
                                             tabPanel("Cancer pathway plot", align="center",
                                                      uiOutput('oncopathplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "oncopathplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("oncopathpltwidth","Figure width",value = 10),
                                                                                   numericInput("oncopathpltheight","Figure height",value = 10)),
                                                                                 downloadButton('saveoncopathness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             ),
                                             tabPanel("Cancer pathway table", DT::dataTableOutput('oncopathtable'),align="center",
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "oncopathtable_wrapper",
                                                                                 downloadButton('saveoncopathtable', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             )
                                             
                                      )
                      ))
    ),
    column(3,box(
      title = "Correlation with cancer pathway",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      
      selectizeInput(
        "oncopathgene",
        label = "Official gene symbol",
        choices = NULL,
        multiple = F
      ),
      bsTooltip("oncopathgene", "Input one gene (biomarker candidate) with official gene symbol that you want to analyze the correlation between the gene and the cancer pathway","left"),
      selectizeInput(
        "oncopathmethod",
        label = "Correlation method",
        choices = c("pearson", "spearman"),
        multiple = F,
        selected="spearman"
      ),
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("oncopathwidth", "Plot Width (%)", min = 0, max = 100, value = 80),
          sliderInput("oncopathheight", "Plot Height (px)", min = 0, max = 1000, value = 800)
      ),
      actionButton("oncopathbt",
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
        actionButton(inputId = 'page_before_oncopath',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Correlation with cytolytic activity</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Correlation with hallmark signature</i>'),
        actionButton(inputId = 'page_after_oncopath',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)