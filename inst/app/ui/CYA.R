tabItem(
  tabName = "CYA",
  fluidRow(
    column(9,
           bsAlert("CYAmess"),
           bsCollapse(id = "collapseCYA", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("CYA.html"), style = "default"),
                      bsCollapsePanel("Correlation with cytotoxic activity",style = "default",
                                      tabBox(width = NULL,
                                             tabPanel("Cytotoxic activity plot", align="center",
                                                      uiOutput('CYAplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "CYAplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("CYApltwidth","Figure width",value = 10),
                                                                                   numericInput("CYApltheight","Figure height",value = 10)),
                                                                                 downloadButton('saveCYAness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             ),
                                             tabPanel("Cytotoxic activity table", DT::dataTableOutput('CYAtable'),align="center",
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "CYAtable_wrapper",
                                                                                 downloadButton('saveCYAtable', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             )
                                             
                                      )
                      ))
    ),
    column(3,box(
      title = "Correlation with cytotoxic activity",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      
      selectizeInput(
        "CYAgene",
        label = "Official gene symbol",
        choices = NULL,
        multiple = F
      ),
      bsTooltip("CYAgene", "Input one gene (biomarker candidate) with official gene symbol that you want to analyze the correlation between the gene and the cytotoxic activity","left"),
      selectizeInput(
        "CYAmethod",
        label = "Correlation method",
        choices = c("pearson", "spearman"),
        multiple = F,
        selected="spearman"
      ),
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("CYAwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
          sliderInput("CYAheight", "Plot Height (px)", min = 0, max = 1000, value = 430)
      ),
      actionButton("CYAbt",
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
        actionButton(inputId = 'page_before_CYA',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Correlation with IFN-gamma score</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Correlation with cancer pathway</i>'),
        actionButton(inputId = 'page_after_CYA',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)