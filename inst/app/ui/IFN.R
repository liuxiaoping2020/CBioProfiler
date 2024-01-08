tabItem(
  tabName = "IFN",
  fluidRow(
    column(9,
           bsAlert("IFNmess"),
           bsCollapse(id = "collapseIFN", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("IFN.html"), style = "default"),
                      bsCollapsePanel("Correlation with interferon-gamma score",style = "default",
                                      tabBox(width = NULL,
                                             tabPanel("Interferon-gamma score plot", align="center",
                                                      uiOutput('IFNplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "IFNplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("IFNpltwidth","Figure width",value = 10),
                                                                                   numericInput("IFNpltheight","Figure height",value = 10)),
                                                                                 downloadButton('saveIFNness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             ),
                                             tabPanel("Interferon-gamma score table", DT::dataTableOutput('IFNtable'),align="center",
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "IFNtable_wrapper",
                                                                                 downloadButton('saveIFNtable', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             )
                                             
                                      )
                      ))
    ),
    column(3,box(
      title = "Correlation with interferon-gamma score",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      
      selectizeInput(
        "IFNgene",
        label = "Official gene symbol",
        choices = NULL,
        multiple = F
      ),
      bsTooltip("IFNgene", "Input one gene (biomarker candidate) with official gene symbol that you want to analyze the correlation between the gene and the IFN score","left"),
      selectizeInput(
        "IFNmethod",
        label = "Correlation method",
        choices = c("pearson", "spearman"),
        multiple = F,
        selected="spearman"
      ),
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("IFNwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
          sliderInput("IFNheight", "Plot Height (px)", min = 0, max = 1000, value = 430)
      ),
      actionButton("IFNbt",
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
        actionButton(inputId = 'page_before_IFN',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Correlation with stemness score</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Biological annotationf</i>'),
        actionButton(inputId = 'page_after_IFN',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)