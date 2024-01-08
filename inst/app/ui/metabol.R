tabItem(
  tabName = "metapath",
  fluidRow(
    column(9,
           bsAlert("metapathmess"),
           bsCollapse(id = "collapsemetapath", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("metapath.html"), style = "default"),
                      bsCollapsePanel("Correlation with metabolism pathway",style = "default",
                                      tabBox(width = NULL,
                                             tabPanel("Metabolism pathway plot", align="center",
                                                      uiOutput('metapathplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "metapathplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("metapathpltwidth","Figure width",value = 10),
                                                                                   numericInput("metapathpltheight","Figure height",value = 15)),
                                                                                 downloadButton('savemetapathness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             ),
                                             tabPanel("Metabolism pathway table", DT::dataTableOutput('metapathtable'),align="center",
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "metapathtable_wrapper",
                                                                                 downloadButton('savemetapathtable', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             )
                                             
                                      )
                      ))
    ),
    column(3,box(
      title = "Correlation with metabolism pathway",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      
      selectizeInput(
        "metapathgene",
        label = "Official gene symbol",
        choices = NULL,
        multiple = F
      ),
      bsTooltip("metapathgene", "Input one gene (biomarker candidate) with official gene symbol that you want to analyze the correlation between the gene and the metabolism pathway","left"),
      selectizeInput(
        "metapathmethod",
        label = "Correlation method",
        choices = c("pearson", "spearman"),
        multiple = F,
        selected="spearman"
      ),
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("metapathwidth", "Plot Width (%)", min = 0, max = 100, value = 80),
          sliderInput("metapathheight", "Plot Height (px)", min = 0, max = 1000, value = 1000)
      ),
      actionButton("metapathbt",
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
        actionButton(inputId = 'page_before_metapath',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Correlation with stemness score</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Biological annotationf</i>'),
        actionButton(inputId = 'page_after_metapath',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)