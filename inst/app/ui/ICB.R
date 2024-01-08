tabItem(
  tabName = "ICB",
  fluidRow(
    column(9,
           bsAlert("ICBmess"),
           bsCollapse(id = "collapseICB", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("ICB.html"), style = "default"),
                      bsCollapsePanel("Correlation with immune checkpoint",style = "default",
                                      tabBox(width = NULL,
                                             tabPanel("Immune checkpoint plot", align="center",
                                                      uiOutput('ICBplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "ICBplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("ICBpltwidth","Figure width",value = 15),
                                                                                   numericInput("ICBpltheight","Figure height",value = 15)),
                                                                                 downloadButton('saveICBness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             ),
                                             tabPanel("Immune checkpoint table", DT::dataTableOutput('ICBtable'),align="center",
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "ICBtable_wrapper",
                                                                                 downloadButton('saveICBtable', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             )
                                             
                                      )
                      ))
    ),
    column(3,box(
      title = "Correlation with immune checkpoint",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,

      selectizeInput(
        "ICBgene",
        label = "Official gene symbol",
        choices = NULL,
        multiple = F
      ),
      bsTooltip("ICBgene", "Input one gene (biomarker candidate) with official gene symbol that you want to analyze the correlation between the gene and the ICB score","left"),
      selectizeInput(
        "ICBmethod",
        label = "Correlation method",
        choices = c("pearson", "spearman"),
        multiple = F,
        selected="spearman"
      ),

      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("ICBwidth", "Plot Width (%)", min = 0, max = 100, value = 90),
          sliderInput("ICBheight", "Plot Height (px)", min = 0, max = 1000, value = 900)
      ),
      actionButton("ICBbt",
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
        actionButton(inputId = 'page_before_ICB',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Correlation with ESTIMATE score</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Correlation with IFN-gamma score</i>'),
        actionButton(inputId = 'page_after_ICB',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)