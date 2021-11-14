tabItem(
  tabName = "stemness",
  fluidRow(
    column(9,box(
      title="Correlation with stemness score",width = NULL,
      solidHeader=T,collapsible=T,status="success",collapsed = F,
      bsAlert("stemnessmess"),
      tabBox(width = NULL,
             tabPanel("Correlation plot", align="center",
                      uiOutput('stemnessplot'
                      )
                      ,
                      useShinyjs(),
                      fluidRow(column(4,align="left",
                                      hidden(div(id = "stemnessplot_wrapper",
                                                 splitLayout(
                                                   numericInput("stempltwidth","Figure width",value = 10),
                                                   numericInput("stempltheight","Figure height",value = 10)),
                                                 downloadButton('savestemness', 'Download figure', class = "butt2")
                                      ))
                      ))
             ),

             tabPanel("Stemness table", dataTableOutput('stemtable'),
                      useShinyjs(),
                      fluidRow(column(4,
                                      hidden(div(id = "stemtable_wrapper",
                                                 downloadButton('savestemtable', 'Download table', class = "butt2")
                                      ))
                      ))
             )

      )
    )
    ),
    column(3,box(
      title = "Correlation with Stemness score",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
     selectizeInput(
        "stemgene",
        label = "Official gene symbol",
        choices = NULL,
        multiple = F
      ),
      bsTooltip("stemgene", "Input one gene (biomarker candidate) with official gene symbol that you want to analyze the correlation between the gene and tummor microenvironment ","left"),
      selectizeInput(
        "stemmethod",
        label = "Correlation method",
        choices = c("pearson", "spearman"),
        multiple = F,
        selected="spearman"
      ),
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("stemwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
          sliderInput("stemheight", "Plot Height (px)", min = 0, max = 1000, value = 430)
      ),
      actionButton("stembt",
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
        actionButton(inputId = 'page_before_stemness',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Correlation with Immune infiltration</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Biological annotationf</i>'),
        actionButton(inputId = 'page_after_stemness',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )




)

