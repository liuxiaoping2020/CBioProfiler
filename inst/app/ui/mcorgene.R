tabItem(
  tabName = "mcorgene",
  fluidRow(column(9,
# box(
#     title = "Most correlated genes",
#     width = NULL,
#     solidHeader = T,
#     collapsible = T,
#     status = "success",
    bsAlert("mcorgenemess"),
    
bsCollapse(id = "collapsemcorgene", open = "Descriptions and parameters",
           bsCollapsePanel("Descriptions and parameters",  includeHTML("mcorgene.html"), style = "default"),
           bsCollapsePanel("Most correlated genes",style = "default",
    tabBox(width = NULL,
      tabPanel("Bar plot", align="center",uiOutput('mcorbarplot'),
               useShinyjs(),
               fluidRow(column(4,align="left",
                               hidden(div(id = "mcorgenebarplot_wrapper",

                                          splitLayout(

                                            numericInput("mcorgenebarwidthdl","Figure width",value = 10),
                                            numericInput("mcorgenebarheightdl","Figure width",value = 10)),
                                          downloadButton('downmcorgenebarplot', 'Download figure', class = "butt2")
                               ))
               ))
               ),
      tabPanel("Bubble plot",align="center", uiOutput('mcorbubplot') ,
               useShinyjs(),
               fluidRow(column(4,align="left",
                               hidden(div(id = "mcorgenebubbleplot_wrapper",
                                          splitLayout(
                                           numericInput("mcorgenebubblewidthdl","Figure width",value = 10),
                                           numericInput("mcorgenebubbleheightdl","Figure height",value = 10)),
                                          downloadButton('downmcorgenebubbleplot', 'Download figure', class = "butt2")
                               ))
               ))


               ),
      tabPanel("Table",align="center", DT::dataTableOutput('mcortable'),
               useShinyjs(),
               fluidRow(column(4,align="left",
                               hidden(div(id = "mcorgenetable_wrapper",
                                          downloadButton('downmcorgenetable', 'Download table', class = "butt2")
                               ))
               ))


               )
    )
)
  )),
  column(3,(box(
    title = "Most correlated gene analysis",
    width = NULL,
    status = "danger",
    solidHeader = T,
    collapsible = T,
    selectInput('mcorgene', 'Official gene symbol', choices = NULL, selected = NULL),
    bsTooltip("mcorgene", "Select an official gene symbol you interested in ","left"),
    selectizeInput('mcormethod', "Correlation method", choices =c("pearson","spearman","kendall") , selected = "pearson"),
    bsTooltip("mcormethod", "Specify the correction method","left"),
    checkboxInput("padjust", "Adjust P value ?", value = TRUE, width = NULL),
    conditionalPanel(
      condition = "input.padjust == true ",
      selectInput('mcorrectmethod', "Correction method", choices =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr") , selected = "fdr"),
      bsTooltip("mcormethod", "Specify the correction method for P values. For more details, please refer to https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust","left")
     ),
    numericInput("corsigcut","Significance criterion", value=0.05,max =0.1, step =0.001),
    bsTooltip("corsigcut", "Specify the criterion for significant correlation ","left"),
    box(title = "Siz control",width = NULL,#status = "info",
          solidHeader = F,collapsible = TRUE,collapsed=TRUE,
          box(title = "Bar plot",width = NULL,collapsible = TRUE,collapsed=TRUE,
            sliderInput("mcorbarwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
          sliderInput("mcorbarheight", "Plot Height (px)", min = 0, max = 1000, value = 430)),
        box(title = "Bubble plot",width = NULL,collapsible = TRUE,collapsed=TRUE,
            sliderInput("mcorbubwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
            sliderInput("mcorbubheight", "Plot Height (px)", min = 0, max = 1000, value = 430))),
    useShinyjs(),

    actionButton("mcorgenebt",
                 "Submit",
                 style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                 icon = icon("picture-o"))#,

    )
  )

   )
  ),

  fluidRow(
    div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
        actionButton(inputId = 'page_before_mcorgene',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Time-dependent ROC</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Correlation with specific gene</i>'),
        actionButton(inputId = 'page_after_mcorgene',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
  )



