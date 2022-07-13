tabItem(
  tabName = "estimate",
  fluidRow(
    column(9,
           bsAlert("estimatemess"),
           bsCollapse(id = "collapseestimate", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("Estimate.html"), style = "default"),
                      bsCollapsePanel("Correlation with ESTIMATE score",style = "default",
                                      tabBox(width = NULL,
                                             tabPanel("ESTIMATE plot", align="center",
                                                      uiOutput('estimateplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "estimateplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("estimatepltwidth","Figure width",value = 10),
                                                                                   numericInput("estimatepltheight","Figure height",value = 10)),
                                                                                 downloadButton('saveestimateness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             ),
                                             tabPanel("ESTIMATE table", DT::dataTableOutput('estimatetable'),align="center",
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "estimatetable_wrapper",
                                                                                 downloadButton('saveestimatetable', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             )
                                             
                                      )
                      ))
    ),
    column(3,box(
      title = "Correlation with ESTIMATE score",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      # selectizeInput(
      #   "cbaestimatecor",
      #   label = "Select the cohort",
      #   choices = c("Training set","Validation set"),
      #   selected="Training set",
      #   multiple = F
      # ),
      # bsTooltip("cbaestimate", "Select the cohort (training set or validation set) to characterize the subtype","left"),
      selectizeInput(
        "estimategene",
        label = "Official gene symbol",
        choices = NULL,
        multiple = F
      ),
      bsTooltip("estimategene", "Input one gene (biomarker candidate) with official gene symbol that you want to analyze the correlation between the gene and the ESTIMATE score","left"),
      selectizeInput(
        "estimatemethod",
        label = "Correlation method",
        choices = c("pearson", "spearman"),
        multiple = F,
        selected="spearman"
      ),
      
      # selectizeInput(
      #   "estimatemeth",
      #   label = "Group comparison method",
      #   choices = c("Nonparametric test","Parametric test"),
      #   multiple = F,
      #   selected="Nonparametric test"
      # ),
      # bsTooltip("cbaestimatemeth", "Whether to use a parametric or nonparametric test for comparisons between groups","left"),
      # 
      # checkboxInput("cbaestimatepair", "Perform paire-wised comparision ?", value = F, width = NULL),
      # selectizeInput(
      #   "cbaestimateplot",
      #   label = "Plot type",
      #   choices = c("Bar plot", "Box plot", "Violin plot"),
      #   multiple = F,
      #   selected="Box plot"
      # ),
      # 
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("estimatewidth", "Plot Width (%)", min = 0, max = 100, value = 50),
          sliderInput("estimateheight", "Plot Height (px)", min = 0, max = 1000, value = 430)
      ),
      actionButton("estimatebt",
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
        actionButton(inputId = 'page_before_estimate',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Correlation with stemness score</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Correlation with immune checkpoint</i>'),
        actionButton(inputId = 'page_after_estimate',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)