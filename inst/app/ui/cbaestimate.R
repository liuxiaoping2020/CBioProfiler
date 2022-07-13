tabItem(
  tabName = "cbaestimate",
  fluidRow(
    column(9,
           bsAlert("cbaestimatemess"),
           bsCollapse(id = "cbacollapseestimate", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("cbaestimate.html"), style = "default"),
                      bsCollapsePanel("Correlation with ESTIMATE score",style = "default",
                                      tabBox(width = NULL,
                                      tabPanel("ESTIMATE plot", align="center",
                                      uiOutput('cbaestimateplot',align = "center"),
                                      useShinyjs(),
                                      fluidRow(column(4,align="left",
                                                      hidden(div(id = "cbaestimateplot_wrapper",
                                                                 splitLayout(
                                                                   numericInput("cbaestimatepltwidth","Figure width",value = 10),
                                                                   numericInput("cbaestimatepltheight","Figure height",value = 10)),
                                                                 downloadButton('cbasaveestimateness', 'Download figure', class = "butt2")
                                                      ))
                                      ))
                                      ),
                                      tabPanel("ESTIMATE table", DT::dataTableOutput('cbaestimatetable'),align="center",
                                               useShinyjs(),
                                               fluidRow(column(4,
                                                               hidden(div(id = "cbaestimatetable_wrapper",
                                                                          downloadButton('cbasaveestimatetable', 'Download table', class = "butt2")
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
      selectizeInput(
        "cbaestimatecor",
        label = "Select the cohort",
        choices = c("Training set","Validation set"),
        selected="Training set",
        multiple = F
      ),
      bsTooltip("cbaestimatecor", "Select the cohort (training set or validation set) to characterize the subtype","left"),
      
      selectizeInput(
        "cbaestimatemeth",
        label = "Group comparison method",
        choices = c("Nonparametric test","Parametric test"),
        multiple = F,
        selected="Nonparametric test"
      ),
      bsTooltip("cbaestimatemeth", "Whether to use a parametric or nonparametric test for comparisons between groups","left"),
      
      checkboxInput("cbaestimatepair", "Perform paire-wised comparision ?", value = F, width = NULL),
      selectizeInput(
        "cbaestimateplot",
        label = "Plot type",
        choices = c("Bar plot", "Box plot", "Violin plot"),
        multiple = F,
        selected="Box plot"
      ),
      
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("cbaestimatewidth", "Plot Width (%)", min = 0, max = 100, value = 50),
          sliderInput("cbaestimateheight", "Plot Height (px)", min = 0, max = 1000, value = 430)
      ),
      actionButton("cbaestimatebt",
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
        actionButton(inputId = 'page_before_cbaestimate',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Stemness score among subtypes</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Immune checkpoints among subtypes</i>'),
        actionButton(inputId = 'page_after_cbaestimate',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)