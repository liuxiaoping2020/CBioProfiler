tabItem(
  tabName = "cbaCYA",
  fluidRow(
    column(9,
           bsAlert("cbaCYAmess"),
           bsCollapse(id = "cbacollapseCYA", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("cbaCYA.html"), style = "default"),
                      bsCollapsePanel("Comparision of cytotoxic activity distribtion",style = "default",
                                      tabBox(width = NULL,
                                             tabPanel("Cytotoxic activity plot", align="center",
                                                      uiOutput('cbaCYAplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "cbaCYAplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("cbaCYApltwidth","Figure width",value = 10),
                                                                                   numericInput("cbaCYApltheight","Figure height",value = 10)),
                                                                                 downloadButton('cbasaveCYAness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             ),
                                             tabPanel("Cytotoxic activity table", DT::dataTableOutput('cbaCYAtable'),align="center",
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "cbaCYAtable_wrapper",
                                                                                 downloadButton('cbasaveCYAtable', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             )
                                             
                                      )
                      ))
    ),
    column(3,box(
      title = "Comparision of cytotoxic activity distribtion",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      selectizeInput(
        "cbaCYAcor",
        label = "Select the cohort",
        choices = c("Training set","Validation set"),
        selected="Training set",
        multiple = F
      ),
      bsTooltip("cbaCYAcor", "Select the cohort (training set or validation set) to characterize the subtype","left"),
      
      selectizeInput(
        "cbaCYAmeth",
        label = "Group comparison method",
        choices = c("Nonparametric test","Parametric test"),
        multiple = F,
        selected="Nonparametric test"
      ),
      bsTooltip("cbaCYAmeth", "Whether to use a parametric or nonparametric test for comparisons between groups","left"),
      
      checkboxInput("cbaCYApair", "Perform paire-wised comparision ?", value = F, width = NULL),
      selectizeInput(
        "cbaCYAplot",
        label = "Plot type",
        choices = c("Bar plot", "Box plot", "Violin plot"),
        multiple = F,
        selected="Box plot"
      ),
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("cbaCYAwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
          sliderInput("cbaCYAheight", "Plot Height (px)", min = 0, max = 1000, value = 430)
      ),
      actionButton("cbaCYAbt",
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
        actionButton(inputId = 'page_before_cbaCYA',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Interferon-gamma among subtypes</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Cancer pathway score among subtypes</i>'),
        actionButton(inputId = 'page_after_cbaCYA',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)