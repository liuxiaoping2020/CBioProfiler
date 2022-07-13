tabItem(
  tabName = "cbaIFN",
  fluidRow(
    column(9,
           bsAlert("cbaIFNmess"),
           bsCollapse(id = "cbacollapseIFN", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("cbaIFN.html"), style = "default"),
                      bsCollapsePanel("Comparision of interferon-gamma distribtion",style = "default",
                                      tabBox(width = NULL,
                                             tabPanel("Interferon-gamma score plot", align="center",
                                                      uiOutput('cbaIFNplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "cbaIFNplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("cbaIFNpltwidth","Figure width",value = 10),
                                                                                   numericInput("cbaIFNpltheight","Figure height",value = 10)),
                                                                                 downloadButton('cbasaveIFNness', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             ),
                                             tabPanel("Interferon-gamma score table", DT::dataTableOutput('cbaIFNtable'),align="center",
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "cbaIFNtable_wrapper",
                                                                                 downloadButton('cbasaveIFNtable', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             )
                                             
                                      )
                      ))
    ),
    column(3,box(
      title = "Comparision of interferon-gamma distribtion",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      selectizeInput(
        "cbaIFNcor",
        label = "Select the cohort",
        choices = c("Training set","Validation set"),
        selected="Training set",
        multiple = F
      ),
      bsTooltip("cbaIFNcor", "Select the cohort (training set or validation set) to characterize the subtype","left"),
      
      selectizeInput(
        "cbaIFNmeth",
        label = "Group comparison method",
        choices = c("Nonparametric test","Parametric test"),
        multiple = F,
        selected="Nonparametric test"
      ),
      bsTooltip("cbaIFNmeth", "Whether to use a parametric or nonparametric test for comparisons between groups","left"),
      
      checkboxInput("cbaIFNpair", "Perform paire-wised comparision ?", value = F, width = NULL),
      selectizeInput(
        "cbaIFNplot",
        label = "Plot type",
        choices = c("Bar plot", "Box plot", "Violin plot"),
        multiple = F,
        selected="Box plot"
      ),
      
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("cbaIFNwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
          sliderInput("cbaIFNheight", "Plot Height (px)", min = 0, max = 1000, value = 430)
      ),
      actionButton("cbaIFNbt",
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
        actionButton(inputId = 'page_before_cbaIFN',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Immune checkpoints among subtypes</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Cytolytic activity among subtypes</i>'),
        actionButton(inputId = 'page_after_cbaIFN',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)