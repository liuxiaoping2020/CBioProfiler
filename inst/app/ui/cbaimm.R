tabItem(
  tabName = "cbaimmune",
  fluidRow(
    column(9,
           # box(
           # title="Correlation with Immune cell infiltration",width = NULL,
           # solidHeader=T,collapsible=T,status="success",collapsed = F,
           bsAlert("cbaimmunemess"),
           bsCollapse(id = "cbacollapseimmune", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("cbaimm.html"), style = "default"),
                      bsCollapsePanel("Immune cell infiltration distribution accross different subtypes",style = "default",
                                      # tabBox(width = NULL,
                                             # tabPanel("Distribution plot", align="center",
                                                      uiOutput('cbaimmuneplot',align = "center"),
                                                      useShinyjs(),
                                                      fluidRow(column(4,align="left",
                                                                      hidden(div(id = "cbaimmuneplot_wrapper",
                                                                                 splitLayout(
                                                                                   numericInput("cbaimmpltwidth","Figure width",value = 10),
                                                                                   numericInput("cbaimmpltheight","Figure height",value = 10)),
                                                                                 downloadButton('cbasaveimune', 'Download figure', class = "butt2")
                                                                      ))
                                                      ))
                                             # ),
                                             
                                             # tabPanel("Immune cell composition matrix", dataTableOutput('cbaimmcelltable'),
                                             #          useShinyjs(),
                                             #          fluidRow(column(4,
                                             #                          hidden(div(id = "cbaimmunecell_wrapper",
                                             #                                     downloadButton('cbasaveimmcelltable', 'Download table', class = "butt2")
                                             #                          ))
                                             #          ))
                                             # ),
                                             # tabPanel("Correlation table", dataTableOutput('cbaimmcortable'),
                                             #          useShinyjs(),
                                             #          fluidRow(column(4,
                                             #                          hidden(div(id = "cbaimmcor_wrapper",
                                             #                                     downloadButton('cbasaveimmcortable', 'Download table', class = "butt2")
                                             #                          ))
                                             #          ))
                                             # )
                                      # )
                      ))
    ),
    column(3,box(
      title = "Immune cell infiltration distribution",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      # selectizeInput(
      #   "immgene",
      #   label = "Official gene symbol",
      #   choices = NULL,
      #   multiple = F
      # ),
      # bsTooltip("immgene", "Input one gene (biomarker candidate) with official gene symbol that you want to analyze the correlation between the gene and tummor microenvironment.","left"),
      selectizeInput(
        "cbaimmcor",
        label = "Select the cohort",
        choices = c("Training set","Validation set"),
        selected="Training set",
        multiple = F
      ),
      bsTooltip("cbaimmcor", "Select the cohort (training set or validation set) to characterize the subtype","left"),
      selectizeInput(
        "cbageneset",
        label = "Reference immune geneset",
        choices = c("Bindea","Danaher","Davoli","MCP.Counter","xCell"),
        multiple = F
      ),
      bsTooltip("cbageneset", "Immune gene signatures from Bindea and colleagues, Danaher and colleagues, Davoli and colleagues, MCP-Counter, and xCell, respectively","left"),
      
      selectizeInput(
        "cbaimmeth",
        label = "Group comparison method",
        choices = c("Nonparametric test","Parametric test"),
        multiple = F,
        selected="Nonparametric test"
      ),
      bsTooltip("cbaimmeth", "Whether to use a parametric or nonparametric test for comparisons between groups","left"),
      
      checkboxInput("cbaimmpair", "Perform paired-wised comparision ?", value = F, width = NULL),
      selectizeInput(
        "cbaimmplot",
        label = "Plot type",
        choices = c("Bar plot", "Box plot", "Violin plot"),
        multiple = F,
        selected="Box plot"
      ),
      
      box(title = "Size control",width = NULL,
          solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          sliderInput("cbaimmwidth", "Plot Width (%)", min = 0, max = 100, value = 100),
          sliderInput("cbaimmheight", "Plot Height (px)", min = 0, max = 1000, value = 800)
      ),
      shinyjs::useShinyjs(),
      actionButton("cbaimmunebt",
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
        actionButton(inputId = 'page_before_cbaimmune',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Differentially expressed genes</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Stemness score among subtypes</i>'),
        actionButton(inputId = 'page_after_cbaimmune',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)
