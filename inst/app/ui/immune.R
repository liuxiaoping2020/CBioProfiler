tabItem(
  tabName = "immune",
  fluidRow(
    column(9,
      # box(
      # title="Correlation with Immune cell infiltration",width = NULL,
      # solidHeader=T,collapsible=T,status="success",collapsed = F,
      bsAlert("immunemess"),
      bsCollapse(id = "collapseimmune", open = "Descriptions and parameters",
                 bsCollapsePanel("Descriptions and parameters",  includeHTML("immune.html"), style = "default"),
                 bsCollapsePanel("Correlation with immune cell infiltration",style = "default",
      tabBox(width = NULL,
             tabPanel("Correlation plot", align="center",
                      uiOutput('immuneplot'
                      )
                      ,
                      useShinyjs(),
                      fluidRow(column(4,align="left",
                                      hidden(div(id = "immuneplot_wrapper",
                                                splitLayout(
                                                   numericInput("immpltwidth","Figure width",value = 10),
                                                   numericInput("immpltheight","Figure height",value = 10)),
                                                 downloadButton('saveimune', 'Download figure', class = "butt2")
                                      ))
                      ))
             ),

             tabPanel("Immune cell composition matrix", dataTableOutput('immcelltable'),
                      useShinyjs(),
                      fluidRow(column(4,
                                      hidden(div(id = "immunecell_wrapper",
                                                 downloadButton('saveimmcelltable', 'Download table', class = "butt2")
                                      ))
                      ))
             ),
             tabPanel("Correlation table", dataTableOutput('immcortable'),
                      useShinyjs(),
                      fluidRow(column(4,
                                      hidden(div(id = "immcor_wrapper",
                                                 downloadButton('saveimmcortable', 'Download table', class = "butt2")
                                      ))
                      ))
             )
      )
    ))
           ),
    column(3,box(
      title = "Correlation with Immune cell infiltration",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      selectizeInput(
        "geneset",
        label = "Reference immune geneset",
        choices = c("Bindea","Danaher","Davoli","MCP.Counter","xCell"),
        multiple = F
      ),
      bsTooltip("geneset", "Immune gene signatures from Bindea and colleagues, Danaher and colleagues, Davoli and colleagues, MCP-Counter, and xCell, respectively","left"),
      selectizeInput(
        "immgene",
        label = "Official gene symbol",
        choices = NULL,
        multiple = F
      ),
      bsTooltip("immgene", "Input one gene (biomarker candidate) with official gene symbol that you want to analyze the correlation between the gene and tummor microenvironment.","left"),
      selectizeInput(
        "immmethod",
        label = "Correlation method",
        choices = c("pearson", "spearman"),
      multiple = F,
      selected="spearman"
    ),
      selectizeInput(
        "immnorm",
        label = "Adjust method for P",
        choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),
        multiple = F,
        selected="bonferroni"
      ),
      checkboxInput("corselect", "Select most correlated cells ?", value = T, width = NULL),
      bsTooltip("corselect", "Select the most correlated immune cells based on the cutoff of adjusted p value? ","left"),
    conditionalPanel(
      condition = "input.corselect == true ",
      numericInput("immcutoff","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005)
    ),

      selectizeInput(
      "immpltype",
      label = "Plot type",
      choices = c("bar plot","bubble plot"),
      multiple = F,
      selected="bar plot"
  ),
  box(title = "Size control",width = NULL,
      solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
      sliderInput("immwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
      sliderInput("immheight", "Plot Height (px)", min = 0, max = 1000, value = 430)
  ),
  actionButton("immunebt",
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
        actionButton(inputId = 'page_before_immune',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Gene expression in different groups</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Correlation with stemness score</i>'),
        actionButton(inputId = 'page_after_immune',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )





)
