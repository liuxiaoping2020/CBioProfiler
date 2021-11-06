tabItem(tabName = "gene",
        fluidRow(column(9,
          box(
            title = "The correlation between two genes",
            width = NULL,
            solidHeader = T,
            collapsible = F,
            status = "success",
            align="center",
            tabBox(width = NULL,

              tabPanel(
                "Scatter plot",
                uiOutput('genecorplot'),
                useShinyjs(),
                fluidRow(column(4,align="left",
                                hidden(div(id = "genecor_wrapper",
                                          splitLayout(
                                             numericInput("genecorwidthdl","Figure width",value = 10),
                                             numericInput("genecorheightdl","Figure height",value = 10)),
                                           downloadButton('downgenecor', 'Download figure', class = "butt2")
                                ))
                ))
              ),
              tabPanel(
                "Summary",
                verbatimTextOutput("genecorsummary"),align="left",
                tags$head(
                  tags$style(
                    "#genecorsummary{ font-size:12px; font-style:Arial;width: 1000px; max-width: 215%;background: ghostwhite;}"
                  )
                )
              )

            )
          )
        ),
        column(
          3,
          box(
            title = "Gene-gene correlation analysis",width = NULL,status = "danger",solidHeader=T,
            collapsible = T,
            selectizeInput("gene1",label = "Gene 1",choices = NULL,multiple = F),
            selectizeInput("gene2",label = "Gene 2",choices = NULL,multiple = F),
            radioButtons(
              "cormethods","Correlation methods",
              c("Pearson's correlation" = "pearson",
                "Spearman's correlation" = "spearman",
                "Kendall's correlation" = "kendall"
              )
            ),
            box(title = "Size control",width = NULL,
                solidHeader = F, collapsible = TRUE, collapsed = TRUE,
                sliderInput("genecorwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
                sliderInput("genecorheight", "Plot Height (px)", min = 0, max = 1000, value = 400)),
            actionButton(
              "genecorbt","Submit",
              style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
              icon = icon("picture-o")
            )
          )
        )),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_gene',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>Most correlated genes</i>')
          ),
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>Gene expression in different groups</i>'),
              actionButton(inputId = 'page_after_gene',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )






        )
