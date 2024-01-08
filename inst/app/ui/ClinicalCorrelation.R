tabItem(tabName = "clinical",
        fluidRow(column(
                9,
                # box(
                #         title = "Table for the correlations between gene expression and clinical features",
                #         width = NULL,
                #         solidHeader = T,
                #         collapsible = F,
                #         status = "success",
                bsCollapse(id = "collapseclinical", open = "Descriptions and parameters",
                           bsCollapsePanel("Descriptions and parameters",  includeHTML("table1.html"), style = "default"),
                           bsCollapsePanel("Table for the correlations between gene expression and clinical features",style = "default",
                        # align="center",
                        htmlOutput("table1",align="center")
                )
        )),
        column(
                3,
                box(
                        title = "Correlations between gene expression and clinical features",
                        width = NULL,
                        status = "danger",
                        solidHeader = T,
                        collapsible = F,

                        selectizeInput(
                          "table1bygene",
                          label = "Official gene symbol",
                          choices = NULL,
                          multiple = F,
                          selected = "TP53"
                        ),
                        bsTooltip("table1bygene", "Input a gene with official gene symbol","left"),

                        selectizeInput(
                          "tbgroupby",
                          label = "Group by",
                          choices = c("Percentage","Value"),
                          multiple = F,
                          selected = "Percentage"
                        ),
                        conditionalPanel(
                          condition = "input.tbgroupby == 'Percentage'",
                        sliderInput("grouppercent", "Cutoff percentage", min = 0, max = 1, value = 0.5),
                        bsTooltip("grouppercent", "Input a cutoff percentage (range:0-1) to categorize the samples into low and high expression group regarding to your interested gene. Example if 0.25: 0%-25% = Low, 25%-100% high","left")

                        ),
                        conditionalPanel(
                          condition = "input.tbgroupby == 'Value'",
                          numericInput('tabgpvalue', "Cutoff value", value=1,step = 1),
                          bsTooltip("tabgpvalue", "Input a specific cutoff value to categorize the samples into low and high expression group regarding to your interested gene.","left"),


                        ),
                        selectizeInput(
                                "feature",
                                label = "Select clinical features",
                                choices = NULL,
                                multiple = T
                        ),
                        bsTooltip("feature", "Select clinical features you want to include in the table","left"),
                        useShinyjs(),
                        actionButton(
                                "table1bt",
                                "Draw table1",
                                style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                icon = icon("upload")
                        ))


        )),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_clinical',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>Prediction model: Nomogram</i>')
          ),
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>Kaplan-Meier curve</i>'),
              actionButton(inputId = 'page_after_clinical',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )



        )
