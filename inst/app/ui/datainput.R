tabItem(tabName = "dataset",


         fluidRow(
          column(9,
              box(
              title="Overview of discovery data",width = NULL,
              solidHeader=T,collapsible=T,collapsed=F,status="success",
              bsAlert("inputmess"),
              tabBox(width = NULL,
                     tabPanel("Clinical data", DT::dataTableOutput('viewclinical'),
              useShinyjs(),
              fluidRow(column(4,
                              hidden(div(id = "clinicalinput.table",
                                         downloadButton('saveclinicalinput', 'Download clinical data', class = "butt2")
                              ))
              ))
              ),
              tabPanel("Gene expression data",DT::dataTableOutput('viewexpres'),
              useShinyjs(),
              fluidRow(column(4,
                              hidden(div(id = "express.table",
                                         downloadButton('saveexpresinput', 'Download gene expression data', class = "butt2")
                              ))
              ))


              )))
              # ,
              # box(
              #   title="Overview of gene expression data",width = NULL, solidHeader=T,collapsible=T,collapsed=T,status="success",
              # # bsCollapsePanel("Overview of gene expression data",style = "info",
              #   DT::dataTableOutput('viewexpres'),
              #   useShinyjs(),
              #   fluidRow(column(4,
              #                   hidden(div(id = "express.table",
              #                              downloadButton('saveexpresinput', 'Download gene expression data', class = "butt2")
              #                   ))
              #   ))
              # )
              # )
              ),
          column(3,
                 box(width = NULL,status = "danger", solidHeader=T,title = "Input discovery data",
                   selectInput("inputType", "Dataset type",
                               choices = c("Public dataset", "Customized dataset"),
                               selected = "Public dataset"),
                   conditionalPanel(
                     condition = "input.inputType == 'Customized dataset'",
                     # box(title = "Upload a gene expression matrix", width = NULL,
                     # h3('Upload a gene expression matrix'),
                         # solidHeader = TRUE, collapsible = TRUE,
                         # fileInput("expfile" ,"Gene expression matrix (.csv) with gene in row and sample in column.", accept = ".csv")
                     # )
                     # ,
                     # bsTooltip("expfile", "A matrix (.csv) with gene in row and sample in column.","left"),
                     div(
                       div(
                         # edit1
                         style="width:85%; display:inline-block; vertical-align: middle;",
                         fileInput("expfile", label = h4("Gene expression matrix (.csv)"),
                                   accept = ".csv")
                       ),
                       div(
                         # edit2
                         style="display:inline-block; vertical-align: middle;",
                         bsButton("q1", label = "", icon = icon("question"),
                                  style = "default"
                                  ),
                         bsPopover(id = "q1", title = NULL,
                                   content = paste0("A matrix (.csv) with gene in row and sample in column."),
                                   placement = "left",
                                   trigger = "hover",
                                   options = list(container = "body")
                         )
                       )
                     ),
                     div(
                       div(
                         # edit1
                         style="width:85%; display:inline-block; vertical-align: middle;",
                         fileInput("clin", label = h4("Clinical data matrix (.csv)"),
                                   accept = ".csv")
                       ),
                       div(
                         # edit2
                         style="display:inline-block; vertical-align: middle;",
                         bsButton("q2", label = "", icon = icon("question"),
                                  style = "default"
                         ),
                         bsPopover(id = "q2", title = NULL,
                                   content = paste0("A matrix (.csv) with sample in row and clinical feature in column."),
                                   placement = "left",
                                   trigger = "hover",
                                   options = list(container = "body")
                         )
                       )
                     )

                     # box(title = "Upload clinical data",width = NULL,
                     # h4('Upload clinical data'),
                         # solidHeader = TRUE, collapsible = TRUE,
                         # fileInput("clin", "Clinical data matrix", accept = ".csv"),
                     # bsTooltip("clin", "A matrix (.csv) with sample in row and clinical feature in columns.","left")
                     # )
                   ),


                   conditionalPanel(
                     condition = "input.inputType == 'Public dataset'",
                       selectInput("cancer", "Cancer", choices = unique(dataset$Cancer),selected=NULL),
                       selectInput("accession", "Accession", choices = NULL,selected=NULL)
                   ),
                 actionButton("data",
                              "Submit dataset",
                              style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                              icon = icon("picture-o"))
                 ))
        ),

###################################################################################################################

        fluidRow(
          column(9,
                 box(
                   title="Overview of external validation data",width = NULL,
                   solidHeader=T,collapsible=T,collapsed=F,status="success",

                   tabBox(width = NULL,
                          # tabPanel("Clinical data", DT::dataTableOutput('viewclinical'),
                          #          useShinyjs(),
                          #          fluidRow(column(4,
                          #                          hidden(div(id = "clinicalinput.table",
                          #                                     downloadButton('saveclinicalinput', 'Download clinical data', class = "butt2")
                          #                          ))
                          #          ))
                          # ),

                          tabPanel("Clinical data of external set", DT::dataTableOutput('viewclinical2'),
                                   fluidRow(column(4,align="left",
                                                   hidden(div(id = "viewclinical2_wrapper",
                                                              downloadButton('saveviewclinical2', 'Download table', class = "butt2")
                                                   ))
                                   ))
                          ),

                          # tabPanel("Gene expression data",DT::dataTableOutput('viewexpres'),
                          #          useShinyjs(),
                          #          fluidRow(column(4,
                          #                          hidden(div(id = "express.table",
                          #                                     downloadButton('saveexpresinput', 'Download gene expression data', class = "butt2")
                          #                          ))
                          #          ))
                          #
                          #
                          # )
                          tabPanel("Gene expression data of external set", DT::dataTableOutput('viewexpress2'),
                                   fluidRow(column(4,align="left",
                                                   hidden(div(id = "viewexpress2_wrapper",
                                                              downloadButton('saveviewexpress2', 'Download table', class = "butt2")
                                                   ))
                                   ))
                          )
                          ))
                 # ,
                 # box(
                 #   title="Overview of gene expression data",width = NULL, solidHeader=T,collapsible=T,collapsed=T,status="success",
                 # # bsCollapsePanel("Overview of gene expression data",style = "info",
                 #   DT::dataTableOutput('viewexpres'),
                 #   useShinyjs(),
                 #   fluidRow(column(4,
                 #                   hidden(div(id = "express.table",
                 #                              downloadButton('saveexpresinput', 'Download gene expression data', class = "butt2")
                 #                   ))
                 #   ))
                 # )
                 # )
          ),
          column(3,
                 box(
                   title = "Input validation set",
                   width = NULL,
                   status = "danger",
                   solidHeader = T,
                   collapsible = T,
                   collapsed = F,

                   selectInput("inputType2", "Dataset type",
                               choices = c("Public dataset", "Customized dataset"),
                               selected = "Public dataset"),
                   conditionalPanel(
                     condition = "input.inputType2 == 'Customized dataset'",

                     div(
                       div(
                         # edit1
                         style="width:85%; display:inline-block; vertical-align: middle;",
                         fileInput("expfile2", label = h4("Gene expression matrix (.csv)"),
                                   accept = ".csv")
                       ),
                       div(
                         # edit2
                         style="display:inline-block; vertical-align: middle;",
                         bsButton("q5", label = "", icon = icon("question"),
                                  style = "default"
                         ),
                         bsPopover(id = "q5", title = NULL,
                                   content = paste0("A matrix (.csv) with gene in row and sample in column."),
                                   placement = "left",
                                   trigger = "hover",
                                   options = list(container = "body")
                         )
                       )
                     ),
                     div(
                       div(
                         style="width:85%; display:inline-block; vertical-align: middle;",
                         fileInput("clin2", label = h4("Clinical data matrix (.csv)"),
                                   accept = ".csv")
                       ),
                       div(
                         style="display:inline-block; vertical-align: middle;",
                         bsButton("q6", label = "", icon = icon("question"),
                                  style = "default"
                         ),
                         bsPopover(id = "q6", title = NULL,
                                   content = paste0("A matrix (.csv) with sample in row and clinical feature in column."),
                                   placement = "left",
                                   trigger = "hover",
                                   options = list(container = "body")
                         )
                       )
                     )


                   ),
                   conditionalPanel(
                     condition = "input.inputType2 == 'Public dataset'",
                     selectInput("cancer2", "Cancer", choices = unique(dataset$Cancer),selected=NULL),
                     selectInput("accession2", "Accession", choices = NULL,selected=NULL)
                   ),
                   actionButton("data2",
                                "Submit dataset",
                                style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                icon = icon("picture-o"))
                 ),


                 # box(width = NULL,status = "danger", solidHeader=T,title = "Input data",
                 #
                 #
                 #     selectInput("inputType", "Dataset type",
                 #                 choices = c("Public dataset", "Customized dataset"),
                 #                 selected = "Public dataset"),
                 #     conditionalPanel(
                 #       condition = "input.inputType == 'Customized dataset'",
                 #       # box(title = "Upload a gene expression matrix", width = NULL,
                 #       # h3('Upload a gene expression matrix'),
                 #       # solidHeader = TRUE, collapsible = TRUE,
                 #       # fileInput("expfile" ,"Gene expression matrix (.csv) with gene in row and sample in column.", accept = ".csv")
                 #       # )
                 #       # ,
                 #       # bsTooltip("expfile", "A matrix (.csv) with gene in row and sample in column.","left"),
                 #       div(
                 #         div(
                 #           # edit1
                 #           style="width:85%; display:inline-block; vertical-align: middle;",
                 #           fileInput("expfile", label = h4("Gene expression matrix (.csv)"),
                 #                     accept = ".csv")
                 #         ),
                 #         div(
                 #           # edit2
                 #           style="display:inline-block; vertical-align: middle;",
                 #           bsButton("q1", label = "", icon = icon("question"),
                 #                    style = "default"
                 #           ),
                 #           bsPopover(id = "q1", title = NULL,
                 #                     content = paste0("A matrix (.csv) with gene in row and sample in column."),
                 #                     placement = "left",
                 #                     trigger = "hover",
                 #                     options = list(container = "body")
                 #           )
                 #         )
                 #       ),
                 #       div(
                 #         div(
                 #           # edit1
                 #           style="width:85%; display:inline-block; vertical-align: middle;",
                 #           fileInput("clin", label = h4("Clinical data matrix (.csv)"),
                 #                     accept = ".csv")
                 #         ),
                 #         div(
                 #           # edit2
                 #           style="display:inline-block; vertical-align: middle;",
                 #           bsButton("q2", label = "", icon = icon("question"),
                 #                    style = "default"
                 #           ),
                 #           bsPopover(id = "q2", title = NULL,
                 #                     content = paste0("A matrix (.csv) with sample in row and clinical feature in column."),
                 #                     placement = "left",
                 #                     trigger = "hover",
                 #                     options = list(container = "body")
                 #           )
                 #         )
                 #       )
                 #
                 #       # box(title = "Upload clinical data",width = NULL,
                 #       # h4('Upload clinical data'),
                 #       # solidHeader = TRUE, collapsible = TRUE,
                 #       # fileInput("clin", "Clinical data matrix", accept = ".csv"),
                 #       # bsTooltip("clin", "A matrix (.csv) with sample in row and clinical feature in columns.","left")
                 #       # )
                 #     ),
                 #
                 #
                 #     conditionalPanel(
                 #       condition = "input.inputType == 'Public dataset'",
                 #       selectInput("cancer", "Cancer", choices = unique(dataset$Cancer),selected=NULL),
                 #       selectInput("accession", "Accession", choices = NULL,selected=NULL)
                 #     ),
                 #     actionButton("data",
                 #                  "Submit dataset",
                 #                  style = "background-color: #000080;
                 #                            color: #FFFFFF;
                 #                            margin-left: auto;
                 #                            margin-right: auto;
                 #                            width: 100%",
                 #                  icon = icon("picture-o"))
                 # )

#########################################################################################################################

                 )
),







        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_datainput',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>Introduction</i>')
          ),
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>Dimensionality reduction: WGCNA</i>'),
              actionButton(inputId = 'page_after_datainput',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )


)


