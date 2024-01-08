tabItem(tabName = "meta",
        fluidRow(
          column(9,
                 bsAlert("metamess"),
                 bsCollapse(id = "collapsedatainputmeta", open = "Select dataset for meta-analysis",
                            bsCollapsePanel("Descriptions and parameters",  includeHTML("meta.html"), style = "default"),
                            bsCollapsePanel("Select dataset for meta-analysis",DT::dataTableOutput('dataset'), style = "default"),
                            
                            bsCollapsePanel("Meta-analysis of biomarker",style = "default",
                                            uiOutput('metaplot',align = "center"), 
                                            fluidRow(column(4,align="left",
                                                            useShinyjs(),
                                                            hidden(div(id = "meta_wrapper",
                                                                       splitLayout(
                                                                         numericInput("metawidthdl","Figure width",value = 10),
                                                                         numericInput("metaheightdl","Figure height",value = 10)),
                                                                       downloadButton('downloadmeta', 'Download figure', class = "butt2")
                                                             )
                                                            )
                                              )
                                             )
                                            )
                 )
          ),
          column(3,
                 box(width = NULL,status = "danger", solidHeader=T,title = "Input datasets",
                     # selectInput("metacancer", "Cancer", choices = unique(dataset$Cancer),selected=NULL),
                     # selectInput("metaaccession", "Accession", choices = NULL,selected=NULL,multiple = T),
                     # # actionButton("metadatabt",
                     #              "Submit dataset",
                     #              style = "background-color: #000080;
                     #                        color: #FFFFFF;
                     #                        margin-left: auto;
                     #                        margin-right: auto;
                     #                        width: 100%",
                     #              icon = icon("upload")),
                     textInput("metamarker", "Official gene symbol", placeholder="TP53"),
                     bsTooltip("metamarker", "Input one gene (biomarker candidate) with official gene symbol that you want to analyze the correlation between the gene and the prognsis information of the patients based on meta-analysis","left"),
                     
                     selectizeInput("metatime",label = "Survival time",choices = c("OS.time", "RFS.time", "PFS.time","DMFS.time","CSS.time","DFS.time"),multiple = F,selected = "OS.time"),
                     bsTooltip("metatime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),
                     
                     selectizeInput("metastatus",label = "Survival status",choices = c("OS", "RFS", "PFS","DMFS","CSS","DFS"), multiple = F,selected = "OS"),
                     bsTooltip("metastatus", "Select survival status column. Example: OS, RFS, PFS","left"),
                     box(title = "Size control",width = NULL,
                         solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                         sliderInput("metawidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                         sliderInput("metaheight", "Plot Height (px)", min = 0, max = 1000, value = 400)
                     ),

                     actionButton("metabt",
                                  "Submit dataset",
                                  style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                  icon = icon("upload"))
                 ))
        ),
        
        ###################################################################################################################
        # 
        # fluidRow(
        #   column(9,
        #          bsCollapse(id = "collapsedatainput1", open = "Overview of external validation set",
        #                     bsCollapsePanel("Overview of external validation set",   style = "default",
        #                                     
        #                                     tabBox(width = NULL,
        #                                 
        #                                            tabPanel("Clinical data of external set", DT::dataTableOutput('viewclinical2'),
        #                                                     fluidRow(column(4,align="left",
        #                                                                     hidden(div(id = "viewclinical2_wrapper",
        #                                                                                downloadButton('saveviewclinical2', 'Download clinical data', class = "butt2")
        #                                                                     ))
        #                                                     ))
        #                                            ),
        #                                            
        #                
        #                                            tabPanel("Gene expression data of external set", DT::dataTableOutput('viewexpress2'),
        #                                                     fluidRow(column(4,align="left",
        #                                                                     hidden(div(id = "viewexpress2_wrapper",
        #                                                                                downloadButton('saveviewexpress2', 'Download gene expression data', class = "butt2")
        #                                                                     ))
        #                                                     ))
        #                                            )
        #                                     ))
        # 
        #          )
        #   ),
        #   column(3,
        #          box(
        #            title = "Input validation set",
        #            width = NULL,
        #            status = "danger",
        #            solidHeader = T,
        #            collapsible = T,
        #            collapsed = F,
        #            
        #            selectInput("inputType2", "Dataset type",
        #                        choices = c("Public dataset", "Customized dataset"),
        #                        selected = "Public dataset"),
        #            conditionalPanel(
        #              condition = "input.inputType2 == 'Customized dataset'",
        #              
        #              div(
        #                div(
        #                  # edit1
        #                  style="width:85%; display:inline-block; vertical-align: middle;",
        #                  fileInput("expfile2", label = h4("Gene expression matrix (.csv)"),
        #                            accept = ".csv")
        #                ),
        #                div(
        #                  # edit2
        #                  style="display:inline-block; vertical-align: middle;",
        #                  bsButton("q5", label = "", icon = icon("question"),
        #                           style = "default"
        #                  ),
        #                  bsPopover(id = "q5", title = NULL,
        #                            content = paste0("A matrix (.csv) with gene in row and sample in column."),
        #                            placement = "left",
        #                            trigger = "hover",
        #                            options = list(container = "body")
        #                  )
        #                )
        #              ),
        #              div(
        #                div(
        #                  style="width:85%; display:inline-block; vertical-align: middle;",
        #                  fileInput("clin2", label = h4("Clinical data matrix (.csv)"),
        #                            accept = ".csv")
        #                ),
        #                div(
        #                  style="display:inline-block; vertical-align: middle;",
        #                  bsButton("q6", label = "", icon = icon("question"),
        #                           style = "default"
        #                  ),
        #                  bsPopover(id = "q6", title = NULL,
        #                            content = paste0("A matrix (.csv) with sample in row and clinical feature in column."),
        #                            placement = "left",
        #                            trigger = "hover",
        #                            options = list(container = "body")
        #                  )
        #                )
        #              )
        #            ),
        #            conditionalPanel(
        #              condition = "input.inputType2 == 'Public dataset'",
        #              selectInput("cancer2", "Cancer", choices = unique(dataset$Cancer),selected=NULL),
        #              selectInput("accession2", "Accession", choices = NULL,selected=NULL)
        #            ),
        #            actionButton("data2",
        #                         "Submit dataset",
        #                         style = "background-color: #000080;
        #                                     color: #FFFFFF;
        #                                     margin-left: auto;
        #                                     margin-right: auto;
        #                                     width: 100%",
        #                         icon = icon("upload"))
        #          ),
        #   )
        # ),
        
        br(),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_meta',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>Biology annotation</i>')
          ),

          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>Subtype identification</i>'),
              actionButton(inputId = 'page_after_meta',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )
      )

