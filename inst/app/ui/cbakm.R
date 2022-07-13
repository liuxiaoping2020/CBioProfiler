tabItem(tabName = "cbaKM",
        fluidRow(column(9,
                        # box(
                        #      title = "Kaplan-Meier plot",
                        #      width = NULL,
                        #      solidHeader = T,
                        #      collapsible = T,
                        #      collapsed = F,
                        #      bsAlert("kmmess"),
                        #      status = "success",
                        #      align = "center",
                        bsAlert("cbakmmess"),
                        bsCollapse(id = "collapsecbaKM", open = "Descriptions and parameters",
                                   bsCollapsePanel("Descriptions and parameters",  includeHTML("cbakm.html"), style = "default"),
                                   bsCollapsePanel("Kaplan-Meier plot",style = "default",
                                                   uiOutput('cbaKMplot',align = "center")#,
                                                   
                                                   ,
                                                   fluidRow(column(4,align="left",
                                                                   useShinyjs(),
                                                                   hidden(div(id = "cbakm_wrapper",
                                                                              
                                                                              splitLayout(
                                                                                
                                                                                
                                                                                numericInput("cbakmwidth","Figure width",value = 10),
                                                                                numericInput("cbakmheight","Figure height",value = 10)),
                                                                              
                                                                              downloadButton('downloadcbaKM', 'Download figure', class = "butt2")
                                                                   ))
                                                   ))
                                                   
                                   ))
        ),
        column(3,box(
          title = "Kaplan-Meier plot",
          width = NULL,
          status = "danger",
          solidHeader = T,
          collapsible = T,
          
          # selectizeInput(
          #   "KMgene",
          #   label = "Official gene symbol",
          #   choices = NULL,
          #   multiple = F,
          #   selected = "TP53"
          # ),
          # bsTooltip("KMgene", "Input a gene with official gene symbol","left"),
          # selectizeInput(
          #   "kmgroupby",
          #   label = "Group by",
          #   choices = c("Percentage","Value"),
          #   multiple = F,
          #   selected = "Percentage"
          # ),
          # 
          # conditionalPanel(
          #   condition = "input.kmgroupby == 'Percentage'",
          #   sliderInput("genecut", "Cutoff percentage", min = 0, max = 1, value = 0.5),
          #   bsTooltip("genecut", "Input a cutoff (range:0-1) to categorize the samples into low and high expression group regarding to your interested gene. Example if 0.25: 0%-25% = Low, 25%-100% high","left"),
          #   
          # ),
          # conditionalPanel(
          #   condition = "input.kmgroupby == 'Value'",
          #   numericInput('kmgpvalue', "Cutoff value", value=1,step = 1),
          #   bsTooltip("kmgpvalue", "Input a specific cutoff value to categorize the samples into low and high expression group regarding to your interested gene.","left")
          #   
          #   
          # ),
          selectizeInput(
            "cbakmcor",
            label = "Select the cohort",
            choices = c("Training set","Validation set"),
            selected="Training set",
            multiple = F
          ),
          bsTooltip("cbakmcor", "Select the cohort (training set or validation set) to characterize the subtype","left"),
          
          selectizeInput(
            "cbasurvivaltime",
            label = "Survival time",
            choices = NULL,
            multiple = F,
            selected = "OS.time"
          ),
          bsTooltip("cbasurvivaltime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),
          
          selectizeInput(
            "cbasurvivalstatus",
            label = "Survival status",
            choices = NULL,
            multiple = F,
            selected = "OS"
          ),
          bsTooltip("cbasurvivalstatus", "Select survival time column. Example: OS, RFS, PFS","left"),
          box(title="Customized setting",width=NULL,solidHeader=F,collapsible=T,status=NULL,collapsed =T,
              textInput("cbasurvxlab","Label of X axis",placeholder="Overall survival (months)",value="Overall survival (months)"),
              bsTooltip("cbasurvxlab", "Define the label of X axis","left"),
              checkboxInput("cbasurvP", "Show p-value?", value = TRUE, width = NULL),
              checkboxInput("cbasurvRT", "Show risk table?", value = TRUE, width = NULL),
              checkboxInput("cbasurvCI", "Show confidence interval?", value = F, width = NULL)
              
          ),
          box(title="Size control",width=NULL,solidHeader=F,collapsible=T,status=NULL,collapsed =T,
              sliderInput("cbasurvwidth", "Plot Width (%)", min = 0, max = 100, value =50),
              sliderInput("cbasurvheight", "Plot Height (px)", min = 0, max = 1200, value = 500),
              sliderInput("cbakmlg", "Figure legend position", value = c(0.8, 0.90), min = 0, max = 1)
          ),
          useShinyjs(),
          actionButton(
            "cbaKMplotbt",
            "Draw KM curve",
            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
            icon = icon("picture-o")
          )
          
        ))
        ),
        
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_cbaKM',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>Correlation with clinical features</i>')
          ),
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>CoxPH model</i>'),
              actionButton(inputId = 'page_after_cbaKM',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )
        
        
        
)
