tabItem(
  tabName = "wgcna",
  fluidPage(
    column(9,
           box(
             title = "Data Preprocessing",
             width = NULL,
             height=NULL,
             solidHeader = T,
             collapsible = T,
             collapsed = F,
             status = "success",
             align="center",
             bsAlert("wgcnapremess"),

                   uiOutput('sampletree'),
                   textOutput(outputId = "datacleandescrip"),
             useShinyjs(),
             fluidRow(column(4,align="left",
                             hidden(div(id = "sampletree_wrapper",
                                        splitLayout(
                                          numericInput("sampletreewidthdl","Figure width",value = 10),
                                          numericInput("sampletreeheightdl","Figure height",value = 10)),
                                        downloadButton('downloadsampletree', 'Download figure', class = "butt2")
                             ))
             ))
                   )
  ),
  column(3,
         box(title = "Data preprocess",
             solidHeader = T,
             collapsible = T, width = NULL, collapsed = F, status = "danger",
             checkboxInput("setopvar", "Select most variable genes ?", value = T, width = NULL),
             conditionalPanel(
               condition = "input.setopvar == true ",
               numericInput("topvar","Most variable genes to include",2000,min = 1,max = 15000,step =100),
               bsTooltip("topvar", "Selection top n variable genes to be included for network constuction based on the variances of genes cross all samples","left")
             ),
             numericInput("ZK","Z.K",-2.5,min = -10,max = 10,step =1),
             bsTooltip("ZK", "Cutoff in sample network for outlier detection, please refer to Horvath S (2011) Weighted Network Analysis. Applications in Genomics and Systems Biology. Springer Book. ISBN: 978-1-4419-8818-8 for more details","left"),
             selectizeInput(
               'wgvariable',
               "Numeric clinical variable",
               choices = NULL,
               selected = NULL,
               multiple=T
             ),
             bsTooltip("wgvariable", "Select 2 or more numeric clinical variables to be included for correlation analysis","left"),
             box(title="Size control",solidHeader=F,collapsible=T,width=NULL,collapsed = T,
               sliderInput("sampwidth", "Heatmap Width (%)", min = 0, max = 100, value = 100),
               sliderInput("samheight", "Heatmap Height (px)", min = 0, max = 2000, value = 400)
             ),

             actionButton("datacleanbt",
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
hidden(div(id = "wgcna1_wrapper",
fluidPage(column(9,
                 box(
                   title = "Network construction and module detection",
                   width = NULL,
                   solidHeader = T,
                   collapsible = T,
                   collapsed = F,
                   align="center",
                   status = "success",
                   tabBox(
                     width=NULL,
                     tabPanel(
                       title = "Soft threshold selection",
                       value = "sft",
                       uiOutput('sftdistribution'),
                       textOutput(outputId = "sftdiscrp"),
                       useShinyjs(),
                       fluidRow(column(4,align="left",
                                       hidden(div(id = "sftdistribution_wrapper",
                                                  splitLayout(
                                                    numericInput("sftdistributionwidthdl","Figure width",value = 10),
                                                    numericInput("sftdistributionheightdl","Figure height",value = 10)),
                                                  downloadButton('downloadsftdistribution', 'Download figure', class = "butt2")
                                       ))
                       ))
                       ),
                     tabPanel(title="Module dendrogram",value="modeng",uiOutput('modengram'),
                     fluidRow(column(4,align="left",
                                     hidden(div(id = "modengram_wrapper",
                                                splitLayout(
                                                  numericInput("modengramwidthdl","Figure width",value = 10),
                                                  numericInput("modengramheightdl","Figure height",value = 10)),
                                                downloadButton('downloadmodengram', 'Download figure', class = "butt2")
                                     ))
                     ))
                   )
                 )
  )
),
column(3,box(title="Network & module detection",
             solidHeader=T,
             collapsible=T,width=NULL,collapsed = F,status = "danger",
              bsCollapse(id = "Networkmoduledetection",open="Network construction",
              bsCollapsePanel("Network construction",style = "info",


               selectizeInput(
               'networkType',
               "NetworkType",
               choices = c("unsigned", "signed", "signed hybrid"),
               selected = "unsigned",
               multiple=F),
             bsTooltip("networkType", "network type. Allowed values are (unique abbreviations of) 'signed','unsigned','signed hybrid'. More detail can be found at https://cran.r-project.org/web/packages/WGCNA/WGCNA.pdf","left"),
             numericInput("RsquaredCut","Rsquared Cutoff",0.9,min = 0,max = 1,step =0.1),
             bsTooltip("RsquaredCut", "Desired minimum scale free topology fitting index R^2.","left"),
               selectizeInput(
               'corType',
               "Correlation Type",
               choices = c("pearson", "bicor"),
               selected = "pearson",
               multiple=F
             )
             ,
             bsTooltip("corType", "Specifying the correlation type, 'pearson' means 'Pearson's Correlation','bicor' means 'Bidweight midcorrelation' ","left"),
               selectizeInput(
               'TOMType',
               "TOM Type",
               choices = c("none", "unsigned", "signed", "signed Nowick", "unsigned 2", "signed 2","signed Nowick 2"),
               selected ="unsigned",
               multiple=F
             ),
             bsTooltip("TOMType", "Specifying the  Topology Overlap Matrix (TOM) type' ","left")

             ),


             bsCollapsePanel("Module detection",style = "info",
             numericInput("deepSplit","Deep split",2,min = 0,max = 4,step =1),
             bsTooltip("deepSplit", "Provides a simplified control over how sensitive module detection should be to module splitting, with 0 least and 4 most sensitive' ","left"),

             numericInput("detectCutHeight","Maximum joining heights",0.995,min = 0,max = 1,step =0.001),
             bsTooltip("detectCutHeight", "For method=='tree' it defaults to 0.99. For method=='hybrid' it defaults to 99% of the range between the 5th percentile and the maximum of the joining heights on the dendrogram' ","left"),
             numericInput("minModuleSize","Minimum Module size",30,min = 1,max = 1000,step =5),
             bsTooltip("minModuleSize", "Minimum module size for module detection' ","left"),
             numericInput("reassignThreshold","P value  threshold", 0,min = 0,max = 1,step =0.0001),
             bsTooltip("reassignThreshold", "p-value ratio threshold for reassigning genes between modules' ","left"),
             numericInput("mergeCutHeight","Threshold to Merge Modules", 0.25,min = 0,max = 1,step =0.0001),
             bsTooltip("mergeCutHeight", "Dendrogram cut height for module merging.","left"),
             checkboxInput("numericLabels", "Should modules be labeled numbers ?", value = T, width = NULL),
             bsTooltip("numericLabels", "Should the returned modules be labeled by colors (FALSE), or by numbers (TRUE)?","left"),
             checkboxInput("pamStage", "PAM-like stage ?", value = F, width = NULL),
             bsTooltip("pamStage", "Only used for method 'hybrid'. If TRUE, the second (PAM-like) stage will be performed.","left"),
             checkboxInput("pamRespectsDendro", "PAM stage will respect the dendrogram ?", value = F, width = NULL),
             bsTooltip("pamRespectsDendro", "Only used for method 'hybrid'. If TRUE, the PAM stage will respect the dendrogram in the sense that objects and small clusters will only be assigned to clusters that belong to the same branch that the objects or small clusters being assigned belong to.","left")
             ),

             bsCollapsePanel("Size control",style = "info",
               bsCollapse(id = "WGCNASizecontrol",open=NULL,
              bsCollapsePanel("Soft threshold selection",style = "info",
                 sliderInput("sftwidth", "Heatmap Width (%)", min = 0, max = 100, value = 50),
                 sliderInput("sftheight", "Heatmap Height (px)", min = 0, max = 1000, value = 600)),
              bsCollapsePanel("Module dendrogram",style = "info",
                sliderInput("MDwidth", "Heatmap Width (%)", min = 0, max = 100, value = 50),
                sliderInput("MDheight", "Heatmap Height (px)", min = 0, max = 1000, value = 450))


             )))
             ,

             actionButton("WGCNAbt",
                          "Submit",
                          style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                          icon = icon("picture-o"))

   )
  )

)
)),
 hidden(div(id = "wgcna2_wrapper",
fluidPage(column(9,
                   box(
                     title = "Module-trait relationships",
                     width = NULL,
                     solidHeader = T,
                     collapsible = T,
                     collapsed = F,
                     status = "success",
                     align="center",
                     tabBox(
                       width=NULL,
                       tabPanel(title="Module - trait relationships",value="MTR",uiOutput('MTRplot'),
                                fluidRow(column(4,align="left",
                                                hidden(div(id = "MTR_wrapper",
                                                           splitLayout(
                                                             numericInput("MTRwidthdl","Figure width",value = 10),
                                                             numericInput("MTRheightdl","Figure height",value = 10)),
                                                           downloadButton('downloadMTR', 'Download figure', class = "butt2")
                                                ))
                                ))
                                ),
                       tabPanel(title="GS vs MM",value="GM",uiOutput('GSplot'),
                                fluidRow(column(4,align="left",
                                                hidden(div(id = "GM_wrapper",
                                                           splitLayout(
                                                             numericInput("GMwidthdl","Figure width",value = 10),
                                                             numericInput("GMheightdl","Figure height",value = 10)),
                                                           downloadButton('downloadGM', 'Download figure', class = "butt2")
                                                ))
                                ))
                                ),
                       tabPanel(title="Output genes",value="netsum",DT::dataTableOutput('netsum'),
                                fluidRow(column(4,align="left",
                                                hidden(div(id = "netsum_wrapper",
                                                           downloadButton('downloadnetsum', 'Download table', class = "butt2")
                                                ))
                                ))
                                )
                     )
                   )
                   ),
            column(3,
            box(title="Module-trait relationships",
                solidHeader=T,
                collapsible=T,width=NULL,collapsed =F,status = "danger",
                selectizeInput(
                  'trait',
                  "Clinical trait",
                  choices = NULL,
                  selected = NULL,
                  multiple=F
                ),
                bsTooltip("trait", "Select the clinical trait you are interested and to identify correlated modules and genes","left"),
                selectizeInput(
                  'module',
                  "Module",
                  choices = NULL,
                  selected = NULL,
                  multiple=F
                ),
                bsTooltip("module", "Select the module you are interested","left"),
                selectizeInput(
                  'output',
                  "Output ",
                  choices = NULL,
                  selected = NULL,
                  multiple=F
                ),
                bsTooltip("output", "Output genes in modules for downstream analysis","left"),
                box(title="Size control",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                  box(title="Module trait relationships",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                      sliderInput("MTRwidth", "Heatmap Width (%)", min = 0, max = 100, value = 50),
                      sliderInput("MTRheight", "Heatmap Height (px)", min = 0, max = 1000, value = 750),
                      sliderInput("BMar", "Bottom margin", min = 0, max = 20, value = 5),
                      sliderInput("LMar", "Left margin", min = 0, max = 20, value = 5),
                      sliderInput("TMar", "Top margin", min = 0, max = 20, value = 2),
                      sliderInput("RMar", "Right margin", min = 0, max = 20, value = 2)
                  ),
                  box(title="GS vs MM",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                      sliderInput("GVMwidth", "Heatmap Width (%)", min = 0, max = 100, value = 50),
                      sliderInput("GVMheight", "Heatmap Height (px)", min = 0, max = 1000, value = 600)
                  )
                ),
                actionButton("MTRbt",
                             "Submit",
                             style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                             icon = icon("picture-o"))
             )
            )
          )
)
),

fluidRow(
  div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
      actionButton(inputId = 'page_before_WGCNA',label = '',icon = icon('arrow-left'),
                   style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
      HTML('<i>Data input</i>')
  ),
  div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
      HTML('<i>Dimensionality redction: Survival related genes</i>'),
      actionButton(inputId = 'page_after_WGCNA',label = '',icon = icon('arrow-right'),
                   style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
  )
)

)

