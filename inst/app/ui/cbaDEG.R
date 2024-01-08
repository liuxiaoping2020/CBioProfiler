tabItem(
  tabName = "cbaDEG",
  fluidRow(
    column(9,

           bsAlert("cbaDEGmsg"),
           bsCollapse(id = "cbacollapseDEG", open = "Descriptions and parameters for differentially expressed gene analysis",
                      bsCollapsePanel("Descriptions and parameters for differentially expressed gene analysis",  includeHTML("cbaDEG.html"), style = "default"),
                      bsCollapsePanel("Differentially expressed gene table",style = "default",
                                      DT::dataTableOutput('cbaDEGtable'),
                                      fluidRow(column(4,align="left",
                                                      hidden(div(id = "cbaDEG_wrapper",
                                                                 
                                                                 downloadButton('cbadownloadDEGtable', 'Download table', class = "butt2")
                                                      ))
                                      ))
                      ))
    ),
    
    column(3,box(title="Differentially expressed gene",width=NULL,
                 solidHeader=T,
                 collapsible=T,status="danger",collapsed = F,
                 selectizeInput(
                   "cbaDEGcor",
                   label = "Select the cohort",
                   choices = c("Training set","Validation set"),
                   selected="Training set",
                   multiple = F
                 ),
                 bsTooltip("cbaDEGcor", "Select the cohort (training set or validation set) to characterize the subtype","left"),
                 

                 selectizeInput('cbaDEGmethod', "Method", choices ="limma", selected = "limma"),
                 bsTooltip("cbaDEGmethod", "'limma' means conducting moderated contrast t-test for each gene in limma","left"),

                 # selectizeInput('cbaDEGfactor', "Factor/Condition", choices =NULL , selected = "Sex"),
                 # bsTooltip("cbaDEGfactor", "A categorical variable to define the group for DEG analysis","left"),
                 
                 selectizeInput('cbaDEGpadj', "P adjust methods", choices =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none") , selected = "fdr"),
                 bsTooltip("cbaDEGpadj", "Specify the correction method for P values. For more details, please refer to https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust","left"),

                 actionButton("cbaDEGbt",
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
  
  hidden(div(id = "cbaDEG1_wrapper",
             fluidRow(column(9,
                             bsAlert("cbaDEGvismess"),

                             bsCollapse(id = "cbacollapseDEGvis", open = "Descriptions and parameters for DEG output",
                                        bsCollapsePanel("Descriptions and parameters for DEG output",  includeHTML("cbaDEGoutput.html"), style = "default"),
                                        bsCollapsePanel("DEG output",style = "default",
                                                        tabBox(
                                                          width=NULL,
                                                          tabPanel(title="Output genes",
                                                                   value="cbaDEGog",
                                                                   align="center",
                                                                   DT::dataTableOutput('cbaDEGog'),
                                                                   fluidRow(column(4,align="left",
                                                                                   hidden(div(id = "cbaDEGog_wrapper",
                                                                                              
                                                                                              downloadButton('cbadownloadopgene', 'Download table', class = "butt2")
                                                                                   ))
                                                                   ))
                                                          ),
                                                          tabPanel(title="DEG Visualization",
                                                                   value="cbaDEGvis",
                                                                   align="center",
                                                                   uiOutput('cbaDEGvis'),
                                                                   fluidRow(column(4,align="left",
                                                                                   hidden(div(id = "cbaDEGvis_wrapper",
                                                                                              splitLayout(
                                                                                                numericInput("cbaDEGviswidthdl","Figure width",value = 5),
                                                                                                numericInput("cbaDEGvisheightdl","Figure height",value = 5)),
                                                                                              downloadButton('cbadownloadDEGvis', 'Download figure', class = "butt2")
                                                                                   ))
                                                                   ))
                                                                   
                                                          )
                                                          
                                                        )
                                                        
                                                        
                                        )
                             ))
                      ,column(3,
                              box(title="DEG output",solidHeader=T,collapsible=T,width=NULL,collapsed = F,status = "danger",
                                  selectizeInput(
                                    'cbamlDEGsel',
                                    "DEG cutoff method",
                                    choices = c("Adjusted P", "LogFC", "Adjusted P & LogFC"),
                                    selected = "LogFC"
                                  ),
                                  conditionalPanel(
                                    condition = "input.cbamlDEGsel == 'Adjusted P'",
                                    numericInput("cbamlgeneselp","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005)),
                                  conditionalPanel(
                                    condition = "input.cbamlDEGsel == 'LogFC'",
                                    numericInput("cbamlgeneselogFC","LogFC cutoff",2,min = 0,max = 100,step =0.5)
                                  ),
                                  conditionalPanel(
                                    condition = "input.cbamlDEGsel == 'Adjusted P & LogFC'",
                                    numericInput("cbamlgeneselp1","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005),
                                    numericInput("cbamlgeneselogFC1","LogFC cutoff",2,min = 0,max = 100,step =0.5)
                                  ),
                                  selectizeInput(
                                    'cbaDEGvismeth',
                                    "Visualization method",
                                    choices = c("Heatmap", "Volcano plot", "MA plot","Adjusted P plot"),
                                    selected = "Heatmap"
                                  ),
                                  bsTooltip("cbaDEGvismeth", "Choose visualization method for DEGs, which includes Heatmap, Volcano plot, MA plot,Adjusted P plot.","left"),
                                  
                                  conditionalPanel(
                                    condition= "input.cbaDEGvismeth == 'Heatmap'",
                                    
                                    bsCollapse(id = "cbacollheat",open=NULL,
                                               bsCollapsePanel("Basic setting",style = "info",
                                                               textInput("cbaheatname", "Heatmap name", value = "Expression", width = NULL, placeholder = NULL),
                                                               colourpicker::colourInput("cbaheatcol1", "Min colour", "blue",closeOnClick = TRUE,allowTransparent = TRUE),
                                                               colourpicker::colourInput("cbaheatcol2", "Median colour", "white",closeOnClick = TRUE,allowTransparent = TRUE),
                                                               colourpicker::colourInput("cbaheatcol3", "Max colour", "red",closeOnClick = TRUE,allowTransparent = TRUE),
                                                               bsTooltip("cbaheatname", "Name of the heatmap. By default the heatmap name is used as the title of the heatmap legend.","left"),
                                                               checkboxInput("cbaheatnorm", "Normalize the heatmap ?", value = TRUE, width = NULL),
                                                               conditionalPanel(
                                                                 condition = "input.cbaheatnorm == true ",
                                                                 selectizeInput('cbaDEGscale', "Normalization method", choices =c("Scale","Center","Log","Z-score","0-1 normalization") , selected = "Scale"),
                                                                 bsTooltip("cbaDEGscale", "Select a normalzation method for the heatmap ","left")
                                                               )),
                                               bsCollapsePanel("Row setting",style = "info",
                                                               checkboxInput("cbaClusterR", "Cluster on rows ?", value = TRUE, width = NULL),
                                                               conditionalPanel(
                                                                 condition = "input.cbaClusterR == true",
                                                                 checkboxInput("cbacluster_row_slices", "Cluster on row slice ?", value = TRUE, width = NULL),
                                                                 bsTooltip("cbacluster_row_slices", "If rows are split into slices, whether perform clustering on the slice means?","left"),
                                                                 selectizeInput('cbacludistanrow', "Cluster distance on rows", choices =c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall") , selected ="euclidean"),
                                                                 selectizeInput('cbaclumethodrow', "Cluster method on rows", choices =c("ward.D", "ward.D2", "single", "complete", "average","mcquitty","median","centroid") , selected ="complete"),
                                                                 selectizeInput('cbaRdend_side', "Row dendrogram side", choices =c("right","left") , selected ="left"),
                                                               ),
                                                               checkboxInput("cbashowRname", "Show row names ?", value = TRUE, width = NULL),
                                                               conditionalPanel(
                                                                 condition = "input.cbashowRname == true",
                                                                 selectizeInput('cbaRnameside', "Row name side", choices =c("right","left") , selected ="right"),
                                                               ),
                                                               checkboxInput("cbashowFDR", "Show adjusted P value ?", value = F, width = NULL),
                                                               checkboxInput("cbashowFC", "Show logFC ?", value = F, width = NULL)
                                               ),
                                               bsCollapsePanel("column title setting",style = "info",
                                                               # uiOutput('manySliders'),
                                                               numericInput("cbacoltitsize","Column title font size",20,min = 1,max = 50,step = 5)
                                               ),
                                               bsCollapsePanel("Column  setting",style = "info",
                                                               checkboxInput("cbaClusterC", "Cluster on columns ?", value = TRUE, width = NULL),
                                                               conditionalPanel(
                                                                 condition = "input.cbaClusterC == true",
                                                                 checkboxInput("cbacluster_column_slices", "Cluster on column slices?", value = TRUE, width = NULL),
                                                                 bsTooltip("cbacluster_column_slices", "If columns are split into slices, whether perform clustering on the slice means?","left"),
                                                                 selectizeInput('cbacludistancol', "Cluster distance on rows", choices =c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall") , selected ="euclidean"),
                                                                 selectizeInput('cbaclumethodcol', "Cluster method on rows", choices =c("ward.D", "ward.D2", "single", "complete", "average","mcquitty","median","centroid") , selected ="complete"),
                                                                 selectizeInput('cbaCdend_side', "Row dendrogram side", choices =c("top","bottom") , selected ="top")
                                                               ),
                                                               checkboxInput("cbashowCname", "Show column names ?", value = F, width = NULL),
                                                               conditionalPanel(
                                                                 condition = "input.cbashowCname == true",
                                                                 selectizeInput('cbaCnameside', "Column name side", choices =c("top","bottom") , selected ="bottom"),
                                                               )
                                               )
                                               # ,
                                               #                             bsCollapsePanel("Size control",style = "info",
                                               #                                             sliderInput("heatwidth", "Heatmap Width (%)", min = 0, max = 100, value = 50),
                                               #                                             sliderInput("heatheight", "Heatmap Height (px)", min = 0, max = 1000, value = 430)
                                    )
                                    # )
                                  ),
                                  conditionalPanel(
                                    condition =  "input.cbaDEGvismeth == 'Volcano plot'",
                                    bsCollapse(id = "cbacollvolca",open=NULL,
                                               bsCollapsePanel("Cutoff for significance",style = "info",
                                                               numericInput("cbavocacutp","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005),
                                                               numericInput("cbavocacutfc","LogFC cutoff",2,min = 0,max = 1,step =0.005)
                                               ),
                                               bsCollapsePanel("Axis label",style = "info",
                                                               textInput("cbavocaXlab", label = "X-axis label", value =  "Default"),
                                                               textInput("cbavocaYlab", label = "Y-axis label", value = "Default"),
                                                               selectInput("cbalegendPosition","Legend position", choices =c("top","bottom","left", "right") , selected ="top")
                                               )
                                               
                                               # bsCollapsePanel("Size control",style = "info",
                                               #                 sliderInput("Volcanowidth", "Heatmap Width", min = 0, max = 100, value = 450),
                                               #                 sliderInput("Volcanoheight", "Heatmap Height (px)", min = 0, max = 1000, value = 430)
                                    )
                                    
                                  ),
                                  conditionalPanel(
                                    condition= "input.cbaDEGvismeth == 'MA plot'",
                                    bsCollapse(id = "cbacollMA",open=NULL,
                                               bsCollapsePanel("Cutoff for significance",style = "info",
                                                               numericInput("cbaMAcutp","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005),
                                                               numericInput("cbaMAfc","LogFC cutoff",2,min = 0,max = 1,step =0.005)
                                               ),
                                               bsCollapsePanel("Gene selection",style = "info",
                                                               selectInput("cbaMASelmeth","Selection methods", choices =c("Top genes","Cutoff","Gene symbol"), selected ="Top genes"),
                                                               numericInput("cbaTopgene","Top gene",20,min = 0,max = Inf,step =1),
                                                               selectInput("cbaMAstopmeth","Select top method", choices =c("Adjusted P","LogFC"), selected ="Adjusted P"),
                                                               selectizeInput("cbaMAgenesym","Gene symbol", choices =NULL , selected =NULL,multiple = T)
                                               ),
                                               bsCollapsePanel("Axis label",style = "info",
                                                               textInput("cbaMAXlab", label = "X-axis label", value =  "Log2 mean expression"),
                                                               textInput("cbaMAYlab", label = "Y-axis label", value = "Log2 fold change"),
                                                               selectInput("cbaMAlegendPosition","Legend position", choices =c("top","bottom","left", "right") , selected ="top")
                                               )
                                               # ,
                                               # bsCollapsePanel("Size control",style = "info",
                                               #                 sliderInput("MAwidth", "MA-plot Width", min = 0, max = 1000, value = 450),
                                               #                 sliderInput("MAheight", "MA-plot Height (px)", min = 0, max = 1000, value = 430)
                                    )
                                    
                                    
                                  ),
                                  conditionalPanel(
                                    condition= "input.cbaDEGvismeth == 'Adjusted P plot'",
                                    colourpicker::colourInput("cbapadjcol", "Color of the histogram", value = "gray")
                                  ),
                                  box(title="Size control",width=NULL,solidHeader=F,collapsible=T,collapsed =T,
                                      sliderInput("cbapadjwidth", "Width (%)", min = 0, max = 100, value =50),
                                      sliderInput("cbapadjheight", "Height (px)", min = 0, max = 1000, value = 430))
                                  ,
                                  actionButton("cbaDEGvisbt",
                                               "Submit",
                                               style = "background-color: #000080;
              color: #FFFFFF;
              margin-left: auto;
              margin-right: auto;
              width: 100%",
                                               icon = icon("upload"))
                              )
                      )
             ))),
  
  fluidRow(
    div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
        actionButton(inputId = 'page_before_cbaDEG',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Time-dependent ROC</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Immune infiltration among subtypes</i>'),
        actionButton(inputId = 'page_after_cbaDEG',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
        
    )
  )
)

