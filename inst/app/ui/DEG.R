tabItem(
  tabName = "DEG",
  fluidRow(
    column(9,
           # box(
           #   title = "Differentially expressed gene table",width = NULL,align="center",
           #   bsAlert("DEGmsg"),
           #   solidHeader=T,collapsible=T,status="success", 
           bsAlert("DEGmsg"),
           bsCollapse(id = "collapseDEG", open = "Descriptions and parameters for differentially expressed gene analysis",
                      bsCollapsePanel("Descriptions and parameters for differentially expressed gene analysis",  includeHTML("DEG.html"), style = "default"),
                      bsCollapsePanel("Differentially expressed gene table",style = "default",
             DT::dataTableOutput('DEGtable'),
             fluidRow(column(4,align="left",
                             hidden(div(id = "DEG_wrapper",

                                        downloadButton('downloadDEGtable', 'Download table', class = "butt2")
                             ))
             ))
             ))
           ),

    column(3,box(title="Differentially expressed gene",width=NULL,
                 solidHeader=T,
                 collapsible=T,status="danger",collapsed = F,
                 # checkboxInput("filterDEG", "Filter low expression ?", value = F, width = NULL),
                 # conditionalPanel(condition = "input.filterDEG == true",
                 #                  numericInput(inputId = "filnumb",label = "Minimun CPM reads",value=1,min=0,max=50,step = 1),
                 #                  bsTooltip("filnumb", "Set the minimun CPM reads of genes to be retained for subsequent analysis, CPM reads row sums <n","left"),
                 #                  numericInput(inputId = "leasamp",label="Least number of samples",value=5,min=1,step=1),
                 #                  bsTooltip("leasamp", "Set the least samples with minimun CPM counts","left")
                 # ),
                 #
                 selectizeInput('DEGmethod', "Method", choices ="limma", selected = "limma"),
                 bsTooltip("DEGmethod", "'limma' means conducting moderated contrast t-test for each gene in limma","left"),
                 # conditionalPanel(
                 #   condition = "input.DEGmethod == 'edgeR-LRT' ",
                 #   selectizeInput('DEGnorm', "Normalization method", choices =c("TMM","TMMwsp","RLE","upperquartile","none"), selected = "TMM")
                 #
                 # ),
                 # conditionalPanel(
                 #   condition = "input.DEGmethod == 'edgeR-QLF'",
                 #   selectizeInput('DEGnorm', "Normalization method", choices =c("TMM","TMMwsp","RLE","upperquartile","none"), selected = "TMM")
                 #
                 # ),
                 selectizeInput('DEGfactor', "Factor/Condition", choices =NULL , selected = "Sex"),
                 bsTooltip("DEGfactor", "A categorical variable to define the group for DEG analysis","left"),

                 selectizeInput('DEGpadj', "P adjust methods", choices =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none") , selected = "fdr"),
                 bsTooltip("DEGpadj", "Specify the correction method for P values. For more details, please refer to https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust","left"),
                 # box(
                  #  selectizeInput(
                  #    'mlDEGsel',
                  #    "DEG cutoff method",
                  #    choices = c("Adjusted P", "LogFC", "Adjusted P & LogFC"),
                  #    selected = "LogFC"
                  #  ),
                  #  conditionalPanel(
                  #    condition = "input.mlDEGsel == 'Adjusted P'",
                  #    numericInput("mlgeneselp","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005)),
                  #  conditionalPanel(
                  #    condition = "input.mlDEGsel == 'LogFC'",
                  #    numericInput("mlgeneselogFC","LogFC cutoff",2,min = 0,max = 1,step =0.005)
                  #  ),
                  #  conditionalPanel(
                  #    condition = "input.mlDEGsel == 'Adjusted P & LogFC'",
                  #    numericInput("mlgeneselp1","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005),
                  #    numericInput("mlgeneselogFC1","LogFC cutoff",0.05,min = 0,max = 1,step =0.005)
                  #  # )
                  # ),
                 actionButton("DEGbt",
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

  hidden(div(id = "DEG1_wrapper",
  fluidRow(column(9,
                  bsAlert("DEGvismess"),
                  # box(
                  #   title = "DEG output",
                  #   width = NULL,
                  #   solidHeader = T,
                  #   collapsible = T,
                  #   collapsed = F,
                  #   status = "success",
                  #   bsAlert("DEGvismess"),
                  #   align="center",
                  bsCollapse(id = "collapseDEGvis", open = "Descriptions and parameters for DEG output",
                             bsCollapsePanel("Descriptions and parameters for DEG output",  includeHTML("DEGoutput.html"), style = "default"),
                             bsCollapsePanel("DEG output",style = "default",
                    tabBox(
                      width=NULL,
                      tabPanel(title="Output genes",
                               value="DEGog",
                               align="center",
                               DT::dataTableOutput('DEGog'),
                               fluidRow(column(4,align="left",
                                               hidden(div(id = "DEGog_wrapper",

                                                          downloadButton('downloadopgene', 'Download table', class = "butt2")
                                               ))
                               ))


                               ),
                      tabPanel(title="DEG Visualization",
                               value="DEGvis",
                               align="center",
                               uiOutput('DEGvis'),
                               fluidRow(column(4,align="left",
                                               hidden(div(id = "DEGvis_wrapper",
                                                          splitLayout(
                                                            numericInput("DEGviswidthdl","Figure width",value = 10),
                                                            numericInput("DEGvisheightdl","Figure height",value = 10)),
                                                          downloadButton('downloadDEGvis', 'Download figure', class = "butt2")
                                               ))
                               ))

                               )

                    )


                  )
  ))
  ,column(3,
           box(title="DEG output",solidHeader=T,collapsible=T,width=NULL,collapsed = F,status = "danger",
               selectizeInput(
                 'mlDEGsel',
                 "DEG cutoff method",
                 choices = c("Adjusted P", "LogFC", "Adjusted P & LogFC"),
                 selected = "LogFC"
               ),
               conditionalPanel(
                 condition = "input.mlDEGsel == 'Adjusted P'",
                 numericInput("mlgeneselp","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005)),
               conditionalPanel(
                 condition = "input.mlDEGsel == 'LogFC'",
                 numericInput("mlgeneselogFC","LogFC cutoff",2,min = 0,max = 100,step =0.5)
               ),
               conditionalPanel(
                 condition = "input.mlDEGsel == 'Adjusted P & LogFC'",
                 numericInput("mlgeneselp1","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005),
                 numericInput("mlgeneselogFC1","LogFC cutoff",2,min = 0,max = 100,step =0.5)
               ),
               selectizeInput(
                 'DEGvismeth',
                 "Visualization method",
                 choices = c("Heatmap", "Volcano plot", "MA plot","Adjusted P plot"),
                 selected = "Heatmap"
               ),
               bsTooltip("DEGvismeth", "Choose visualization method for DEGs, which includes Heatmap, Volcano plot, MA plot,Adjusted P plot.","left"),

               conditionalPanel(
                 condition= "input.DEGvismeth == 'Heatmap'",

                bsCollapse(id = "collheat",open=NULL,
                          bsCollapsePanel("Basic setting",style = "info",
                                                                      textInput("heatname", "Heatmap name", value = "Expression", width = NULL, placeholder = NULL),
                                                                      colourpicker::colourInput("heatcol1", "Min colour", "blue",closeOnClick = TRUE,allowTransparent = TRUE),
                                                                      colourpicker::colourInput("heatcol2", "Median colour", "white",closeOnClick = TRUE,allowTransparent = TRUE),
                                                                      colourpicker::colourInput("heatcol3", "Max colour", "red",closeOnClick = TRUE,allowTransparent = TRUE),
                                                                      bsTooltip("heatname", "Name of the heatmap. By default the heatmap name is used as the title of the heatmap legend.","left"),
                                                                      checkboxInput("heatnorm", "Normalize the heatmap ?", value = TRUE, width = NULL),
                                                                      conditionalPanel(
                                                                        condition = "input.heatnorm == true ",
                                                                        selectizeInput('DEGscale', "Normalization method", choices =c("Scale","Center","Log","Z-score","0-1 normalization") , selected = "Scale"),
                                                                        bsTooltip("DEGscale", "Select a normalzation method for the heatmap ","left")
                                                                      )),
                                                      bsCollapsePanel("Row setting",style = "info",
                                                                      checkboxInput("ClusterR", "Cluster on rows ?", value = TRUE, width = NULL),
                                                                      conditionalPanel(
                                                                        condition = "input.ClusterR == true",
                                                                        checkboxInput("cluster_row_slices", "Cluster on row slice ?", value = TRUE, width = NULL),
                                                                        bsTooltip("cluster_row_slices", "If rows are split into slices, whether perform clustering on the slice means?","left"),
                                                                        selectizeInput('cludistanrow', "Cluster distance on rows", choices =c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall") , selected ="euclidean"),
                                                                        selectizeInput('clumethodrow', "Cluster method on rows", choices =c("ward.D", "ward.D2", "single", "complete", "average","mcquitty","median","centroid") , selected ="complete"),
                                                                        selectizeInput('Rdend_side', "Row dendrogram side", choices =c("right","left") , selected ="left"),
                                                                      ),
                                                                      checkboxInput("showRname", "Show row names ?", value = TRUE, width = NULL),
                                                                      conditionalPanel(
                                                                        condition = "input.showRname == true",
                                                                        selectizeInput('Rnameside', "Row name side", choices =c("right","left") , selected ="right"),
                                                                      ),
                                                                      checkboxInput("showFDR", "Show adjusted P value ?", value = F, width = NULL),
                                                                      checkboxInput("showFC", "Show logFC ?", value = F, width = NULL)
                                                      ),
                                                      bsCollapsePanel("column title setting",style = "info",
                                                                      # uiOutput('manySliders'),
                                                                      numericInput("coltitsize","Column title font size",20,min = 1,max = 50,step = 5)
                                                      ),
                                                      bsCollapsePanel("Column  setting",style = "info",
                                                                      checkboxInput("ClusterC", "Cluster on columns ?", value = TRUE, width = NULL),
                                                                      conditionalPanel(
                                                                        condition = "input.ClusterC == true",
                                                                        checkboxInput("cluster_column_slices", "Cluster on column slices?", value = TRUE, width = NULL),
                                                                        bsTooltip("cluster_column_slices", "If columns are split into slices, whether perform clustering on the slice means?","left"),
                                                                        selectizeInput('cludistancol', "Cluster distance on rows", choices =c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall") , selected ="euclidean"),
                                                                        selectizeInput('clumethodcol', "Cluster method on rows", choices =c("ward.D", "ward.D2", "single", "complete", "average","mcquitty","median","centroid") , selected ="complete"),
                                                                        selectizeInput('Cdend_side', "Row dendrogram side", choices =c("top","bottom") , selected ="top")
                                                                      ),
                                                                      checkboxInput("showCname", "Show column names ?", value = F, width = NULL),
                                                                      conditionalPanel(
                                                                        condition = "input.showCname == true",
                                                                        selectizeInput('Cnameside', "Column name side", choices =c("top","bottom") , selected ="bottom"),
                                                                      )
                                                      )
                          # ,
                          #                             bsCollapsePanel("Size control",style = "info",
                          #                                             sliderInput("heatwidth", "Heatmap Width (%)", min = 0, max = 1000, value = 840),
                          #                                             sliderInput("heatheight", "Heatmap Height (px)", min = 0, max = 1000, value = 430)
                                                      )
                                           # )
               ),
               conditionalPanel(
                condition =  "input.DEGvismeth == 'Volcano plot'",
                                           bsCollapse(id = "collvolca",open=NULL,
                                                      bsCollapsePanel("Cutoff for significance",style = "info",
                                                                      numericInput("vocacutp","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005),
                                                                      numericInput("vocacutfc","LogFC cutoff",2,min = 0,max = 1,step =0.005)
                                                      ),
                                                      bsCollapsePanel("Axis label",style = "info",
                                                                      textInput("vocaXlab", label = "X-axis label", value =  "Default"),
                                                                      textInput("vocaYlab", label = "Y-axis label", value = "Default"),
                                                                      selectInput("legendPosition","Legend position", choices =c("top","bottom","left", "right") , selected ="top")
                                                      )

                                                      # bsCollapsePanel("Size control",style = "info",
                                                      #                 sliderInput("Volcanowidth", "Heatmap Width", min = 0, max = 1000, value = 450),
                                                      #                 sliderInput("Volcanoheight", "Heatmap Height (px)", min = 0, max = 1000, value = 430)
                                                      )

                 ),
               conditionalPanel(
                 condition= "input.DEGvismeth == 'MA plot'",
                                           bsCollapse(id = "collMA",open=NULL,
                                                      bsCollapsePanel("Cutoff for significance",style = "info",
                                                                      numericInput("MAcutp","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005),
                                                                      numericInput("MAfc","LogFC cutoff",2,min = 0,max = 1,step =0.005)
                                                      ),
                                                      bsCollapsePanel("Gene selection",style = "info",
                                                                      selectInput("MASelmeth","Selection methods", choices =c("Top genes","Cutoff","Gene symbol"), selected ="Top genes"),
                                                                      numericInput("Topgene","Top gene",20,min = 0,max = Inf,step =1),
                                                                      selectInput("MAstopmeth","Select top method", choices =c("Adjusted P","LogFC"), selected ="Adjusted P"),
                                                                      selectizeInput("MAgenesym","Gene symbol", choices =NULL , selected =NULL,multiple = T)
                                                      ),
                                                      bsCollapsePanel("Axis label",style = "info",
                                                                      textInput("MAXlab", label = "X-axis label", value =  "Log2 mean expression"),
                                                                      textInput("MAYlab", label = "Y-axis label", value = "Log2 fold change"),
                                                                      selectInput("MAlegendPosition","Legend position", choices =c("top","bottom","left", "right") , selected ="top")
                                                      )
                                                      # ,
                                                      # bsCollapsePanel("Size control",style = "info",
                                                      #                 sliderInput("MAwidth", "MA-plot Width", min = 0, max = 1000, value = 450),
                                                      #                 sliderInput("MAheight", "MA-plot Height (px)", min = 0, max = 1000, value = 430)
                                                      )


                 ),
             conditionalPanel(
               condition= "input.DEGvismeth == 'Adjusted P plot'",
               colourpicker::colourInput("padjcol", "Color of the histogram", value = "gray")
                           ),
      box(title="Size control",width=NULL,solidHeader=F,collapsible=T,collapsed =T,
          sliderInput("padjwidth", "Width (px)", min = 0, max = 100, value =50),
          sliderInput("padjheight", "Height (px)", min = 0, max = 1000, value = 430))
          ,
             actionButton("DEGvisbt",
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
        actionButton(inputId = 'page_before_DEG',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Dimensionality reduction: Survival related genes</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Benchmark experiment</i>'),
        actionButton(inputId = 'page_after_DEG',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")

    )
  )




)

