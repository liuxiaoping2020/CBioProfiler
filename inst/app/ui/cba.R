tabItem(
  tabName = "cba",
  fluidRow(
    column(9,
           bsAlert("cbamess"),
           bsCollapse(id = "collapsecba", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("cba.html"), style = "default"),
                      bsCollapsePanel("Data preprocessing result",style = "default",
                                       # tabBox(width = NULL,
                                             # tabPanel("Correlation plot", align="center",
                                             #          uiOutput('immuneplot'),
                                             #          useShinyjs(),
                                             #          fluidRow(column(4,align="left",
                                             #                          hidden(div(id = "immuneplot_wrapper",
                                             #                                     splitLayout(
                                             #                                       numericInput("immpltwidth","Figure width",value = 10),
                                             #                                       numericInput("immpltheight","Figure height",value = 10)),
                                             #                                     downloadButton('saveimune', 'Download figure', class = "butt2")
                                             #                          ))
                                             #          ))
                                             # ),
                                             # 
                                             # tabPanel("Immune cell composition matrix", dataTableOutput('immcelltable'),
                                             #          useShinyjs(),
                                             #          fluidRow(column(4,
                                             #                          hidden(div(id = "immunecell_wrapper",
                                             #                                     downloadButton('saveimmcelltable', 'Download table', class = "butt2")
                                             #                          ))
                                             #          ))
                                             # ),
                                             # tabPanel("Correlation table", 
                                                      DT::dataTableOutput('cbapreproc'),
                                                      useShinyjs(),
                                                      fluidRow(column(4,
                                                                      hidden(div(id = "cbapreproc_wrapper",
                                                                                 downloadButton('savecbapreproc', 'Download table', class = "butt2")
                                                                      ))
                                                      ))
                                             # )
                                      # )
                      ))
    ),
    column(3,box(
      title = "Data preprocessing",
      width = NULL,
      status = "danger",
      solidHeader = T,
      collapsible = T,
      collapsed = F,
      selectizeInput(
        "cbapreproctime",
        label = "Survival time",
        choices = NULL,
        multiple = F,
        selected = "OS.time"
      ),
      bsTooltip("cbapreproctime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),
      
      selectizeInput(
        "cbapreprocstatus",
        label = "Survival status",
        choices = NULL,
        multiple = F,
        selected = "OS"
      ),
      bsTooltip("cbapreprocstatus", "Select survival time column. Example: OS, RFS, PFS","left"),      
      selectizeInput(
        "fsm",
        label = "Feature selection method",
        choices = c("Variance","MAD","PCA","CoxPH"),
        multiple = F
      ),
      bsTooltip("fsm", "Specify the method for feature selection. 'Variance' means performing feature selection based on the variace of gene expression across the samples; 'MAD' means performing feature based on median absolute deviation; 'PCA' means performing feature selection based on principal component analysis; 'CoxPH' means performing feature selection based on univariate Cox proportional hazards regression model.","left"),
      conditionalPanel(
        condition = "input.fsm == 'Variance' || input.fsm ==  'MAD'",
      selectizeInput(
        "FScutmethod",
        label = "Cut method",
        choices = c("TopK" ,"Cutoff"),
        multiple = F
      ),
      bsTooltip("FScutmethod", "Method for select the most significant features based on 'Variance' or 'MAD'.","left"),
      conditionalPanel(
        condition = "input.FScutmethod == 'TopK' ",
        numericInput("FSvaluet","Top K value",1000, min = 1,max = 70000,step =100),
        bsTooltip("FSvaluet", "Define the K value for the 'TopK' method to retain K features","left")
      ),
      conditionalPanel(
        condition = "input.FScutmethod == 'Cutoff' ",
        numericInput("FSvaluec","Cutoff value",0.05,min = 0,max = 1,step =0.05),
        bsTooltip("FSvaluec", "Define the cutoff value for the 'Cutoff' method to certain features","left")
      )
      ),
      conditionalPanel(
        condition = "input.fsm == 'PCA'",
        numericInput("PC_percent","Ratio of principal component",0.95,min = 0,max = 1,step =0.1),
        bsTooltip("PC_percent", "Define the ratio (0-1) of principal component, with which ratio the features will be retain","left"),
        checkboxInput("PCAscale", "Scale the data ?", value = T, width = NULL),
        bsTooltip("PCAscale", "Whether normalize the data before performing PCA analysis.","left")
      ),
      conditionalPanel(
        condition = "input.fsm == 'CoxPH'",
        numericInput("Coxphcutoff","Cutoff value",0.05,min = 0,max = 1,step =0.05),
        bsTooltip("Coxphcutoff", "Define the cutoff value for CoxPH model to select the most survival related features","left")
      ),
      actionButton("cbeprebt",
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
  hidden(div(id = "cba_wrapper",
  fluidRow(
    column(9,
           bsAlert("cbamess1"),
           bsCollapse(id = "collapsecba1", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("cba-ide.html"), style = "default"),
                      bsCollapsePanel("Subtype identification",style = "default",
                                      tabBox(width = NULL,
                                      tabPanel("Consensus matrix plot", align="center",
                                               uiOutput('cmp'),
                                               useShinyjs(),
                                               fluidRow(column(4,align="left",
                                                               hidden(div(id = "cmp_wrapper",
                                                                          splitLayout(
                                                                            numericInput("cmpwidth","Figure width",value = 10),
                                                                            numericInput("cmpheight","Figure height",value = 10)),
                                                                          downloadButton('savecmp', 'Download figure', class = "butt2")
                                                               ))
                                               ))
                                      ),
                                      tabPanel("CDF plot", align="center",
                                               # uiOutput('silhouette'),
                                               uiOutput('cdf'),
                                               useShinyjs(),
                                               fluidRow(column(4,align="left",
                                                               hidden(div(id = "cdf_wrapper",
                                                                          splitLayout(
                                                                            numericInput("cdfwidth","Figure width",value = 10),
                                                                            numericInput("cdfheight","Figure height",value = 10)),
                                                                          downloadButton('savecdf', 'Download figure', class = "butt2")
                                                               ))
                                               ))
                                      )
                                      ,
                                      tabPanel("Silhouette plot", align="center",
                                               uiOutput('silhouette'),
                                               # uiOutput('cdf'),
                                               useShinyjs(),
                                               fluidRow(column(4,align="left",
                                                               hidden(div(id = "silhouette_wrapper",
                                                                          splitLayout(
                                                                            numericInput("Silhouettewidth","Figure width",value = 10),
                                                                            numericInput("Silhouetteheight","Figure height",value = 10)),
                                                                          downloadButton('saveSilhouette', 'Download figure', class = "butt2")
                                                               ))
                                               ))
                                      ),
                                      tabPanel("Subtype assignment",
                                      DT::dataTableOutput('cbaassign'),
                                      useShinyjs(),
                                      fluidRow(column(4,
                                                      hidden(div(id = "cbaassign_wrapper",
                                                                 downloadButton('savecbaassign', 'Download table', class = "butt2")
                                                      ))
                                      ))
                                      )
                                      )
                      ))
           
           ),
    column(3,
           box(
             title = "Subtype identification",
             width = NULL,
             status = "danger",
             solidHeader = T,
             collapsible = T,
             collapsed = F,
             selectizeInput(
               "clustermethod",
               label = "Clustering strategy",
               choices = c("Monti consensus clustering","M3C"),
               selected="Monti consensus clustering",
               multiple = F
             ),
             bsTooltip("clustermethod", "Define the clustering strategy for clustering analysis. 'Monti consensus clustering' means performing consensus clustering analysis based on the approach introduced by Monti et al.(Monti et al., 2003); 'M3C' means performing clustering based on Monte Carlo Reference-based Consensus Clustering.","left"),
             
             selectizeInput(
               "clusterAlg",
               label = "Cluster algorithm",
               choices = c("hc","pam","km"),
               multiple = F
             ),
             bsTooltip("clusterAlg", "Define the clustering algorithm consensus clustering. 'hc' hierarchical (hclust), 'pam' for paritioning around medoids, 'km' for k-means upon data matrix, or a function that returns a clustering.","left"),
             
             numericInput("ccseed","Random seed",123, min = 1,max = 10000,step =1),
             bsTooltip("ccseed", "Set the random seed to reproduce the result.","left"),
             
             numericInput("maxK","Max K",5, min = 1,max = 10,step =1),
             bsTooltip("maxK", "Maximum cluster number to evaluate","left"),
             numericInput("pItem","pItem",0.8, min = 0,max = 1,step =0.1),
             bsTooltip("pItem", "Proportion of items to sample","left"),
             numericInput("ccreps","Number of subsamples",100, min = 5,max = 1500,step =10),
             bsTooltip("ccreps", "Integer value. Number of subsamples","left"),
             
             #CC
             conditionalPanel(
               condition = "input.clustermethod == 'Monti consensus clustering'",
               selectizeInput("distance",label = "Distance", choices = c( 'pearson','spearman', 'euclidean','binary', 'maximum', 'canberra', 'minkowski' ),
                              multiple = F,selected="pearson"
               ),
             bsTooltip("distance", "character value. 'pearson': (1 - Pearson correlation), 'spearman' (1 - Spearman correlation), 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski' or custom distance function.","left"),
             numericInput("pFeature","pFeature",1, min = 0,max = 1,step =0.1),
             bsTooltip("pFeature", "Numerical value. Proportion of features to sample.","left"),
             selectizeInput("corUse",label = "Correlation use", choices = c( 'everything','pairwise.complete.obs', 'complete.obs' ),
               multiple = F
             ),
             bsTooltip("corUse", "Optional character value. specifies how to handle missing data in correlation distances 'everything','pairwise.complete.obs', 'complete.obs' see cor() for description.","left")
             ),
             #M3C
             # iters
             conditionalPanel(
               condition = "input.clustermethod == 'M3C'",
             numericInput("m3citers","iteration",25, min = 5,max = 100,step =5),
             bsTooltip("m3citers", "Numerical value: how many Monte Carlo iterations to perform (default: 25, recommended: 5-100)","left"),
             
             selectizeInput("ref_method",label = "Reference method", choices = c("reverse-pca", "chol"),multiple = F,selected = "reverse-pca"),
             bsTooltip("ref_method", "Character string: refers to which reference method to use","left"),
             
             numericInput("repsref","Reference resamples",100, min = 50,max = 250,step =5),
             bsTooltip("repsref", "Numerical value: how many resampling reps to use for reference (default: 100, recommended: 100-250)","left"),
             numericInput("repsreal","Real resamples",100, min = 50,max = 250,step =5),
             bsTooltip("repsreal", "Numerical value: how many resampling reps to use for real data (default: 100, recommended: 100-250)","left"),
             numericInput("pacx1","Pacx1",0.1, min = 0.0001,max = 1,step =0.1),
             bsTooltip("pacx1", "Numerical value: The 1st x co-ordinate for calculating the pac score from the CDF (default: 0.1)","left"),
             numericInput("pacx2","Pacx2",0.9, min = 0.0001,max = 1,step =0.1),
             bsTooltip("pacx2", "Numerical value: The 2nd x co-ordinate for calculating the pac score from the CDF (default: 0.1)","left"),
             selectizeInput("objective",label = "Objective function", choices = c('PAC','entropy'),multiple = F,selected = "entropy"),
             bsTooltip("objective", "Character string: whether to use 'PAC' or 'entropy' objective function (default = entropy)","left"),
             # selectizeInput("objective",label = "Objective function", choices = c('PAC','entropy'),multiple = F,selected = "entropy"),
             selectizeInput("M3Cmethod",label = "Simulation method", choices = c(1,2),multiple = F,selected = 1),
             bsTooltip("M3Cmethod", "1 refers to the Monte Carlo simulation method, 2 to regularised consensus clustering","left"),
             
             # lambdadefault
             checkboxInput("tunelambda", "Tune lambda ?", value = T, width = NULL),
             bsTooltip("tunelambda", "Logical flag: whether to tune lambda or not.","left"),
             numericInput("lambdadefault","Default lambda value",0.1, min = 0.0001,max = 1,step =0.1),
             bsTooltip("lambdadefault", "Numerical value: if not tuning fixes the default (default: 0.1)","left")
             ),
             box(title="Size control",width=NULL,solidHeader=F,collapsible=T,status=NULL,collapsed =T,
                 sliderInput("comwidth", "Consensus matrix plot Width (%)", min = 0, max = 100, value =50),
                 sliderInput("comheight", "Consensus matrix plot Height (px)", min = 0, max = 1200, value = 500),
                 sliderInput("cdfpwidth", "CDF plot Width (%)", min = 0, max = 100, value =50),
                 sliderInput("cdfpheight", "CDF plot Height (px)", min = 0, max = 1200, value = 500),
                 sliderInput("silhouettepwidth", "CDF plot Width (%)", min = 0, max = 100, value =50),
                 sliderInput("silhouetteppheight", "CDF plot Height (px)", min = 0, max = 1200, value = 500)
             ),
             actionButton("cbabt",
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
  hidden(div(id = "cbavalidation_wrapper",
  fluidRow(
    column(9,
           bsAlert("cbamess2"),
           bsCollapse(id = "collapsecbavalidate", open = "Descriptions and parameters",
                      bsCollapsePanel("Descriptions and parameters",  includeHTML("cba-valid.html"), style = "default"),
                      bsCollapsePanel("Subtype Validation",style = "default",
                                      
                                      DT::dataTableOutput('cbavalid'),
                                      useShinyjs(),
                                      fluidRow(column(4,
                                                      hidden(div(id = "cbavalid_wrapper",
                                                                 downloadButton('savecbavalid', 'Download table', class = "butt2")
                                                      ))
                                      ))
                      ))
           
           ),
    column(3,
           box(
             title = "Subtype validation",
             width = NULL,
             status = "danger",
             solidHeader = T,
             collapsible = T,
             collapsed = F,
             selectizeInput(
               "cbavalimethod",
               label = "Method",
               choices = c("pearson", "spearman"), 
               selected="pearson",
               multiple = F
             ),
             bsTooltip("cbavalimethod", "Assign subtype on new data based on 'Spearman's' distance or 'Pearson's' distance","left"),
             actionButton("cbavalidationbt",
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
  
  fluidRow(
    div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
        actionButton(inputId = 'page_before_cba',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Biological annotation</i>')
    )
    ,
      div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
          HTML('<i>Correlation with clinical features</i>'),
          actionButton(inputId = 'page_after_cba',label = '',icon = icon('arrow-right'),
                       style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
      )

  )
)
