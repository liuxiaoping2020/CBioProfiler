tabItem(
  tabName = "nestr",
  fluidPage(
    column(9,
           box(
    title = "Benchmark experiment output",
    width = NULL,
    height=NULL,
    solidHeader = T,
    collapsible = T,
    collapsed = F,
    align="center",
    status = "success",
    bsAlert("nestrmess"),
    tabBox(width = NULL,
    tabPanel("Cindex comparison", align="center",
     uiOutput('modelcomp'),
     textOutput(outputId = "modeldescrip"),
    fluidRow(column(4,align="left",
                    hidden(div(id = "modelcomp_wrapper",
                               splitLayout(
                                 numericInput("modelcompwidthdl","Figure width",value = 10),
                                 numericInput("modelcompheightdl","Figure height",value = 10)),
                               downloadButton('downloadmodelcomp', 'Download figure', class = "butt2")
                    ))
    ))
  ),
    tabPanel("Biomarker output", dataTableOutput('biomarkout'),
             useShinyjs(),
             fluidRow(column(4,align="left",
                             hidden(div(id = "biomarkout_wrapper.table",
                                        downloadButton('savebiomarkout', 'Download table', class = "butt2")
                             ))
             ))
    )
    
           )
  )
                   ),
  column(3,
        box(title="Benchmark experiment setting",
            solidHeader=T,
            collapsible=T,
            width=NULL,
            collapsed = F,
            status = "danger",
        selectizeInput(
          'Learner',
          "Survival learner",
          choices = c("Lasso","ElasticNet","Ridge","Glmboost","Coxboost","RandomForestSRC"),
          selected ="Lasso",
          multiple=T
        ),
    bsTooltip("Learner", "Select the survivals you want to train you data","left"),
    selectizeInput(
      'validat',
      "Method",
      choices = c("Cross-validation","Nested cross-validation"),
      selected ="Cross-validation",
      multiple=F
    ),
    bsTooltip("validat", "Select a method for model and feature selection","left"),
    textInput("seed","Random Seed",value="123",placeholder = "1234"),
    bsTooltip("seed", "Set the random seed for benchmark experiment","left"),
    selectizeInput("nrtime",label = "Survival time",choices = NULL,multiple = F,selected = "OS.time"),
    bsTooltip("nrtime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),
    selectizeInput("nrstatus",label = "Survival status",choices = NULL,multiple = F,selected = "OS"),
    bsTooltip("nrstatus", "Select survival time column. Example: OS, RFS, PFS","left"),
    selectizeInput(
      'optimization',
      "Optimization algorithm",
      choices = c("Grid search","Random search"),
      selected ="Grid search",
      multiple=F
    ),
    bsTooltip("optimization", "Choose an optimization algorithm to search an appropriate set of parameters from a given search space for given learners.","left"),
    conditionalPanel(
      condition = "input.optimization == 'Grid search'",
      numericInput("resolution","Resolution",10,min = 5,max = 50,step =1),
      bsTooltip("resolution", "Resolution of the grid for each numeric/integer parameter in par.set. For vector parameters, it is the resolution per dimension. Either pass one resolution for all parameters, or a named vector. See ParamHelpers::generateGridDesign. Default is 10","left")
    ),
    conditionalPanel(
      condition = "input.optimization == 'Random search'",
      numericInput("maxit","Maxit",100,min = 10,max = 300,step =5),
      bsTooltip("maxit", "Number of iterations for random search. Default is 100.","left")
    ),
    numericInput("ifold","Inner fold",5,min = 3,max = 20,step =1),
    bsTooltip("ifold", "Set the fold number for inner K-fold cross-validation","left"),
    conditionalPanel(
      condition = "input.validat == 'Nested cross-validation'",
    numericInput("ofold","Outer fold",5,min = 3,max = 20,step =1),
    bsTooltip("ofold", "Set the fold number for outer K-fold cross-validation","left")
    ),
    conditionalPanel(
      condition = "input.validat == 'Cross-validation'",
      numericInput("biter","Bootstrap iterations",100,min = 50,max = 1000,step =100),
      bsTooltip("biter", "Set the iteration for bootstrap validation","left")
    ),

    checkboxInput("split", "Split the dataset ?", value = T, width = NULL),
    bsTooltip("split", "Whether split the whole dataset into internal training and test set based on stratified sampling.","left"),
    conditionalPanel(
      condition = "input.split == true",
      numericInput("sratio","Split ratio",0.6,min = 0,max = 1,step =0.1),
      bsTooltip("sratio", "0-1. The percentage of samples that goes to training","left")
    ),
    box(title="Size control",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
        sliderInput("nrwidth", "Width (%)", min = 0, max = 100, value = 50),
        sliderInput("nrheight", "Height (px)", min = 0, max = 1000, value = 450)
    ),
    actionButton("nrbt",
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
        actionButton(inputId = 'page_before_nestr',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Dimensionality reduction: Differentially expressed genes</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Prediction model: Construct model</i>'),
        actionButton(inputId = 'page_after_nestr',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")

    )
  )


)
