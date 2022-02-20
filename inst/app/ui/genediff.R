tabItem(
  tabName="genediff",
  fluidRow(column(9,
  # box(
  #   title="Gene expression in different groups",
  #   width = NULL,
  #   solidHeader = T,
  #   status = "success",
  #   align="center",
  bsAlert("genediffmess"),
  bsCollapse(id = "collapsegenediff", open = "Descriptions and parameters",
             bsCollapsePanel("Descriptions and parameters",  includeHTML("genediff.html"), style = "default"),
             bsCollapsePanel("Gene expression in different groups",style = "default",

    uiOutput('genediffplot',align="center"),
    useShinyjs(),
    fluidRow(column(4,align="left",
                    hidden(div(id = "genediff_wrapper",
                               splitLayout(
                                 numericInput("genediffwidthdl","Figure width",value = 10),
                                 numericInput("genediffheightdl","Figure width",value = 10)),
                               downloadButton('downloadgenediff', 'Download figure', class = "butt2")
                    ))
    ))
     )
)
),
  column(3,box(
    title="Gene expression in different groups",
    width = NULL,
    status = "danger",
    solidHeader = T,
    collapsible = T,
    selectizeInput('genediff', 'Official gene symbol', choices = NULL, selected = "A1BG"),
    bsTooltip("genediff", "Select an official gene symbol you are interested in ","left"),
    selectizeInput('genediffgroup', "Group",choices = NULL, selected = "Grade"),
    bsTooltip("genediffgroup", "Select a clinical variable to divide the gene expression into different groups","left"),
    checkboxInput("gdsubgroup", "Perform subgroup analysis ?", value = F, width = NULL),
    conditionalPanel(
      condition = "input.gdsubgroup == true ",
      selectizeInput('genediffsub', "Subgroup", choices =NULL , selected = "Sex"),
      bsTooltip("genediffsub", "Select a clinical variable as a subgroup indicator to compare the expression of the interested gene in subgroups","left")
    ),

    checkboxInput("genediffp", "Show P value ?", value = TRUE, width = NULL),
    conditionalPanel(
      condition = "input.genediffp == true",
      selectizeInput('genediffmethod', "Statistical test", choices=c("T", "Wilcoxon"), selected = NULL),
      bsTooltip("genediffmethod", "Select a statistical test for the comparisons, 'T' means 'T test', and 'wilcoxon' means 'wilcoxon test' ","left"),
      selectizeInput('comparisons', "Comparisons", choices=c("Pairwise comparison", "Comparisons to a reference"), selected = NULL),
      conditionalPanel(
        condition = "input.comparisons == 'Comparisons to a reference'",
        selectizeInput('comparefer', "Comparison reference", choices=NULL, selected = NULL),
      )
    ),
    selectizeInput('genediffplottype', "Plot type", choices=c("Box plot", "Bar plot", "Violin plot"), selected = "Box plot",multiple=F),

    box(title = "Size control",width = NULL,status = NULL,
      solidHeader = T,collapsible = TRUE,collapsed=TRUE,
          sliderInput("genediffbarwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
          sliderInput("genediffbarheight", "Plot Height (px)", min = 0, max = 1000, value = 430)
    ),

    actionButton("genediffbt",
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
      actionButton(inputId = 'page_before_genediff',label = '',icon = icon('arrow-left'),
                   style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
      HTML('<i>Correlation with specific gene</i>')
  ),
  div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
      HTML('<i>Correlation with Immune infiltration</i>'),
      actionButton(inputId = 'page_after_genediff',label = '',icon = icon('arrow-right'),
                   style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
  )
)


)
