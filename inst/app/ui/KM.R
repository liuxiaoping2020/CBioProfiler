tabItem(tabName = "KM",
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
     bsAlert("kmmess"),
     bsCollapse(id = "collapseKM", open = "Descriptions and parameters",
                bsCollapsePanel("Descriptions and parameters",  includeHTML("KM.html"), style = "default"),
                bsCollapsePanel("Kaplan-Meier plot",style = "default",
                                
                                
          uiOutput('KMplot',align = "center")#,

        ,
     fluidRow(column(4,align="left",
                     useShinyjs(),
                     hidden(div(id = "km_wrapper",

                                 splitLayout(


                                            numericInput("kmwidth","Figure width",value = 10),
                                            numericInput("kmheight","Figure height",value = 10)),

                                downloadButton('downloadKM', 'Download figure', class = "butt2")
                      )
                     )
                    )
                   )
     

  )
  
   )
    ),
        column(3,box(
          title = "Kaplan-Meier plot",
          width = NULL,
          status = "danger",
          solidHeader = T,
          collapsible = T,

          selectizeInput(
                  "KMgene",
                  label = "Official gene symbol",
                  choices = NULL,
                  multiple = F,
                  selected = "TP53"
          ),
          bsTooltip("KMgene", "Input a gene with official gene symbol","left"),
          selectizeInput(
                  "kmgroupby",
                  label = "Group by",
                  choices = c("Percentage","Value"),
                  multiple = F,
                  selected = "Percentage"
          ),

          conditionalPanel(
                  condition = "input.kmgroupby == 'Percentage'",
                  sliderInput("genecut", "Cutoff percentage", min = 0, max = 1, value = 0.5),
                  bsTooltip("genecut", "Input a cutoff (range:0-1) to categorize the samples into low and high expression group regarding to your interested gene. Example if 0.25: 0%-25% = Low, 25%-100% high","left"),

          ),
          conditionalPanel(
                  condition = "input.kmgroupby == 'Value'",
                  numericInput('kmgpvalue', "Cutoff value", value=1,step = 1),
                  bsTooltip("kmgpvalue", "Input a specific cutoff value to categorize the samples into low and high expression group regarding to your interested gene.","left")


          ),
          selectizeInput(
                  "survivaltime",
                  label = "Survival time",
                  choices = NULL,
                  multiple = F,
                  selected = "OS.time"
          ),
          bsTooltip("survivaltime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),

          selectizeInput(
                  "survivalstatus",
                  label = "Survival status",
                  choices = NULL,
                  multiple = F,
                  selected = "OS"
          ),
          bsTooltip("survivalstatus", "Select survival time column. Example: OS, RFS, PFS","left"),
          box(title="Customized setting",width=NULL,solidHeader=F,collapsible=T,status=NULL,collapsed =T,
          textInput("survxlab","Label of X axis",placeholder="Overall survival (months)",value="Overall survival (months)"),
          bsTooltip("survxlab", "Define the label of X axis","left"),
          checkboxInput("survP", "Show p-value?", value = TRUE, width = NULL),
          checkboxInput("survRT", "Show risk table?", value = TRUE, width = NULL),
          checkboxInput("survCI", "Show confidence interval?", value = F, width = NULL),
          colourpicker::colourInput("kmcolor1", "Color of group 1", value = "#FC4E07", showColour = "background",closeOnClick = TRUE),
          colourpicker::colourInput("kmcolor2", "Color of group 2", value = "#00AFBB", showColour = "background",closeOnClick = TRUE),
          ),
          box(title="Size control",width=NULL,solidHeader=F,collapsible=T,status=NULL,collapsed =T,
          sliderInput("survwidth", "Plot Width (%)", min = 0, max = 100, value =50),
          sliderInput("survheight", "Plot Height (px)", min = 0, max = 1200, value = 500),
          sliderInput("kmlg", "Figure legend position", value = c(0.8, 0.90), min = 0, max = 1)
          ),
          useShinyjs(),
          actionButton(
            "KMplotbt",
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
                 actionButton(inputId = 'page_before_KM',label = '',icon = icon('arrow-left'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
                 HTML('<i>Correlation with clinical features</i>')
             ),
             div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
                 HTML('<i>CoxPH model</i>'),
                 actionButton(inputId = 'page_after_KM',label = '',icon = icon('arrow-right'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
             )
     )



     )
