tabItem(
  tabName = "nomo",
  fluidRow(
    column(9,
           box(
             title = "Nomogram",
             width = NULL,height=NULL,
             solidHeader = T,
             collapsible = T,
             bsAlert("nomomess"),
             status = "success",
             align="center",
             uiOutput('nomoplot'),
             fluidRow(column(4,align="left",
                             hidden(div(id = "nomogram_wrapper",
                                        splitLayout(
                                          numericInput("nomogramwidthdl","Figure width",value = 10),
                                          numericInput("nomogramheightdl","Figure heigth",value = 6)),
                                        downloadButton('savenomogram', 'Download figure', class = "butt2")
                             ))
             ))
   )),
   column(3,
          box(
            title = "Nomogram",
            width = NULL,
            status = "danger",
            solidHeader = T,
            collapsible = T,
            selectizeInput(
              "nomotime",
              label = "Survival time",
              choices = NULL,
              multiple = F,
              selected = "OS.time"
            ),
            bsTooltip("nomotime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),

            selectizeInput(
              "nomostatus",
              label = "Survival status",
              choices = NULL,
              multiple = F,
              selected = "OS"
            ),
            bsTooltip("nomostatus", "Select survival status column. Example: OS, RFS, PFS","left"),

            selectizeInput(
              "nomovar",
              label = "Clinical variable",
              choices = NULL,
              multiple = T,
              selected = "Risk"
            ),
            bsTooltip("nomovar", "Select survival variable you want to included in the nomogram","left"),
            textInput("nomovarlab", "Variable names", placeholder="Age|Gender|Stage|Grade"),
            bsTooltip("nomovarlab", "Define the variable names you interested,the number of variable names should be identical with the number of variables you included in the CoxPH model to contrust the nomogram and use '|' to separate multiple variable names","left"),
            selectizeInput(
              "nomoypoint",
              label = "Prediction years",
              choices = NULL,
              multiple = T,
              options = list(maxItems = 3),
              selected="1-year"
            ),
            bsTooltip(
              "nomoypoint",
              "Define the time points in years you want to predict based on the nomogram. The longest time point should not exceed the max of survival (relapse) duration",
              "left"
            ),
            box(title = "Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                sliderInput("nomowidth", "Plot Width (%)", min = 0, max =100, value = 50),
                sliderInput("nomoheight", "Plot Height (px)", min = 0, max = 2500, value = 430)
            ),
            actionButton(
              "nomogrambt",
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
  hidden(div(id = "invalcal_wrapper",
             fluidRow(
    column(9,
           box(
             title = "Internal validation and calibration",
             width = NULL,height=NULL,
             solidHeader = T,
             collapsible = T,
             bsAlert("innomomess"),
             status = "success",
             align="center",
             tabBox(width = NULL,
             tabPanel("Internal validation", uiOutput('invalid'),
                      textOutput(outputId = "invaliddescrip"),
             fluidRow(column(4,align="left",
                             hidden(div(id = "invalid_wrapper",
                                        splitLayout(
                                          numericInput("invalidwidthdl","Figure width",value = 10),
                                          numericInput("invalidheightdl","Figure height",value = 10)),
                                        downloadButton('saveinvalid', 'Download figure', class = "butt2")
                             ))
             ))),
             tabPanel("Internal calibration", uiOutput('incalib'),
                      fluidRow(column(4,align="left",
                                      hidden(div(id = "incalib_wrapper",
                                                 splitLayout(
                                                   numericInput("incalibwidthdl","Figure width",value = 10),
                                                   numericInput("incalibheightdl","Figure heigth",value = 10)),
                                                 downloadButton('saveincalib', 'Download figure', class = "butt2")
                                      ))
                      )))
           )

           )
  ),column(3,
           box(
             title = "Internal validation and calibration",
             width = NULL,
             status = "danger",
             solidHeader = T,
             collapsible = T,
             selectizeInput(
               "invalidypoint",
               label = "Prediction years",
               choices = NULL,
               multiple = T,
               options = list(maxItems = 3),
               selected="1-year"
             ),
             bsTooltip(
               "invalidypoint",
               "Define the time points in years you want to predict based on the nomogram. The longest time point should not exceed the max of survival (relapse) duration",
               "left"
             ),
             numericInput("inreps","Bootstrap iterations",1000,min = 50,max=3000,step = 100),
             bsTooltip("inreps", "Number of bootstrap resamplings","left"),
             numericInput("inratio","Ratio",0.8,min = 0.1,max=1,step = 0.1),
             bsTooltip("inratio", "Ratio of resamplings","left"),
             box(title = "C-index Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                 sliderInput("invalidwidth", "Plot Width (%)", min = 0, max =100, value = 50),
                 sliderInput("invalidheight", "Plot Height (px)", min = 0, max = 2500, value = 430)
             ),
             box(title = "Calibration plot Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                 sliderInput("incalibwidth", "Plot Width (%)", min = 0, max =100, value = 50),
                 sliderInput("incalibheight", "Plot Height (px)", min = 0, max = 2500, value = 430),
                 checkboxInput("incalibsubtitle","Show subtitle ?",value = F)
             ),
             actionButton(
               "invalidbt",
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
  hidden(div(id = "exvalcal_wrapper",
             fluidRow(
               column(9,
                      box(
                        title = "External validation and calibration",
                        width = NULL,height=NULL,
                        solidHeader = T,
                        collapsible = T,
                        bsAlert("exnomomess"),
                        status = "success",
                        align="center",
                        tabBox(width = NULL,
                               tabPanel("External validation", uiOutput('exvalid'),
                                        textOutput(outputId = "exvaliddescrip"),
                                        fluidRow(column(4,align="left",
                                                        hidden(div(id = "exvalid_wrapper",
                                                                   splitLayout(
                                                                     numericInput("exvalidwidthdl","Figure width",value = 10),
                                                                     numericInput("exvalidheightdl","Figure height",value = 10)),
                                                                   downloadButton('saveexvalid', 'Download figure', class = "butt2")
                                                        ))
                                        ))),
                               tabPanel("External calibration", uiOutput('extcalib'),
                                        fluidRow(column(4,align="left",
                                                        hidden(div(id = "excalib_wrapper",
                                                                   splitLayout(
                                                                     numericInput("excalibwidthdl","Figure width",value = 10),
                                                                     numericInput("excalibheightdl","Figure height",value = 10)),
                                                                   downloadButton('saveexcalib', 'Download figure', class = "butt2")
                                                        ))
                                        )))


                        )
                      )
               ),column(3,
                        box(
                          title = "External validation and calibration",
                          width = NULL,
                          status = "danger",
                          solidHeader = T,
                          collapsible = T,

   
                          splitLayout(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                      cellWidths = c("0%","50%", "50%"),
                            selectizeInput(
                            "extvaltime",
                            label = "Survival time",
                            choices = NULL,
                            multiple = F,
                            selected = "OS.time"
                          ),
                          selectizeInput(
                            "extvalstatus",
                            label = "Survival status",
                            choices = NULL,
                            multiple = F,
                            selected = "OS"
                          )),
                          bsTooltip("extvaltime", "Select survival time column. Example: OS.time, RFS.time, PFS.time, etc.","left"),
                          bsTooltip("extvalstatus", "Select survival status column. Example: OS, RFS, PFS, etc.","left"),

                         splitLayout(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                     cellWidths = c("0%","50%", "50%"),
                            selectizeInput(
                              "exvalidypoint",
                              label = "Prediction years",
                              choices = NULL,
                              multiple = T,
                              options = list(maxItems = 3),
                              selected="1-year"
                            ),
                            selectizeInput(
                            "extvalvar",
                            label = "Clinical variable",
                            choices = NULL,
                            multiple = T,
                            selected = "Risk"
                          )
                          ),
                          bsTooltip("extvalvar", "Select survival variable you want to included","left"),
                          bsTooltip(
                            "exvalidypoint",
                            "Define the time points in years you want to predict based on the nomogram. The longest time point should not exceed the max of survival (relapse) duration",
                            "left"
                          ),
                          splitLayout(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                      cellWidths = c("0%","50%", "50%"),
                          numericInput("exreps","Boostrap iterations",1000,min = 50,max=3000,step = 100),
                          numericInput("exratio","Ratio",0.8,min = 0.1,max=1,step = 0.1)),
                          bsTooltip("exreps", "Number of resamplings","left"),
                          bsTooltip("exratio", "Ratio of resamplings","left"),

                          box(title = "C-index Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                              sliderInput("exvalidwidth", "Plot Width (%)", min = 0, max =100, value = 50),
                              sliderInput("exvalidheight", "Plot Height (px)", min = 0, max = 2500, value = 430)
                          ),
                          box(title = "Calibration plot Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                              sliderInput("excalibwidth", "Plot Width (%)", min = 0, max =100, value = 50),
                              sliderInput("excalibheight", "Plot Height (px)", min = 0, max = 2500, value = 430),
                              checkboxInput("excalibsubtitle","Show subtitle ?",value = F)
                          ),
                          actionButton(
                            "exvalidbt",
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
        actionButton(inputId = 'page_before_nomo',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Prediction model: Validate model</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Correlation with clinical features</i>'),
        actionButton(inputId = 'page_after_nomo',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )



)


