tabItem(
  tabName = "msurv",
  fluidRow(
    column(9,
           bsAlert("msurvmess"),
    #         box(
    # title="Univariate CoxPH table",
    # width = NULL,
    # solidHeader = T,
    # status = "success",
    # align="center",
    # bsAlert("msurvmess"),
    bsCollapse(id = "collapsemsurv", open = "Descriptions and parameters for Univariate CoxPH analysis",
               bsCollapsePanel("Descriptions and parameters for Univariate CoxPH analysis",  includeHTML("unicox.html"), style = "default"),
               bsCollapsePanel("Univariate CoxPH table",style = "default",
               DT::dataTableOutput('unicoxtable'),
    fluidRow(column(4,align="left",
                    hidden(div(id = "unicoxtable_wrapper",
                               downloadButton('downloadunicoxtable', 'Download table', class = "butt2")
                    ))
    ))
   )
  )
  )
  ,
  column(3,box(title="Univariate CoxPH",
               width = NULL,
               status = "danger",
               solidHeader = T,
               collapsible = T,
               selectizeInput('msurvtime', "Survival time",choices = NULL, selected = "OS.time"),
               bsTooltip("msurvtime", "Select the survival time column for univariate CoxPH analysis","left"),
               selectizeInput('msurvstatus', "Survival status",choices = NULL, selected = "OS"),
               bsTooltip("msurvstatus", "Select the survival Status column for univariate CoxPH analysis","left"),
               # selectizeInput('msurvcor', "Adjusting method", choices =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr") , selected = "fdr"),
               # bsTooltip("msurvcor", "Specify the correction method","left"),
               actionButton("msurvbt",
                            "Submit",
                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                            icon = icon("picture-o"))
               ))),
  hidden(div(id = "unicoxanalysis_wrapper",
  fluidRow(column(9,
                  # box(
                  #   title="Significant univariate CoxPH result",
                  #   width = NULL,
                  #   solidHeader = T,
                  #   status = "success",
                  #   align="center",
                  bsCollapse(id = "collapsemsurvsig", open = "Descriptions and parameters",
                             bsCollapsePanel("Descriptions and parameters",  includeHTML("unicoxsig.html"), style = "default"),
                             bsCollapsePanel("Significant univariate CoxPH result",style = "default",
                    tabBox(width=NULL,
                    tabPanel(title="Forestplot",
                             value="sigfor",
                             align="center",
                             uiOutput('sigfor'),
                             fluidRow(column(4,align="left",
                                             hidden(div(id = "unicoxForestplot_wrapper",
                                                        splitLayout(
                                                          numericInput("unicoxFPwidthdl","Figure width",value = 10),
                                                          numericInput("unicoxFPheightdl","Figure height",value = 10)),
                                                        downloadButton('downloadunicoxForestplot', 'Download figure', class = "butt2")
                                             ))
                             ))
                             ),
                    tabPanel(title="Significant CoxPH table",
                             value="sigunt",
                             align="center",
                             DT::dataTableOutput('sigcoxtable'),
                             fluidRow(column(4,align="left",
                                             hidden(div(id = "sigcoxtable_wrapper",
                                                        downloadButton('downloadsigcoxtable', 'Download table', class = "butt2")
                                             ))
                             ))
                             )
                    )
                  )


                  )
                  ),
           column(3,box(
             title="Significant univariate CoxPH",
             width = NULL,
             status = "danger",
             solidHeader = T,
             collapsible = T,
             selectizeInput('msurvcor', "Adjusting method", choices =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr") , selected = "fdr"),
             bsTooltip("msurvcor", "Specify the correction method","left"),
             selectizeInput('msurvsel', "Significance cutoff method", choices =c("P value", "Adjusted P value") , selected = "Adjusted P value"),
             bsTooltip("msurvsel", "Specify the method to select to most signifcant survival related genes.","left"),
             conditionalPanel(
               condition = "input.msurvsel == 'Adjusted P value'",
               numericInput("msurvap","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005)),
             conditionalPanel(
               condition = "input.msurvsel == 'P value'",
               numericInput("msurvp","P value cutoff",0.05,min = 0,max = 1,step =0.005)
             ),
             box(title="Size control",width=NULL,solidHeader=F,collapsible=T,collapsed =T,
                 sliderInput("msurvwidth", "Width (px)", min = 0, max = 100, value = 50),
                 sliderInput("msurvheight", "Height (px)", min = 0, max = 3500, value = 700),
                 sliderInput("msurvxtick", "xticks", min = 0, max = 20, value = 5),
                 bsTooltip("msurvxtick", "Specified x-axis largest tick mark","left")
             ),

             actionButton("msurvopbt",
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
        actionButton(inputId = 'page_before_unicox',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Dimensionality reduction: WGCNA</i>')
    ),
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Dimensionality reduction: Differentially expressed genes</i>'),
        actionButton(inputId = 'page_after_unicox',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")

    )
  )



)



