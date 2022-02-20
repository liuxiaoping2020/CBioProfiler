

tabItem(
  tabName = "welcome1",
  fluidRow(column(
    12,
    HTML('<img style="width: 25%; display: block; margin-left: auto; margin-right: auto;" src="logo.png"/>'),
    # br(),
    HTML('<div align="justify">
                          <h3>Introduction </h3>
                          CBioExplorer (Cancer Biomarker Explorer) was developed to facilitate researchers and clinicians to screen, characterize, annotate and translate cancer biomarkers from molecular level to clinical settings more comfortably with graphical user interfaces (GUI). The whole pipeline of CBioExplorer includes data collection, data curation, dimensionality reduction using three methods of WGCNA, univariate Cox proportional hazards regression model, differentially expressed gene analysis, benchmark experiment with 6 machine learning learners (Lasso, Ridge, Elastic net, Glmboost, Coxboost, Randomforest) using cross validation and nested cross validation based on R package mlr, prediction model construction using Cox proportional hazards regression model and nomogram, clinical annotation using a variety of clinical approaches, and biological annotation using over-representation analysis (ORA) and gene set enrichment analysis (GSEA). The overview of CBioExplorer is summarized below:
                          </div>'),
    br(),
    HTML('<img style="width: 85%; display: block; margin-left: auto; margin-right: auto;" src="overview.png"/>'),
    HTML('<div align="justify">
                          <h3>Notes</h3>
                          <ul>
                          <li>Thanks for considering CBioExplorer for your study. Due to the limited computing power of the server, when multiple users use CBioExplorer at the same time, the response of the program may become slower, please be patient and only click the button once and wait until one step done.</li>
                          <li>If you plan to use CBioExplorer for high-iterative nested cross validation calculations, we strongly recommend that you download the CBioExplorer source code to your R software for corresponding calculations. This is caused by the limited computing power of the server. We apologize to you for this.</li>
                          <li>The App will be disconnected from our server after an hour if there is no mouse action on the web browser.</li>
                          </ul>
                          </div>'),


    h3("Citation"),
    p("To be added"),
    h3("Licence"),
    p("Open source under GPLV3.0. Both CBioExplorer software and curated cancer gene expression data are free for academic and non-commercial use."),

    h3("Issue report & feedback"),
    p("The user guide for CBioExplorer can be found "
      , a("here.", href="https://github.com/liuxiaoping2020/CBioExplorer_tutorial/blob/main/CBioExplorer_tutorial.pdf")
      , style="padding-left: 0em"),
    p("If you have questions to raise or are experiencing difficulties using the CBioExplorer, please report it at "
      , a("CBioExplorer Github issues.", href="https://github.com/liuxiaoping2020/CBioExplorer/issues")
      , style="padding-left: 0em"),

    h3("Contact"),
    p("This web app is developed and maintained by Xiao-Ping Liu at Wang's Lab, Department of Urology, Zhongnan hospital of WUhan University"),
    p("If you have any questions, comments, or suggestions, please feel free to contact the developer at liuxiaoping@whu.edu.cn"),
    br(),

    p("All rights reserved.",style = "display: flex; flex-direction: column; align-items: center; width: 100%")





  )
  ),
  fluidRow(
    div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
        HTML('<i>Data input</i>'),
        actionButton(inputId = 'page_after_introduction',label = '',icon = icon('arrow-right'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
    )
  )
)




