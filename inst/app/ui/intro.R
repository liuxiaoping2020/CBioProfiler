tabItem(
  tabName = "welcome1",
  fluidRow(column(
    12,
    HTML('<img style="width: 20%; display: block; margin-left: auto; margin-right: auto;" src="lg.png"/>'),
    # br(),
#     HTML('<div align="justify">
#                           <h3>Introduction </h3>
#                           CBIoProfiler (Cancer Biomarker and subtype Explorer) was developed to facilitate researchers and clinicians to screen, characterize, annotate and translate cancer biomarkers and subtypes from molecular level to clinical settings more comfortably with graphical user interfaces (GUI), which will help implement targeted clinical diagnosis and treatment measures for different patients to achieve precision medicine.
# 
# CBIoProfiler integrated a novel R package CuratedCancerPrognosisData that reviewed, curated and integrated the gene expression data and corresponding clinical data of 47,210 clinical samples from 268 gene expression studies of 43 common blood and solid tumors.
# 
# The whole pipeline of CBIoProfiler includes two main pipelines: cancer biomarker pipeline and cancer subtype pipeline. The cancer biomarker pipeline includes 5 modules: (1) dimensionality reduction using three methods of weighted gene co-expression network analysis(WGCNA), univariate Cox proportional hazards regression model(CoxPH), differentially expressed gene (DEG) analysis, (2) benchmark experiment with 6 machine learning learners (Lasso, Ridge, Elastic net, Glmboost, Coxboost, Randomforest) using cross validation (CV) and nested cross validation (nCV) based on R package mlr, (3) prediction model construction using Cox proportional hazards regression model and nomogram, (4) clinical annotation using a variety of clinical approaches (correlation with clinical features, Kaplan-Meier curve, CoxPH model, time-dependent ROC, most correlated genes, correlation with specific gene, gene expression in different groups, correlation with immune infiltration, correlation with stemness score, correlation with ESTIMATE score, correlation with immune checkpoint, correlation with IFN-gamma score, correlation with cytolytic activity, correlation with cancer pathway, correlation with metabolism pathway, correlation with hallmark signature, correlation with drug response
#  ), and (5) biological annotation using over-representation analysis (ORA) and gene set enrichment analysis (GSEA).
#  
# The subtype pipeline includes 3 modules: (1) data preprocessing (feature selection based on variance, median absolute deviation(MAD),CoxPH model, and principal component analysis (PCA), (2) subtype identification (integration of multiple unsupervised machine learning methods (K-means clustering (K-means), hierarchical clustering, partitioning around medoids (PAM) clustering, etc.) using two popular consensus clustering methods (ConsensusClusterPlus and M3C), (3) subtype evaluation and validation.
#                           </div>'),
h3("Introduction"),
p("CBioProfiler (Cancer Biomarker and Subtype Profiler) was developed to facilitate researchers and clinicians to screen, characterize, annotate and translate cancer biomarkers and subtypes from molecular level to clinical settings more comfortably with graphical user interfaces (GUI), which will help implement targeted clinical diagnosis and treatment measures for different patients to achieve precision medicine."),
p("CBioProfiler integrated a novel R package CuratedCancerPrognosisData that reviewed, curated and integrated the gene expression data and corresponding clinical data of 47,210 clinical samples from 268 gene expression studies of 43 common blood and solid tumors. The whole pipeline of CBioProfiler includes two main pipelines: cancer biomarker pipeline and cancer subtype pipeline. "),
p("The cancer biomarker pipeline includes 5 modules: (1) dimensionality reduction using three methods of weighted gene co-expression network analysis (WGCNA), univariate Cox proportional hazards regression model(CoxPH), differentially expressed gene (DEG) analysis, (2) benchmark experiment with 6 machine learning learners (Lasso, Ridge, Elastic net, Glmboost, Coxboost, Randomforest) using cross validation (CV) and nested cross validation (nCV) based on R package mlr, (3) prediction model construction using Cox proportional hazards regression model and nomogram, (4) clinical annotation using a variety of clinical approaches (correlation with clinical features, Kaplan-Meier curve, CoxPH model, time-dependent ROC, most correlated genes, correlation with specific gene, gene expression in different groups, correlation with immune infiltration, correlation with stemness score, correlation with ESTIMATE score, correlation with immune checkpoint, correlation with IFN-gamma score, correlation with cytolytic activity, correlation with cancer pathway, correlation with metabolism pathway, correlation with hallmark signature, correlation with drug response
 ), and (5) biological annotation using over-representation analysis (ORA) and gene set enrichment analysis (GSEA)."),
p("The subtype pipeline includes 3 modules: (1) data preprocessing (feature selection based on variance, median absolute deviation (MAD),CoxPH model, and principal component analysis (PCA), (2) subtype identification (integration of multiple unsupervised machine learning methods (K-means clustering (K-means), hierarchical clustering, partitioning around medoids (PAM) clustering, etc.) using two popular consensus clustering methods (ConsensusClusterPlus and M3C), (3) subtype evaluation and validation."),

    HTML('<img style="width: 65%; display: block; margin-left: auto; margin-right: auto;" src="ov-1.png"/>'),
    HTML('<div align="justify">
                          <h3>Notes</h3>
                          <ul>
                          <li>Thanks for considering CBioProfiler for your study. Due to the limited computing power of the server, when multiple users use CBioProfiler at the same time, the response of the program may become slower, please be patient and only click the button once and wait until one step done.</li>
                          <li>If you plan to use CBioProfiler for high-iterative nested cross validation calculations, we strongly recommend that you download the CBioProfiler source code to your R software for corresponding calculations. This is caused by the limited computing power of the server. We apologize to you for this.</li>
                          <li>The App will be disconnected from our server after an hour if there is no mouse action on the web browser.</li>
                          </ul>
                          </div>'),


    h3("Citation"),
    p("To be added"),
    h3("Licence"),
    p("Open source under GPLV3.0. Both CBioProfiler software and curated cancer gene expression data are free for academic and non-commercial use."),

    h3("Issue report & feedback"),
    p("The user guide for CBioProfiler can be found "
      , a("here.", href="https://github.com/liuxiaoping2020/CBioProfiler_tutorial/blob/main/CBioProfiler_tutorial.pdf")
      , style="padding-left: 0em"),
    p("If you have questions to raise or are experiencing difficulties using the CBioProfiler, please report it at "
      , a("CBioProfiler Github issues.", href="https://github.com/liuxiaoping2020/CBioProfiler/issues")
      , style="padding-left: 0em"),

    h3("Contact"),
    p("This web app is developed and maintained by Xiaoping Liu at Wang's Lab, Department of Urology, Zhongnan Hospital of Wuhan University"),
    p("If you have any questions, comments, or suggestions, please feel free to contact the developer via liuxiaoping@whu.edu.cn"),
    br(),

   # p("All Rights Reserved.",style = "display: flex; flex-direction: column; align-items: center; width: 100%"),

# HTML('<div style="width:300px; margin:0 auto; padding:20px 0;">
#   <a target="_blank" href="http://www.beian.gov.cn/portal/registerSystemInfo?recordcode=42010602004397" style="display:inline-block; text-decoration:none; height:20px; line-height:20px;">
#   <img src="ba.png" style="float:left;"/>
#   <p style="float:left; height:20px; line-height:20px; margin:0 0 0 5px; color:#939393;">鄂公网安备 42010602004397号</p>
#   </a>
#   </div>')
# ,
# HTML('<div class="content">
#         <img src="ba.png" alt="BA Icon" class="icon">
#         <p style= "align-items: center; width: 100%;">鄂公网安备 42010602004397号 | All Right Reserved.</p>
#     </div>'),
# HTML('<div style="display: flex; justify-content: center; align-items: center;">
#         <img src="ba.png" alt="BA Icon" class="icon">
#         <p>鄂公网安备 42010602004397号 | All Right Reserved.</p>
#     </div>'),
# HTML('<div style="display: flex; justify-content: center; align-items: center;">
#         <img src="ba.png" alt="BA Icon" class="icon" style="vertical-align: middle;">
#         <p style="margin: 0;">鄂ICP备2023018207号。Copyright © Xiao-Ping Liu All Right Reserved.</p>
# 
#     </div>')
# ,
HTML('<div style="display: flex; justify-content: center; align-items: center;">
        <p style="margin: 0;">Copyright © Xiaoping Liu. &nbsp;All Right Reserved.&nbsp;</p>
        <img src="ba.png" alt="BA Icon" class="icon" style="vertical-align: middle;">
        <p style="margin: 0;">鄂ICP备2023018207号</p>
    </div>')

#   
# HTML('<img style="width: 15%; display: block; margin-left: auto; margin-right: auto;" src="lg.png"/> <br> 鄂公网安备 42010602004397号 <br> All Right Reserved.')

# HTML('<div style="width:300px; margin:0 auto; padding:20px 0;">
#        <img src="ba.png" style="float:left;"/>
#        <p style="display: flex; flex-direction: column; align-items: center; width: 100%;">鄂公网安备 42010602004397号 | All Rights Reserved.</p>
#        </a>
#        </div>')

# ,
# HTML('<div style="width:300px;margin:0 auto; padding:20px 0;">
#        <a target="_blank" href="http://www.beian.gov.cn/portal/registerSystemInfo?recordcode=42010602004397" style="display:inline-block;text-decoration:none;height:20px;line-height:20px;"><img src="" style="float:left;"/><p style="float:left;height:20px;line-height:20px;margin: 0px 0px 0px 5px; color:#939393;">鄂公网安备 42010602004397号</p></a>
#        </div>')

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





#<p style="float:left; height:20px; line-height:20px; margin:0 0 0 5px; color:#939393;">鄂公网安备 42010602004397号 | All Rights Reserved.</p>
  




