tabItem(
  tabName = "bioan",
               fluidRow(column(9,
                               box(
                                 title = "Functional Enrichment analysis",
                                 width = NULL,
                                 solidHeader = T,
                                 collapsible = T,
                                 collapsed = F,
                                 align="center",
                                 status = "success",
                                 bsAlert("FEAmess"),
                                 tabBox(width = NULL,
                                        tabPanel("Enrichment table", dataTableOutput('enrichtable'),
                                                 fluidRow(column(4,align="left",
                                                                 hidden(div(id = "enrichtable_wrapper",
                                                                            downloadButton('downloadenrichtable', 'Download table', class = "butt2")
                                                                 ))
                                                 ))
                                        ),
                                        tabPanel(title="Enrichment plot",value="DEGen",uiOutput('DEGenplot'),
                                                 fluidRow(column(4,align="left",
                                                                 hidden(div(id = "DEGenplot_wrapper",
                                                                            splitLayout(
                                                                              numericInput("DEGenplotwidthdl","Figure width",value = 10),
                                                                              numericInput("DEGenplotheightdl","Figure height",value = 10)),
                                                                            downloadButton('downloadDEGenplot', 'Download figure', class = "butt2")
                                                                 ))
                                                 ))
                                        )
                                 )
                               )),

                        column(3,
                               box(title="Enrichment analysis",
                                   solidHeader=T,
                                   collapsible=T,
                                   width=NULL,
                                   collapsed = F,
                                   status = "danger",
                                   bsCollapse(id = "colEnrichment",open="Enrichment method",
                                              bsCollapsePanel("Enrichment method",style = "info",

                                                              selectizeInput(
                                                                'biomakers',
                                                                "Biomarkers",
                                                                choices = NULL,
                                                                multiple = F
                                                                ),
                                                              bsTooltip("biomakers", "'Significant DEGs' means significantly differential experssion genes at the cutoff you specified in the 'Differentially expressed genes' module; 'Significant SRGs' means significant survival-related genes at the cutoff you specified in the 'Survival related genes' module; 'Network hub genes' means genes in one or all non-grey modules you selected in the 'WGCNA' module; 'Genes from benchmark experiment' means genes derived from benchmark experiment based on cross-validation or nested cross-validation.","left"),
                                                              selectizeInput(
                                                                'FEAana',
                                                                "Functional analysis",
                                                                choices = c("GO", "KEGG", "MSigDb","Reactome Pathway"),
                                                                selected = "GO"
                                                              ),
                                                              bsTooltip("FEAana", "Define the functional enrichment analysis methods,GO, gene ontology analysis; KEGG, Kyoto Encyclopedia of Genes and Genomes analysis; MsigDb, Molecular Signatures Database analysis; Reactome Pathway, Reactome Pathway analysis","left"),

                                                              conditionalPanel(
                                                                condition = "input.FEAana == 'GO'",
                                                                selectizeInput(
                                                                  'ont',
                                                                  "GO categories",
                                                                  choices = c("BP", "CC", "MF","ALL"),
                                                                  selected = "BP"
                                                                ),
                                                                bsTooltip("ont", "Define the subcategory of gene ontology.BP, biological process; CC, cellular component;MF, molecular function; 'All' means performing functional GO analysis based on the all three subcategories","left"),
                                                                checkboxInput("readable", "Map gene ID to gene Name ?", value = F, width = NULL)
                                                              ),

                                                              selectizeInput(
                                                                'FEAmethod',
                                                                "Functional analysis method",
                                                                choices = c("ORA", "GSEA"),
                                                                selected = "ORA"
                                                              ),

                                                              conditionalPanel(
                                                                condition = "input.FEAmethod == 'GSEA'",
                                                                selectizeInput(
                                                                  'gseamethod',
                                                                  "GSEA algorithm",
                                                                  choices = c("DOSE", "fgsea"),
                                                                  selected = "fgsea"
                                                                )
                                                              ),

                                                              numericInput("minGSSize","Minimal size of genes",10,min = 1,max = 50,step =1),
                                                              bsTooltip("minGSSize", "Minimal size of genes annotated by functional term for testing","left"),
                                                              numericInput("maxGSSize","Maximal size of genes",500,min = 10,max = 100,step =1),
                                                              bsTooltip("maxGSSize", "Maximal size of genes annotated by functional term for testing","left"),
                                                              numericInput("pvalueCutoff","P value Cutoff",0.05,min = 0,max = 1,step =0.05),
                                                              bsTooltip("pvalueCutoff", "Adjusted pvalue cutoff on enrichment tests to report","left"),
                                                              selectizeInput(
                                                                'pAdjustMethod',
                                                                "P Adjust Method",
                                                                choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                                                selected = "BH"
                                                              ),
                                                              numericInput("qvalueCutoff","q value Cutoff",0.2,min = 0,max = 1,step =0.05)

                                              ),
                                              bsCollapsePanel("Visualization",style = "info",
                                                              selectizeInput(
                                                                'enplottype',
                                                                "Plot type",
                                                                choices = c("Bar plot", "Dot plot", "Gene-concept network", "Heatmap", "Enrichment Map", "Ridgeline plot", "Geneset enrichment plot1", "Geneset enrichment plot2"),
                                                                selected = "Dot plot"
                                                              ),
                                                              bsTooltip("enplottype", "Please note! Only 'Dot plot', 'Gene-concept network',  'Enrichment Map' are suitable for visualizing the results of ORA analysis, while only 'Dot plot', 'Gene-concept network', 'Heatmap', 'Enrichment Map', 'Ridgeline plot', 'Geneset enrichment plot1', 'Geneset enrichment plot2 are suitable for visualizing the result of GSEA","left"),
                                                              conditionalPanel(

                                                                condition = "input.enplottype == 'Bar plot'",

                                                                selectizeInput(
                                                                  'barcolor',
                                                                  "Color",
                                                                  choices = c('pvalue', 'p.adjust', 'qvalue'),
                                                                  selected = "p.adjust"
                                                                ),
                                                                bsTooltip("barcolor", "The value that the color of the bar map to","left"),
                                                                numericInput("barshowCategory","Show Category",5,min = 1,max = Inf,step =1),
                                                                bsTooltip("barshowCategory", "Number of categories to show","left"),
                                                                numericInput("barfont.size","Font size",12,min = 1,max = Inf,step =1),
                                                                bsTooltip("barfont.size", "Font size of the text","left"),

                                                                numericInput("barlabel_format","Label format",30,min = 1,max = Inf,step =1)

                                                              ),
                                                              conditionalPanel(

                                                                condition = "input.enplottype == 'Dot plot'",
                                                                selectizeInput(
                                                                  'dotcolor',
                                                                  "Color",
                                                                  choices = c('pvalue', 'p.adjust', 'qvalue'),
                                                                  selected = "p.adjust"
                                                                ),
                                                                selectizeInput(
                                                                  'dotxvar',
                                                                  "variable for x-axis",
                                                                  choices = c("GeneRatio", 'Count'),
                                                                  selected = "GeneRatio"
                                                                ),
                                                                bsTooltip("dotcolor", "The value that the color of the dot map to","left"),
                                                                numericInput("dotshowCategory","Show Category",5,min = 1,max = Inf,step =1),
                                                                bsTooltip("dotshowCategory", "Number of categories to show","left"),

                                                                numericInput("dotfont.size","Font size",12,min = 1,max = Inf,step =1),
                                                                bsTooltip("dotfont.size", "Font size of the text","left"),

                                                                numericInput("dotlabel_format","Label format",30,min = 1,max = Inf,step =1)


                                                              ),

                                                              conditionalPanel(

                                                                condition = "input.enplottype == 'Enrichment Map'",

                                                                numericInput("emshowCategory","Show Category",30,min = 1,max = Inf,step =1),
                                                                bsTooltip("emshowCategory", "Number of categories to show","left"),
                                                                selectizeInput(
                                                                  'emcolor',
                                                                  "Color",
                                                                  choices = c("pvalue", "p.adjust" , "qvalue"),
                                                                  selected = "p.adjust"
                                                                ),
                                                                bsTooltip("emcolor", "Variable that used to color enriched terms", "left"),
                                                                selectizeInput(
                                                                  'emlayout',
                                                                  "Layout",
                                                                  choices = c('nicely','star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' , 'lgl'),
                                                                  selected = "nicely"
                                                                ),
                                                                bsTooltip("emlayout", "Set the layout of the plot", "left"),
                                                                numericInput("emmin_edge","Min edge",0.2,min = 0,max = 1,step =0.1),
                                                                bsTooltip("emmin_edge", "Minimum percentage of overlap genes to display the edge, should between 0 and 1, default value is 0.2", "left"),
                                                                numericInput("emcex_label_category","Cex label category",1,min = 1,max = Inf,step =1),
                                                                bsTooltip("emcex_label_category", "Scale of category node label size", "left"),
                                                                numericInput("emcex_category","Cex category",1,min = 1,max = Inf,step =1),
                                                                bsTooltip("emcex_category", "Number indicating the amount by which plotting category nodes should be scaled relative to the default", "left"),
                                                                numericInput("emcex_line","cex line",1,min = 1,max = Inf,step =1),
                                                                bsTooltip("emcex_line", "scale of line width", "left")

                                                              ),

                                                              conditionalPanel(

                                                                condition = "input.enplottype == 'Gene-concept network'",

                                                                numericInput("cneshowCategory","Show Category",5,min = 1,max = Inf,step =1),
                                                                bsTooltip("cneshowCategory", "Number of categories to show","left"),
                                                                checkboxInput("cnecolorEdge", "Coloring Edge ?", value = T, width = NULL),
                                                                bsTooltip("cnecolorEdge", "whether coloring edge by enriched terms","left"),
                                                                checkboxInput("cnecircular", "Circular layout ?", value = F, width = NULL),
                                                                bsTooltip("cnecircular", "whether using circular layout ?","left"),
                                                                selectizeInput(
                                                                  'cnenode_label',
                                                                  "Node label",
                                                                  choices = c('category', 'gene', 'all' , 'none'),
                                                                  selected = "all"
                                                                ),
                                                                bsTooltip("cnenode_label", "select which labels to be displayed", "left"),
                                                                numericInput("cnecex_category","Cex category",1,min = 1,max = Inf,step =1),
                                                                bsTooltip("cnecex_category", "Number indicating the amount by which plotting category nodes should be scaled relative to the default.", "left"),
                                                                numericInput("cnecex_gene","Cex gene",1,min = 1,max = Inf,step =1),
                                                                bsTooltip("cnecex_gene", "Number indicating the amount by which plotting gene nodes should be scaled relative to the default", "left"),
                                                                numericInput("cnecex_label_category","Cex label category",1,min = 1,max = Inf,step =1),
                                                                bsTooltip("cnecex_label_category", "Scale of category node label size", "left"),
                                                                numericInput("cnecex_label_gene","Cex label gene",1,min = 1,max = Inf,step =1),
                                                                bsTooltip("cnecex_label_gene", "scale of gene node label size", "left")

                                                              )
                                                              ,
                                                              conditionalPanel(

                                                                condition = "input.enplottype == 'Heatmap'",

                                                                numericInput("heatshowCategory","Show Category",5,min = 1,max = Inf,step =1),
                                                                bsTooltip("heatshowCategory", "Number of categories to show","left")

                                                              ),
                                                              conditionalPanel(

                                                                condition = "input.enplottype == 'Ridgeline plot'",

                                                                numericInput("regshowCategory","Show Category",30,min = 1,max = Inf,step =1),
                                                                bsTooltip("regshowCategory", "Number of categories to show","left"),
                                                                selectizeInput(
                                                                  'regfill',
                                                                  "Fill",
                                                                  choices = c("pvalue", "p.adjust" , "qvalue"),
                                                                  selected = "p.adjust"
                                                                ),
                                                                bsTooltip("regfill", "Variable that used to color enriched terms", "left"),
                                                                checkboxInput("regcore_enrichment", "Coloring Edge ?", value = T, width = NULL),
                                                                bsTooltip("regcore_enrichment", "Whether only using core_enriched genes", "left"),
                                                                numericInput("reglabel_format","Label format",30,min = 1,max = Inf,step =1),
                                                                bsTooltip("reglabel_format", "A numeric value sets wrap length, alternatively a custom function to format axis labels", "left")

                                                              ),
                                                              conditionalPanel(

                                                                condition = "input.enplottype == 'Geneset enrichment plot1'",

                                                                numericInput("gsgeneSetID1","GeneSet ID",1,min = 1,max = 500,step =1),
                                                                bsTooltip("gsgeneSetID1", "The numeric geneSet ID", "left"),
                                                                selectizeInput(
                                                                  'gsby',
                                                                  "By",
                                                                  choices = c("runningScore", "preranked", "all"),
                                                                  selected ="all"
                                                                ),
                                                                bsTooltip("gsby", "plot the geneset enrichment plot by runningScore, preranked or both", "left"),

                                                                colourpicker::colourInput("gscolor", "Color", value = "black", showColour = "background",closeOnClick = TRUE),
                                                                bsTooltip("gscolor", "Color of line segments", "left"),
                                                                colourpicker::colourInput("gscolor.line", "Score Line color", value = "green", showColour = "background",closeOnClick = TRUE),
                                                                bsTooltip("gscolor.line", "Color of running enrichment score line", "left"),
                                                                colourpicker::colourInput("gscolor.vline", "Vertical line color", value = "#FA5860", showColour = "background",closeOnClick = TRUE),
                                                                bsTooltip("gscolor.vline", "color of vertical line which indicating the maximum/minimal running enrichment score", "left")

                                                              ),
                                                              conditionalPanel(

                                                                condition= "input.enplottype == 'Geneset enrichment plot2'",

                                                                numericInput("gsegeneSetID2","GeneSet ID",1,min = 1,max = Inf,step =1),
                                                                bsTooltip("gsegeneSetID2", "The numeric geneSet ID", "left"),
                                                                colourpicker::colourInput("gsecolor", "Color", value = "black", showColour = "background",closeOnClick = TRUE),
                                                                bsTooltip("gsecolor", "Color of line segments", "left"),
                                                                numericInput("gsebase_size","Font size",11,min = 1,max = Inf,step =1),
                                                                bsTooltip("gsebase_size", "Base font size", "left"),
                                                                checkboxInput("gsepvalue_table", "Add pvalue table ?", value = F, width = NULL),
                                                                selectizeInput(
                                                                  'gseES_geom',
                                                                  "Enrichment score geom",
                                                                  choices = c("line", "dot"),
                                                                  selected ="line"
                                                                )

                                                              ),
                                                              box(title="Size control",solidHeader=F,collapsible=T,width=NULL,collapsed = T,
                                                                  sliderInput("enriwidth", "Heatmap Width (%)", min = 0, max = 100, value = 50),
                                                                  sliderInput("enriheight", "Heatmap Height (px)", min = 0, max = 1000, value = 430)

                                                              )

                                              )
                                   ),
                                   actionButton("enrichbt",
                                                "Submit",
                                                style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                                icon = icon("picture-o"))
                               ))),
  fluidRow(
    div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
        actionButton(inputId = 'page_before_bioan',label = '',icon = icon('arrow-left'),
                     style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
        HTML('<i>Correlation with stemness score</i>')
    )




    )
)
