observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'estimategene',
                       choices = {
                         if (!is.null(data)) {
                           if (class(data) ==  class(data.frame()))
                             choices <- as.character(row.names(data))
                           if (class(data) !=  class(data.frame()))
                             choices <-
                               as.character(row.names(as.data.frame(data)))
                           choices
                         }
                       },
                       server = TRUE,
                       selected = "TP53")
})


estimaterun<-eventReactive(input$estimatebt,{
  input$estimatebt
  data <- isolate({
    rawdata()
  })
  withProgress(message = "Correlating with ESTIMATE score",
               detail = "This may take a while...",
               value = 3,
               {
  require(estimate)
  expres<-data$expres

  accession<-isolate({input$accession})

  method<-isolate({input$estimatemethod})
  gene<-isolate({input$estimategene})

  # withProgress(message = "Correlating with ESTIMATE score",
  #              detail = "This may take a while...",
  #              value = 3,
  #              {
                 if(gene %in% row.names(expres)==F){
                   createAlert(
                     session,
                     "estimatemess",
                     "exampleAlert",
                     title = "Please note!",
                     style =  "danger",
                     content = paste(
                       gene, "is not found in the gene expression profile you specific. Please check!!!"
                     ),
                     append = T,
                     dismiss=T
                   )
                   return(NULL)
                 }else{
                   closeAlert(session, "estimatemess")
                 aff<-c(
                   "E_MTAB_863",
                   "E_MTAB_864",
                   "GSE58644",
                   "E_MTAB_6389",
                   "E_MEXP_3488",
                   "E_MTAB_1038",
                   "E_MTAB_3267",
                   "E_MTAB_1216",
                   "E_TABM_117",
                   "E_TABM_346",
                   "E_MTAB_1205",
                   "E_MTAB_951",
                   "E_MTAB_3892",
                   "E_MTAB_1719",
                   "E_MEXP_2780",
                   "E_TABM_1202",
                   "E_MTAB_923",
                   "E_MTAB_2435",
                   "GSE16011",
                   "E_MTAB_2768",
                   "GSE16581",
                   "GSE12276",
                   "GSE1456_GPL96",
                   "GSE3143",
                   "E_TABM_158",
                   "GSE14520",
                   "GSE16125_GPL5175",
                   "GSE30378",
                   "GSE24549_GPL5175",
                   "GSE24550_GPL5175",
                   "GSE68571",
                   "GSE30074",
                   "E_MTAB_4032",
                   "GSE62452",
                   "GSE78229",
                   "GSE42669",
                   "GSE37751",
                   "GSE28735",
                   "E_MTAB_6877",
                   "GSE12417_GPL570",
                   "GSE33371",
                   "GSE19750",
                   "GSE31684",
                   "GSE13041_GPL570",
                   "GSE7696",
                   "GSE2817",
                   "GSE37418",
                   "GSE146558",
                   "GSE16446",
                   "GSE19615",
                   "GSE20711",
                   "GSE21653",
                   "GSE42568",
                   "GSE48390",
                   "GSE58812",
                   "GSE6532_GPL570",
                   "GSE9195",
                   "GSE22762_GPL570",
                   "GSE17536",
                   "GSE17537",
                   "GSE29621",
                   "GSE31595",
                   "GSE38832",
                   "GSE39582",
                   "GSE14333",
                   "GSE10846",
                   "GSE23501",
                   "GSE15459",
                   "GSE34942",
                   "GSE62254",
                   "GSE10300",
                   "E_MTAB_1328",
                   "GSE3141",
                   "GSE30219",
                   "GSE31210",
                   "GSE8894",
                   "GSE157009",
                   "GSE157010",
                   "GSE19188",
                   "GSE37745",
                   "GSE50081",
                   "GSE24080",
                   "GSE57317",
                   "GSE18520",
                   "GSE19829_GPL570",
                   "GSE30161",
                   "GSE63885",
                   "GSE9891",
                   "GSE22138",
                   "GSE136337",
                   "GSE12417_GPL96",
                   "GSE5287",
                   "GSE13041_GPL96",
                   "GSE11121",
                   "GSE12093",
                   "GSE158309",
                   "GSE17705",
                   "GSE2603",
                   "GSE2990",
                   "GSE3494_GPL96",
                   "GSE45255",
                   "GSE4922_GPL96",
                   "GSE6532_GPL96",
                   "GSE4475",
                   "GSE22762_GPL96",
                   "GSE16131_GPL96",
                   "GSE27020",
                   "GSE31547",
                   "GSE68465",
                   "GSE4573",
                   "GSE14814",
                   "GSE9782_GPL96",
                   "GSE14764",
                   "GSE23554",
                   "GSE12417_GPL97",
                   "GSE4271_GPL96",
                   "GSE4271_GPL97",
                   "GSE22762_GPL97",
                   "GSE16131_GPL97",
                   "GSE9782_GPL97",
                   "GSE4412_GPL96",
                   "GSE2034",
                   "GSE5327",
                   "GSE12945",
                   "GSE41258",
                   "GSE26712",
                   "GSE4412_GPL97",
                   "GSE7390",
                   "GSE25055",
                   "GSE25065",
                   "GSE53031",
                   "E_MTAB_6134",
                   "GSE85916",
                   "E_MTAB_3218",
                   "GSE13041_GPL8300",
                   "GSE19829_GPL8300",
                   "GSE31245",
                   "ICGC_BRCA_FR"
                 )
                 if (accession %in% aff) {
                   platform <- "affymetrix"
                   name<-c("Stromal score", "Immune score", "Estimate score", "Tumor purity")
                 } else {
                   platform <- "other"
                   name<-c("Stromal score", "Immune score","Estimate score")
                 }
                 merged.df <- merge(common_genes, expres, by.x = "GeneSymbol", by.y = "row.names")
                 rownames(merged.df) <- merged.df$GeneSymbol
                 merged.df <- merged.df[, -1:-ncol(common_genes)]
                 print(sprintf("Merged dataset includes %d genes (%d mismatched).", nrow(merged.df), nrow(common_genes) - nrow(merged.df)))
                 gct.df<-outputgct(merged.df)
                 
                 ES<-EstimateScore(gct.df,platform=platform)
                 names(ES)<-name
                 index<-intersect(row.names(ES),names(expres))
                 ES<-ES[index,]
                 expres<-expres[,index]
                 ES$gene<-as.numeric(expres[gene,])
                 
                 scatter<-list()
                 for(i in 1:(ncol(ES)-1)){
                   scatter[[i]]<-ggscatter(data=ES, y = paste(names(ES)[i]),x = "gene",
                                           size = 1,add = "reg.line", conf.int = TRUE,
                                           add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                                           cor.coef = TRUE, cor.method = method,cor.coef.size=4.5,
                                           cor.coeff.args = list(label.y=max(ES[,i])+se(ES[,i]),label.sep = ","),
                                           xlab =paste("Normalized expression of",gene,sep=" ") , ylab = paste(names(ES)[i]))
                   
                 }
                 
                 list(estimate=ES,plot=scatter)
                 }
                 
               })
})



observe({
  if(
    is.null(input$estimatemethod) ||
    input$estimatemethod == ""||
    is.null(input$estimategene)||
    input$estimategene == ""  
  ){
    disable("estimatebt")
  }
  else{
    enable("estimatebt")
  }
})

observeEvent(input$estimatebt, {
  updateCollapse(session, "collapseestimate", open = "Correlation with ESTIMATE score", close = "Descriptions and parameters")
})



observeEvent(estimaterun(), {
  shinyjs::show("estimateplot_wrapper")
  shinyjs::show("estimatetable_wrapper")
})

observeEvent(input$estimatebt, {
  output$estimatetable <-  DT::renderDT({
    estimaterun()$estimate
  })
})

output$saveestimatetable <- downloadHandler(
  filename = function() {
    paste0("ESTIMATE-score-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(estimaterun()$estimate, file)
  }
)


observeEvent(input$estimatebt,{
  output$estimateploting  <- renderPlot({
    closeAlert(session, "estimatemess")
    ggarrange(plotlist=estimaterun()$plot,common.legend = T,align = "hv",labels="AUTO" )
  })
})

observeEvent(input$estimatebt, {
  output$estimateplot<- renderUI({
    plotOutput("estimateploting",
               width = paste0(isolate({input$estimatewidth}), "%"),
               height = isolate({input$estimateheight}))
  })})

output$saveestimateness <- downloadHandler(
  filename = function() {
    paste0("ESTIMATE-score-correlation-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$estimatepltwidth}),height = isolate({input$estimatepltheight}))
    print(ggarrange(plotlist=estimaterun()$plot,common.legend = T,align = "hv",labels="AUTO" ))
    dev.off()
  }
)

observeEvent(input$page_before_estimate, {
  newtab <- switch(input$tabs, "estimate" = "stemness","stemness" = "estimate")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_estimate, {
  newtab <- switch(input$tabs, "ICB" = "estimate","estimate" = "ICB")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









