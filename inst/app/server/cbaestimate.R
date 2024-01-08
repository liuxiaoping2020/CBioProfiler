
value <- reactiveValues(data = NULL)
observeEvent(input$cbavalidationbt, {
  value$data <- isolate({validst()})
})

cbaestimaterun<-eventReactive(input$cbaestimatebt,{
  input$cbaestimatebt
  
  res<-isolate({validst()})
  cohort<-isolate({input$cbaestimatecor})
  accession<-isolate({input$accession})
  accession2<-isolate({input$accession2})
  if(cohort=="Training set"){
    data<-res$traindata
  }else{
    data<-res$validata
  }
  method<-isolate({input$cbaestimatemeth})
  paired<-isolate({input$cbaestimatepair})
  plottype<-isolate({input$cbaestimateplot})
  
  # expres<-data$expres
  # clinical<-data$clinical
  # data<-clinical
  withProgress(message = "Correlating with ESTIMATE score",
               detail = "This may take a while...",
               value = 3,
               {
                 if(cohort=="Training set"){
                   expres<-as.data.frame(res$traindata$expres)
                   clinical<-res$traindata$clinical
                 }else{
                   expres<-as.data.frame(res$validata$expres)
                   clinical<-res$validata$clinical
                 }
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
                 
                 if(cohort=="Training set") {
                   if (accession %in% aff) {
                     platform <- "affymetrix"
                   } else {
                     platform <- "other"
                   }
                 } else {
                   if (accession2 %in% aff){
                     platform <- "affymetrix"
                   } else {
                     platform <- "other"
                   }
                 }
                 
                 merged.df <- merge(common_genes, expres, by.x = "GeneSymbol", by.y = "row.names")
                 rownames(merged.df) <- merged.df$GeneSymbol
                 merged.df <- merged.df[, -1:-ncol(common_genes)]
                 print(sprintf("Merged dataset includes %d genes (%d mismatched).", nrow(merged.df), nrow(common_genes) - nrow(merged.df)))
                 gct.df<-outputgct(merged.df)
                 ES<-EstimateScore(gct.df,platform=platform)
                 index<-intersect(row.names(clinical),row.names(ES))
                 clinical<-clinical[index,]
                 ES<-ES[index,]
                 ES$Subtype<-clinical$Subtype
                 
                 if(platform!="affymetrix"){
                   name<-c("Stromal score", "Immune score","Estimate score")
                 }else{
                   name<-c("Stromal score", "Immune score", "Estimate score", "Tumor purity")
                 }
                 compare<-as.data.frame(combn(unique(as.vector(ES$Subtype)),2))
                 comp<-list()
                 for(i in 1:length(compare)){
                   comp[[i]]<-compare[[i]]
                 }
                pp<-barp(ES,paired=paired,plottype=plottype,method=method,comp=comp,name=name)
                list(estimate=ES,plot=pp)
               })
})

# observe({
#   if(is.null(value$data)){
#     # print(value1$data) cbaestimatebt
#     disable("cbaestimatebt")
#   }
#   else{
#     enable("cbaestimatebt")
#   }
# })

observe({
  if(
    is.null(input$cbaestimatemeth) ||
    input$cbaestimatemeth == ""||
    is.null(input$cbaestimatecor)||
    input$cbaestimatecor == "" ||
    is.null(input$cbaestimatepair)||
    input$cbaestimatepair == "" ||
    is.null(input$cbaestimateplot)||
    input$cbaestimateplot == "" 
  ){
    disable("cbaestimatebt")
  }
  else{
    enable("cbaestimatebt")
  }
})


observeEvent(input$cbaestimatebt, {
  disable("cbaestimatebt")
})



observeEvent(input$cbaestimatebt, {
  updateCollapse(session, "cbacollapseestimate", open = "Correlation with ESTIMATE score", close = "Descriptions and parameters")
})

observe({
  if(is.null(value$data)){
    # print(value1$data) cbaestimatebt
    disable("cbaestimatebt")
  }
  else{
    enable("cbaestimatebt")
  }
})

observeEvent(cbaestimaterun(), {
  shinyjs::show("cbaestimateplot_wrapper")
  shinyjs::show("cbaestimatetable_wrapper")
  enable("cbaestimatebt")
})

observeEvent(input$cbaestimatebt, {
  output$cbaestimatetable <-  DT::renderDT({
    cbaestimaterun()$estimate
  })
})

output$cbasaveestimatetable <- downloadHandler(
  filename = function() {
    paste0("ESTIMATE-score-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(cbaestimaterun()$estimate, file)
  }
)


observeEvent(input$cbaestimatebt,{
  output$cbaestimateploting  <- renderPlot({
    closeAlert(session, "cbaestimatemess")
    ggarrange(plotlist=cbaestimaterun()$plot,common.legend = T,align = "hv",labels="AUTO" )
  })
})

observeEvent(input$cbaestimatebt, {
  output$cbaestimateplot<- renderUI({
    plotOutput("cbaestimateploting",
               width = paste0(isolate({input$cbaestimatewidth}), "%"),
               height = isolate({input$cbaestimateheight}))
  })})

output$cbasaveestimateness <- downloadHandler(
  filename = function() {
    paste0("ESTIMATE-score-comparision-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cbaestimatepltwidth}),height = isolate({input$cbaestimatepltheight}))
    print(ggarrange(plotlist=cbaestimaterun()$plot,common.legend = T,align = "hv",labels="AUTO" ))
    dev.off()
  }
)

observeEvent(input$page_before_cbaestimate, {
  newtab <- switch(input$tabs, "cbaestimate" = "cbastemness","cbastemness" = "cbaestimate")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbaestimate, {
  newtab <- switch(input$tabs, "cbaICB" = "cbaestimate","cbaestimate" = "cbaICB")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









