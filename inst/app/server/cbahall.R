
value <- reactiveValues(data = NULL)
observeEvent(input$cbavalidationbt, {
  value$data <- isolate({validst()})
})

cbahallpathrun<-eventReactive(input$cbahallpathbt,{
  input$cbahallpathbt
  
  withProgress(message = "Calculating the hallmark signature",
               detail = "This may take a while...",
               value = 3,
               {
                 createAlert(
                   session,
                   "cbahallpathmess",
                   "exampleAlert",
                   title = "Please note!",
                   style =  "info",
                   content = "CBioExplorer is drawing the plot, which will take some time. Please be patient!" ,
                   append = T,
                   dismiss=T
                 )
                 res<-isolate({validst()})
                 cohort<-isolate({input$cbahallpathcor})
                 method<-isolate({input$cbahallpathmeth})
                 paired<-isolate({input$cbahallpathpair})
                 plottype<-isolate({input$cbahallpathplot})
                 
                 require(GSVA)
                 require(msigdbr)
                 
                 h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
                 geneset<-subset(h_gene_sets,select=c(gene_symbol))
                 geneset<-split(geneset,h_gene_sets$gs_name)
                 fun<-function(x){
                   x<-as.vector(x$gene_symbol)
                   return(x)
                 }
                 geneset<-lapply(geneset,fun)
                 
                 name<-gsub("HALLMARK_","", names(geneset))
                 compgeneset(res=res,cohort=cohort,plottype=plottype,method=method,geneset=geneset,name)
               })
})

observe({
  if(
    is.null(input$cbahallpathmeth) ||
    input$cbahallpathmeth == ""||
    is.null(input$cbahallpathcor)||
    input$cbahallpathcor == "" ||
    is.null(input$cbahallpathpair)||
    input$cbahallpathpair == "" ||
    is.null(input$cbahallpathplot)||
    input$cbahallpathplot == "" 
  ){
    disable("cbahallpathbt")
  }
  else{
    enable("cbahallpathbt")
  }
})


observeEvent(input$cbahallpathbt, {
  disable("cbahallpathbt")
})


observeEvent(input$cbahallpathbt, {
  updateCollapse(session, "cbacollapsehallpath", open = "Comparision of hallmark signature", close = "Descriptions and parameters")
})

observe({
  if(is.null(value$data)){
    disable("cbahallpathbt")
  }
  else{
    enable("cbahallpathbt")
  }
})

observeEvent(cbahallpathrun(), {
  shinyjs::show("cbahallpathplot_wrapper")
  shinyjs::show("cbahallpathtable_wrapper")
  enable("cbahallpathbt")
})

observeEvent(input$cbahallpathbt, {
  output$cbahallpathtable <- DT::renderDT({
    DT::datatable(cbahallpathrun()$table,options=list(scrollX=TRUE))
  })
})

output$cbasavehallpathtable <- downloadHandler(
  filename = function() {
    paste0("Comparision-of-hallmark-signature-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(cbahallpathrun()$table, file)
  }
)
observeEvent(input$cbahallpathbt,{
  output$cbahallpathploting  <- renderPlot({
    closeAlert(session, "cbahallpathmess")
    ggarrange(plotlist=cbahallpathrun()$plot,common.legend = T,align = "hv",nrow=10,ncol=5)
  })
})

observeEvent(input$cbahallpathbt, {
  output$cbahallpathplot<- renderUI({
    plotOutput("cbahallpathploting",
               width = paste0(isolate({input$cbahallpathwidth}), "%"),
               height = isolate({input$cbahallpathheight}))
  })})

output$cbasavehallpathness <- downloadHandler(
  filename = function() {
    paste0("Hallmark-signature-comparision-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cbahallpathpltwidth}),height = isolate({input$cbahallpathpltheight}))
    print(ggarrange(plotlist=cbahallpathrun()$plot,common.legend = T,align = "hv" ,nrow=10,ncol=5))
    dev.off()
  }
)

observeEvent(input$page_before_cbahallpath, {
  newtab <- switch(input$tabs, "cbahall" = "cbametapath","cbametapath" = "cbahall")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbahallpath, {
  newtab <- switch(input$tabs, "cbadrug" = "cbahall","cbahall" = "cbadrug")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









