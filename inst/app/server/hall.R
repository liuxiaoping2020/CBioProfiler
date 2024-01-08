observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'hallpathgene',
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

hallpathrun<-eventReactive(input$hallpathbt,{
  input$hallpathbt
  
  withProgress(message = "Correlating with hallmark signature",
               detail = "This may take a while...",
               value = 3,
               {
                 
  data <- isolate({
    rawdata()
  })
  expres<-data$expres
  
  method<-isolate({input$hallpathmethod})
  gene<-isolate({input$hallpathgene})
  
  require(msigdbr)
  require(GSVA)
  h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
  geneset<-subset(h_gene_sets,select=c(gene_symbol))
  geneset<-split(geneset,h_gene_sets$gs_name)
  fun<-function(x){
    x<-as.vector(x$gene_symbol)
    return(x)
  }
  geneset<-lapply(geneset,fun)
  
  name<-gsub("HALLMARK_","", names(geneset))
  
  # withProgress(message = "Correlating with hallmark signature",
  #              detail = "This may take a while...",
  #              value = 3,
  #              {

                 coroncopath(expres=expres,gene=gene,geneset=geneset,method=method,name=name)
               }
  )
})


observe({
  if(
    is.null(input$hallpathmethod) ||
    input$hallpathmethod == ""||
    is.null(input$hallpathgene)||
    input$hallpathgene == ""  
  ){
    disable("hallpathbt")
  }
  else{
    enable("hallpathbt")
  }
})


observeEvent(input$hallpathbt, {
  disable("hallpathbt")
})

observeEvent(input$hallpathbt, {
  updateCollapse(session, "collapsehallpath", open = "Correlation with hallmark signature", close = "Descriptions and parameters")
})



observeEvent(hallpathrun(), {
  shinyjs::show("hallpathplot_wrapper")
  shinyjs::show("hallpathtable_wrapper")
  enable('hallpathbt')
})

observeEvent(input$hallpathbt, {
  output$hallpathtable <-  DT::renderDT({
    DT::datatable(hallpathrun()$table,options=list(scrollX=TRUE))
  })
})

output$savehallpathtable <- downloadHandler(
  filename = function() {
    paste0("Hallmark-signature-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(hallpathrun()$table, file)
  }
)


observeEvent(input$hallpathbt,{
  output$hallpathploting  <- renderPlot({
    closeAlert(session, "hallpathmess")
    ggarrange(plotlist=hallpathrun()$plot,common.legend = T,align = "hv",nrow = 10,ncol=5)
    
  })
})

observeEvent(input$hallpathbt, {
  output$hallpathplot<- renderUI({
    plotOutput("hallpathploting",
               width = paste0(isolate({input$hallpathwidth}), "%"),
               height = isolate({input$hallpathheight}))
  })})

output$savehallpathness <- downloadHandler(
  filename = function() {
    paste0("Hallmark-signature-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$hallpathpltwidth}),height = isolate({input$hallpathpltheight}))
    print(ggarrange(plotlist=hallpathrun()$plot,common.legend = T,align = "hv",nrow = 10,ncol=5 ))
    dev.off()
  }
)

observeEvent(input$page_before_hallpath, {
  newtab <- switch(input$tabs, "hallpath" = "metapath","metapath" = "hallpath")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_hallpath, {
  newtab <- switch(input$tabs, "drug" = "hallpath","hallpath" = "drug")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









