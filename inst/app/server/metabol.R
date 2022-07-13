observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'metapathgene',
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

metapathrun<-eventReactive(input$metapathbt,{
  input$metapathbt
  data <- isolate({
    rawdata()
  })
  expres<-data$expres
  
  method<-isolate({input$metapathmethod})
  gene<-isolate({input$metapathgene})
  geneset<-metabolism
  name<-names(geneset)
  
  withProgress(message = "Correlating with metabolism pathway",
               detail = "This may take a while...",
               value = 3,
               {
                 # gmgene<-c("GZMA","PRF1")
                 # if(length(setdiff(gmgene,row.names(expres)))!=0){
                 #   createAlert(
                 #     session,
                 #     "metapathmess",
                 #     "exampleAlert",
                 #     title = "Please note!",
                 #     style =  "danger",
                 #     content = paste(setdiff(gmgene),names(expres),"is not found in the gene expression profile"),
                 #     append = T,
                 #     dismiss=T
                 #   )
                 #   return(NULL)
                 # }else{
                 #   cormetapath(expres=expres,gene=gene, method=method)
                 # }
                 coroncopath(expres=expres,gene=gene,geneset=geneset,method=method,name=name)
               }
  )
})


observe({
  if(
    is.null(input$metapathmethod) ||
    input$metapathmethod == ""||
    is.null(input$metapathgene)||
    input$metapathgene == ""  
  ){
    disable("metapathbt")
  }
  else{
    enable("metapathbt")
  }
})

observeEvent(input$metapathbt, {
  updateCollapse(session, "collapsemetapath", open = "Correlation with metabolism pathway", close = "Descriptions and parameters")
})



observeEvent(metapathrun(), {
  shinyjs::show("metapathplot_wrapper")
  shinyjs::show("metapathtable_wrapper")
})

observeEvent(input$metapathbt, {
  output$metapathtable <-  DT::renderDT({
    DT::datatable(metapathrun()$table,options=list(scrollX=TRUE))
  })
})

output$savemetapathtable <- downloadHandler(
  filename = function() {
    paste0("Metabolism-pathway-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(metapathrun()$table, file)
  }
)


observeEvent(input$metapathbt,{
  output$metapathploting  <- renderPlot({
    closeAlert(session, "metapathmess")
    ggarrange(plotlist=metapathrun()$plot,common.legend = T,align = "hv",labels="AUTO" ,nrow = 5,ncol=3)
    
  })
})

observeEvent(input$metapathbt, {
  output$metapathplot<- renderUI({
    plotOutput("metapathploting",
               width = paste0(isolate({input$metapathwidth}), "%"),
               height = isolate({input$metapathheight}))
  })})

output$savemetapathness <- downloadHandler(
  filename = function() {
    paste0("Metabolism-pathway-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$metapathpltwidth}),height = isolate({input$metapathpltheight}))
    print(ggarrange(plotlist=metapathrun()$plot,common.legend = T,align = "hv",labels="AUTO",nrow = 5,ncol=3 ))
    dev.off()
  }
)

observeEvent(input$page_before_metapath, {
  newtab <- switch(input$tabs, "metapath" = "oncopath","oncopath" = "metapath")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_metapath, {
  newtab <- switch(input$tabs, "hallpath" = "metapath","metapath" = "hallpath")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









