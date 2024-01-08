observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'oncopathgene',
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

oncopathrun<-eventReactive(input$oncopathbt,{
  input$oncopathbt

  
  withProgress(message = "Correlating with cancer pathway",
               detail = "This may take a while...",
               value = 3,
               {
                 # gmgene<-c("GZMA","PRF1")
                 # if(length(setdiff(gmgene,row.names(expres)))!=0){
                 #   createAlert(
                 #     session,
                 #     "oncopathmess",
                 #     "exampleAlert",
                 #     title = "Please note!",
                 #     style =  "danger",
                 #     content = paste(setdiff(gmgene),names(expres),"is not found in the gene expression profile"),
                 #     append = T,
                 #     dismiss=T
                 #   )
                 #   return(NULL)
                 # }else{
                 #   coroncopath(expres=expres,gene=gene, method=method)
                 # }
                 data <- isolate({
                   rawdata()
                 })
                 expres<-data$expres
                 require(GSVA)
                 method<-isolate({input$oncopathmethod})
                 gene<-isolate({input$oncopathgene})
                 geneset<-oncopath
                 name<-names(geneset)
                 coroncopath(expres=expres,gene=gene,geneset=geneset,method=method,name=name)
               }
               )
})


observe({
  if(
    is.null(input$oncopathmethod) ||
    input$oncopathmethod == ""||
    is.null(input$oncopathgene)||
    input$oncopathgene == ""  
  ){
    disable("oncopathbt")
  }
  else{
    enable("oncopathbt")
  }
})
observeEvent(input$oncopathbt, {
  disable("oncopathbt")
})

observeEvent(input$oncopathbt, {
  updateCollapse(session, "collapseoncopath", open = "Correlation with cancer pathway", close = "Descriptions and parameters")
})



observeEvent(oncopathrun(), {
  shinyjs::show("oncopathplot_wrapper")
  shinyjs::show("oncopathtable_wrapper")
  enable("oncopathbt")
})

observeEvent(input$oncopathbt, {
  output$oncopathtable <-  DT::renderDT({
    DT::datatable(oncopathrun()$table,options=list(scrollX=TRUE))
  })
})

output$saveoncopathtable <- downloadHandler(
  filename = function() {
    paste0("Cancer-pathway-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(oncopathrun()$table, file)
  }
)


observeEvent(input$oncopathbt,{
  output$oncopathploting  <- renderPlot({
    closeAlert(session, "oncopathmess")
    ggarrange(plotlist=oncopathrun()$plot,common.legend = T,align = "hv",labels="AUTO" ,nrow = 5,ncol=3)
    
  })
})

observeEvent(input$oncopathbt, {
  output$oncopathplot<- renderUI({
    plotOutput("oncopathploting",
               width = paste0(isolate({input$oncopathwidth}), "%"),
               height = isolate({input$oncopathheight}))
  })})

output$saveoncopathness <- downloadHandler(
  filename = function() {
    paste0("Cancer-pathway-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$oncopathpltwidth}),height = isolate({input$oncopathpltheight}))
    print(ggarrange(plotlist=oncopathrun()$plot,common.legend = T,align = "hv",labels="AUTO" ,nrow = 5,ncol=3 ))
    dev.off()
  }
)

observeEvent(input$page_before_oncopath, {
  newtab <- switch(input$tabs, "oncopath" = "CYA","CYA" = "oncopath")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_oncopath, {
  newtab <- switch(input$tabs, "metapath" = "oncopath","oncopath" = "metapath")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









