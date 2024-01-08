observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'CYAgene',
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

CYArun<-eventReactive(input$CYAbt,{
  input$CYAbt
  data <- isolate({
    rawdata()
  })
  expres<-data$expres
  
  method<-isolate({input$CYAmethod})
  gene<-isolate({input$CYAgene})

  withProgress(message = "Correlating with cytotoxic activity",
               detail = "This may take a while...",
               value = 3,
               {
                 gmgene<-c("GZMA","PRF1")
                 if(length(setdiff(gmgene,row.names(expres)))!=0){
                   createAlert(
                     session,
                     "CYAmess",
                     "exampleAlert",
                     title = "Please note!",
                     style =  "danger",
                     content = paste(setdiff(gmgene),names(expres),"is not found in the gene expression profile"),
                     append = T,
                     dismiss=T
                   )
                   return(NULL)
                 }else{
                   corCYA(expres=expres,gene=gene, method=method)
                 }
               })
})


observe({
  if(
    is.null(input$CYAmethod) ||
    input$CYAmethod == ""||
    is.null(input$CYAgene)||
    input$CYAgene == ""  
  ){
    disable("CYAbt")
  }
  else{
    enable("CYAbt")
  }
})

observeEvent(input$CYAbt, {
  updateCollapse(session, "collapseCYA", open = "Correlation with cytotoxic activity", close = "Descriptions and parameters")
})



observeEvent(CYArun(), {
  shinyjs::show("CYAplot_wrapper")
  shinyjs::show("CYAtable_wrapper")
})

observeEvent(input$CYAbt, {
  output$CYAtable <-  DT::renderDT({
    DT::datatable(CYArun()$table,options=list(scrollX=TRUE))
  })
})

output$saveCYAtable <- downloadHandler(
  filename = function() {
    paste0("Cytotoxic-activity-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(CYArun()$table, file)
  }
)


observeEvent(input$CYAbt,{
  output$CYAploting  <- renderPlot({
    closeAlert(session, "CYAmess")
    CYArun()$plot
  })
})

observeEvent(input$CYAbt, {
  output$CYAplot<- renderUI({
    plotOutput("CYAploting",
               width = paste0(isolate({input$CYAwidth}), "%"),
               height = isolate({input$CYAheight}))
  })})

output$saveCYAness <- downloadHandler(
  filename = function() {
    paste0("Cytotoxic-activity-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$CYApltwidth}),height = isolate({input$CYApltheight}))
    print(CYArun()$plot)
    dev.off()
  }
)

observeEvent(input$page_before_CYA, {
  newtab <- switch(input$tabs, "CYA" = "IFN","IFN" = "CYA")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_CYA, {
  newtab <- switch(input$tabs, "oncopath" = "CYA","CYA" = "oncopath")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









