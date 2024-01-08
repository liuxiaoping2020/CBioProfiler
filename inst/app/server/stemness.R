observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'stemgene',
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

stemrun<-eventReactive(input$stembt,{
  input$stembt

  data <- isolate({
    rawdata()
  })
  expres<-data$expres
  withProgress(message = "Correlating with stemness score",
               detail = "This may take a while...",
               value = 3,
               {
                 require(TCGAbiolinks)
                 stem(exp=expres,gene=isolate({input$stemgene}),method=isolate({input$stemmethod}))
               })
})

observe({
  if(
    is.null(input$stemgene) ||
    input$stemgene == "" ||
    is.null(input$stemmethod) ||
    input$stemmethod == ""
  ){
    disable("stembt")
  }
  else{
    enable("stembt")
  }
})

observeEvent(input$stembt, {
  disable("stembt")
})

observeEvent(input$stembt, {
  updateCollapse(session, "collapsestemness", open = "Correlation with stemness score", close = "Descriptions and parameters")
})


observeEvent(input$stembt, {
  output$stemtable <-  DT::renderDT({
    stemrun()$stemnesstable
  })
})

output$savestemtable <- downloadHandler(
  filename = function() {
    paste0("Stemness-score-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(stemrun()$stemnesstable, file)
  }
)


observeEvent(input$stembt,{
  output$stemnessploting  <- renderPlot({
    closeAlert(session, "stemnessmess")
    stemrun()$scatter
  })
})

observeEvent(input$stembt, {
  output$stemnessplot<- renderUI({
    plotOutput("stemnessploting",
               width = paste0(isolate({input$stemwidth}), "%"),
               height = isolate({input$stemheight}))
  })})

observeEvent(stemrun(), {
  shinyjs::show("stemnessplot_wrapper")
  shinyjs::show("stemtable_wrapper")
  enable("stembt")
})


output$savestemness <- downloadHandler(
  filename = function() {
    paste0("Stemness-score-correlation-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$stempltwidth}),height = isolate({input$stempltheight}))
    print(stemrun()$scatter)
    dev.off()
  }
)

observeEvent(input$page_before_stemness, {
  newtab <- switch(input$tabs, "stemness" = "immune","immune" = "stemness")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_stemness, {
  newtab <- switch(input$tabs, "estimate" = "stemness","stemness" = "estimate")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









