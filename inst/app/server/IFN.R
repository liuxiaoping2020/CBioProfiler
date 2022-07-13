observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'IFNgene',
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

IFNrun<-eventReactive(input$IFNbt,{
  
  input$IFNbt
  data <- isolate({
    rawdata()
  })
  expres<-data$expres
  require(GSVA)
  method<-isolate({input$IFNmethod})
  gene<-isolate({input$IFNgene})
  IFN<-c('IFNG',"IFNGR1","IFNGR2","IRF1","JAK1","JAK2","STAT1")
  
  withProgress(message = "Correlating with interferon-gamma score",
               detail = "This may take a while...",
               value = 3,
               {
                 if(length(intersect(IFN,row.names(expres)))==0){
                   createAlert(
                     session,
                     "IFNmess",
                     "exampleAlert",
                     title = "Please note!",
                     style =  "danger",
                     content = paste(
                       "No immune checkpoints were not found in the gene expression profile you specified. Please check!!!"
                     ),
                     append = T,
                     dismiss=T
                   )
                   return(NULL)
                 }else{
                   corIFN(expres=expres,gene=gene, method=method,IFN=IFN)
                 }
               })
})


observe({
  if(
    is.null(input$IFNmethod) ||
    input$IFNmethod == ""||
    is.null(input$IFNgene)||
    input$IFNgene == ""  
  ){
    disable("IFNbt")
  }
  else{
    enable("IFNbt")
  }
})

observeEvent(input$IFNbt, {
  updateCollapse(session, "collapseIFN", open = "Correlation with interferon-gamma score", close = "Descriptions and parameters")
})



observeEvent(IFNrun(), {
  shinyjs::show("IFNplot_wrapper")
  shinyjs::show("IFNtable_wrapper")
})

observeEvent(input$IFNbt, {
  output$IFNtable <-  DT::renderDT({
    DT::datatable(IFNrun()$IFNscore,options=list(scrollX=TRUE))
  })
})

output$saveIFNtable <- downloadHandler(
  filename = function() {
    paste0("Interferon-gamma-score-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(IFNrun()$IFNscore, file)
  }
)


observeEvent(input$IFNbt,{
  output$IFNploting  <- renderPlot({
    closeAlert(session, "IFNmess")
     IFNrun()$plot
  })
})

observeEvent(input$IFNbt, {
  output$IFNplot<- renderUI({
    plotOutput("IFNploting",
               width = paste0(isolate({input$IFNwidth}), "%"),
               height = isolate({input$IFNheight}))
  })})

output$saveIFNness <- downloadHandler(
  filename = function() {
    paste0("Interferon-gamma-score-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$IFNpltwidth}),height = isolate({input$IFNpltheight}))
    print(IFNrun()$plot)
    dev.off()
  }
)

observeEvent(input$page_before_IFN, {
  newtab <- switch(input$tabs, "IFN" = "ICB","ICB" = "IFN")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_IFN, {
  newtab <- switch(input$tabs, "CYA" = "IFN","IFN" = "CYA")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









