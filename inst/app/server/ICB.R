observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'ICBgene',
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

ICBrun<-eventReactive(input$ICBbt,{
  input$ICBbt
  data <- isolate({
    rawdata()
  })
  expres<-data$expres
  expres<-as.data.frame(t(expres))
  method<-isolate({input$ICBmethod})
  gene<-isolate({input$ICBgene})
  ICB<-c("PDCD1","CD274","PDCD1LG2","CTLA4","PVR","LAG3","TIGIT","HAVCR2","VTCN1","CD86","CD28", "CD80","IDO1",
         "CD27","CD40","IL2RB","TNFRSF9","TNFRSF4","TNFRSF18","ICOS","CD276","BTLA","KIR3DL1","CYBB","VSIR","SIGLEC7")
  withProgress(message = "Correlating with immune checkpoint",
               detail = "This may take a while...",
               value = 3,
               {
                 if(length(intersect(ICB,names(expres)))==0){
                   createAlert(
                     session,
                     "ICBmess",
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
                 }else if(length(setdiff(ICB,names(expres)))!=0){
                   createAlert(
                     session,
                     "ICBmess",
                     "exampleAlert",
                     title = "Please note!",
                     style =  "info",
                     content = paste(paste(setdiff(ICB,names(expres)),collapse = ', '), "were not found in the gene expression data") ,
                     append = T,
                     dismiss=T
                   )
                   ICB<-setdiff(ICB,setdiff(ICB,names(expres)))
                   corICB(expres=expres,ICB=ICB,gene=gene,method=method)
                 }else{
                   corICB(expres=expres,ICB=ICB,gene=gene,method=method)
                 }
               })
})


observe({
  if(
    is.null(input$ICBmethod) ||
    input$ICBmethod == ""||
    is.null(input$ICBgene)||
    input$ICBgene == ""  
  ){
    disable("ICBbt")
  }
  else{
    enable("ICBbt")
  }
})

observeEvent(input$ICBbt, {
  updateCollapse(session, "collapseICB", open = "Correlation with immune checkpoint", close = "Descriptions and parameters")
})



observeEvent(ICBrun(), {
  shinyjs::show("ICBplot_wrapper")
  shinyjs::show("ICBtable_wrapper")
})

observeEvent(input$ICBbt, {
  output$ICBtable <-  DT::renderDT({
    DT::datatable(ICBrun()$expres,options=list(scrollX=TRUE))
  })
})

output$saveICBtable <- downloadHandler(
  filename = function() {
    paste0("Immune-checkpoint-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(ICBrun()$expres, file)
  }
)


observeEvent(input$ICBbt,{
  output$ICBploting  <- renderPlot({
    closeAlert(session, "ICBmess")
    ggarrange(plotlist=ICBrun()$plot,common.legend = T,align = "hv",labels="AUTO" )
  })
})

observeEvent(input$ICBbt, {
  output$ICBplot<- renderUI({
    plotOutput("ICBploting",
               width = paste0(isolate({input$ICBwidth}), "%"),
               height = isolate({input$ICBheight}))
  })})

output$saveICBness <- downloadHandler(
  filename = function() {
    paste0("Immune-checkpoint-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$ICBpltwidth}),height = isolate({input$ICBpltheight}))
    print(ggarrange(plotlist=ICBrun()$plot,common.legend = T,align = "hv",labels="AUTO" ))
    dev.off()
  }
)

observeEvent(input$page_before_ICB, {
  newtab <- switch(input$tabs, "ICB" = "estimate","estimate" = "ICB")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_ICB, {
  newtab <- switch(input$tabs, "IFN" = "ICB","ICB" = "IFN")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









