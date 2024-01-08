observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'druggene',
                       choices = {
                         if (!is.null(data)) {
                           if (class(data) ==  class(data.frame()))
                             choices <- as.character(row.names(data))
                           if (class(data) !=  class(data.frame()))
                             choices <- as.character(row.names(as.data.frame(data)))
                           choices
                         }
                       },
                       server = TRUE,
                       selected = "TP53")
})

drugrun<-eventReactive(input$drugbt,{
  input$drugbt
  data <- isolate({rawdata()})
  expres<-data$expres
  method<-isolate({input$drugmethod})
  gene<-isolate({input$druggene})
  drug<-isolate({input$drugs})
  tissueType = isolate({input$drugtissueType})
  batchCorrect<-  isolate({input$batchCorrect})

  withProgress(message = "Correlating with cytotoxic activity",
               detail = "This may take a while...",
               value = 3,
               {
                 # gmgene<-c("GZMA","PRF1")
                 # if(length(setdiff(gmgene,row.names(expres)))!=0){
                 #   createAlert(
                 #     session,
                 #     "drugmess",
                 #     "exampleAlert",
                 #     title = "Please note!",
                 #     style =  "danger",
                 #     content = paste(setdiff(gmgene),names(expres),"is not found in the gene expression profile"),
                 #     append = T,
                 #     dismiss=T
                 #   )
                 #   return(NULL)
                 # }else
                 {
                   library(pRRophetic)
                   library(car)
                   library(ridge)
                   response <- pRRopheticPredict1(testMatrix=as.matrix(expres), 
                                                 drug=drug,
                                                 tissueType = tissueType, 
                                                 batchCorrect = batchCorrect,
                                                 selection=1,
                                                 dataset = "cgp2016")
                   response<-as.data.frame(response)
                   if(identical(names(expres),row.names(response))){
                     response$gene<-as.numeric(expres[gene,])
                   }
                   names(response)<-c(paste("Estimated IC50 of",drug),"gene")
                   scatter<-ggscatter(data=response, y = paste(names(response)[1]),x = "gene",
                                      size = 1,
                                      add = "reg.line", conf.int = TRUE,
                                      add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                                      cor.coef = TRUE, cor.method = method,cor.coef.size=4.5,
                                      cor.coeff.args = list(label.y=max(response[,1])+se(response[,1]),label.sep = ","),
                                      xlab =paste("Normalized expression of",gene,sep=" ") , ylab = paste(names(response)[1]))
                   
                   list(plot=scatter,table=response)
                 }
               })
})


observe({
  if(
    is.null(input$drugmethod) ||
    input$drugmethod == ""||
    is.null(input$druggene)||
    input$druggene == ""  
  ){
    disable("drugbt")
  }
  else{
    enable("drugbt")
  }
})

observeEvent(input$drugbt, {
  disable("drugbt")
})

observeEvent(input$drugbt, {
  updateCollapse(session, "collapsedrug", open = "Correlation with drug response", close = "Descriptions and parameters")
})



observeEvent(drugrun(), {
  shinyjs::show("drugplot_wrapper")
  shinyjs::show("drugtable_wrapper")
  enable("drugbt")
})

observeEvent(input$drugbt, {
  output$drugtable <-  DT::renderDT({
    DT::datatable(drugrun()$table,options=list(scrollX=TRUE))
  })
})

output$savedrugtable <- downloadHandler(
  filename = function() {
    paste0("drug-response-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(drugrun()$table, file)
  }
)


observeEvent(input$drugbt,{
  output$drugploting  <- renderPlot({
    closeAlert(session, "drugmess")
    drugrun()$plot
  })
})

observeEvent(input$drugbt, {
  output$drugplot<- renderUI({
    plotOutput("drugploting",
               width = paste0(isolate({input$drugwidth}), "%"),
               height = isolate({input$drugheight}))
  })})

output$savedrugness <- downloadHandler(
  filename = function() {
    paste0("drug-response-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$drugpltwidth}),height = isolate({input$drugpltheight}))
    print(drugrun()$plot)
    dev.off()
  }
)

observeEvent(input$page_before_drug, {
  newtab <- switch(input$tabs, "drug" = "hallpath","hallpath" = "drug")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_drug, {
  newtab <- switch(input$tabs, "bioan" = "drug","drug" = "bioan")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









