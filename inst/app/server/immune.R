observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'immgene',
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

immrun<-eventReactive(input$immunebt,{
  input$immunebt
  data <- isolate({
    rawdata()
  })
  expres<-data$expres
  withProgress(message = "Calculating immune cell composition",
               detail = "This may take a while...",
               value = 3,
               {
  geneset<-isolate({input$geneset})
  require(ConsensusTME)
  if(geneset=="Bindea"){
    idx<-intersect(row.names(expres),as.vector(unlist(methodSignatures$Bindea)))
  }else if(geneset=="Danaher"){
    idx<-intersect(row.names(expres),as.vector(unlist(methodSignatures$Danaher)))
  }else if(geneset=="Davoli"){
    idx<-intersect(row.names(expres),as.vector(unlist(methodSignatures$Davoli)))
  }else if(geneset=="MCP.Counter"){
    idx<-intersect(row.names(expres),as.vector(unlist(methodSignatures$MCP.Counter)))
  }else if(geneset=="xCell"){
    idx<-intersect(row.names(expres),as.vector(unlist(methodSignatures$xCell)))
  }
  # print(idx)
  if(length(idx)==0){
    createAlert(
      session,
      "immunemess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = paste(
        "No genes in the reference immune geneset could be matched to the identifiers in the expression data you selected, please Check!!"
      ),
      append = T,
      dismiss=T
    )
    return(NULL)
  }else
     {
  score<-ims(exp=as.matrix(expres),geneset=geneset)
  res <-
    plotcor(score = score,
            exp = as.data.frame(expres),
            gene=isolate({input$immgene}),
            select=isolate({input$corselect}),
            normmethod = isolate({input$immnorm}),
            fdrcutoff=isolate({input$immcutoff}),
            type=isolate({input$immmethod}),
            plottype=isolate({input$immpltype}))

  if(nrow(res$cordata)==0){
    createAlert(
      session,
      "immunemess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = paste(
        "No immune cell infiltration score meet the significance cutoff (",
        isolate({input$immcutoff}),
        ") at the correction method",
        paste("'", isolate({input$immnorm}), "'"),
        "you specified",
        sep = " "
      ),
      append = T,
      dismiss=T
    )
    return(NULL)
  }else{
    res
  }
     }
       })
})

observe({
  if(
    is.null(input$geneset) ||
    input$geneset == "" ||
    is.null(input$immgene) ||
    input$immgene == "" ||
    is.null(input$immmethod) ||
    input$immmethod == "" ||
    is.null(input$immnorm) ||
    input$immnorm == "" ||
    is.null(input$immpltype)||
    input$immpltype==""
  ){
    disable("immunebt")
  }
  else{
    enable("immunebt")
  }
})

observeEvent(input$immunebt, {
  disable("immunebt")

})

observeEvent(input$immunebt, {
  updateCollapse(session, "collapseimmune", open = "Correlation with immune cell infiltration", close = "Descriptions and parameters")
})

observeEvent(immrun(), {
  shinyjs::show("immuneplot_wrapper")
  shinyjs::show("immunecell_wrapper")
  shinyjs::show("immcor_wrapper")
  enable("immunebt")
})

observeEvent(is.null(immrun()), {
  enable("immunebt")
})


observeEvent(input$immunebt, {
  output$immcelltable <-  DT::renderDT({
    DT::datatable(immrun()$score,options=list(scrollX=TRUE))
  })
})

output$saveimmcelltable <- downloadHandler(
  filename = function() {
    paste0("Immune-cell-inflitration-composition-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(immrun()$score, file)
  }
)


observeEvent(input$immunebt, {
  output$immcortable <-  DT::renderDT({
    DT::datatable( immrun()$cordata,options=list(scrollX=TRUE))
   
  })
})

output$saveimmcortable <- downloadHandler(
  filename = function() {
    paste0("Immune-cell-inflitration-composition-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(immrun()$cordata, file)
  }
)

observeEvent(input$immunebt,{
  output$immuneploting  <- renderPlot({
     closeAlert(session, "immunemess")
    immrun()$corplot
  })
})

observeEvent(input$immunebt, {
  output$immuneplot<- renderUI({
    plotOutput("immuneploting",
               width = paste0(isolate({input$immwidth}), "%"),
               height = isolate({input$immheight}))
  })})



output$saveimune <- downloadHandler(

  filename = function() {
    paste0("Immune-correlation-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$immpltwidth}),height = isolate({input$immpltheight}))
    print(immrun()$corplot)
    dev.off()

  }
)

observeEvent(input$page_before_immune, {
  newtab <- switch(input$tabs, "genediff" = "immune","immune" = "genediff")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_immune, {
  newtab <- switch(input$tabs, "immune" = "stemness","stemness" = "immune")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})








