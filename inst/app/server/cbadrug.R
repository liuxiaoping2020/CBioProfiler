
value <- reactiveValues(data = NULL)
observeEvent(input$cbavalidationbt, {
  value$data <- isolate({validst()})
})

cbadrugrun<-eventReactive(input$cbadrugbt,{
  input$cbadrugbt
  
  res<-isolate({validst()})
  cohort<-isolate({input$cbadrugcor})
  method<-isolate({input$cbadrugmeth})
  paired<-isolate({input$cbadrugpair})
  plottype<-isolate({input$cbadrugplot})
  drug<-isolate({input$cbadrugs})
  tissueType = isolate({input$cbadrugtissueType})
  batchCorrect<- isolate({input$cbabatchCorrect})

  
  
  
  withProgress(message = "Calculating the drug response",
               detail = "This may take a while...",
               value = 3,
               {
                 require(pRRophetic)
                 compdrug(cohort=cohort,res=res,drug=drug,tissueType=tissueType,batchCorrect=batchCorrect,paired=paired,plottype=plottype,method=method)
                 
               })
})

observe({
  if(
    is.null(input$cbadrugmeth) ||
    input$cbadrugmeth == ""||
    is.null(input$cbadrugcor)||
    input$cbadrugcor == "" ||
    is.null(input$cbadrugpair)||
    input$cbadrugpair == "" ||
    is.null(input$cbadrugplot)||
    input$cbadrugplot == "" 
  ){
    disable("cbadrugbt")
  }
  else{
    enable("cbadrugbt")
  }
})

observeEvent(input$cbadrugbt, {
  disable("cbadrugbt")
})

observeEvent(input$cbadrugbt, {
  updateCollapse(session, "cbacollapsedrug", open = "Comparision of drug response", close = "Descriptions and parameters")
})

observe({
  if(is.null(value$data)){
    disable("cbadrugbt")
  }
  else{
    enable("cbadrugbt")
  }
})

observeEvent(cbadrugrun(), {
  shinyjs::show("cbadrugplot_wrapper")
  shinyjs::show("cbadrugtable_wrapper")
})

observeEvent(input$cbadrugbt, {
  output$cbadrugtable <- DT::renderDT({
    DT::datatable(cbadrugrun()$table,options=list(scrollX=TRUE))
  })
})

output$cbasavedrugtable <- downloadHandler(
  filename = function() {
    paste0("Drug-response-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(cbadrugrun()$table, file)
  }
)
observeEvent(input$cbadrugbt,{
  output$cbadrugploting  <- renderPlot({
    closeAlert(session, "cbadrugmess")
    cbadrugrun()$plot
  })
})

observeEvent(input$cbadrugbt, {
  output$cbadrugplot<- renderUI({
    plotOutput("cbadrugploting",
               width = paste0(isolate({input$cbadrugwidth}), "%"),
               height = isolate({input$cbadrugheight}))
  })})

output$cbasavedrugness <- downloadHandler(
  filename = function() {
    paste0("Drug-response-comparision-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cbadrugpltwidth}),height = isolate({input$cbadrugpltheight}))
    print(cbadrugrun()$plot)
    dev.off()
  }
)

observeEvent(input$page_before_cbadrug, {
  newtab <- switch(input$tabs, "cbadrug" = "cbahall","cbahall" = "cbadrug")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})
# 
# observeEvent(input$page_after_cbadrugness, {
#   newtab <- switch(input$tabs, "bioan" = "cbadrug","cbadrug" = "bioan")
#   updateTabItems(session, "tabs", newtab)
#   shinyjs::runjs("window.scrollTo(0, 50)")
# })









