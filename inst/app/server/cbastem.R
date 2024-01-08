
value <- reactiveValues(data = NULL)
observeEvent(input$cbavalidationbt, {
  value$data <- isolate({validst()})
})

cbastemrun<-eventReactive(input$cbastembt,{
  input$cbastembt
  
  withProgress(message = "Correlating with stemness score",
               detail = "This may take a while...",
               value = 3,
               {
                 res<-isolate({validst()})
                 cohort<-isolate({input$cbastemcor})
                 
                 if(cohort=="Training set"){
                   data<-res$traindata
                 }else{
                   data<-res$validata
                 }
                 method<-isolate({input$cbastemmeth})
                 paired<-isolate({input$cbastempair})
                 plottype<-isolate({input$cbastemplot})
                 
                 expres<-data$expres
                 clinical<-data$clinical
                 data<-clinical
                 require(TCGAbiolinks)
                 stemgroup(exp=expres,clinical=clinical,method=method,plottype=plottype,paired=paired) 
               })
})

observe({
  if(is.null(value$data)){
    disable("cbastembt")
  }
  else{
    enable("cbastembt")
  }
})

observeEvent(input$cbastembt, {
  disable("cbastembt")
})

observe({
  if(
    is.null(input$cbastemmethod) ||
    input$cbastemmethod == ""||
    is.null(input$cbastemcor)||
    input$cbastemcor == "" ||
    is.null(input$cbastempair)||
    input$cbastempair == "" ||
    is.null(input$cbastemplot)||
    input$cbastemplot == "" 
    
  ){
    disable("cbastembt")
  }
  else{
    enable("cbastembt")
  }
})


observeEvent(cbastemrun(), {
  shinyjs::show("cbastemnessplot_wrapper")
  # shinyjs::show("cbastemtable_wrapper")
  enable("cbastembt")
})

observeEvent(input$cbastembt, {
  updateCollapse(session, "cbacollapsestemness", open = "Correlation with stemness score", close = "Descriptions and parameters")
})


# observeEvent(input$cbastembt, {
#   output$stemtable <-  DT::renderDT({
#     stemrun()$stemnesstable
#   })
# })

# output$savestemtable <- downloadHandler(
#   filename = function() {
#     paste0("Stemness-score-table-",Sys.Date(), ".csv")
#   },
#   content = function(file) {
#     write.csv(stemrun()$stemnesstable, file)
#   }
# )


observeEvent(input$cbastembt,{
  output$cbastemnessploting  <- renderPlot({
    closeAlert(session, "cbastemnessmess")
    cbastemrun()
  })
})

observeEvent(input$cbastembt, {
  output$cbastemnessplot<- renderUI({
    plotOutput("cbastemnessploting",
               width = paste0(isolate({input$cbastemwidth}), "%"),
               height = isolate({input$cbastemheight}))
  })})

output$cbasavestemness <- downloadHandler(
  filename = function() {
    paste0("Stemness-score-correlation-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cbastempltwidth}),height = isolate({input$cbastempltheight}))
    print(cbastemrun())
    dev.off()
  }
)

observeEvent(input$page_before_cbastemness, {
  newtab <- switch(input$tabs, "cbastemness" = "cbaimmune","cbaimmune" = "cbastemness")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbastemness, {
  newtab <- switch(input$tabs, "cbaestimate" = "cbastemness","cbastemness" = "cbaestimate")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









