
value <- reactiveValues(data = NULL)
observeEvent(input$cbavalidationbt, {
  value$data <- isolate({validst()})
})

cbaoncopathrun<-eventReactive(input$cbaoncopathbt,{
  input$cbaoncopathbt
  require(GSVA)
  res<-isolate({validst()})
  cohort<-isolate({input$cbaoncopathcor})
  method<-isolate({input$cbaoncopathmeth})
  paired<-isolate({input$cbaoncopathpair})
  plottype<-isolate({input$cbaoncopathplot})

  withProgress(message = "Calculating the cancer pathway score",
               detail = "This may take a while...",
               value = 3,
               {
                 compgeneset(res=res,cohort=cohort,plottype=plottype,method=method,geneset=oncopath,name=names(oncopath))
                 })
})

observe({
  if(
    is.null(input$cbaoncopathmeth) ||
    input$cbaoncopathmeth == ""||
    is.null(input$cbaoncopathcor)||
    input$cbaoncopathcor == "" ||
    is.null(input$cbaoncopathpair)||
    input$cbaoncopathpair == "" ||
    is.null(input$cbaoncopathplot)||
    input$cbaoncopathplot == "" 
  ){
    disable("cbaoncopathbt")
  }
  else{
    enable("cbaoncopathbt")
  }
})

observeEvent(input$cbaoncopathbt, {
  disable("cbaoncopathbt")
})

observeEvent(input$cbaoncopathbt, {
  updateCollapse(session, "cbacollapseoncopath", open = "Comparision of cancer pathway", close = "Descriptions and parameters")
})

observe({
  if(is.null(value$data)){
    disable("cbaoncopathbt")
  }
  else{
    enable("cbaoncopathbt")
  }
})

observeEvent(cbaoncopathrun(), {
  shinyjs::show("cbaoncopathplot_wrapper")
  shinyjs::show("cbaoncopathtable_wrapper")
  enable("cbaoncopathbt")
})

observeEvent(input$cbaoncopathbt, {
  output$cbaoncopathtable <- DT::renderDT({
    DT::datatable(cbaoncopathrun()$table,options=list(scrollX=TRUE))
  })
})

output$cbasaveoncopathtable <- downloadHandler(
  filename = function() {
    paste0("Cancer-pathway-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(cbaoncopathrun()$table, file)
  }
)
observeEvent(input$cbaoncopathbt,{
  output$cbaoncopathploting  <- renderPlot({
    closeAlert(session, "cbaoncopathmess")
    ggarrange(plotlist=cbaoncopathrun()$plot,common.legend = T,align = "hv",labels="AUTO" )
  })
})

observeEvent(input$cbaoncopathbt, {
  output$cbaoncopathplot<- renderUI({
    plotOutput("cbaoncopathploting",
               width = paste0(isolate({input$cbaoncopathwidth}), "%"),
               height = isolate({input$cbaoncopathheight}))
  })})

output$cbasaveoncopathness <- downloadHandler(
  filename = function() {
    paste0("Cancer-pathway-comparision-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cbaoncopathpltwidth}),height = isolate({input$cbaoncopathpltheight}))
    print(ggarrange(plotlist=cbaoncopathrun()$plot,common.legend = T,align = "hv",labels="AUTO" ))
    dev.off()
  }
)

observeEvent(input$page_before_cbaoncopath, {
  newtab <- switch(input$tabs, "cbaoncopath" = "cbaCYA","cbaCYA" = "cbaoncopath")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbaoncopath, {
  newtab <- switch(input$tabs, "cbametapath" = "cbaoncopath","cbaoncopath" = "cbametapath")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









