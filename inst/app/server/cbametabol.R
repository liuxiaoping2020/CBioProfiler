
value <- reactiveValues(data = NULL)
observeEvent(input$cbavalidationbt, {
  value$data <- isolate({validst()})
})

cbametapathrun<-eventReactive(input$cbametapathbt,{
  input$cbametapathbt
  
  res<-isolate({validst()})
  cohort<-isolate({input$cbametapathcor})
  method<-isolate({input$cbametapathmeth})
  paired<-isolate({input$cbametapathpair})
  plottype<-isolate({input$cbametapathplot})
  
  withProgress(message = "Calculating the metabolism pathway score",
               detail = "This may take a while...",
               value = 3,
               {
                 compgeneset(res=res,cohort=cohort,plottype=plottype,method=method,geneset=metabolism,name=names(metabolism))
               })
})

observe({
  if(
    is.null(input$cbametapathmeth) ||
    input$cbametapathmeth == ""||
    is.null(input$cbametapathcor)||
    input$cbametapathcor == "" ||
    is.null(input$cbametapathpair)||
    input$cbametapathpair == "" ||
    is.null(input$cbametapathplot)||
    input$cbametapathplot == "" 
  ){
    disable("cbametapathbt")
  }
  else{
    enable("cbametapathbt")
  }
})

observeEvent(input$cbametapathbt, {
  updateCollapse(session, "cbacollapsemetapath", open = "Comparision of metabolisim pathway", close = "Descriptions and parameters")
})

observe({
  if(is.null(value$data)){
    disable("cbametapathbt")
  }
  else{
    enable("cbametapathbt")
  }
})

observeEvent(cbametapathrun(), {
  shinyjs::show("cbametapathplot_wrapper")
  shinyjs::show("cbametapathtable_wrapper")
})

observeEvent(input$cbametapathbt, {
  output$cbametapathtable <- DT::renderDT({
    DT::datatable(cbametapathrun()$table,options=list(scrollX=TRUE))
  })
})

output$cbasavemetapathtable <- downloadHandler(
  filename = function() {
    paste0("Metabolism-pathway-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(cbametapathrun()$table, file)
  }
)
observeEvent(input$cbametapathbt,{
  output$cbametapathploting  <- renderPlot({
    closeAlert(session, "cbametapathmess")
    ggarrange(plotlist=cbametapathrun()$plot,common.legend = T,align = "hv",labels="AUTO",nrow=5,ncol=3 )
  })
})

observeEvent(input$cbametapathbt, {
  output$cbametapathplot<- renderUI({
    plotOutput("cbametapathploting",
               width = paste0(isolate({input$cbametapathwidth}), "%"),
               height = isolate({input$cbametapathheight}))
  })})

output$cbasavemetapathness <- downloadHandler(
  filename = function() {
    paste0("Metabolism-comparision-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cbametapathpltwidth}),height = isolate({input$cbametapathpltheight}))
    print(ggarrange(plotlist=cbametapathrun()$plot,common.legend = T,align = "hv",labels="AUTO" ,nrow=5,ncol=3))
    dev.off()
  }
)

observeEvent(input$page_before_cbametapath, {
  newtab <- switch(input$tabs, "cbametapath" = "cbaoncopath","cbaoncopath" = "cbametapath")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbametapath, {
  newtab <- switch(input$tabs, "cbahall" = "cbametapath","cbametapath" = "cbahall")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









