
value <- reactiveValues(data = NULL)
observeEvent(input$cbavalidationbt, {
  value$data <- isolate({validst()})
})

cbaCYArun<-eventReactive(input$cbaCYAbt,{
  input$cbaCYAbt
  res<-isolate({validst()})
  cohort<-isolate({input$cbaCYAcor})
  method<-isolate({input$cbaCYAmeth})
  paired<-isolate({input$cbaCYApair})
  plottype<-isolate({input$cbaCYAplot})
  name<-"Cytotoxic activity"
  
  
  if(cohort=="Training set"){
    expres<-as.data.frame(res$traindata$expres)
    clinical<-res$traindata$clinical
  }else{
    expres<-as.data.frame(res$validata$expres)
    clinical<-res$validata$clinical
  }
  
  index<-intersect(row.names(clinical),names(expres))
  clinical<-clinical[index,]
  expres<-expres[,index]
  expres<-as.data.frame(t(expres))
  expres$Subtype<-clinical$Subtype
  withProgress(message = "Calculating the CYA score",
               detail = "This may take a while...",
               value = 3,
               {
                 if(length(setdiff(c("GZMA","PRF1"),names(expres)))!=0){
                   createAlert(
                     session,
                     "cbaCYAmess",
                     "exampleAlert",
                     title = "Please note!",
                     style =  "danger",
                     content = paste(setdiff(c("GZMA","PRF1"),names(expres))), "was not found in the expression profile",
                     append = T,
                     dismiss=T
                   )
                   return(NULL)
                 } else {
                   closeAlert(session, "cbaCYAmess")
                   require(GSVA)
                   compca(expres=expres,clinical=clinical,paired=paired,plottype=plottype,method=method,name=name)
                 }
               })
})

observe({
  if(
    is.null(input$cbaCYAmeth) ||
    input$cbaCYAmeth == ""||
    is.null(input$cbaCYAcor)||
    input$cbaCYAcor == "" ||
    is.null(input$cbaCYApair)||
    input$cbaCYApair == "" ||
    is.null(input$cbaCYAplot)||
    input$cbaCYAplot == "" 
  ){
    disable("cbaCYAbt")
  }
  else{
    enable("cbaCYAbt")
  }
})

observeEvent(input$cbaCYAbt, {
  disable("cbaCYAbt")
})

observeEvent(input$cbaCYAbt, {
  updateCollapse(session, "cbacollapseCYA", open = "Comparision of cytotoxic activity distribtion", close = "Descriptions and parameters")
})

observe({
  if(is.null(value$data)){
    disable("cbaCYAbt")
  }
  else{
    enable("cbaCYAbt")
  }
})

observeEvent(cbaCYArun(), {
  shinyjs::show("cbaCYAplot_wrapper")
  shinyjs::show("cbaCYAtable_wrapper")
  enable("cbaCYAbt")
})

observeEvent(input$cbaCYAbt, {
  output$cbaCYAtable <- DT::renderDT({
    DT::datatable(cbaCYArun()$table,options=list(scrollX=TRUE))
  })
})

output$cbasaveCYAtable <- downloadHandler(
  filename = function() {
    paste0("Cytotoxic-activity-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(cbaCYArun()$table, file)
  }
)
observeEvent(input$cbaCYAbt,{
  output$cbaCYAploting  <- renderPlot({
    closeAlert(session, "cbaCYAmess")
    ggarrange(plotlist=cbaCYArun()$plot,common.legend = T,align = "hv",labels="AUTO" )
  })
})

observeEvent(input$cbaCYAbt, {
  output$cbaCYAplot<- renderUI({
    plotOutput("cbaCYAploting",
               width = paste0(isolate({input$cbaCYAwidth}), "%"),
               height = isolate({input$cbaCYAheight}))
  })})

output$cbasaveCYAness <- downloadHandler(
  filename = function() {
    paste0("Cytotoxic-activity-comparision-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cbaCYApltwidth}),height = isolate({input$cbaCYApltheight}))
    print(ggarrange(plotlist=cbaCYArun()$plot,common.legend = T,align = "hv",labels="AUTO" ))
    dev.off()
  }
)

observeEvent(input$page_before_cbaCYA, {
  newtab <- switch(input$tabs, "cbaCYA" = "cbaIFN","cbaIFN" = "cbaCYA")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbaCYA, {
  newtab <- switch(input$tabs, "cbaoncopath" = "cbaCYA","cbaCYA" = "cbaoncopath")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









