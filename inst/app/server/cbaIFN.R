
value <- reactiveValues(data = NULL)
observeEvent(input$cbavalidationbt, {
  value$data <- isolate({validst()})
})

cbaIFNrun<-eventReactive(input$cbaIFNbt,{
  input$cbaIFNbt
  
  res<-isolate({validst()})
  cohort<-isolate({input$cbaIFNcor})
  method<-isolate({input$cbaIFNmeth})
  paired<-isolate({input$cbaIFNpair})
  plottype<-isolate({input$cbaIFNplot})
  IFN<-c('IFNG',"IFNGR1","IFNGR2","IRF1","JAK1","JAK2","STAT1")
  name<-"Interferon-gamma"
  
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
  
  withProgress(message = "Calculating the IFN score",
               detail = "This may take a while...",
               value = 3,
               {
                 if(length(intersect(IFN,row.names(expres)))==0){
                   createAlert(
                     session,
                     "cbaIFNmess",
                     "exampleAlert",
                     title = "Please note!",
                     style =  "danger",
                     content = paste(
                       "No Interferon-gamma genes were not found in the gene expression profile you specified. Please check!!!"
                     ),
                     append = T,
                     dismiss=T
                   )
                   return(NULL)
                 } else {
                   closeAlert(session, "cbaIFNmess")
                   require(GSVA)
                   compIFN(expres=expres,clinical=clinical,IFN=IFN,plottype=plottype,method=method,name=name)
                 }
               })
})

observe({
  if(
    is.null(input$cbaIFNmeth) ||
    input$cbaIFNmeth == ""||
    is.null(input$cbaIFNcor)||
    input$cbaIFNcor == "" ||
    is.null(input$cbaIFNpair)||
    input$cbaIFNpair == "" ||
    is.null(input$cbaIFNplot)||
    input$cbaIFNplot == "" 
  ){
    disable("cbaIFNbt")
  }
  else{
    enable("cbaIFNbt")
  }
})

observeEvent(input$cbaIFNbt, {
  updateCollapse(session, "cbacollapseIFN", open = "Comparision of interferon-gamma distribtion", close = "Descriptions and parameters")
})

observe({
  if(is.null(value$data)){
    disable("cbaIFNbt")
  }
  else{
    enable("cbaIFNbt")
  }
})

observeEvent(cbaIFNrun(), {
  shinyjs::show("cbaIFNplot_wrapper")
  shinyjs::show("cbaIFNtable_wrapper")
})

observeEvent(input$cbaIFNbt, {
  output$cbaIFNtable <- DT::renderDT({
    DT::datatable(cbaIFNrun()$IFN,options=list(scrollX=TRUE))
  })
})

output$cbasaveIFNtable <- downloadHandler(
  filename = function() {
    paste0("IFN-expression-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(cbaIFNrun()$IFN, file)
  }
)
observeEvent(input$cbaIFNbt,{
  output$cbaIFNploting  <- renderPlot({
    closeAlert(session, "cbaIFNmess")
    ggarrange(plotlist=cbaIFNrun()$plot,common.legend = T,align = "hv",labels="AUTO" )
  })
})

observeEvent(input$cbaIFNbt, {
  output$cbaIFNplot<- renderUI({
    plotOutput("cbaIFNploting",
               width = paste0(isolate({input$cbaIFNwidth}), "%"),
               height = isolate({input$cbaIFNheight}))
  })})

output$cbasaveIFNness <- downloadHandler(
  filename = function() {
    paste0("Interferon-gamma-score-comparision-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cbaIFNpltwidth}),height = isolate({input$cbaIFNpltheight}))
    print(ggarrange(plotlist=cbaIFNrun()$plot,common.legend = T,align = "hv",labels="AUTO" ))
    dev.off()
  }
)

observeEvent(input$page_before_cbaIFN, {
  newtab <- switch(input$tabs, "cbaIFN" = "cbaICB","cbaICB" = "cbaIFN")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbaIFN, {
  newtab <- switch(input$tabs, "cbaCYA" = "cbaIFN","cbaIFN" = "cbaCYA")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









