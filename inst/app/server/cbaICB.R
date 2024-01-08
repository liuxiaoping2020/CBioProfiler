
value <- reactiveValues(data = NULL)
observeEvent(input$cbavalidationbt, {
  value$data <- isolate({validst()})
})

cbaICBrun<-eventReactive(input$cbaICBbt,{
  input$cbaICBbt

  res<-isolate({validst()})
  cohort<-isolate({input$cbaICBcor})

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
  
  ICB<-c("PDCD1","CD274","PDCD1LG2","CTLA4","PVR","LAG3","TIGIT","HAVCR2","VTCN1","CD86","CD28", "CD80","IDO1",
         "CD27","CD40","IL2RB","TNFRSF9","TNFRSF4","TNFRSF18","ICOS","CD276","BTLA","KIR3DL1","CYBB","VSIR","SIGLEC7")
  
  method<-isolate({input$cbaICBmeth})
  paired<-isolate({input$cbaICBpair})
  plottype<-isolate({input$cbaICBplot})
  
  withProgress(message = "Calculating the ICB distribution",
               detail = "This may take a while...",
               value = 3,
               {
                 if(length(intersect(ICB,names(expres)))==0){
                   createAlert(
                     session,
                     "cbaICBmess",
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
                     "cbaICBmess",
                     "exampleAlert",
                     title = "Please note!",
                     style =  "info",
                     content = paste(paste(setdiff(ICB,names(expres)),collapse = ', '), "were not found in the gene expression data") ,
                     append = T,
                     dismiss=T
                   )
                   ICB<-setdiff(ICB,setdiff(ICB,names(expres)))
                   compicb(ICB=ICB,expres=expres,clinical=clinical,plottype=plottype,method=method, paired=paired)
                   
                 } else {
                   closeAlert(session, "cbaICBmess")
                   compicb(ICB=ICB,expres=expres,clinical=clinical,plottype=plottype,method=method, paired=paired)
                   }
               })
})

observe({
  if(
    is.null(input$cbaICBmeth) ||
    input$cbaICBmeth == ""||
    is.null(input$cbaICBcor)||
    input$cbaICBcor == "" ||
    is.null(input$cbaICBpair)||
    input$cbaICBpair == "" ||
    is.null(input$cbaICBplot)||
    input$cbaICBplot == "" 
  ){
    disable("cbaICBbt")
  }
  else{
    enable("cbaICBbt")
  }
})

observeEvent(input$cbaICBbt, {
  updateCollapse(session, "cbacollapseICB", open = "Comparision of immune checkpoint distribtion", close = "Descriptions and parameters")
})

observe({
  if(is.null(value$data)){
    disable("cbaICBbt")
  }
  else{
    enable("cbaICBbt")
  }
})

observeEvent(cbaICBrun(), {
  shinyjs::show("cbaICBplot_wrapper")
  shinyjs::show("cbaICBtable_wrapper")
})

observeEvent(input$cbaICBbt, {
  output$cbaICBtable <- DT::renderDT({
    DT::datatable(cbaICBrun()$expres,options=list(scrollX=TRUE))
  })
})

output$cbasaveICBtable <- downloadHandler(
  filename = function() {
    paste0("ICB-expression-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(cbaICBrun()$expres, file)
  }
)
observeEvent(input$cbaICBbt,{
  output$cbaICBploting  <- renderPlot({
    closeAlert(session, "cbaICBmess")
    ggarrange(plotlist=cbaICBrun()$plot,common.legend = T,align = "hv",labels="AUTO" )
  })
})

observeEvent(input$cbaICBbt, {
  output$cbaICBplot<- renderUI({
    plotOutput("cbaICBploting",
               width = paste0(isolate({input$cbaICBwidth}), "%"),
               height = isolate({input$cbaICBheight}))
  })})

output$cbasaveICBness <- downloadHandler(
  filename = function() {
    paste0("ICB-expression-comparision-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cbaICBpltwidth}),height = isolate({input$cbaICBpltheight}))
    print(ggarrange(plotlist=cbaICBrun()$plot,common.legend = T,align = "hv",labels="AUTO" ))
    dev.off()
  }
)

observeEvent(input$page_before_cbaICB, {
  newtab <- switch(input$tabs, "cbaICB" = "cbaestimate","cbaestimate" = "cbaICB")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbaICB, {
  newtab <- switch(input$tabs, "cbaIFN" = "cbaICB","cbaICB" = "cbaIFN")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









