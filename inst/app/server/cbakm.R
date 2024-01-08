
cbakmfeat <- reactive({
  req(input$cbakmcor)
  req(validst())
  res<-isolate({validst()})
  cohort<-isolate({input$cbakmcor})
  if(cohort=="Training set"){
    clinical<-res$traindata$clinical
  }else{
    clinical<-res$validata$clinical
  }
  colnames(clinical)
})

observeEvent(cbakmfeat(), {
  choices <- cbakmfeat()
  updateSelectInput(session, "cbasurvivaltime", choices = choices
,

selected = "OS.time"
  )
 }
)

observeEvent(cbakmfeat(), {
  choices <- cbakmfeat()
  updateSelectInput(session, "cbasurvivalstatus", choices = choices
,

selected = "OS"
)
}
)




plotcbaKM <- eventReactive(input$cbaKMplotbt,{
  
  input$cbaKMplotbt
  res<-isolate({validst()})
  cohort<-isolate({input$cbakmcor})

  time <- isolate({ input$cbasurvivaltime})
  status <- isolate({input$cbasurvivalstatus})

  if(cohort=="Training set"){
    clinical<-res$traindata$clinical
  }else{
    clinical<-res$validata$clinical
  }

  clinical<-clinical[complete.cases(clinical[,time]) & clinical[,time]>0 & complete.cases(clinical$Subtype),]
  
  OS.time<-clinical[,time]
  OS<-clinical[,status]
  if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
    createAlert(
      session,
      "cbakmmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "Either the survival time or survival status column is not numeric data, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
    
  }else if(any(levels(factor(OS))== c("0", "1"))!=T){
    createAlert(
      session,
      "cbakmmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival status column is not numerically in 1, and 0, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  } else{
    
    survxlab <- isolate({input$cbasurvxlab})
    survP <- isolate({input$cbasurvP})
    survRT <- isolate({input$cbasurvRT})
    survCI <- isolate({input$cbasurvCI})
    legend<-isolate({input$cbakmlg})
    require(survival)
    require(survminer)
    clinical$time<-clinical[,time]
    clinical$status<-clinical[,status]
    clinical$Subtype<-factor(gsub("Subtype ",'',clinical$Subtype,fixed = T))
    fit <-
      survfit(Surv(time,status) ~ Subtype, data = clinical)
  ggsurvplot(
      fit,
      data = clinical,
      risk.table = survRT,
      risk.table.height = 0.3,
      risk.table.y.text = FALSE,
      risk.table.title = "",
      main = "Survival curve",
      pval = survP,
      pval.method = T,
      conf.int = survCI,
      risk.table.y.text.col = T,
      legend = legend,
      legend.title = "",
      xlab = survxlab
    )

  }
})

observe({
  if(is.null(input$cbasurvivaltime) ||input$cbasurvivaltime == "" ||is.null(input$cbasurvivalstatus) ||input$cbasurvivalstatus == ""){
    disable("cbaKMplotbt")
  }
  else{
    enable("cbaKMplotbt")
  }
})


observeEvent(input$cbaKMplotbt, {
  updateCollapse(session, "collapsecbaKM", open = "Kaplan-Meier plot", close = "Descriptions and parameters")
})

observeEvent(input$cbaKMplotbt, {
  output$cbaKMploting <-  renderPlot(
    
    {print(plotcbaKM())}
  )
})

observeEvent(input$cbaKMplotbt, {
  
  output$cbaKMplot <- renderUI({
    plotOutput("cbaKMploting",
               
               width =   paste0(isolate({input$cbasurvwidth}),"%"),
               height = isolate({input$cbasurvheight}))
  })})
observeEvent(input$cbaKMplotbt, {
  shinyjs::show("cbakm_wrapper")
})

output$downloadcbaKM <- downloadHandler(
  
  filename = function() {
    paste0("K-M-curve-",Sys.Date(),'.pdf')
  },
  
  content = function(file) {
    
    pdf(file,width = isolate({input$cbakmwidth}), height =isolate({input$cbakmheight})) # open the pdf device
    print(plotcbaKM())
    dev.off()
  }
)

observeEvent(input$page_before_cbaKM, {
  newtab <- switch(input$tabs, "cbaclinical" = "cbaKM","cbaKM" = "cbaclinical")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbaKM, {
  newtab <- switch(input$tabs, "cbaKM" = "cbaCoxPH","cbaCoxPH" = "cbaKM")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})











