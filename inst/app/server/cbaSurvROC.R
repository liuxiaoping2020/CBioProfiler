options(shiny.fullstacktrace = TRUE)
cbaROCfeat <- reactive({
  req(input$cbasurvROCcor)
  req(validst())
  # res<-isolate({validst()})
  res<-validst()
  cohort<-isolate({input$cbasurvROCcor})
  if(cohort=="Training set"){
    clinical<-res$traindata$clinical
  }else{
    clinical<-res$validata$clinical
  }
  colnames(clinical)
})

observeEvent(cbaROCfeat(), {
  choices <- cbaROCfeat()
  updateSelectInput(session, "cbaSurvROCtime", choices = choices,
                    selected = "OS.time"
    )
  }
)

observeEvent(cbaROCfeat(), {
  choices <- cbaROCfeat()
  updateSelectInput(session, "cbaSurvROCstatus", choices = choices,
                    selected = "OS"
    )
  }
)



cbayearlabelx <- reactive({
  req(input$cbaSurvROCtime)
  req(validst())
  
  req(input$cbasurvROCcor)
  # res<-isolate({validst()})
  res<-validst()
  SurvROCtime<-isolate({input$cbaSurvROCtime})
  cohort<-isolate({input$cbasurvROCcor})
  if(cohort=="Training set"){
    data<-res$traindata$clinical
  }else{
    data<-res$validata$clinical
  }
  
  if(SurvROCtime %in% names(data)){
    time<-data[,SurvROCtime]
    
    if(is.numeric(time)==F){
      createAlert(
        session,
        "cbasurvrocmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival time you selected is not numeric data, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    } else {
      
      if(max(na.omit(data[,SurvROCtime]))>600){
        by<- 365
      } else{
        by<-12
      }
      
      years<-seq(from=0,to=quantile(na.omit(data[,SurvROCtime]),0.95),by=by)
      years<-years[-1]
      yearlabel<-c()
      for(i in 1:length(years)){
        yearlabel[i]<- paste(i,"year",sep="-")
      }
      yearlabel
    }
  }
})

observeEvent(cbayearlabelx(), {
  choices <- cbayearlabelx()
  updateSelectInput(session, "cbapredictyear", choices = choices)
})


cbaSurvROCplot<-eventReactive(input$cbaSurvROCbt,{
  input$cbaSurvROCbt
  # browser()
  
  res<-isolate({validst()})
  cohort<-isolate({input$cbasurvROCcor})
  
  SurvROCtime <- isolate({ input$cbaSurvROCtime})
  SurvROCstatus <- isolate({input$cbaSurvROCstatus})
  
  if(cohort=="Training set"){
    clinical<-res$traindata$clinical
  }else{
    clinical<-res$validata$clinical
  }
  
  clinical<-clinical[complete.cases(clinical[,SurvROCtime]) & clinical[,SurvROCtime]>0 & complete.cases(clinical$Subtype),]

  method<-isolate({input$cbaSurvROCmethod})
  
  OS.time<-clinical[,SurvROCtime]
  OS<-clinical[,SurvROCstatus]
  
  if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
    createAlert(
      session,
      "cbasurvrocmess",
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
      "cbasurvrocmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival status column is not numerically in 1, and 0, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  }else{
    withProgress(
      message = "Constructing the prediction model",
      detail = "It may take a while, please be patient",
      
      value = 5,{
    require(survival)
    clinical$Subtype<-factor(gsub("Subtype ",'',clinical$Subtype,fixed = T))
    Surv<-Surv(clinical[,SurvROCtime],clinical[,SurvROCstatus])
    uni<-coxph(Surv~Subtype,clinical)
    clinical$pred<-predict(uni,clinical)


        survROC(
          data =clinical,
          bmtime = SurvROCtime,
          bmstatus = SurvROCstatus,
          marker = "pred",
          method = method,
          predyear = isolate({input$cbapredictyear}),
          cutpoint=isolate({input$cbacutoff})
        )
      })
    
    
  }}
  
)

observe({
  if (

      is.null(input$cbaSurvROCtime) ||
      input$cbaSurvROCtime == "" ||
      is.null(input$cbaSurvROCstatus) ||
      input$cbaSurvROCstatus == "" ||
      is.null(input$cbapredictyear) || input$cbapredictyear == "") {
    disable("cbaSurvROCbt")
  }
  else{
    enable("cbaSurvROCbt")
  }
})

observeEvent(input$cbaSurvROCbt, {
  updateCollapse(session, "cbacollapseSurvROC", open = "Time-dependent ROC analysis", close = "Descriptions and parameters")
})

observeEvent(input$cbaSurvROCbt, {
  output$cbaSurvROCplotting <- renderPlot({
    cbaSurvROCplot()
  })
})

observeEvent(input$cbaSurvROCbt, {
  output$cbaSurvROCplot <- renderUI({
    plotOutput("cbaSurvROCplotting",
               width = paste0(isolate({input$cbaSurvROCwidth}), "%"),
               height = isolate({input$cbaSurvROCheight}))
  })})

observeEvent(input$cbaSurvROCbt, {
  shinyjs::show("cbasurvROC_wrapper")
})

output$cbadownloadsurvROC <- downloadHandler(
  
  filename = function() {
    paste0("SurvivalROC-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cbasurvROCwidth.dl}),height = isolate({input$cbasurvROCheight.dl})) # open the pdf device
    print(cbaSurvROCplot())
    dev.off()
  }
)

observeEvent(input$page_before_cbaSurvROC, {
  newtab <- switch(input$tabs, "cbaCoxPH" = "cbaSurvROC","cbaSurvROC" = "cbaCoxPH")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbaSurvROC, {
  newtab <- switch(input$tabs, "cbaSurvROC" = "cbaDEG","cbaDEG" = "cbaSurvROC")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

