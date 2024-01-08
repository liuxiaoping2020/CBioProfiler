observeEvent(rawdata(), {
  data <- rawdata()$clinical
  updateSelectInput(session, "bmtime", choices = {
    if (!is.null(data)) {
      if (class(data) ==  class(data.frame()))
        choices <- as.character(colnames(data))
      if (class(data) !=  class(data.frame()))
        choices <-
          as.character(colnames(as.data.frame(data)))
      choices
    }
  },
  selected = "OS.time"
  # server = TRUE
  )
})

observeEvent(rawdata(), {
  data <- rawdata()$clinical
  updateSelectInput(session, "bmstatus", choices = {
    if (!is.null(data)) {
      if (class(data) ==  class(data.frame()))
        choices <- as.character(colnames(data))
      if (class(data) !=  class(data.frame()))
        choices <-
          as.character(colnames(as.data.frame(data)))
      choices
    }
  }
  ,
  selected = "OS"
  )
})

yearlabel <- eventReactive(nestrun(),{
  req(input$bmtime)
  req(nestrun())
  res<-isolate({nestrun()})
  split<-isolate({input$split})
  if(split==F){
    data =res$clinical
  }else{
    data<-res$Testsurv
  }
  bmtime<-isolate({input$bmtime})
  # data <- rawdata()$clinical
  time<-data[,bmtime]
  if(is.numeric(time)==F){
    createAlert(
      session,
      "bmmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival time you selected is not numeric data, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  } else {
    if(max(na.omit(data[,bmtime]))>600){
    by<- 365
  } else{
    by<-12
  }
  years<-seq(from=0,to=quantile(na.omit(data[,bmtime]),0.95),by=by)
  years<-years[-1]
  yearlabel<-c()
  for(i in 1:length(years)){
    yearlabel[i]<- paste(i,"year",sep="-")
  }
  yearlabel
  }
})

observeEvent(yearlabel(), {
  choices <- yearlabel()
  updateSelectInput(session, "bmpredictyear", choices = choices)
})

rocrun<-eventReactive(input$pmbt,{
  input$pmbt
  res<-isolate({nestrun()})
  
predyear = isolate({input$bmpredictyear})
clinical<-res$clinical
OS.time<-clinical[,isolate({input$bmtime})]
OS<-clinical[,isolate({input$bmstatus})]
  if(is.null(predyear)){
    createAlert(
      session,
      "bmmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = paste(
        "Prediction years not provided, please provide the prediction years"),
      append = T,
      dismiss=T
    )
    return(NULL)
  }else if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
      createAlert(
        session,
        "bmmess",
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
        "bmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival status column is not numerically in 1, and 0, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    }else{

  split<-isolate({input$split})
  withProgress(
    message = "Constructing the prediction model",
    detail = "It may take a while, please be patient",
    value = 5,{
  if(split==F){
    survROC(
      data =res$clinical,
      bmtime = isolate({input$bmtime}),
      bmstatus = isolate({input$bmstatus}),
      marker = "Risk",
      method = isolate({input$bmSurvROCmethod}),
      predyear = predyear,
      cutpoint=isolate(input$bmcutoff)
    )

  }else{
    rocp <-
      lapply(
        list(res$Trainsurv, res$Testsurv),
        FUN = survROC,
        bmtime = isolate({input$bmtime}),
        bmstatus = isolate({input$bmstatus}),
        marker = "Risk",
        method = isolate({input$bmSurvROCmethod}),
        predyear = predyear,
        cutpoint = isolate(input$bmcutoff)
      )
    ggarrange(plotlist = rocp,ncol=2,nrow=1,align="hv",labels="AUTO")
  }
})
  }
  })

observe({
  if (is.null(input$bmtime) ||
      input$bmtime == "" ||
      is.null(input$bmstatus) ||
      input$bmstatus == "" ||
      is.null(input$bmpredictyear) ||
      input$bmpredictyear == "" || is.null(input$bmcoxclinvar) ||
      input$bmcoxclinvar == "")
  {
    disable("pmbt")
  }
  else{
    enable("pmbt")
  }
})


observeEvent(input$pmbt,{
  output$bmsurvROCing  <- renderPlot({
    closeAlert(session, "bmmess")
    rocrun()
  })
})

observeEvent(input$pmbt, {
  updateCollapse(session, "collapsebmpm", open = "Prediction model", close = "Descriptions and parameters")
})


observeEvent(input$pmbt, {
  output$bmsurvROC <- renderUI({
    plotOutput("bmsurvROCing",
               width = paste0(isolate({input$bmSurvROCwidth}), "%"),
               height = isolate({input$bmSurvROCheight}))
  })})

observeEvent(input$pmbt, {
  shinyjs::show("bmsurvROC_wrapper")
  shinyjs::show("bmKMplot_wrapper")
  shinyjs::show("bmCoxforest_wrapper")
  shinyjs::show("bmCoxtable_wrapper")
  disable("pmbt")
})
observeEvent(rocrun(), {
  enable("pmbt")
})


output$savebmsurvROC <- downloadHandler(
  filename = function(){
    paste0("Survival-ROC-validation-of-benchmark-experiment-",Sys.Date(),'.pdf')
  },
  content = function(file){
    pdf(file,width = isolate({input$bmsurvROCwidthdl}),height = isolate({input$bmsurvROCheightdl})) # open the pdf device
    print(rocrun())
    dev.off()
  }
)


kmrun<-eventReactive(input$pmbt,{
  input$pmbt
  res<-isolate({nestrun()})
  split<-isolate({input$split})
  clinical<-res$clinical
  OS.time<-clinical[,isolate({input$bmtime})]
  OS<-clinical[,isolate({input$bmstatus})]
  if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
    createAlert(
      session,
      "bmmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "Eihter the survival time or survival status column is not numeric data, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)

  }else if(any(levels(factor(OS))== c("0", "1"))!=T){
    createAlert(
      session,
      "bmmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival status column is not numerically in 1, and 0, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  }else {
  if(split==F){
    km(
      data=res$clinical,
      time=isolate({input$bmtime}),
      status=isolate({input$bmstatus}),
      marker="Risk",
      groupby=isolate({input$bmkmgroupby}),
      ratio= isolate({input$riskcut}),
      value=isolate({input$bmkmgpvalue}),
      high="High risk group",
      low="Low risk group",
      survxlab=isolate({input$bmsurvxlab}),
      survP=isolate({input$bmsurvP}),
      survRT=isolate({input$bmsurvRT}),
      survCI=isolate({input$bmsurvCI}),
      color1=isolate({input$bmkmcolor1}),
      color2=isolate({input$bmkmcolor2})
    )
  } else{
    kmp<-lapply(
      list(res$Trainsurv,res$Testsurv),
      FUN=km,
      time=isolate({input$bmtime}),
      status=isolate({input$bmstatus}),
      marker="Risk",
      groupby=isolate({input$bmkmgroupby}),
      ratio= isolate({input$riskcut}),
      value=isolate({input$bmkmgpvalue}),
      high="High risk group",
      low="Low risk group",
      survxlab=isolate({input$bmsurvxlab}),
      survP=isolate({input$bmsurvP}),
      survRT=isolate({input$bmsurvRT}),
      survCI=isolate({input$bmsurvCI}),
      color1=isolate({input$bmkmcolor1}),
      color2=isolate({input$bmkmcolor2})
    )
    if(isolate({input$bmsurvRT})==T){
    ggarrange(
      ggarrange(
        kmp[[1]]$plot,
        kmp[[1]]$table,
        ncol = 1,
        align = "hv",
        heights = c(1.8, 0.7)
      ),
      ggarrange(
        kmp[[2]]$plot,
        kmp[[2]]$table,
        ncol = 1,
        align = "hv",
        heights = c(1.8, 0.7)
      ),
      ncol = 2,
      align = "hv",
      labels="AUTO"
    ) } else{
      ggarrange(kmp[[1]]$plot,kmp[[2]]$plot,ncol = 2,align = "v",labels="AUTO")
    }
  }
  }
})


observeEvent(input$pmbt,{
  output$bmKMploting  <- renderPlot({
    closeAlert(session, "bmmess")
    kmrun()
  })
})

observeEvent(input$pmbt, {
  output$bmKMplot <- renderUI({
    plotOutput("bmKMploting",
               width = paste0(isolate({input$bmsurvwidth}), "%"),
               height = isolate({input$bmsurvheight}))
  })})

output$savebmKMplot <- downloadHandler(
  filename = function(){
    paste0("KM-curve-of-benchmark-experiment-",Sys.Date(),'.pdf')
  },
  content = function(file){
    pdf(file,width = isolate({input$bmKMplotwidthdl}),height = isolate({input$bmKMplotheightdl})) # open the pdf device
    print( kmrun())
    dev.off()
  }
)



observeEvent(rawdata(), {
  data <- rawdata()$clinical
  updateSelectInput(session, "bmcoxclinvar", choices = {
    if (!is.null(data)) {
      if (class(data) ==  class(data.frame()))
        choices <- as.character(colnames(data))
      if (class(data) !=  class(data.frame()))
        choices <-
          as.character(colnames(as.data.frame(data)))
      choices
    }
  },
  selected = "Age"
  )
})


coxrun<-eventReactive(input$pmbt,{
  input$pmbt
  res<-isolate({nestrun()})
  split<-isolate({input$split})
  varname<-isolate({input$bmcoxfeat})
   # print(varname)
   # print(length(varname))
   # print(mode(varname))
  feature<-isolate({input$bmcoxclinvar})
  clinical<-res$clinical
  clinfeat<-subset(clinical,select=feature)
  name<-names(clinfeat)
  fc<-function(x){
    x<-factor(x)
    x<-length(levels(x))
    return(x)
  }
  clinfeature<-lapply(clinfeat,fc)
  marker<-"Risk"

  # clinical<-res$clinical
  OS.time<-clinical[,isolate({input$bmtime})]
  OS<-clinical[,isolate({input$bmstatus})]
  if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
    createAlert(
      session,
      "bmmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "Eihter the survival time or survival status column is not numeric data, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  }else if(any(levels(factor(OS))== c("0", "1"))!=T){
    createAlert(
      session,
      "bmmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival status column is not numerically in 1, and 0, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  }else if(all(feature=="")!=F){
          createAlert(
            session,
            "bmmess",
            "exampleAlert",
            title = "Please note!",
            style =  "danger",
            content = paste(
              "No clinical covariates are provided, please provide the covarites you want to include in the CoxPH model"),
            append = T,
            dismiss=T
          )
          return(NULL)

        }else if(1 %in% clinfeature){
            createAlert(
              session,
              "bmmess",
              "exampleAlert",
              title = "Please note!",
              style =  "danger",
              content = paste("Variable",paste(name[which(clinfeat==1)]),"you selected has only one level (only one value), which is not accepted by CoxPH model. Please check!"),
              append = T,
              dismiss = T
            )
            return(NULL)
          }else{
          var1 <- strsplit(varname, "|", fixed = T)[[1]]
          if(varname=="" || length(var1)==(length(feature)+1)){
  cox1(
    data = res,
    time = isolate({input$bmtime}),
    status = isolate({input$bmstatus}),
    feature = feature,
    marker = marker,
    maxtick = isolate({input$bmmaxtick}),
    split = split,
    varname = varname,
    legend.pos=isolate({input$bmlegend.pos})
  )
      }else{

            createAlert(
              session,
              "bmmess",
              "exampleAlert",
              title = "Please note!",
              style =  "danger",
              content = "The length of variable names does not equal to the clinical variables included in the CoxPH model",
              append = T,
              dismiss = T
            )
            return(NULL)
          }
        }
})


observeEvent(input$pmbt,{
  output$bmCoxforesting  <- renderPlot({
    closeAlert(session, "bmmess")
    if(isolate({input$split})==F){
      coxrun()$plotcox
    }else{
      ggarrange(
        coxrun()[[1]]$plotcox,
        coxrun()[[2]]$plotcox,
        nrow = 2,
        labels="AUTO",
        align ="hv"
      )
    }

  })
})

observeEvent(input$pmbt, {
  output$bmCoxforest <- renderUI({
    plotOutput("bmCoxforesting",
               width = paste0(isolate({input$bmCoxwidth}), "%"),
               height = isolate({input$bmCoxheight}))
  })})

output$savebmCoxforest <- downloadHandler(
  filename = function(){
    paste0("Cox-forest-plot-of-benchmark-experiment-",Sys.Date(),'.pdf')
  },
  content = function(file){
    pdf(file,width = isolate({input$bmCoxforestwidthdl}),height = isolate({input$bmCoxforestheightdl})) # open the pdf device
    if(isolate({input$split})==F){
      print(coxrun()$plotcox)
    }else{
      print(ggarrange(
        coxrun()[[1]]$plotcox,
        coxrun()[[2]]$plotcox,
        nrow = 2,
        labels="AUTO",
        align ="hv"
      ))
    }
    dev.off()
  }
)

observeEvent(input$pmbt, {
  output$bmCoxtable <-  DT::renderDT(
    {
      if(isolate({input$split})==F){
      coxrun()$tablecox
    }else{
      tab<-list(TrainingSet=coxrun()[[1]]$tablecox, TestSet=coxrun()[[2]]$tablecox)
      do.call(rbind,tab)
    }
    }
  )
})


output$savebmCoxtable <- downloadHandler(
  filename = function() {
    paste0("CoxPH-table-of-benchmark-experiment-",Sys.Date(), ".csv")
  },
  content = function(file) {

    if(isolate({input$split})==F){
      write.csv(coxrun()$tablecox, file)
    }else{
      tab<-list(TrainingSet=coxrun()[[1]]$tablecox, TestSet=coxrun()[[2]]$tablecox)
      do.call(rbind,tab)
      write.csv(do.call(rbind,tab), file)
    }

  }
)


observeEvent(input$page_before_bmpm, {
  newtab <- switch(input$tabs, "bmpm" = "nestr","nestr" = "bmpm")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_bmpm, {
  newtab <- switch(input$tabs, "bmpm" = "valmo","valmo" = "bmpm")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})






























