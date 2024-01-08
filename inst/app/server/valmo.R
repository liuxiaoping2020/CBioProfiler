observe({
  data <- rawdata2()$clinical
  updateSelectizeInput(session,
                       'Valtime',
                       choices = {
                         if (!is.null(data)) {
                           if (class(data) ==  class(data.frame()))
                             choices <- as.character(colnames(data))
                           if (class(data) !=  class(data.frame()))
                             choices <-
                               as.character(colnames(as.data.frame(data)))
                           choices
                         }
                       },
                       server = TRUE,
                       selected = "OS.time")
})

observe({
  data <- rawdata2()$clinical
  updateSelectizeInput(session,
                       'Valstatus',
                       choices = {
                         if (!is.null(data)) {
                           if (class(data) ==  class(data.frame()))
                             choices <- as.character(colnames(data))
                           if (class(data) !=  class(data.frame()))
                             choices <-
                               as.character(colnames(as.data.frame(data)))
                           choices
                         }
                       },
                       server = TRUE,
                       selected = "OS")
})

# yearlabel1 <- eventReactive(rawdata2(),{
yearlabel1 <- reactive({
  req(input$Valtime)
  Valtime<-paste(input$Valtime)
  data <- rawdata2()$clinical
  if(Valtime %in% names(data)){
  time<-data[,Valtime]
  if(is.numeric(time)==F){
    createAlert(
      session,
      "valmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival time you selected is not numeric data, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  } else {
  if(max(na.omit(data[,Valtime]))>600){
    by<- 365
  } else{
    by<-12
  }
  years1<-seq(from=0,to=quantile(na.omit(data[,Valtime]),0.95),by=by)
  years1<-years1[-1]
  yearlabel1<-c()
  for(i in 1:length(years1)){
    yearlabel1[i]<- paste(i,"year",sep="-")
  }
  yearlabel1
  }
  }
})

observeEvent(yearlabel1(), {
  choices <- yearlabel1()
  updateSelectInput(session, "Valpredictyear", choices = choices)
})




predres<-eventReactive(input$Valbt,{
  input$Valbt
  withProgress(
    message = "Starting resample experiment",
    detail = "According to your parameter settings, it may take a long time, please be patient",
    value = 5,{

  data<-isolate({rawdata2()})
  res<-isolate({nestrun()})

  expval<-data$expres
  expval<-as.data.frame(expval)

  selvar=res$`Fitted model`$features
  clinical<-data$clinical
  OS.time<-clinical[,isolate({input$Valtime})]
  
  OS<-clinical[,isolate({input$Valstatus})]
  if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
    createAlert(
      session,
      "valmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "Either survival time or survival status you selected is not numeric data, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)

  }else if(length(setdiff(selvar, row.names(expval)))!=0){
    createAlert(
      session,
      "valmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = paste(
        paste(paste(setdiff(selvar, row.names(expval)),collapse =","), "were not found in the expression file you provided")
      ),
      append = T,
      dismiss=T
    )
    return(NULL)
  } else if(any(levels(factor(OS))== c("0", "1"))!=T){
      createAlert(
        session,
        "valmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival status column is not numerically in 1, and 0, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    }else{

  predm(expclin=data$clinical,
                expval=data$expres,
                time=isolate({input$Valtime}),
                status=isolate({input$Valstatus}),
                lnid=res$lrnid,
                modeltrain=res$`Fitted model`)
  }
    })

})

observe({
  if (
      # is.null(input$data1) ||
      # input$data1 == "" ||
      is.null(input$Valtime) ||
      input$Valtime == "" ||
      is.null(input$Valstatus) ||
      input$Valstatus == "" || is.null(input$Valpredictyear) ||
      input$Valpredictyear == ""||
      is.null(input$Valcoxclinvar)||
      input$Valcoxclinvar==""
      )
  {
    disable("Valbt")
  }
  else{
    enable("Valbt")
  }
})


observeEvent(predres(), {
  shinyjs::show("valsurvROC_wrapper")
  shinyjs::show("valKMplot_wrapper")
  shinyjs::show("valCoxforest_wrapper")
  shinyjs::show("valCoxtable_wrapper")
})

vrocrun<-eventReactive(input$Valbt,{
  input$Valbt
  withProgress(
    message = "Starting resample experiment",
    detail = "According to your parameter settings, it may take a long time, please be patient",
    value = 5,{
  predyear<-isolate({input$Valpredictyear})

  clinical = isolate({predres()})

  OS.time<-clinical[,isolate({input$Valtime})]
  OS<-clinical[,isolate({input$Valstatus})]
  if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
    createAlert(
      session,
      "valmess",
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
      "valmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival status column is not numerically in 1, and 0, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  }else if(is.null(predyear)){
    createAlert(
      session,
      "valmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = paste(
        "Prediction years not provided, please provide the prediction years."),
      append = T,
      dismiss=T
    )
    return(NULL)
  }else{

 survROC(
    data = isolate({predres()}),
    bmtime = isolate({input$Valtime}),
    bmstatus = isolate({input$Valstatus}),
    marker = "Risk",
    method = isolate({input$ValSurvROCmethod}),
    predyear = predyear,
    cutpoint = isolate({input$Valcutoff})
  )
  }
    }
  )
})


observeEvent(input$Valbt,{
  output$valsurvROCing  <- renderPlot({
    closeAlert(session, "valmess")
    vrocrun()
  })
})

observeEvent(input$Valbt, {
  updateCollapse(session, "collapsevalmo", open = "Validation model", close = "Descriptions and parameters")
})

observeEvent(input$Valbt, {
  output$valsurvROC <- renderUI({
    plotOutput("valsurvROCing",
               width = paste0(isolate({input$ValSurvROCwidth}), "%"),
               height = isolate({input$ValSurvROCheight}))
  })})


output$savevalsurvROC <- downloadHandler(
  filename = function(){
    paste0("Survival-ROC-in-validation-set-",Sys.Date(),'.pdf')
  },
  content = function(file){
    pdf(file,width = isolate({input$valsurvROCwidthdl}),height = isolate({input$valsurvROCheightdl})) # open the pdf device
    print(vrocrun())
    dev.off()
  }
)


vkmrun<-eventReactive(input$Valbt,{
  input$Valbt
  clinical = isolate({predres()})

  OS.time<-clinical[,isolate({input$Valtime})]
  OS<-clinical[,isolate({input$Valstatus})]
  if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
    createAlert(
      session,
      "valmess",
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
      "valmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival status column is not numerically in 1, and 0, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  }else{
  km(data=isolate({predres()}),
     time=isolate({input$Valtime}),
     status=isolate({input$Valstatus}),
     marker="Risk",
     groupby=isolate({input$Valkmgroupby}),
     ratio= isolate({input$Valriskcut}),
     value=isolate({input$Valkmgpvalue}),
     high="High risk group",
     low="Low risk group",
     survxlab=isolate({input$Valsurvxlab}),
     survP=isolate({input$ValsurvP}),
     survRT=isolate({input$ValsurvRT}),
     survCI=isolate({input$ValsurvCI}),
     color1=isolate({input$Valkmcolor1}),
     color2=isolate({input$Valkmcolor2})
  )
}}
  )

observeEvent(input$Valbt,{
  output$valKMploting  <- renderPlot({
    closeAlert(session, "valmess")
    vkmrun()
  })
})

observeEvent(input$Valbt, {

  output$valKMplot <- renderUI({
    plotOutput("valKMploting",
               width = paste0(isolate({input$Valsurvwidth}), "%"),
               height = isolate({input$Valsurvheight}))
  })})


output$savevalvalKMplot <- downloadHandler(
  filename = function(){
    paste0("KM-curve-in-validation-set-",Sys.Date(),'.pdf')
  },
  content = function(file){
    pdf(file,width = isolate({input$valKMplotwidthdl}),height = isolate({input$valKMplotheightdl})) # open the pdf device
    print( vkmrun())
    dev.off()
  }
)


observeEvent(rawdata2(), {
  data <- rawdata2()$clinical
  updateSelectInput(session, "Valcoxclinvar", choices = {
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



vcoxrun<-eventReactive(input$Valbt,{
  input$Valbt
  res<-isolate({predres()})

  varname<-isolate({input$Valcoxfeat})
  feature<-isolate({input$Valcoxclinvar})
  marker<-"Risk"
  clinical<-res
  clinfeat<-subset(clinical,select=feature)
  name<-names(clinfeat)
  fc<-function(x){
    x<-factor(x)
    x<-length(levels(x))
    return(x)
  }
  clinfeature<-lapply(clinfeat,fc)

  OS.time<-clinical[,isolate({input$Valtime})]
  OS<-clinical[,isolate({input$Valstatus})]
  if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
    createAlert(
      session,
      "valmess",
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
      "valmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival status column is not numerically in 1, and 0, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  } else if(all(feature=="")!=F){
    createAlert(
      session,
      "valmess",
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
      "valmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = paste("Variable",paste(name[which(clinfeat==1)]),"you selected has only one level (only one value), which is not accepted by CoxPH model. Please check!"),
      append = T,
      dismiss = T
    )
    return(NULL)
  }else{
    feature<-c(feature,marker)

    var1 <- strsplit(varname, "|", fixed = T)[[1]]

    if(varname=="" || length(var1)==length(feature)){

      if(varname==""){
        varname<-NULL
      }else{
        varname <- strsplit(varname, "|", fixed = T)[[1]]
      }

      cox(data=res,
          time=isolate({input$Valtime}),
          status=isolate({input$Valstatus}),
          feature=feature,
          maxtick=isolate({input$Valmaxtick}),
          varname=varname,
          legend.pos=isolate({input$valegend.pos})
      )
      }else{

      createAlert(
        session,
        "valmess",
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



observeEvent(input$Valbt,{
  output$valCoxforesting  <- renderPlot({
    closeAlert(session, "valmess")
    vcoxrun()$plotcox
  })
})

observeEvent(input$Valbt, {

  output$valCoxforest <- renderUI({
    plotOutput("valCoxforesting",
               width = paste0(isolate({input$ValCoxwidth}), "%"),
               height = isolate({input$ValCoxheight}))
  })})


output$savevalvalCoxforest <- downloadHandler(
  filename = function(){
    paste0("Forestplot-in-validation-set-",Sys.Date(),'.pdf')
  },
  content = function(file){
    pdf(file,width = isolate({input$valCoxforestwidthdl}),height = isolate({input$valCoxforestheightdl})) # open the pdf device
    print(vcoxrun()$plotcox)
    dev.off()
  }
)


observeEvent(input$Valbt, {
  output$valCoxtable <-  DT::renderDT(
    {vcoxrun()$tablecox}
  )
})


output$savevalCoxtable <- downloadHandler(
  filename = function() {
    paste0("CoxPH-table-in-the-validation-set-",Sys.Date(), ".csv")
  },
  content = function(file) {

      write.csv(vcoxrun()$tablecox, file)

  }
)


observeEvent(input$page_before_valmo, {
  newtab <- switch(input$tabs, "bmpm" = "valmo","valmo" = "bmpm")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_valmo, {
  newtab <- switch(input$tabs, "nomo" = "valmo","valmo" = "nomo")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})








