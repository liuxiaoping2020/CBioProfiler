observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'KMgene',
                       choices = {
                         if (!is.null(data)) {
                           if (class(data) ==  class(data.frame()))
                             choices <- as.character(row.names(data))
                           if (class(data) !=  class(data.frame()))
                             choices <-
                               as.character(row.names(as.data.frame(data)))
                           choices
                         }
                       },
                       server = TRUE,
                       selected = "TP53")
})


observe({
  data <- rawdata()$clinical
  updateSelectizeInput(session,
                       'survivaltime',
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
  data <- rawdata()$clinical
  updateSelectizeInput(session,
                       'survivalstatus',
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

plotKM <- eventReactive(input$KMplotbt,{
  input$KMplotbt
  data <- isolate({
    rawdata()
  })

  data <- lapply(data, as.data.frame)
  time <- isolate({ input$survivaltime})
  status <- isolate({input$survivalstatus})
  data$clinical <-data$clinical[complete.cases(data$clinical[, time]) & data$clinical[, time] > 0,]
  index <-intersect(colnames(data$expres), row.names(data$clinical))
  clinical <- data$clinical[index,]
  clinical[clinical == ""] <- NA
  expres <- data$expres[, index]
  KMgene <- isolate({
    input$KMgene
  })

  clinical$KMgene <- as.numeric(expres[KMgene,])
  clinical<-clinical[complete.cases(clinical$KMgene),]


  OS.time<-clinical[,time]
  OS<-clinical[,status]
  if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
    createAlert(
      session,
      "kmmess",
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
      "kmmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival status column is not numerically in 1, and 0, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  } else if(isolate({input$kmgroupby})=="Value"&& isolate({input$kmgpvalue})>=max(clinical$KMgene) ){

      createAlert(
        session,
        "kmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste("The cutoff value you set exceeds the range: (", paste(range(clinical$KMgene),collapse = ", "), ") of the expression level of", paste(KMgene), "you specified, please check!"),
        append = T,
        dismiss=T
      )
      return(NULL)

  }else if(isolate({input$kmgroupby})=="Value" && isolate({input$kmgpvalue})<= min(clinical$KMgene)){
    createAlert(
      session,
      "kmmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = paste("The cutoff value you set exceeds the range: (", paste(range(clinical$KMgene),collapse = ", "), ") of the expression level of", paste(KMgene), "you specified, please check!"),
      append = T,
      dismiss=T
    )
    return(NULL)
  }else {
  genecut <- isolate({input$genecut})
  if(isolate({input$kmgroupby})=="Percentage"){
  clinical$KMgroup <-ifelse(
      clinical$KMgene <= quantile(clinical$KMgene, genecut),
      paste(KMgene, "low expression group", sep = " "),
      paste(KMgene, "high expression group", sep = " ")
    )
  }else{

    clinical$KMgroup<-ifelse(clinical$KMgene<=isolate({input$kmgpvalue}),paste(KMgene,"low expression group",sep=" "),paste(KMgene,"high expression group", sep=" "))
    }

  survxlab <- isolate({input$survxlab})
  survP <- isolate({input$survP})
  survRT <- isolate({input$survRT})
  survCI <- isolate({input$survCI})
  KMcurvename <- isolate({input$KMcurvename})
  clinical$time <- clinical[, time]
  clinical$status <- clinical[, status]
  require(survival)
  require(survminer)
  fit <- survfit(Surv(time, status) ~ KMgroup, data = clinical)
  level <- levels(factor(clinical$KMgroup))
   c <-
    ifelse(
      level[1] == paste(KMgene, "low expression group", sep = " "),

      paste(isolate({input$kmcolor2})),

      paste(isolate({input$kmcolor1}))
    )

  d <-
    ifelse(
      level[2] == paste(KMgene, "high expression group", sep = " "),
      paste(isolate({input$kmcolor1})),
      paste(isolate({input$kmcolor2}))
    )

  e <-
    ifelse(
      level[1] == paste(KMgene, "low expression group", sep = " "),
      paste(KMgene, "low expression group", sep = " "),
      paste(KMgene, "high expression group", sep = " ")
    )

  f <-
    ifelse(
      level[2] == paste(KMgene, "low expression group", sep = " "),
      paste(KMgene, "low expression group", sep = " "),
      paste(KMgene, "high expression group", sep = " ")
    )

  res <- ggsurvplot(
    fit,
    data = clinical,
    risk.table = survRT,
    risk.table.height = 0.3,
    risk.table.y.text = FALSE,
    risk.table.title = "",
    main = "Survival curve",
    palette = c(c, d),
    pval = survP,
    pval.method = T,
    conf.int = survCI,
    risk.table.y.text.col = T,
    legend = isolate({input$kmlg}),#c(0.8, 0.90),
    legend.title = "",
    xlab = survxlab,
    legend.labs = c(e, f)
  )
  res
}
})

observe({
  if(is.null(input$KMgene) || input$KMgene== "" ||is.null(input$survivaltime) ||input$survivaltime == "" ||is.null(input$survivalstatus) ||input$survivalstatus == ""){
    disable("KMplotbt")
  }
  else{
    enable("KMplotbt")
  }
})


observeEvent(input$KMplotbt, {
  updateCollapse(session, "collapseKM", open = "Kaplan-Meier plot", close = "Descriptions and parameters")
})

observeEvent(input$KMplotbt, {
  output$KMploting <-  renderPlot(

    {print(plotKM())}
    )
})

observeEvent(input$KMplotbt, {

  output$KMplot <- renderUI({
    plotOutput("KMploting",

               width =   paste0(isolate({input$survwidth}),"%"),
               height = isolate({input$survheight}))
  })})
observeEvent(input$KMplotbt, {
  shinyjs::show("km_wrapper")
})




output$downloadKM <- downloadHandler(

  filename = function() {
    paste0("K-M-curve-",Sys.Date(),'.pdf')
  },

  content = function(file) {

    pdf(file,width = isolate({input$kmwidth}), height =isolate({input$kmheight})) # open the pdf device
    print(plotKM())
    dev.off()
  }
)

observeEvent(input$page_before_KM, {
  newtab <- switch(input$tabs, "clinical" = "KM","KM" = "clinical")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_KM, {
  newtab <- switch(input$tabs, "KM" = "CoxPH","CoxPH" = "KM")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})











