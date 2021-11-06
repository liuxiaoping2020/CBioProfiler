observe({
  data <- rawdata()$clinical
  updateSelectizeInput(session,
                       'nrtime',
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
                       'nrstatus',
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

value <- reactiveValues(data = NULL)
observeEvent(input$MTRbt, {
  value$data <- isolate({netsummary()})

})

observeEvent(input$DEGvisbt, {
  value$data <- isolate({DEGop()})

})

observeEvent(input$msurvopbt, {
  value$data <- isolate({msurvsig()})
})

  nestrun<-eventReactive(input$nrbt,{
  input$nrbt
  req(value$data)

  data <- isolate(value$data)

  withProgress(
    message = "Starting resample experiment",
    detail = "According to your parameter settings, it may take a long time, please be patient",
    value = 5,{

      OS.time<-isolate({input$nrtime})
      OS<-isolate({input$nrstatus})
      clinical <- data$clinical
      expres <- data$expres

      clinical<-clinical[complete.cases(clinical[[OS.time]]) & complete.cases(clinical[[OS]]) & clinical[[OS.time]]>0,]

      os.time<-clinical[,OS.time]
      os<-clinical[,OS]
      if(mode(os.time)!="numeric" || mode(os)!="numeric" ){
        createAlert(
          session,
          "nestrmess",
          "exampleAlert",
          title = "Please note!",
          style =  "danger",
          content = "Either the survival time or survival status column is not numeric data, please check!",
          append = T,
          dismiss=T
        )
        return(NULL)

      }else if(any(levels(factor(os))== c("0", "1"))!=T){
        createAlert(
          session,
          "nestrmess",
          "exampleAlert",
          title = "Please note!",
          style =  "danger",
          content = "The survival status column is not numerically in 1, and 0, please check!",
          append = T,
          dismiss=T
        )

        return(NULL)
      }else{

      index <- intersect(names(expres), row.names(clinical))
      clinical <- clinical[index, ]
      expres <- expres[, index]
      expres<-as.data.frame(t(expres))
      genesel<-row.names(data$Output)
      expres<-expres[,genesel]
      names(expres)<-make.names(names(expres))

      optimization<-isolate({input$optimization})
      if(optimization=="Grid search"){
        ctrl = makeTuneControlGrid(resolution=isolate({input$resolution}))
      } else{
        ctrl = makeTuneControlRandom(maxit=isolate({input$maxit}))
      }

      inner = makeResampleDesc("CV", iters = isolate({input$ifold}))
      outer = makeResampleDesc("CV", iters = isolate({input$ofold}))

      Learner<- isolate({input$Learner})
      if("RandomForestSRC" %in% Learner){

          createAlert(
            session,
            "nestrmess",
            "exampleAlert",
            title = "Please note!",
            style =  "danger",
            content = "Random forest is a machine learning method based on ensemble learning. Its calculation process will consume a lot of time and computing resources, especially when integrated into nested cross validation. Considering the limited computing power of the server, we recommend that users download and install the standalone application of CBioExplorer to perform related calculations locally.",
            append = T,
            dismiss=T
          )
      }

      split<-isolate({input$split})
      ratio<-isolate({input$sratio})
      biter<-isolate({input$biter})

      rdesc<-makeResampleDesc("Bootstrap", iters = biter)

      if(isolate({input$validat})=="Nested cross-validation"){
      res<-nestml(
        expres = expres,
        clinical = clinical,
        endpoint = c(OS.time,OS),
        ratio = ratio,
        ctrl = ctrl,
        inner = inner,
        outer = outer,
        lnid = Learner,
        split=split
      )

      }else {
        res<-cv(expres=expres,
                clinical=clinical,
                split=split,
                endpoint =c(OS.time,OS),
                ratio=ratio,
                lnid=Learner,
                inner=inner,
                ctrl=ctrl,
                rdesc=rdesc
                )

      }
  if(is.null(res)){
    createAlert(
      session,
      "nestrmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = paste(
        "At least one iteration of the model failed to converge, please reset the data or parameters."),
      append = T,
      dismiss=T
    )
    return(NULL)
  }else if(length(res[[2]])==0){
  createAlert(
    session,
    "nestrmess",
    "exampleAlert",
    title = "Please note!",
    style =  "danger",
    content = paste(
      "No biomarkers identified based on the benchmark experiment you specified."),
    append = T,
    dismiss=T
  )
  return(NULL)

}else{
  res
} }
    }
  )
})

observe({
  if(
    is.null(value$data)||
    is.null(input$Learner) ||
    input$Learner == "" ||
    is.null(input$validat) ||
    input$validat == "" ||
    is.null(input$nrtime) ||
    input$nrtime == "" ||
    is.null(input$nrstatus) ||
    input$nrstatus == ""
  ){
    disable("nrbt")
  }
  else{
    enable("nrbt")
  }
})

observeEvent(input$nrbt, {
  disable("nrbt")
})

observeEvent(nestrun(), {
  shinyjs::show("modelcomp_wrapper")
  enable("nrbt")
})

observeEvent(is.null(nestrun()), {
  enable("nrbt")
})


observeEvent(input$nrbt,{
  output$modeldescrip <- renderText(
    {
        paste("Features selected base on", nestrun()$lrnid,":",paste(nestrun()[[2]],collapse = ", "),sep=" ")
    }
  )
})


observeEvent(input$nrbt,{
  output$modelcomp <- renderPlot(
    width = isolate({input$nrwidth}),
    height = isolate({input$nrheight}),
    res = 96,
    {
        df<-as.data.frame(nestrun()$`Benchmark result`)
        df$learner.id<-gsub(".tuned","",df$learner.id)
        ggboxplot(df, "learner.id", "cindex",
                  fill = "learner.id",ylab="C-index",xlab=NULL
        )+theme(axis.title.x = element_blank(),legend.position = 'none')
  }
  )
})


observeEvent(input$nrbt, {
  output$modelcomping  <- renderPlot({
    closeAlert(session, "nestrmess")
      df<-as.data.frame(nestrun()$`Benchmark result`)
      df$learner.id<-gsub(".tuned","",df$learner.id)
      ggboxplot(df, "learner.id", "cindex",
                fill = "learner.id",ylab="C-index",xlab=NULL
      )+theme(axis.title.x = element_blank(),legend.position = 'none')
  })
})

observeEvent(input$nrbt, {
  output$modelcomp <- renderUI({
    plotOutput("modelcomping",
               width = paste0(isolate({input$nrwidth}), "%"),
               height = isolate({input$nrheight}))
  })})


output$downloadmodelcomp <- downloadHandler(
  filename = function() {
    paste0("C-index-comparison",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$modelcompwidthdl}),height = isolate({input$modelcompheightdl})) # open the pdf device

      df<-as.data.frame(nestrun()$`Benchmark result`)
      df$learner.id<-gsub(".tuned","",df$learner.id)
      print(ggboxplot(df, "learner.id", "cindex",
                fill = "learner.id",ylab="C-index",xlab=NULL
      )+theme(axis.title.x = element_blank(),legend.position = 'none'))
    dev.off()
  }
)

observeEvent(input$page_before_nestr, {
  newtab <- switch(input$tabs, "DEG" = "nestr","nestr" = "DEG")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_nestr, {
  newtab <- switch(input$tabs, "bmpm" = "nestr","nestr" = "bmpm")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})




