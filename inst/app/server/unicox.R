observeEvent(rawdata(),{
  data <- rawdata()$clinical

  updateSelectizeInput(session, "msurvtime", choices = {
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
  selected = "OS.time"
  )
 }
)

observeEvent(rawdata(),{
  data <- rawdata()$clinical
  updateSelectizeInput(session, "msurvstatus", choices = {
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
  server = TRUE,
  selected = "OS"
   )
  }
 )

msurvfun<-eventReactive(input$msurvbt,{
  input$msurvbt
  withProgress(message = "Performing Univariate Cox proportional hazards analysis",
               detail = "This may take a while...",
               value =5,
               {
  data <- isolate({
  rawdata()})
  data <- lapply(data, as.data.frame)
  expres<-data$expres
  clinical<-data$clinical
  OS.time<-isolate({input$msurvtime})
  OS<-isolate({input$msurvstatus})

  # method<-isolate({input$msurvcor})
  clinical<-clinical[complete.cases(clinical[,OS.time])&clinical[,OS.time]>0,]
  index<-intersect(row.names(clinical),names(expres))
  clinical<-clinical[index,]
  expres<-expres[,index]
  OS.time<-clinical[,OS.time]
  OS<-clinical[,OS]
  if(mode(OS.time)!="numeric" || mode(OS)!="numeric" ){
    createAlert(
      session,
      "msurvmess",
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
      "msurvmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival status column is not numerically in 1, and 0, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  }else{
  res<-unicox(OS.time,OS,expres)
  res<-as.data.frame(res)
  res<-na.omit(res)
  # res$P.adjusted<-p.adjust(res$PValue, method = method)
  list(res=res,clinical=clinical,expres=expres)
  }
               })
})

observe({
  if(is.null(input$msurvtime) || input$msurvtime == ""||is.null(input$msurvstatus) || input$msurvstatus == "" ){
    disable("msurvbt")
  }
  else{
    enable("msurvbt")
  }
})

observeEvent(input$msurvbt, {
  updateCollapse(session, "collapsemsurv", open = "Univariate CoxPH table", close = "Descriptions and parameters for Univariate CoxPH analysis")
})

observeEvent(input$msurvbt, {
  output$unicoxtable <-  DT::renderDT(
    {
      DT::datatable(msurvfun()$res
      ,options=list(scrollX=TRUE))
      }
  )
})

observeEvent(input$msurvbt, {
  disable("msurvtime")
  disable("msurvstatus")
  disable("msurvcor")
  disable("msurvbt")
  disable("msurvopbt")
})
observeEvent(is.null(msurvfun()), {
  # shinyjs::show("unicoxanalysis_wrapper")
  enable("msurvtime")
  enable("msurvstatus")
  enable("msurvcor")
  enable("msurvbt")
  enable("msurvopbt")
  # shinyjs::show("unicoxtable_wrapper")
})
observeEvent(msurvfun(), {
  shinyjs::show("unicoxanalysis_wrapper")
  enable("msurvtime")
  enable("msurvstatus")
  enable("msurvcor")
  enable("msurvbt")
  enable("msurvopbt")
  shinyjs::show("unicoxtable_wrapper")
})



output$downloadunicoxtable <- downloadHandler(
  filename = function() {
    paste0("Univariate-CoxPH-analysis-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(msurvfun()$res, file)
  }
)
#########
msurvsig<-eventReactive(input$msurvopbt,{
  input$msurvopbt
  res<-isolate({msurvfun()$res})
  method<-isolate({input$msurvcor})
  res$P.adjusted<-p.adjust(res$PValue, method = method)
  if(isolate({input$msurvsel})=="P value"){
    sigres<-res[res$PValue<isolate({input$msurvp}),]
  } else if(isolate({input$msurvsel})=="Adjusted P value"){
    sigres<-res[res$P.adjusted<isolate({input$msurvap}),]
  }else{
    # sigres<-res[res$P.adjusted<isolate({input$msurvap}),]
    res<-res[order(res$PValue,decreasing = FALSE),]
    sigres<-res[1:isolate({input$msurvtp}),,drop=F]

  }

  # sigres$P.adjusted<-NULL

  if(nrow(sigres)==0){
    createAlert(
      session,
      "msurvmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = paste(
        "No gene meet the significance cutoff you specified"),
      append = T,
      dismiss=T
    )
    list(res = res,
         clinical = isolate({
           msurvfun()$clinical
         }),
         expres = isolate({
           msurvfun()$expres
         }))
  } else{

    list(
      res = res,
      clinical = isolate({
        msurvfun()$clinical
      }),
      expres = isolate({
        msurvfun()$expres
      }),
      Output = sigres
    )
  }
})
observeEvent(input$msurvopbt, {
  updateCollapse(session, "collapsemsurvsig", open = "Significant univariate CoxPH result", close = "Descriptions and parameters")
})

observeEvent(input$msurvopbt, {
  output$sigcoxtable <-  DT::renderDT(
    {
      DT::datatable(msurvsig()$Output,
                    options=list(scrollX=TRUE))
      }

  )
})

output$downloadsigcoxtable <- downloadHandler(
  filename = function() {
    paste0("significant-Univariate-CoxPH-analysis-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(msurvsig()$Output, file)
  }
)

forestp<-eventReactive(input$msurvopbt,{
  input$msurvopbt
  withProgress(message = "Drawing forestplot",
               detail = "This may take a while...",
               value =5,
               {
  res<-isolate({msurvsig()$Output})
  unicoxforestp(res,xlim = isolate(input$msurvxtick))
  })
})

observeEvent(input$msurvopbt, {
  disable("msurvopbt")
  disable("msurvbt")
})

observeEvent(msurvsig(), {
  shinyjs::show("unicoxForestplot_wrapper")
  enable("msurvopbt")
  enable("msurvbt")
  shinyjs::show("sigcoxtable_wrapper")
})




observeEvent(input$msurvopbt, {
  output$sigforing  <- renderPlot({
     closeAlert(session, "msurvmess")
    forestp()
  })
})

observeEvent(input$msurvopbt, {
  output$sigfor <- renderUI({
    plotOutput("sigforing",
               width = paste0(isolate({input$msurvwidth}), "%"),
               height = isolate({input$msurvheight}))
  })})


output$downloadunicoxForestplot <- downloadHandler(

  filename = function() {
    paste0("significant-univariate-forestplot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$unicoxFPwidthdl}),height = isolate({input$unicoxFPheightdl})) # open the pdf device
    print(forestp())
    dev.off()
  }
)


observeEvent(input$page_before_unicox, {
  newtab <- switch(input$tabs, "wgcna" = "msurv","msurv" = "wgcna")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_unicox, {
  newtab <- switch(input$tabs, "DEG" = "msurv","msurv" = "DEG")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

