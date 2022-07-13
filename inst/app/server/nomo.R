observeEvent(rawdata(), {
  data <- rawdata()$clinical
  updateSelectInput(session, "nomotime", choices = {
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
  )
})

observeEvent(rawdata(), {
  data <- rawdata()$clinical
  updateSelectInput(session, "nomostatus", choices = {
    if (!is.null(data)) {
      if (class(data) ==  class(data.frame()))
        choices <- as.character(colnames(data))
      if (class(data) !=  class(data.frame()))
        choices <-
          as.character(colnames(as.data.frame(data)))
      choices
    }
  },
  selected = "OS"
  )
})

yearlabel2 <- eventReactive(nestrun(),
  {

  req(input$nomotime)
  nomotime<-isolate({input$nomotime})
  req(nestrun())
  res<-isolate({nestrun()})
  split<-isolate({input$split})
  if(split==F){
    data =res$clinical
  }else{
    data<-res$Trainsurv
  }
  # data <- rawdata()$clinical

  time<-data[,nomotime]
  if(is.numeric(time)==F){
    createAlert(
      session,
      "nomomess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival time you selected is not numeric data, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  } else{

  if(max(na.omit(data[,nomotime]))>600){
    by<- 365
  } else{
    by<-12
  }
  years2<-seq(from=0,to=quantile(na.omit(data[,nomotime]),0.95),by=by)
  years2<-years2[-1]
  yearlabel2<-c()
  for(i in 1:length(years2)){
    yearlabel2[i]<- paste(i,"year",sep="-")
  }
  yearlabel2
  }
})

observeEvent(yearlabel2(), {
  choices <- yearlabel2()
  updateSelectInput(session, "nomoypoint", choices = choices)
})

observeEvent(nestrun(), {
  if(isolate({input$split})==T){
    data<-nestrun()$Trainsurv
  }else{
    data<-nestrun()$clinical
  }
  updateSelectInput(session, "nomovar", choices = {
    if (!is.null(data)) {
      if (class(data) ==  class(data.frame()))
        choices <- as.character(colnames(data))
      if (class(data) !=  class(data.frame()))
        choices <-
          as.character(colnames(as.data.frame(data)))
      choices
    }
  },
  selected = "Risk"
  )
})

observe({
  if(
    is.null(input$nomotime)||
    input$nomotime == "" ||
    is.null(input$nomostatus) ||
    input$nomostatus == "" ||
    is.null(input$nomovar) ||
    input$nomovar == "" ||
    is.null(input$nomoypoint) ||
    input$nomoypoint == ""
  )
  {
    disable("nomogrambt")
  }
  else{
    enable("nomogrambt")
  }
})

nomorun<-eventReactive(input$nomogrambt,{
  input$nomogrambt
  withProgress(
    message = "Constructing the nomogram",
    detail = "It may take a while, please be patient",
    value = 5,{
  split<-isolate({input$split})
  res<-isolate({nestrun()})
  if(split==T){
    train<-res$Trainsurv
  }else{
    train<-res$clinical
  }
  variable <-isolate({input$nomovar})
  varlabels<-isolate({input$nomovarlab})
  yearpoint<-isolate({input$nomoypoint})
  time<-isolate({input$nomotime})
  status<-isolate({input$nomostatus})

  clinfeat<-subset(train,select=variable)
  name<-names(clinfeat)
  fc<-function(x){
    x<-factor(x)
    x<-length(levels(x))
    return(x)
  }
  clinfeature<-lapply(clinfeat,fc)
  OS.time<-train[,time]
  OS<-train[,status]
  if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
    createAlert(
      session,
      "nomomess",
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
      "nomomess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival status column is not numerically in 1, and 0, please check!",
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
  var1 <- strsplit(varlabels, "|", fixed = T)[[1]]
  if(varlabels=="" || length(var1)==length(variable)){
    if(varlabels==""){
      varlabels<-NULL
    }else{
      varlabels <- strsplit(varlabels, "|", fixed = T)[[1]]
    }
    nomo(
      train = train,
      time = time,
      status = status,
      variable = variable,
      varlabels = varlabels,
      yearpoint = yearpoint
    )

  } else{
    createAlert(
      session,
      "nomomess",
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
})

observeEvent(nomorun(), {
  shinyjs::show("nomogram_wrapper")
  shinyjs::show("invalcal_wrapper")
})

observeEvent(input$nomogrambt,{
  updateCollapse(session, "collapsenomo", open = "Nomogram", close = "Descriptions and parameters")
})

observeEvent(input$nomogrambt,{
  output$nomoploting  <- renderPlot({
    print(plot(nomorun()))
  })
})

observeEvent(input$nomogrambt, {
  output$nomoplot <- renderUI({
    plotOutput("nomoploting",
               width = paste0(isolate({input$nomowidth}), "%"),
               height = isolate({input$nomoheight}))
  })})

output$savenomogram <- downloadHandler(
  filename = function(){
    paste0("Nomogram-",Sys.Date(),'.pdf')
  },
  content = function(file){
    pdf(file,width = isolate({input$nomogramwidthdl}),height = isolate({input$nomogramheightdl})) # open the pdf device
    print(plot(nomorun()))
    dev.off()
  }
)

yearlabel3 <- eventReactive(nestrun(),{
  req(input$nomotime)
  nomotime<-isolate({input$nomotime})
  split<-isolate({input$split})
  if(split==F){
    data <- nestrun()$clinical
  }else{
    data<-nestrun()$Testsurv
  }
  
  if(max(na.omit(data[,nomotime]))>600){
    by<- 365
  } else{
    by<-12
  }
  years3<-seq(from=0,to=quantile(na.omit(data[,nomotime]),0.95),by=by)
  years3<-years3[-1]
  yearlabel3<-c()
  for(i in 1:length(years3)){
    yearlabel3[i]<- paste(i,"year",sep="-")
  }
  yearlabel3
})

observeEvent(yearlabel3(), {
  choices <- yearlabel3()[-length(yearlabel3())]
  updateSelectInput(session, "invalidypoint", choices = choices)
})


invalidrun<-eventReactive(input$invalidbt,{
  input$invalidbt
  withProgress(
    message = "Internally validating and calibrating the CoxPH based nomogram with bootstrap",
    detail = "It may take a while, please be patient",
    value = 5,{
      split<-isolate({input$split})
      res<-isolate({nestrun()})
        internalvc(split=split,
                   time=isolate({input$nomotime}),
                   status=isolate({input$nomostatus}),
                   variable=isolate({input$nomovar}),
                   ratio=isolate({input$inratio}),
                   bootstrap=isolate({input$inreps}),
                   yearpoint=isolate({input$invalidypoint}),
                   res=res)
      })
})


observe({
  if(
    is.na(input$invalidypoint)||
    length(input$invalidypoint)==0||
    is.null(input$invalidypoint)||
    input$invalidypoint==""||

    is.na(input$inreps)||
    length(input$inreps)==0||
    is.null(input$inreps)||
    input$inreps == "" ||

    is.na(input$inratio)||
    length(input$inratio)==0||
    is.null(input$inratio) ||
    input$inratio == ""||
    input$inratio >1 ||
    input$inratio <0
  )
  {
    disable("invalidbt")
  }
  else{
    enable("invalidbt")
  }
})


observeEvent(invalidrun(), {
  shinyjs::show("invalid_wrapper")
  shinyjs::show("incalib_wrapper")
  shinyjs::show("exvalcal_wrapper")

})

observeEvent(input$invalidbt, {
  output$invaliddescrip <- renderText({
    if(isolate({input$split==T})){
      paste("Training set",invalidrun()$idex[[1]],"; Test set",invalidrun()$idex[[2]])
    }else{
      paste(invalidrun()$idex[[1]])
    }
  })
})

observeEvent(input$invalidbt,{
  updateCollapse(session, "collapsevalnomo", open = "Internal validation and calibration", close = "Descriptions and parameters")
})

observeEvent(input$invalidbt,{
  output$invaliding  <- renderPlot({
    print(invalidrun()$p)
  })
})

observeEvent(input$invalidbt, {
  output$invalid <- renderUI({
    plotOutput("invaliding",
               width = paste0(isolate({input$invalidwidth}), "%"),
               height = isolate({input$invalidheight}))
  })})

output$saveinvalid <- downloadHandler(
  filename = function(){
    paste0("Internal-validation-",Sys.Date(),'.pdf')
  },
  content = function(file){
    pdf(file,width = isolate({input$invalidwidthdl}),height = isolate({input$invalidheightdl})) # open the pdf device
    print(invalidrun()$p)
    dev.off()
  }
)

observeEvent(input$invalidbt,{
  output$incalibing  <- renderPlot({
    if(isolate({input$split})==T){
    calibtrain<-invalidrun()$validation$validtrain$calibration
    calibtest<-invalidrun()$validation$validtest$calibration
    yearpoint<-isolate({input$invalidypoint})
    calitration<-append(calibtrain,calibtest)
    yp<-c(yearpoint,yearpoint)
    par(mfrow = c(2, length(yearpoint)))
    for(i in 1:length(calitration)){
      plot(calitration[[i]],xlab=paste("Predicted", yp[i], "Survival Probability"),
           ylab=paste("Actual", yp[i], "survival Probability"),  subtitles=isolate({input$incalibsubtitle}))
      my.label <- paste(LETTERS[i])
      mtext(my.label, side = 3, adj = 0,font = 2,cex = 1,line = 1)
    }
    }else{
      calibration<-result1$validation$validtrain$calibration
      yearpoint<-isolate({input$invalidypoint})
      par(mfrow = c(1, length(yearpoint)))
      for(i in 1:length(calitration)){
        plot(calitration[[i]],xlab=paste("Predicted", yearpoint[i], "Survival Probability"),
             ylab=paste("Actual", yearpoint[i], "survival Probability"),  subtitles=isolate({input$incalibsubtitle}))
        my.label <- paste(LETTERS[i])
        mtext(my.label, side = 3, adj = 0,font = 2,cex = 1,line = 1)
      }
    }
  })
})

observeEvent(input$invalidbt, {
  output$incalib<- renderUI({
    plotOutput("incalibing",
               width = paste0(isolate({input$incalibwidth}), "%"),
               height = isolate({input$incalibheight}))
  })})

output$saveincalib <- downloadHandler(
  filename = function(){
    paste0("Internal-calibration-",Sys.Date(),'.pdf')
  },
  content = function(file){
    pdf(file,width = isolate({input$incalibwidthdl}),height = isolate({input$incalibheightdl})) # open the pdf device
    if(isolate({input$split})==T){
      calibtrain<-invalidrun()$validation$validtrain$calibration
      calibtest<-invalidrun()$validation$validtest$calibration
      yearpoint<-isolate({input$invalidypoint})
      calitration<-append(calibtrain,calibtest)
      yp<-c(yearpoint,yearpoint)
      par(mfrow = c(2, length(yearpoint)))
      for(i in 1:length(calitration)){
        plot(calitration[[i]],xlab=paste("Predicted", yp[i], "Survival Probability"),
             ylab=paste("Actual", yp[i], "survival Probability"),  subtitles=isolate({input$incalibsubtitle}))
        my.label <- paste(LETTERS[i])
        mtext(my.label, side = 3, adj = 0,font = 2,cex = 1,line = 1)
      }
    }else{
      calibration<-result1$validation$validtrain$calibration
      yearpoint<-isolate({input$invalidypoint})
      par(mfrow = c(1, length(yearpoint)))
      for(i in 1:length(calitration)){
        plot(calitration[[i]],xlab=paste("Predicted", yearpoint[i], "Survival Probability"),
             ylab=paste("Actual", yearpoint[i], "survival Probability"),  subtitles=isolate({input$incalibsubtitle}))
        my.label <- paste(LETTERS[i])
        mtext(my.label, side = 3, adj = 0,font = 2,cex = 1,line = 1)
      }
    }
    dev.off()
  }
)

observeEvent(rawdata2(), {
  data <- rawdata2()$clinical
  updateSelectInput(session, "extvaltime", choices = {
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
  )
})

observeEvent(rawdata2(), {
  data <- rawdata2()$clinical
  updateSelectInput(session, "extvalstatus", choices = {
    if (!is.null(data)) {
      if (class(data) ==  class(data.frame()))
        choices <- as.character(colnames(data))
      if (class(data) !=  class(data.frame()))
        choices <-
          as.character(colnames(as.data.frame(data)))
      choices
    }
  },
  selected = "OS"
  )
})

observeEvent(rawdata2(), {
  data <- rawdata2()$clinical
  updateSelectInput(session, "extvalvar", choices = {
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
####################################################################################################################

# yearlabel4 <- eventReactive(rawdata2(),{
yearlabel4 <- reactive({
  
  req(input$extvaltime)
  nomotime<-isolate({input$extvaltime})
  data <- rawdata2()$clinical

  time<-data[,nomotime]
  if(is.numeric(time)==F){
    createAlert(
      session,
      "exnomomess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival time you selected is not numeric data, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  } else{
  if(max(na.omit(data[,nomotime]))>600){
    by<- 365
  } else{
    by<-12
  }
  years4<-seq(from=0,to=quantile(na.omit(data[,nomotime]),0.95),by=by)
  years4<-years4[-1]
  yearlabel4<-c()
  for(i in 1:length(years4)){
    yearlabel4[i]<- paste(i,"year",sep="-")
  }
  yearlabel4
  }
})

observeEvent(yearlabel4(), {
  choices <- yearlabel4()[-length(yearlabel4())]
  updateSelectInput(session, "exvalidypoint", choices = choices)
})


extvcrun<-eventReactive(input$exvalidbt,{
  input$exvalidbt
  data<- isolate({rawdata2()})
  time=isolate({input$extvaltime})
  status=isolate({input$extvalstatus})
  variable = isolate({input$extvalvar})

  clinical<-data$clinical
  OS.time<-clinical[,time]
  OS<-clinical[,status]

  clinfeat<-subset(clinical,select=variable)
  name<-names(clinfeat)
  fc<-function(x){
    x<-factor(x)
    x<-length(levels(x))
    return(x)
  }
  clinfeature<-lapply(clinfeat,fc)

  if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
    createAlert(
      session,
      "exnomomess",
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
      "exnomomess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival status column is not numerically in 1, and 0, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  } else if(1 %in% clinfeature){
    createAlert(
      session,
      "exnomomess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = paste("Variable",paste(name[which(clinfeat==1)]),"you selected has only one level (only one value), which is not accepted by CoxPH model. Please check!"),
      append = T,
      dismiss = T
    )
    return(NULL)
  } else{

  withProgress(
    message = "Externally validating and calibrating the CoxPH based nomogram with bootstrap",
    detail = "It may take a while, please be patient",
    value = 5,{
      res<-isolate({nestrun()})
      extvalid(
          extdata = isolate({rawdata2()}),
          variable = variable,
          res = res,
          ratio = isolate({input$exratio}),
          boostrap = isolate({input$exreps}),
          time=time,
          status=status,
          yearpoint=isolate({input$exvalidypoint})
        )
    })

}
  })

observe({
  if(
    is.na(input$extvaltime)||
    is.null(input$extvaltime)||
    input$extvaltime == "" ||

    is.null(input$extvalstatus)||
    is.na(input$extvalstatus)||
    input$extvalstatus== "" ||

    is.na(input$extvalvar)||
    is.null(input$extvalvar)||
    input$extvalvar== "" ||

    is.na(input$exvalidypoint) ||
    is.null(input$exvalidypoint) ||
    input$exvalidypoint == "" ||

    is.null(input$exreps) ||
    is.na(input$exreps) ||
    input$exreps == "" ||
    input$exratio > 1 ||
    input$exratio< 0||
    is.null(input$exratio) ||
    is.na(input$exratio) ||
    input$exratio == ""
  )
  {
    disable("exvalidbt")
  }
  else{
    enable("exvalidbt")
  }
})

observeEvent(input$exvalidbt, {
  disable("exvalidbt")
})

observeEvent(input$exvalidbt, {
    updateCollapse(session, "collapseextnomo", open = "External validation and calibration", close = "Descriptions and parameters")
})


observeEvent(extvcrun(), {
  shinyjs::show("exvalid_wrapper")
  shinyjs::show("excalib_wrapper")
  enable("exvalidbt")
})


observeEvent(input$exvalidbt, {
  output$exvaliddescrip <- renderText({
      paste(extvcrun()$idex[[1]])
  })
})

observeEvent(input$exvalidbt,{
  output$exvaliding  <- renderPlot({
    print(extvcrun()$p)
  })
})

observeEvent(input$exvalidbt, {
  output$exvalid <- renderUI({
    plotOutput("exvaliding",
               width = paste0(isolate({input$exvalidwidth}), "%"),
               height = isolate({input$exvalidheight}))
  })})

output$saveexvalid <- downloadHandler(
  filename = function(){
    paste0("External-validation-",Sys.Date(),'.pdf')
  },
  content = function(file){
    pdf(file,width = isolate({input$exvalidwidthdl}),height = isolate({input$exvalidheightdl})) # open the pdf device
    print(extvcrun()$p)
    dev.off()
  }
)

observeEvent(input$exvalidbt,{
  output$extcalibing  <- renderPlot({

    calibration<-extvcrun()$validation$validexternal$calibration
    yearpoint<-isolate({input$exvalidypoint})
    par(mfrow = c(1, length(yearpoint)))
    for(i in 1:length(calibration)){
      plot(calibration[[i]],xlab=paste("Predicted", yearpoint[i], "Survival Probability"),
           ylab=paste("Actual", yearpoint[i], "survival Probability"),  subtitles=isolate({input$excalibsubtitle}))
      my.label <- paste(LETTERS[i])
      mtext(my.label, side = 3, adj = 0,font = 2,cex = 1,line = 1)
    }
  })
})

observeEvent(input$exvalidbt, {
  output$extcalib<- renderUI({
    plotOutput("extcalibing",
               width = paste0(isolate({input$excalibwidth}), "%"),
               height = isolate({input$excalibheight}))
  })})

output$saveexcalib <- downloadHandler(
  filename = function(){
    paste0("External-calibration-",Sys.Date(),'.pdf')
  },
  content = function(file){
    pdf(file,width = isolate({input$excalibwidthdl}),height = isolate({input$excalibheightdl})) # open the pdf device

    calibration<-extvcrun()$validation$validexternal$calibration
    yearpoint<-isolate({input$exvalidypoint})
    par(mfrow = c(1, length(yearpoint)))
    for(i in 1:length(calibration)){
      plot(calibration[[i]],xlab=paste("Predicted", yearpoint[i], "Survival Probability"),
           ylab=paste("Actual", yearpoint[i], "survival Probability"),  subtitles=isolate({input$excalibsubtitle}))
      my.label <- paste(LETTERS[i])
      mtext(my.label, side = 3, adj = 0,font = 2,cex = 1,line = 1)
    }
    dev.off()
  }
)

observeEvent(input$page_before_nomo, {
  newtab <- switch(input$tabs, "valmo" = "nomo","nomo" = "valmo")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_nomo, {
  newtab <- switch(input$tabs, "nomo" = "clinical","clinical" = "nomo")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})





