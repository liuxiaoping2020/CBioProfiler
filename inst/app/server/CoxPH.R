observe({
  data <- rawdata()$clinical
  updateSelectizeInput(session,
                       'CoxPHtime',
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
                       'CoxPHstatus',
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

observe({
  data <- rawdata()$clinical

  updateSelectizeInput(session, 'coxclinvar', choices = {
    if (!is.null(data)) {
      if (class(data) ==  class(data.frame())) choices <- as.character(colnames(data))
      if (class(data) !=  class(data.frame())) choices <- as.character(colnames(as.data.frame(data)))
      choices
    }
  }, server = TRUE)
}
)

observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'coxgene',
                       choices = {
                         if (!is.null(data)) {
                           if (class(data) ==  class(data.frame()))

                             choices <- as.character(make.names(row.names(data)))
                           if (class(data) !=  class(data.frame()))

                             choices <- as.character(make.names(row.names(as.data.frame(data))))
                           choices
                         }
                       },
                       server = TRUE,
                       selected = "TP53")
})


coxPHmodel <- eventReactive(input$CoxPHbt,{
  input$CoxPHbt
  data <- isolate({
    rawdata()
  })

  data <- lapply(data, as.data.frame)
  time <- isolate({input$CoxPHtime})
  status <- isolate({input$CoxPHstatus})
  data$clinical <-data$clinical[complete.cases(data$clinical[, time]) & data$clinical[, time] > 0, ]
  index <- intersect(colnames(data$expres), row.names(data$clinical))
  clinical <- data$clinical[index, ]
  clinical[clinical == ""] <- NA
  expres <- data$expres[, index]
  row.names(expres)<-gsub("-","_",row.names(expres))
  coxclinvar <- isolate({input$coxclinvar})

  coxgene <- isolate({input$coxgene})

  maxtick<-isolate(input$maxtick)

  OS.time<-clinical[,time]
  OS<-clinical[,status]

  clinfeat<-subset(clinical,select=coxclinvar)
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
      "CoxPHmess",
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
      "CoxPHmess",
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
      "CoxPHmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = paste("Variable",paste(name[which(clinfeat==1)]),"you selected has only one level (only one value), which is not accepted by CoxPH model. Please check!"),
      append = T,
      dismiss = T
    )
    return(NULL)
  }else{
  clinical$time<-clinical[,time]

  clinical$status<-clinical[,status]
  if(!is.null(coxgene)){

    feature<-c(coxclinvar,coxgene)
    coxgene<-subset(as.data.frame(t(expres)),select=coxgene)
    clinical<-merge(clinical,coxgene,by=0)

  }else{
    feature<-coxclinvar
    clinical<-clinical
  }

    Surv<-Surv(clinical$time, clinical$status)
    feature1<-paste(feature,collapse ="+")
    fomu<-as.formula(paste("Surv","~",feature1,sep=''))
    fit<-coxph(fomu,data=clinical)
    sumfit<-summary(fit)
  
    a<-cbind(sumfit$conf.int[,1:4],sumfit$coefficients[,5])
    colnames(a)<-c("HR","exp(-coef)","LCI","UCI","P value")
    multicox <- as.data.frame(a)

    unicox <- c()
    for (i in 1:length(feature)) {
      fomu <- as.formula(paste("Surv~", feature[i], sep = ""))
      cox <- coxph(fomu, data = clinical)
      cox <- summary(cox)
      conf.int <- cox$conf.int
      coef1 <- cox$coef
      a <- cbind(conf.int, coef1[, 5])
      colnames(a) <- c("HR", "exp(-coef)", "LCI", "UCI", "P value")
      unicox <- rbind(unicox, a)
    }

    unicox<-as.data.frame(unicox)
    index<-intersect(rownames(multicox),row.names(unicox))

    varname<-isolate({input$coxPHvarname})
    var1 <- strsplit(varname, "|", fixed = T)[[1]]

    if(length(var1)==length(index)){
      unicox<-unicox[index,]
      multicox<-multicox[index,]
      row.names(unicox)<-var1
      row.names(multicox)<-var1
    }else if(varname==""){
      unicox<-unicox[index,]
      multicox<-multicox[index,]
    }else{
      createAlert(
        session,
        "CoxPHmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The length of variable names does not equal to the clinical variables included in the CoxPH model",
        append = T,
        dismiss = T
      )
      return(NULL)
    }
    
    plotcox<-coxforestp(unicox=unicox,multicox=multicox,legend.pos=isolate({input$coxfiglg}) ,xlim=isolate({input$maxtick}),varname=NULL)
    sumcox<-sumfit
    unicox<-round(unicox,3)
    unicox<-as.data.frame(unicox)
    unicox$Pvalue<- format.pval(unicox$`P value`, digits = 3, eps = 0.001)

    multicox<-round(multicox,3)

    multicox$Pvalue<- format.pval(multicox$`P value`, digits = 3, eps = 0.001)
    unicox$summaryx<-paste(unicox$HR,paste("(",unicox$LCI,"-",paste(paste(unicox$UCI,",",sep=""),unicox$Pvalue,sep=" "),")",sep=""),sep=" ")
    multicox$summaryy<-paste(multicox$HR,paste("(",multicox$LCI,"-",paste(paste(multicox$UCI,",",sep=""),multicox$Pvalue,sep=" "),")",sep=""),sep=" ")

    tab<-merge(unicox,multicox,by=0)
    rownames(tab)<-row.names(multicox)

    tablecox<-subset(tab,select=c(summaryx,summaryy))
    names(tablecox)<-c("Univariate HR(95% CI,Pvalue)","Multivariable HR(95% CI,Pvalue)")
    res<-list(tablecox,plotcox,sumcox)
    names(res)<-c("tablecox","plotcox","sumcox")

    res
  }
})



observe({
  if(is.null(input$coxclinvar) || input$coxclinvar == ""||is.null(input$CoxPHtime) || input$CoxPHtime == "" || is.null(input$CoxPHstatus) || input$CoxPHstatus == ""){
    disable("CoxPHbt")
  }
  else{
    enable("CoxPHbt")
  }
})

observeEvent(input$CoxPHbt, {
  output$Coxtable <-  DT::renderDT({
    coxPHmodel()$tablecox
  })
})


observeEvent(input$CoxPHbt, {
  shinyjs::show("mybox_wrapper")
})



observeEvent(input$CoxPHbt, {
  output$Coxforestploting  <- renderPlot({
    coxPHmodel()$plotcox
  })
})

observeEvent(input$CoxPHbt, {
  output$Coxforestplot<- renderUI({
    plotOutput("Coxforestploting",
               width = paste0(isolate({input$Coxwidth}), "%"),
               height = isolate({input$Coxheight}))
  })})


observeEvent(input$CoxPHbt, {
  shinyjs::show("mybox_wrapper.table")
})


observeEvent(input$CoxPHbt, {
  output$Coxsummary <-   renderPrint({
   coxPHmodel()$sumcox
  })
})

output$saveforest <- downloadHandler(

  filename = function() {
    paste0("CoxPH-forestplot-",Sys.Date(),'.pdf')
  },
  content = function(file) {

    pdf(file,width = isolate({input$FPwidth}),height = isolate({input$FPheight})) # open the pdf device
    print(coxPHmodel()$plotcox)
    dev.off()

  }
)

output$saveforesttable <- downloadHandler(
  filename = function() {
    paste0("CoxPH-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(coxPHmodel()$tablecox, file)
  }
)

observeEvent(input$page_before_CoxPH, {
  newtab <- switch(input$tabs, "CoxPH" = "KM","KM" = "CoxPH")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_CoxPH, {
  newtab <- switch(input$tabs, "SurvROC" = "CoxPH","CoxPH" = "SurvROC")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})
