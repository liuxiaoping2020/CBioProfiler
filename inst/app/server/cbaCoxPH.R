
cbakmfeat1 <- reactive({
  req(input$cbaCoxphcor)
  req(validst())
  res<-isolate({validst()})
  cohort<-isolate({input$cbaCoxphcor})
  if(cohort=="Training set"){
    clinical<-res$traindata$clinical
  }else{
    clinical<-res$validata$clinical
  }
  colnames(clinical)
})



observeEvent(cbakmfeat1(), {
  choices <- cbakmfeat1()
  updateSelectInput(session, "cbaCoxPHtime", choices = choices
                    ,
                    selected = "OS.time"
  )
}
)

observeEvent(cbakmfeat1(), {
  choices <- cbakmfeat1()
  updateSelectInput(session, "cbaCoxPHstatus", choices = choices
                    ,
                    selected = "OS"
  )
}
)


observeEvent(cbakmfeat1(), {
  choices <- cbakmfeat1()
  updateSelectInput(session, "cbacoxclinvar", choices = choices
                    ,
                    selected = "Age"
  )
}
)


cbacoxPHmodel <- eventReactive(input$cbaCoxPHbt,{
  input$cbaCoxPHbt
  res<-isolate({validst()})
  cohort<-isolate({input$cbaCoxphcor})
  
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
  
  
  # data <- lapply(data, as.data.frame)
  # time <- isolate({input$CoxPHtime})
  # status <- isolate({input$CoxPHstatus})
  # data$clinical <-data$clinical[complete.cases(data$clinical[, time]) & data$clinical[, time] > 0, ]
  # index <- intersect(colnames(data$expres), row.names(data$clinical))
  # clinical <- data$clinical[index, ]
  # clinical[clinical == ""] <- NA
  # expres <- data$expres[, index]
  # row.names(expres)<-gsub("-","_",row.names(expres))
  coxclinvar <- isolate({input$cbacoxclinvar})
  
  # coxgene <- isolate({input$coxgene})
  
  maxtick<-isolate(input$cbamaxtick)
  
  # OS.time<-clinical[,time]
  # OS<-clinical[,status]
  
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
      "cbaCoxPHmess",
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
      "cbaCoxPHmess",
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
      "cbaCoxPHmess",
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
    # if(!is.null(coxgene)){
    #   
    #   feature<-c(coxclinvar,coxgene)
    #   coxgene<-subset(as.data.frame(t(expres)),select=coxgene)
    #   clinical<-merge(clinical,coxgene,by=0)
    #   
    # }else{
      feature<-coxclinvar
      clinical<-clinical
    # }
    require(survival)
    clinical$Subtype<-factor(gsub("Subtype ",'',clinical$Subtype,fixed = T))
    # print(clinical$Subtype)
    Surv<-Surv(clinical$time, clinical$status)
    feature<-c("Subtype",feature)
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
    
    varname<-isolate({input$cbacoxPHvarname})
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
        "cbaCoxPHmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The length of variable names does not equal to the clinical variables included in the CoxPH model",
        append = T,
        dismiss = T
      )
      return(NULL)
    }
    
    plotcox<-coxforestp(unicox=unicox,multicox=multicox,legend.pos=isolate({input$cbacoxfiglg}) ,xlim=isolate({input$cbamaxtick}),varname=NULL)
    sumcox<-sumfit
    unicox<-round(unicox,3)
    unicox<-as.data.frame(unicox)
    unicox$Pvalue<- format.pval(unicox$`P value`, digits = 3, eps = 0.001)
    
    multicox<-round(multicox,3)
    
    multicox$Pvalue<- format.pval(multicox$`P value`, digits = 3, eps = 0.001)
    unicox$summaryx<-paste(unicox$HR,paste("(",unicox$LCI,"-",paste(paste(unicox$UCI,",",sep=""),unicox$Pvalue,sep=" "),")",sep=""),sep=" ")
    multicox$summaryy<-paste(multicox$HR,paste("(",multicox$LCI,"-",paste(paste(multicox$UCI,",",sep=""),multicox$Pvalue,sep=" "),")",sep=""),sep=" ")

    tab<-merge(unicox,multicox,by=0)
    rownames(tab)<-tab$Row.names

    tablecox<-subset(tab,select=c(summaryx,summaryy))
    names(tablecox)<-c("Univariate HR(95% CI,Pvalue)","Multivariable HR(95% CI,Pvalue)")
    res<-list(tablecox,plotcox,sumcox)
    names(res)<-c("tablecox","plotcox","sumcox")
    
    res
  }
})



observe({
  if(is.null(input$cbacoxclinvar) || input$cbacoxclinvar == ""||is.null(input$cbaCoxPHtime) || input$cbaCoxPHtime == "" || is.null(input$cbaCoxPHstatus) || input$cbaCoxPHstatus == ""){
    disable("cbaCoxPHbt")
  }
  else{
    enable("cbaCoxPHbt")
  }
})


observeEvent(input$cbaCoxPHbt, {
  updateCollapse(session, "cbacollapseCoxPH", open = "Cox proportional hazard regression model", close = "Descriptions and parameters")
})


observeEvent(input$cbaCoxPHbt, {
  output$cbaCoxtable <-  DT::renderDT({
    cbacoxPHmodel()$tablecox
  })
})


observeEvent(input$cbaCoxPHbt, {
  shinyjs::show("cbamybox_wrapper")
})


observeEvent(input$cbaCoxPHbt, {
  output$cbaCoxforestploting  <- renderPlot({
    cbacoxPHmodel()$plotcox
  })
})

observeEvent(input$cbaCoxPHbt, {
  output$cbaCoxforestplot<- renderUI({
    plotOutput("cbaCoxforestploting",
               width = paste0(isolate({input$cbaCoxwidth}), "%"),
               height = isolate({input$cbaCoxheight}))
  })})


observeEvent(input$cbaCoxPHbt, {
  shinyjs::show("cbamybox_wrapper.table")
})


observeEvent(input$cbaCoxPHbt, {
  output$cbaCoxsummary <-   renderPrint({
    cbacoxPHmodel()$sumcox
  })
})

output$cbasaveforest <- downloadHandler(
  
  filename = function() {
    paste0("CoxPH-forestplot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    
    pdf(file,width = isolate({input$cbaFPwidth}),height = isolate({input$cbaFPheight})) # open the pdf device
    print(cbacoxPHmodel()$plotcox)
    dev.off()
    
  }
)

output$cbasaveforesttable <- downloadHandler(
  filename = function() {
    paste0("CoxPH-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(cbacoxPHmodel()$tablecox, file)
  }
)

observeEvent(input$page_before_cbaCoxPH, {
  newtab <- switch(input$tabs, "cbaCoxPH" = "cbaKM","cbaKM" = "cbaCoxPH")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbaCoxPH, {
  newtab <- switch(input$tabs, "cbaSurvROC" = "cbaCoxPH","cbaCoxPH" = "cbaSurvROC")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})
