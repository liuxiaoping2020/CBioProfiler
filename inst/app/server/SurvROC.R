
observeEvent(rawdata(), {
  data <- rawdata()$clinical
  updateSelectInput(session, "SurvROCtime", choices = {
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
  updateSelectInput(session, "SurvROCstatus", choices = {
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

observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'SurvROCgene',
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

yearlabelx <- reactive({
  req(input$SurvROCtime)
  SurvROCtime<-isolate({input$SurvROCtime})
  data <- rawdata()$clinical
  if(SurvROCtime %in% names(data)){
  time<-data[,SurvROCtime]

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

observeEvent(yearlabelx(), {
  choices <- yearlabelx()
  updateSelectInput(session, "predictyear", choices = choices)
})


SurvROCplot<-reactive({
  input$SurvROCbt
  data <- isolate({
    rawdata()
  })
  SurvROCgene<-isolate({input$SurvROCgene})
  SurvROCtime<-isolate({input$SurvROCtime})

  SurvROCstatus<-isolate({input$SurvROCstatus})
  method<-isolate({input$SurvROCmethod})
  clinical<- data$clinical

  OS.time<-clinical[,SurvROCtime]
  OS<-clinical[,SurvROCstatus]

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
  }else{

  expres<- data$expres
  clinical<-clinical[complete.cases(clinical[,SurvROCtime])& clinical[,SurvROCtime]>0,]
  index<-intersect(row.names(clinical),colnames(expres))
  clinical<-clinical[index,]
  expres<-expres[,index]
  clinical$SurvROCgene<-as.numeric(expres[SurvROCgene,])
  withProgress(
    message = "Constructing the prediction model",
    detail = "It may take a while, please be patient",
    value = 5,{
  survROC(
    data =clinical,
    bmtime = SurvROCtime,
    bmstatus = SurvROCstatus,
    marker = "SurvROCgene",
    method = method,
    predyear = isolate({input$predictyear}),
    cutpoint=isolate({input$cutoff})
  )
    })
  
  # if(max(clinical[,SurvROCtime])>600){
  #   by<- 365
  # } else{
  #   by<-12
  # }
  # years<-seq(from=0,to=quantile(na.omit(clinical[,SurvROCtime]),0.95),by=by)
  # years<-years[-1]
  # yearlabel<-c()
  # for(i in 1:length(years)){
  #   yearlabel[i]<- paste(i,"year",sep="-")
  # }
  # 
  # sumROC<-list()
  # for (i in 1:length(years)){
  #   sumROC[[i]] <- survivalROC(Stime = clinical[,SurvROCtime],status = clinical[,SurvROCstatus],marker = clinical$SurvROCgene,
  #                              predict.time =years[i],method = method,span = 0.25*nrow(clinical)^(-0.20))
  # }
  # sumAUC<-list()
  # for (i in 1:length(sumROC)){
  #   sumAUC[[i]]<-sumROC[[i]]$AUC
  # }
  # 
  # ROCdata<-c()
  # for(i in 1:length(sumROC)){
  #   predict.time<-sumROC[[i]]$predict.time
  #   TP<-sumROC[[1]]$TP
  #   FP<-sumROC[[i]]$FP
  #   auc<-sumROC[[i]]$AUC
  #   tmp<-c(predict.time,TP,FP,auc)
  #   ROCdata<-rbind(ROCdata,tmp)
  # }
  # 
  # survivalROC_helper <- function(t) {
  #   survivalROC(Stime = clinical[,SurvROCtime],status = clinical[,SurvROCstatus],marker = clinical$SurvROCgene,
  #               predict.time =t,method = method,span = 0.25*nrow(clinical)^(-0.20))
  # }
  # 
  # timeponit<-isolate({input$predictyear})
  # yearid<-as.numeric(substr(timeponit,1,1))
  # time<-years[yearid]
  # 
  # survivalROC_data <- data_frame(t = time) %>%
  #   mutate(survivalROC = map(t, survivalROC_helper),
  #          auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
  #          df_survivalROC = map(survivalROC, function(obj) {
  #            as_data_frame(obj[c("cut.values","TP","FP")])
  #          })) %>%
  #   dplyr::select(-survivalROC) %>%
  #   unnest() %>%
  #   arrange(t, FP, TP)
  # 
  # survivalROC_data1<-mutate(survivalROC_data,auc =sprintf("%.3f",auc))
  # survivalROC_data1$years<-survivalROC_data1$t/by
  # survivalROC_data1<-unite(survivalROC_data1,year, years,auc,sep = " year AUC: " )
  # AUC =factor(survivalROC_data1$year)
  # sumAUC1<-list()
  # for(i in 1:length(yearid)){
  #   sumAUC1[[i]]<-sumAUC[[yearid[i]]]
  # }
  # 
  # sumROC1<-list()
  # for(i in 1:length(yearid)){
  #   sumROC1[[i]]<-sumROC[[yearid[i]]]
  # }
  # 
  # 
  # ROC.1<-sumROC1[[which.max(sumAUC1)]]
  # 
  # dot <- data.frame(TP = ROC.1$TP[with(ROC.1, which.min((1-TP)^2+ FP^2))],
  #                   FP = ROC.1$FP[with(ROC.1, which.min((1-TP)^2+ FP^2))])
  # dot <- rbind(c(1,0),dot)
  # 
  # if(isolate({input$cutoff})==T){
  #   cutoff.imm <- ROC.1$cut.values[with(ROC.1, which.min((1-TP)^2+ FP^2))]
  #   ROC.plot<-ggplot(survivalROC_data1,mapping = aes(x = FP, y = TP)) +
  #     geom_path(aes(color= AUC))+
  #     geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  #     theme_bw() +scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+xlab("False positive rate")+
  #     ylab("True positive rate") +
  #     theme(legend.position = c(0.7,0.2))+
  #     geom_path(mapping = aes(x = FP,y = TP),data = dot)+
  #     annotate("text",x = dot$FP[2] + 0.05,y = dot$TP[2],label = paste0("Cutoff: ",round(cutoff.imm,3)))
  # } else{
  #   ROC.plot<-ggplot(survivalROC_data1,mapping = aes(x = FP, y = TP)) +
  #     geom_path(aes(color= AUC))+
  #     geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  #     theme_bw() +scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+xlab("False positive rate")+
  #     ylab("True positive rate") +
  #     theme(legend.position = c(0.7,0.2))
  # }
  # ROC.plot
}}

)

observe({
  if (is.null(input$SurvROCgene) ||
      input$SurvROCgene == "" ||
      is.null(input$SurvROCtime) ||
      input$SurvROCtime == "" ||
      is.null(input$SurvROCstatus) ||
      input$SurvROCstatus == "" ||
      is.null(input$predictyear) || input$predictyear == "") {
    disable("SurvROCbt")
  }
  else{
    enable("SurvROCbt")
  }
})

observeEvent(input$SurvROCbt, {
  output$SurvROCplotting <- renderPlot({
    SurvROCplot()
  })
})

observeEvent(input$SurvROCbt, {
  output$SurvROCplot <- renderUI({
    plotOutput("SurvROCplotting",
               width = paste0(isolate({input$SurvROCwidth}), "%"),
               height = isolate({input$SurvROCheight}))
  })})

observeEvent(input$SurvROCbt, {
  shinyjs::show("survROC_wrapper")
})

output$downloadsurvROC <- downloadHandler(

  filename = function() {
    paste0("SurvivalROC-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$survROCwidth.dl}),height = isolate({input$survROCheight.dl})) # open the pdf device
    print(SurvROCplot())
    dev.off()
   }
)

observeEvent(input$page_before_SurvROC, {
  newtab <- switch(input$tabs, "CoxPH" = "SurvROC","SurvROC" = "CoxPH")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_SurvROC, {
  newtab <- switch(input$tabs, "SurvROC" = "mcorgene","mcorgene" = "SurvROC")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

