observe({
  data <- rawdata()$clinical
  updateSelectizeInput(session, 'feature', choices = {
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
                       'table1bygene',
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


percent<-function(x,per){
  n<-length(x)
  half <- (n + 1L) %/% (1/per)
  if(n %% (1/per) == 0) mean(sort(x, partial = half + 0L:1L)[half + 0L:1L])
  else   sort(x, partial = half)[half]
}


drawtable1<-eventReactive(input$table1bt,{

  input$table1bt
  data<-isolate({rawdata()})
  data<-lapply(data,as.data.frame)
  index<-intersect(colnames(data$expres),row.names(data$clinical))
  expres<-data$expres[,index]
  clinical<-data$clinical[index,]
  clinical[clinical == ""]<- NA
  gene<-isolate({input$table1bygene})
  clinical$gene<-as.numeric(expres[gene,])
  grouppercent<-isolate({input$grouppercent})
  if(isolate({input$tbgroupby=="Percentage"})){
  clinical$Group<-ifelse(clinical$gene<=quantile(clinical$gene,grouppercent),paste(gene,"low expression group",sep=" "),paste(gene,"high expression group", sep=" "))
  }else{
    clinical$Group<-ifelse(clinical$gene<=isolate({input$tabgpvalue}),paste(gene,"low expression group",sep=" "),paste(gene,"high expression group", sep=" "))
  }
  if("Age" %in% names(clinical)){
    units(clinical$Age) <- "years"
  }
if("OS.time" %in% names(clinical)){
  if(max(na.omit(clinical$OS.time))>600){
    units(clinical$OS.time)="Day"
    } else{
    units(clinical$OS.time)="Month"

    }
}
  if("RFS.time" %in% names(clinical)){
    if(max(na.omit(clinical$RFS.time))>600){
      units(clinical$RFS.time)="Day"
    } else{
      units(clinical$RFS.time)="Month"
    }
  }

  if("PFS.time" %in% names(clinical)){
    if(max(na.omit(clinical$PFS.time))>600){
      units(clinical$PFS.time)="Day"
    } else{
      units(clinical$PFS.time)="Month"
    }
  }

  if("DFS.time" %in% names(clinical)){
    if(max(na.omit(clinical$DFS.time))>600){
      units(clinical$DFS.time)="Day"
    } else{
      units(clinical$DFS.time)="Month"
    }
  }
  if("DMFS.time" %in% names(clinical)){
    if(max(na.omit(clinical$DMFS.time))>600){
      units(clinical$DMFS.time)="Day"
    } else{
      units(clinical$DMFS.time)="Month"
    }
  }
  if("DRFS.time" %in% names(clinical)){
    if(max(na.omit(clinical$DRFS.time))>600){
      units(clinical$DRFS.time)="Day"
    } else{
      units(clinical$DRFS.time)="Month"
    }
  }

  clinical$Group<-as.factor(clinical$Group)
  clinical$Group<-factor(clinical$Group,levels=c(levels(clinical$Group),"P value"))


  rndr <- function(x, name, ...) {
    if (length(x) == 0) {
      y <- clinical[[name]]
      s <- rep("", length(render.default(x=y, name=name, ...)))
      if (is.numeric(y)) {
        p <- t.test(y ~ clinical$Group)$p.value
      } else {
        p <- chisq.test(table(y, droplevels(clinical$Group)))$p.value
      }
      s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
      s
    } else {
      render.default(x=x, name=name, ...)
    }
  }
  rndr.strat <- function(label, n, ...) {
    ifelse(n==0, label, render.strat.default(label, n, ...))
  }

  feature<-isolate({input$feature})
  feature<-paste(feature,collapse ="+")
  fomu<-as.formula(paste("~",feature,"|Group",sep=''))
  table1(fomu,data=clinical,droplevels=F, render=rndr, render.strat=rndr.strat,overall = F)

    }
)

observe({
  if(is.null(input$table1bygene) || input$table1bygene == "" ||is.null(input$feature) ||input$feature == "" ){
    disable("table1bt")
  }
  else{
    enable("table1bt")
  }
})



observeEvent(input$table1bt , {

  output$table1 <- renderText({

      drawtable1()
  })
})
observeEvent(input$table1bt , {
  shinyjs::show("clincor_wrapper")
})

output$downloadtable1 <- downloadHandler(

  filename = function() {
    paste0("Table1-clinical-correlation-",Sys.Date(),'.html')
  },
  content = function(file) {
    saveWidget(print(drawtable1()), file,selfcontained = T)
  }
)


observeEvent(input$page_before_clinical, {
  newtab <- switch(input$tabs, "clinical" = "nomo","nomo" = "clinical")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_clinical, {
  newtab <- switch(input$tabs, "KM" = "clinical","clinical" = "KM")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

