# observe({
#   data <- rawdata()$clinical
#   updateSelectizeInput(session, 'feature', choices = {
#     if (!is.null(data)) {
#       if (class(data) ==  class(data.frame())) choices <- as.character(colnames(data))
#       if (class(data) !=  class(data.frame())) choices <- as.character(colnames(as.data.frame(data)))
#       choices
#     }
#   }, server = TRUE)
# }
# )
# 
# observe({
#   data <- rawdata()$expres
#   updateSelectizeInput(session,
#                        'table1bygene',
#                        choices = {
#                          if (!is.null(data)) {
#                            if (class(data) ==  class(data.frame()))
#                              choices <- as.character(row.names(data))
#                            if (class(data) !=  class(data.frame()))
#                              choices <-
#                                as.character(row.names(as.data.frame(data)))
#                            choices
#                          }
#                        },
#                        server = TRUE,
#                        selected = "TP53")
# })
# 
# 
# percent<-function(x,per){
#   n<-length(x)
#   half <- (n + 1L) %/% (1/per)
#   if(n %% (1/per) == 0) mean(sort(x, partial = half + 0L:1L)[half + 0L:1L])
#   else   sort(x, partial = half)[half]
# }

cbafeat <- reactive({
  req(input$cbaclinicalfeaturecohort)
  req(validst())
  res<-isolate({validst()})
  cohort<-isolate({input$cbaclinicalfeaturecohort})
  if(cohort=="Training set"){
    clinical<-res$traindata$clinical
  }else{
    clinical<-res$validata$clinical
  }
  colnames(clinical)
})

observeEvent(cbafeat(), {
  choices <- cbafeat()
  updateSelectInput(session, "cbaclinicalfeature", choices = choices)
})


cbaclin<-eventReactive(input$cbaclinicaltable1bt,{

  require(table1)
  input$cbaclinicaltable1bt
  res<-isolate({validst()})
  cohort<-isolate({input$cbaclinicalfeaturecohort})
  if(cohort=="Training set"){
    clinical<-res$traindata$clinical
  }else{
    clinical<-res$validata$clinical
  }

  clinical[clinical == ""]<- NA

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
  # print(clinical$Group)
  # print(names(clinical))
  clinical$Group<-as.factor(clinical$Subtype)
  clinical$Group<-factor(clinical$Group,levels=c(levels(clinical$Group),"P value"))
  
  rndr <- function(x, name, ...) {
    if (length(x) == 0) {
      y <- clinical[[name]]
      s <- rep("", length(render.default(x=y, name=name, ...)))
      if (is.numeric(y)) {
        if(length(levels(factor(clinical$Group)))==2){
          p <- t.test(y ~ clinical$Group)$p.value
        }else{
          aov <- aov(y~factor(clinical$Group))
          aov<-summary(aov)
          aov<-aov[[1]]
          p<-aov[1,5]
        }
        
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
  
  feature<-isolate({input$cbaclinicalfeature})
  feature<-paste(feature,collapse ="+")
  fomu<-as.formula(paste("~",feature,"|Group",sep=''))
  table1(fomu,data=clinical,droplevels=F, render=rndr, render.strat=rndr.strat,overall = F)
}
)


observe({
  if(is.null(input$cbaclinicalfeature) ||input$cbaclinicalfeature == "" ||is.null(isolate({validst()}))){
    disable("cbaclinicaltable1bt")
  }
  else{
    enable("cbaclinicaltable1bt")
  }
})

observeEvent(input$cbaclinicaltable1bt, {
  
  updateCollapse(session, "collapsecbaclinical", open = "Table for the correlations between cancer subtype and clinical features", close = "Descriptions and parameters")
  
})

observeEvent(input$cbaclinicaltable1bt , {
  
  output$cbaclinical1 <- renderText({
    
    cbaclin()
  })
})

# observeEvent(input$table1bt , {
#   shinyjs::show("clincor_wrapper")
# })

# output$downloadtable1 <- downloadHandler(
#   
#   filename = function() {
#     paste0("Table1-clinical-correlation-",Sys.Date(),'.html')
#   },
#   content = function(file) {
#     saveWidget(print(drawtable1()), file,selfcontained = T)
#   }
# )


observeEvent(input$page_before_cbaclinical, {
  newtab <- switch(input$tabs, "cbaclinical" = "cba","cba" = "cbaclinical")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbaclinical, {
  newtab <- switch(input$tabs, "cbaKM" = "cbaclinical","cbaclinical" = "cbaKM")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})
