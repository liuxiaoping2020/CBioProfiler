
value <- reactiveValues(data = NULL)
observeEvent(input$cbavalidationbt, {
  value$data <- isolate({validst()})
})

cbaDEGana <- eventReactive(input$cbaDEGbt,{
  require(limma)
  input$cbaDEGbt
  # data <- isolate({ rawdata()})
  res<-isolate({validst()})
  cohort<-isolate({input$cbaDEGcor})
  
  if(cohort=="Training set"){
   data<-res$traindata
  }else{
   data<-res$validata
  }
  
  DEGmethod <- isolate({input$cbaDEGmethod})
  DEGFDR <- isolate({input$cbaDEGFDR})
  DEGFC <- isolate({input$cbaDEGFC})
  # DEGfactor <- isolate({input$cbaDEGfactor})
  
  DEGpadj <- isolate({input$cbaDEGpadj})
  mlDEGsel<-isolate({input$cbamlDEGsel})
  mlgeneselp<-isolate({input$cbamlgeneselp})
  mlgeneselogFC<-isolate({input$cbamlgeneselogFC})
  mlgeneselp1<-isolate({input$cbamlgeneselp1})
  mlgeneselogFC1<-isolate(input$cbamlgeneselogFC1)
  
  clinical <- data$clinical
  expres <- data$expres
  
  if("" %in% row.names(expres)){
    expres<-expres[-which(row.names(expres)==""),]
  }
  
  clinical[clinical == ""] <- NA
  clinical<-clinical[complete.cases(clinical$Subtype),]
  
  # if(is.numeric(clinical[,DEGfactor])==F){
  #   clinical[,DEGfactor]<-gsub(' ','',clinical[,DEGfactor])
  #   clinical[,DEGfactor]<-gsub('-','_',clinical[,DEGfactor])
  # }
  clinical$Subtype<-make.names(clinical$Subtype)
  
  index <- intersect(names(expres), row.names(clinical))
  clinical <- clinical[index, ]
  expres <- expres[, index]
  
  expres<-as.data.frame(rmz(expres,per=0.9))
  
  group <- factor(clinical$Subtype)
  
  if(length(levels(group))==1 || length(levels(group))>10){
    createAlert(
      session,
      "cbaDEGmsg",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The levels group varibale you selected equals to 1 or greater than 10, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  }else{
    withProgress(message = 'Performing differential expression analysis using Limma',
                 detail = 'This may take a while...',
                 value = 3,
                 {
                   require(limma)
                   expres <- as.matrix(expres)
                   quant <- as.numeric(quantile(expres, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
                   
                   val <- (quant[5] > 100) ||
                     (quant[6] - quant[1] > 50 &&
                        quant[2] > 0)
                   
                   if (val) {
                     expres[which(expres <= 0)] <- NaN
                     expres <- log2(expres)
                   }
                   
                   expres <- as.data.frame(expres)
                   design <- model.matrix( ~ group + 0, expres)
                   colnames(design) <- levels(group)
                   fit <- lmFit(expres, design)
                   Grp <- levels(group)
                   
                   if (length(Grp) == 2) {
                     comp <- paste(Grp[1], Grp[2], sep = "-")
                   } else if (length(Grp) > 2) {
                     comp <- paste(Grp, c(tail(Grp,-1), head(Grp, 1)), sep = "-")
                   }
                   
                   cont.matrix <-makeContrasts(contrasts = comp, levels = design)
                   fit2 <- contrasts.fit(fit, cont.matrix)
                   fit2 <- eBayes(fit2)
                   if(length(Grp) == 2){
                     DEG <- topTable(fit2, adjust.method = DEGpadj,number =(nrow(expres)))
                     names(DEG) <- c("LogFC", "AveExpr", "t", "PValue", "AdjustedP", "B")
                   }else if (length(Grp) > 2){
                     DEG <- topTable(fit2, adjust.method = DEGpadj,number =(nrow(expres)))
                     DEG1<-subset(DEG,select=c("AveExpr"  ,       "F"  ,    "P.Value" ,   "adj.P.Val"))
                     names(DEG1)<-c("AveExpr", "F" , "PValue" ,"AdjustedP")
                     DEG1$LogFC<-as.numeric(DEG[,1])
                     DEG<-DEG1
                   }
                   res<-list(DEG=DEG, group=group, expres=expres,clinical=clinical)
                 })
    
    return(res)
  }
})

observeEvent(input$cbaDEGbt, {
  output$cbaDEGtable <-  DT::renderDT({
    signif(cbaDEGana()$DEG, 3)
  })
})
observeEvent(input$cbaDEGbt, {
  updateCollapse(session, "cbacollapseDEG", open = "Differentially expressed gene table", close = "Descriptions and parameters for differentially expressed gene analysis")
})

observe({
  if(is.null(value$data)){
    disable("cbaDEGbt")
  }
  else{
    enable("cbaDEGbt")
  }
})


observeEvent(input$cbaDEGbt, {
  disable("cbaDEGbt")
})

observeEvent(cbaDEGana(), {
  shinyjs::show("cbaDEG_wrapper")
  enable("cbaDEGbt")
  shinyjs::show("cbaDEG1_wrapper")
})

output$cbadownloadDEGtable <- downloadHandler(
  filename = function() {
    paste0("DEG-analysis-result-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(cbaDEGana()$DEG, file)
  }
)


cbaDEGop <- eventReactive(input$cbaDEGvisbt,{
  input$cbaDEGvisbt
  DEG<-isolate({cbaDEGana()$DEG})
  # browser()
  mlDEGsel<-isolate({input$cbamlDEGsel})
  mlgeneselp<-isolate({input$cbamlgeneselp})
  mlgeneselogFC<-isolate({input$cbamlgeneselogFC})
  mlgeneselp1<-isolate({input$cbamlgeneselp1})
  mlgeneselogFC1<-isolate(input$cbamlgeneselogFC1)
  if(mlDEGsel=="Adjusted P"){
    sigDEG<-DEG[DEG$AdjustedP<mlgeneselp,]
  } else if(mlDEGsel=="LogFC"){
    sigDEG<-DEG[abs(DEG$LogFC)>mlgeneselogFC,]
  } else{
    sigDEG<-DEG[DEG$AdjustedP<mlgeneselp1 & abs(DEG$LogFC)>mlgeneselogFC1,]
  }
  sigDEG<-sigDEG[complete.cases(sigDEG$AveExpr),]
  if(nrow(sigDEG)==0){
    createAlert(
      session,
      "cbaDEGvismess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = paste(
        "No genes meet the DEG significance cutoff you specified"),
      append = T,
      dismiss=T
    )
    return(NULL)
  } else{
    list(Output=sigDEG, DEG=DEG, group=isolate({cbaDEGana()$group}), expres=isolate({cbaDEGana()$expres}),clinical=isolate({cbaDEGana()$clinical}))
    }
 }
)

observeEvent(input$cbaDEGvisbt, {
  updateCollapse(session, "cbacollapseDEGvis", open = "DEG output", close = "Descriptions and parameters for DEG output")
})

observeEvent(input$cbaDEGvisbt, {
  output$cbaDEGog <-  DT::renderDT({
    signif(cbaDEGop()$Output, 3)
  })
})

observeEvent(cbaDEGop(), {
  shinyjs::show("cbaDEGog_wrapper")
  shinyjs::show("cbaDEGvis_wrapper")
  shinyjs::show("cbaDEG2_wrapper")
})
output$cbadownloadopgene <- downloadHandler(
  filename = function() {
    paste0("Significantly-different-expression-gene-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(cbaDEGop()$Output, file)
  }
)


cbaheatmap<- eventReactive(input$cbaDEGvisbt,{
  input$cbaDEGvisbt
  withProgress(message = "Drawing heatmap",
               detail = "This may take a while...",
               value = 3,
               {
                 DEG<-isolate({cbaDEGop()$DEG})
                 group<-isolate({cbaDEGop()$group})
                 expres<-isolate({cbaDEGop()$expres})
                 sigDEG<-isolate({cbaDEGop()$Output})
                 # print(2)
                 color <- c()
                 level <- levels(group)
                 htm(
                   DEG=DEG,
                   sigDEG=sigDEG,
                   group=group,
                   color=color,
                   DEGscale=isolate({input$cbaDEGscale}),
                   expres=expres,
                   coltitsize=isolate({input$cbacoltitsize}),
                   ClusterR=isolate({input$cbaClusterR}),
                   cluster_row_slices=isolate({input$cbacluster_row_slices}),
                   cludistanrow=isolate({input$cbacludistanrow}),
                   clumethodrow=isolate({input$cbaclumethodrow}),
                   Rdend_side=isolate({input$cbaRdend_side}),
                   showRname=isolate({input$cbashowRname}),
                   Rnameside=isolate({input$cbaRnameside}),
                   showFDR=isolate({input$cbashowFDR}),
                   showFC=isolate({input$cbashowFC}),
                   ClusterC=isolate({input$cbaClusterC}),
                   heatname=isolate({input$cbaheatname}),
                   cluster_column_slices=isolate({input$cbacluster_column_slices}),
                   cludistancol=isolate({input$cbacludistancol}),
                   clumethodcol=isolate({input$cbaclumethodcol}),
                   Cdend_side=isolate({input$cbaCdend_side}),
                   showCname=isolate({input$cbashowCname}),
                   Cnameside=isolate({input$cbaCnameside}),
                   heatColors=isolate({input$cbaheatColors})
                 )
               })
})


cbavocanaplt<-eventReactive(input$cbaDEGvisbt,{
  input$cbaDEGvisbt
  DEG<-isolate({cbaDEGana()$DEG})
  
  withProgress(message = "Drawing volcano plot",
               detail = "This may take a while...",
               value = 3,
               {
                 voca(
                   DEG = DEG,
                   vocaXlab = isolate({input$cbavocaXlab}),
                   vocaYlab = isolate({input$cbavocaYlab}),
                   vocacutp = isolate({input$cbavocacutp}),
                   vocacutfc = isolate({input$cbavocacutfc}),
                   legendPosition = isolate({input$cbalegendPosition})
                 )
               })
  
})


cbaMAplot<-eventReactive(input$cbaDEGvisbt,{
  input$cbaDEGvisbt
  DEG<-isolate({cbaDEGana()$DEG})
  withProgress(message = "Drawing MA plot",
               detail = "This may take a while...",
               value = 3,
               {
                 MAplt(
                   DEG=DEG,
                   DEGmethod=isolate({input$cbaDEGmethod}),
                   MAstopmeth=isolate({input$cbaMAstopmeth}),
                   MAcutp=isolate({input$cbaMAcutp}),
                   MAfc=isolate({input$cbaMAfc}),
                   Topgene=isolate({input$cbaTopgene}),
                   MAgenesym=isolate({input$cbaMAgenesym}),
                   MAXlab=isolate({input$cbaMAXlab}),
                   MAYlab=isolate({input$cbaMAYlab}),
                   MAlegendPosition=isolate({input$cbaMAlegendPosition})
                 )
               })
})

cbaadjplot<-eventReactive(input$cbaDEGvisbt,{
  input$cbaDEGvisbt
  DEG<-isolate({cbaDEGana()$DEG})
  col<-isolate({input$cbapadjcol})
  withProgress(message = "Drawing adjusted P plot",
               detail = "This may take a while...",
               value = 3,
               {
                 padjplt(DEG,col)
               })
}
)

observeEvent(cbaDEGana(), {
  data <- cbaDEGana()[[1]]
  updateSelectizeInput(
    session,
    "cbaMAgenesym",
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
    selected = "TP53"
  )
})

observeEvent(input$cbaDEGvisbt, {
  output$cbaDEGvising  <- renderPlot({
    closeAlert(session, "cbaDEGvismess")
    if(isolate({input$cbaDEGvismeth})=="Heatmap"){
      print(cbaheatmap())
      
    } else if(isolate({input$cbaDEGvismeth})=="Volcano plot"){
      print(cbavocanaplt())
      
    } else if(isolate({input$cbaDEGvismeth})=="MA plot"){
      cbaMAplot()
      
    } else {
      cbaadjplot()
    }
  })
})

observeEvent(input$cbaDEGvisbt, {
  # updateCollapse(session, "collapsesurvivalplot", open = "Survival plot")
  output$cbaDEGvis<- renderUI({
    plotOutput("cbaDEGvising",
               width = paste0(isolate({input$cbapadjwidth}), "%"),
               height = isolate({input$cbapadjheight}))
  })})


output$cbadownloadDEGvis <- downloadHandler(
  filename = function() {
    paste0("DEG-visualization-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cbaDEGviswidthdl}),height = isolate({input$cbaDEGvisheightdl})) # open the pdf device
    if(isolate({input$cbaDEGvismeth})=="Heatmap"){
      print(cbaheatmap())
      
    } else if(isolate({input$cbaDEGvismeth})=="Volcano plot"){
      print(cbavocanaplt())
      
    } else if(isolate({input$cbaDEGvismeth})=="MA plot"){
      print(cbaMAplot())
      
    } else {
      print(cbaadjplot())
    }
    dev.off()
  }
)

observeEvent(input$page_before_cbaDEG, {
  newtab <- switch(input$tabs, "cbaDEG" = "cbaSurvROC","cbaSurvROC" = "cbaDEG")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbaDEG, {
  newtab <- switch(input$tabs, "cbaDEG" = "cbaimmune","cbaimmune" = "cbaDEG")
  updateTabItems(session, "ta  bs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})
