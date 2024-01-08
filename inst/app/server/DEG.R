observeEvent(rawdata(), {
  data <- rawdata()$clinical
 
  updateSelectizeInput(
    session,
    "DEGfactor",
    choices = {
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
    selected = "Sex"
  )
})

DEGana <- eventReactive(input$DEGbt,{
  require(limma)
  input$DEGbt
  data <- isolate({ rawdata()})
  DEGmethod <- isolate({input$DEGmethod})

  DEGFDR <- isolate({input$DEGFDR})
  DEGFC <- isolate({input$DEGFC})
  DEGfactor <- isolate({input$DEGfactor})

  DEGpadj <- isolate({input$DEGpadj})
  mlDEGsel<-isolate({input$mlDEGsel})
  mlgeneselp<-isolate({input$mlgeneselp})
  mlgeneselogFC<-isolate({input$mlgeneselogFC})
  mlgeneselp1<-isolate({input$mlgeneselp1})
  mlgeneselogFC1<-isolate(input$mlgeneselogFC1)

  clinical <- data$clinical
  expres <- data$expres

  if("" %in% row.names(expres)){
    expres<-expres[-which(row.names(expres)==""),]
  }

  clinical[clinical == ""] <- NA
  clinical<-clinical[complete.cases(clinical[,DEGfactor]),]

  if(is.numeric(clinical[,DEGfactor])==F){
    clinical[,DEGfactor]<-gsub(' ','',clinical[,DEGfactor])
    clinical[,DEGfactor]<-gsub('-','_',clinical[,DEGfactor])
  }

  index <- intersect(names(expres), row.names(clinical))
  clinical <- clinical[index, ]
  expres <- expres[, index]

  expres<-as.data.frame(rmz(expres,per=0.9))

  group <- factor(clinical[, DEGfactor])

if(is.numeric(clinical[,DEGfactor])){
  createAlert(
    session,
    "DEGmsg",
    "exampleAlert",
    title = "Please note!",
    style =  "danger",
    content = "The group varibale you selected is not categorical variable, please check!",
    append = T,
    dismiss=T
  )
  return(NULL)
}else if(length(levels(group))==1 || length(levels(group))>10){

    createAlert(
      session,
      "DEGmsg",
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

observeEvent(input$DEGbt, {
  output$DEGtable <-  DT::renderDT({
    signif(DEGana()$DEG, 3)
  })
})
observeEvent(input$DEGbt, {
  updateCollapse(session, "collapseDEG", open = "Differentially expressed gene table", close = "Descriptions and parameters for differentially expressed gene analysis")
})

observe({
  if(is.null(input$DEGfactor) || input$DEGfactor == ""){
    disable("DEGbt")
  }
  else{
    enable("DEGbt")
  }
})

observeEvent(input$DEGbt, {
  disable("DEGbt")
})

observeEvent(DEGana(), {
  shinyjs::show("DEG_wrapper")
  enable("DEGbt")
  shinyjs::show("DEG1_wrapper")
})

output$downloadDEGtable <- downloadHandler(
  filename = function() {
    paste0("DEG-analysis-result-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(DEGana()$DEG, file)
  }
)


DEGop <- eventReactive(input$DEGvisbt,{
  input$DEGvisbt
  DEG<-isolate({DEGana()$DEG})

  mlDEGsel<-isolate({input$mlDEGsel})
  mlgeneselp<-isolate({input$mlgeneselp})
  mlgeneselogFC<-isolate({input$mlgeneselogFC})
  mlgeneselp1<-isolate({input$mlgeneselp1})
  mlgeneselogFC1<-isolate(input$mlgeneselogFC1)
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
      "DEGvismess",
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
    list(Output=sigDEG, DEG=DEG, group=isolate({DEGana()$group}), expres=isolate({DEGana()$expres}),clinical=isolate({DEGana()$clinical}))

  }
})

observeEvent(input$DEGvisbt, {
  updateCollapse(session, "collapseDEGvis", open = "DEG output", close = "Descriptions and parameters for DEG output")
})

observeEvent(input$DEGvisbt, {
  output$DEGog <-  DT::renderDT({
    signif(DEGop()$Output, 3)
  })
})

observeEvent(DEGop(), {
  shinyjs::show("DEGog_wrapper")

  shinyjs::show("DEGvis_wrapper")
  shinyjs::show("DEG2_wrapper")
})
output$downloadopgene <- downloadHandler(
  filename = function() {
    paste0("Significantly-different-expression-gene-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(DEGop()$Output, file)
  }
)


heatmap<- eventReactive(input$DEGvisbt,{
  input$DEGvisbt
  withProgress(message = "Drawing heatmap",
                                detail = "This may take a while...",
                                value = 3,
                                {
  DEG<-isolate({DEGop()$DEG})
  group<-isolate({DEGop()$group})
  expres<-isolate({DEGop()$expres})
  sigDEG<-isolate({DEGop()$Output})

  color <- c()
  level <- levels(group)
      htm(
        DEG=DEG,
        sigDEG=sigDEG,
        group=group,
        color=color,
        DEGscale=isolate({input$DEGscale}),
        expres=expres,
        coltitsize=isolate({input$coltitsize}),
        ClusterR=isolate({input$ClusterR}),
        cluster_row_slices=isolate({input$cluster_row_slices}),
        cludistanrow=isolate({input$cludistanrow}),
        clumethodrow=isolate({input$clumethodrow}),
        Rdend_side=isolate({input$Rdend_side}),
        showRname=isolate({input$showRname}),
        Rnameside=isolate({input$Rnameside}),
        showFDR=isolate({input$showFDR}),
        showFC=isolate({input$showFC}),
        ClusterC=isolate({input$ClusterC}),
        heatname=isolate({input$heatname}),
        cluster_column_slices=isolate({input$cluster_column_slices}),
        cludistancol=isolate({input$cludistancol}),
        clumethodcol=isolate({input$clumethodcol}),
        Cdend_side=isolate({input$Cdend_side}),
        showCname=isolate({input$showCname}),
        Cnameside=isolate({input$Cnameside}),
        heatColors=isolate({input$heatColors})
      )
                   })
})


vocanaplt<-eventReactive(input$DEGvisbt,{
  input$DEGvisbt
  DEG<-isolate({DEGana()$DEG})

  withProgress(message = "Drawing volcano plot",
               detail = "This may take a while...",
               value = 3,
               {

  voca(
    DEG = DEG,
    vocaXlab = isolate({input$vocaXlab}),
    vocaYlab = isolate({input$vocaYlab}),
    vocacutp = isolate({input$vocacutp}),
    vocacutfc = isolate({input$vocacutfc}),
    legendPosition = isolate({input$legendPosition})
  )
               })

})


MAplot<-eventReactive(input$DEGvisbt,{
  input$DEGvisbt

  DEG<-isolate({DEGana()$DEG})
  withProgress(message = "Drawing MA plot",
               detail = "This may take a while...",
               value = 3,
               {
  MAplt(
    DEG=DEG,
    DEGmethod=isolate({input$DEGmethod}),
    MAstopmeth=isolate({input$MAstopmeth}),
    MAcutp=isolate({input$MAcutp}),
    MAfc=isolate({input$MAfc}),
    Topgene=isolate({input$Topgene}),
    MAgenesym=isolate({input$MAgenesym}),
    MAXlab=isolate({input$MAXlab}),
    MAYlab=isolate({input$MAYlab}),
    MAlegendPosition=isolate({input$MAlegendPosition})
     )
  })
})

adjplot<-eventReactive(input$DEGvisbt,{
  input$DEGvisbt
  DEG<-isolate({DEGana()$DEG})
  col<-isolate({input$padjcol})
  withProgress(message = "Drawing adjusted P plot",
               detail = "This may take a while...",
               value = 3,
               {
  padjplt(DEG,col)
               })
  }
)

observeEvent(DEGana(), {
  data <- DEGana()[[1]]
  updateSelectizeInput(
    session,
    "MAgenesym",
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




observeEvent(input$DEGvisbt, {
  output$DEGvising  <- renderPlot({
     closeAlert(session, "DEGvismess")
    if(isolate({input$DEGvismeth})=="Heatmap"){
      print(heatmap())

    } else if(isolate({input$DEGvismeth})=="Volcano plot"){
      print(vocanaplt())

    } else if(isolate({input$DEGvismeth})=="MA plot"){
      MAplot()

    } else {
      adjplot()
    }
  })
})

observeEvent(input$DEGvisbt, {
  # updateCollapse(session, "collapsesurvivalplot", open = "Survival plot")
  output$DEGvis<- renderUI({
    plotOutput("DEGvising",
               width = paste0(isolate({input$padjwidth}), "%"),
               height = isolate({input$padjheight}))
  })})


output$downloadDEGvis <- downloadHandler(

  filename = function() {
    paste0("DEG-visualization-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$DEGviswidthdl}),height = isolate({input$DEGvisheightdl})) # open the pdf device
    if(isolate({input$DEGvismeth})=="Heatmap"){
      print(heatmap())

    } else if(isolate({input$DEGvismeth})=="Volcano plot"){
      print(vocanaplt())

    } else if(isolate({input$DEGvismeth})=="MA plot"){
      print(MAplot())

    } else {
      print(adjplot())
    }
    dev.off()
  }
)

observeEvent(input$page_before_DEG, {
  newtab <- switch(input$tabs, "DEG" = "msurv","msurv" = "DEG")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_DEG, {
  newtab <- switch(input$tabs, "DEG" = "nestr","nestr" = "DEG")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})
