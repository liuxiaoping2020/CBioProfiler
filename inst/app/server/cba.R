observe({
  data <- rawdata()$clinical
  updateSelectizeInput(session,
                       'cbapreproctime',
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
                       'cbapreprocstatus',
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

cbaprerun<-eventReactive(input$cbeprebt,{
  input$cbeprebt
  data <- isolate({
    rawdata()
  })
  require(survival)
  data <- lapply(data, as.data.frame)
  time <- isolate({ input$cbapreproctime})
  status <- isolate({input$cbapreprocstatus})
  fsm<-isolate({input$fsm})
  FScutmethod<-isolate({input$FScutmethod})
  FSvaluet<-isolate({input$FSvaluet})
  FSvaluec<-isolate({input$FSvaluec})
  PC_percent<-isolate({input$PC_percent})
  PCAscale<-isolate({input$PCAscale})
  Coxphcutoff<-isolate({input$Coxphcutoff})
  exp<-data$expres
  clin<-data$clinical
  OS.time<-clin[,time]
  OS<-clin[,status]
  if(is.numeric(OS.time)==F || is.numeric(OS)==F){
    createAlert(
      session,
      "cbamess",
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
      "cbamess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The survival status column is not numerically in 1, and 0, please check!",
      append = T,
      dismiss=T
    )
    return(NULL)
  }else{
  withProgress(message = "Preprocessing the input data",
               detail = "This may take a while...",
               value = 3,
               {
  res<-SBpreproc(
    df=exp,
    fsm=fsm,
    FScutmethod=FScutmethod,
    FSvaluet=FSvaluet,
    PC_percent=PC_percent, 
    scale=PCAscale, 
    clin=clin,
    time=time,
    status=status,
    cutoff=Coxphcutoff
  )
               })
  }

  list(expres=res,clinical=clin,data=data)
  }
)

observe({
  if(is.null(input$cbapreproctime) || input$cbapreproctime== "" ||is.null(input$cbapreprocstatus) ||input$cbapreprocstatus == ""){
    disable("cbeprebt")
  }else{
    enable("cbeprebt")
  }
})

observeEvent(input$cbeprebt, {
  disable("cbeprebt")
})

observeEvent(input$cbeprebt, {
  updateCollapse(session, "collapsecba", open = "Data preprocessing result", close = "Descriptions and parameters")
})

observeEvent(cbaprerun(), {
  shinyjs::show("cbapreproc_wrapper")
  # shinyjs::show("biomarkout_wrapper.table")
  enable("cbeprebt")
})

observeEvent(is.null(cbaprerun()), {
  shinyjs::enable("cbeprebt")
  shinyjs::show("cba_wrapper")
})

observeEvent(input$cbeprebt, {

  output$cbapreproc <-  DT::renderDT({
    DT::datatable(cbaprerun()$expres,options=list(scrollX=TRUE))
    
  })

})

output$savecbapreproc <- downloadHandler(
  filename = function() {
    paste0("Preprocessed-gene-expression-data",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(cbaprerun()$expres, file)
  }
)

# observeEvent(cbaprerun(), {
#   shinyjs::show("sampletree_wrapper")
#   enable("setopvar")
#   enable("topvar")
#   enable("ZK")
#   enable("wgvariable")
#   shinyjs::show("wgcna1_wrapper")
# })


cbarun<-eventReactive(input$cbabt,{
  input$cbabt
  data <- isolate({
    cbaprerun()
  })
  # browser()
  clustermethod<-isolate({input$clustermethod})
  clusterAlg<-isolate({input$clusterAlg})
  ccseed<-isolate({input$ccseed})
  maxK<-isolate({input$maxK})
  pItem<-isolate({input$pItem})
  ccreps<-isolate({input$ccreps})
  pFeature<-isolate({input$pFeature})
  corUse<-isolate({input$corUse})
  m3citers<-isolate({input$m3citers})
  ref_method<-isolate({input$ref_method})
  repsref<-isolate({input$repsref})
  repsreal<-isolate({input$repsreal})
  pacx1<-isolate({input$pacx1})
  pacx2<-isolate({input$pacx2})
  objective<-isolate({input$objective})
  M3Cmethod<-isolate({input$M3Cmethod})
  tunelambda<-isolate({input$tunelambda})
  lambdadefault<-isolate({input$lambdadefault})
  distance<-isolate({input$distance})
  
  expres<-data$expres
  clinical<-data$clinical
  withProgress(
    message = "Starting clustering analysis",
    detail = "According to your parameter settings, it may take a long time, please be patient.",
    value = 5,{
  clusterun(clustermethod=clustermethod,
                 data=expres,
                 maxK=maxK,
                 pItem=pItem,
                 clusterAlg=clusterAlg,
                 seed=ccseed,
                 reps=ccreps,
                 pFeature=pFeature,
                 distance=distance,
                 corUse=corUse,
                 iters=m3citers,
                 ref_method=ref_method,
                 repsref=repsref,
                 repsreal=repsreal,
                 pacx1=pacx1,
                 pacx2=pacx2,
                 objective=objective,
                 # silent,
                 method=M3Cmethod,
                 lambdadefault=lambdadefault,
                 tunelambda=tunelambda 
  )
  })
  
})

observeEvent(input$cbabt, {
  disable("cbabt")
})

observeEvent(input$cbabt, {
  updateCollapse(session, "collapsecba1", open = "Subtype identification", close = "Descriptions and parameters")
})

observeEvent(cbarun(), {
  shinyjs::show("cmp_wrapper")
  shinyjs::show("cdf_wrapper")
  shinyjs::show("silhouette_wrapper")
  shinyjs::show("cbaassign_wrapper")
  shinyjs::show("cbavalidation_wrapper")
  enable("cbabt")
})

observeEvent(is.null(cbarun()), {
  enable("cbabt")
})

observeEvent(cbarun(), {
  output$cmping <-  renderPlot(
    
    {
      res<-cbarun()
      clustermethod<-isolate({input$clustermethod})
      if(clustermethod =="M3C"){
        K<-length(unique(res$assignments))
        fm<-res$realdataresults[[K]]$consensus_matrix
        finalLinkage<-"average"
        hc = hclust(as.dist(1 - fm), method = finalLinkage)
        ct = cutree(hc, K)
        names(ct) = colnames(data)
        if (any(class(data) == "dist")) {
          names(ct) = colnames(as.matrix(data))
        }
        colorList = list()
        colorM = rbind()
        thisPal <- c(
          "#A6CEE3",
          "#1F78B4",
          "#B2DF8A",
          "#33A02C",
          "#FB9A99",
          "#E31A1C",
          "#FDBF6F",
          "#FF7F00",
          "#CAB2D6",
          "#6A3D9A",
          "#FFFF99",
          "#B15928",
          "#bd18ea",
          "#2ef4ca",
          "#f4cced",
          "#f4cc03",
          "#05188a",
          "#e5a25a",
          "#06f106",
          "#85848f",
          "#000000",
          "#076f25",
          "#93cd7f",
          "#4d0776",
          "#ffffff"
        )
        colBreaks = NA
        tmyPal=NULL
        if(is.null(tmyPal) == TRUE) {
          colBreaks = 10
          tmyPal = myPal(colBreaks)
        } else {
          colBreaks = length(tmyPal)
        }
        pc = fm
        pc = pc[hc$order,]
        pc = rbind(pc, 0)
        colorList = setClusterColors(res[[K - 1]][[3]], ct, thisPal, colorList)
        stats::heatmap(pc,
                Colv = as.dendrogram(hc), 
                Rowv = NA, 
                symm = FALSE, 
                scale = "none",
                col = tmyPal,
                na.rm = TRUE, 
                labRow = F,
                labCol = F, 
                mar = c(5,5),
                main = paste("Consensus matrix k=", K, sep = ""),
                ColSideColors=colorList[[1]]
        ) 
        legend("topright", legend = unique(ct), fill = unique(colorList[[1]]),
               horiz = FALSE)  
      } else{
        pc<-res$res$pc
        ct<-res$res$consensusClass
        hc<-res$res$consensusTree
        tmyPal<-res$res$tmyPal
        clrs<-res$res$clrs
        sc=res$res$sc
        K<-length(unique(ct))
         stats::heatmap(x=pc,
                Colv = as.dendrogram(hc),
                Rowv = NA,
                symm = FALSE,
                scale = "none",
                col = tmyPal,
                na.rm = TRUE,
                labRow = F,
                labCol = F,
                mar = c(5,5),
                main = paste("Consensus matrix k=", K, sep = ""),
                ColSideColors=clrs[[1]]
        )
        legend("topright", legend = unique(ct), fill = unique(clrs[[1]]), horiz = FALSE)
      }
      }
  )
})

observeEvent(input$cbabt, {
  output$cmp <- renderUI({
    plotOutput("cmping",
               width =   paste0(isolate({input$comwidth}),"%"),
               height = isolate({input$comheight}))
  })})

output$savecmp <- downloadHandler(
  filename = function() {
    paste0("Consensus-matrix-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cmpwidth}),height = isolate({input$cmpheight})) # open the pdf device
    res<-cbarun()
    clustermethod<-isolate({input$clustermethod})
    if(clustermethod =="M3C"){
      K<-length(unique(res$assignments))
      fm<-res$realdataresults[[K]]$consensus_matrix
      finalLinkage<-"average"
      hc = hclust(as.dist(1 - fm), method = finalLinkage)
      ct = cutree(hc, K)
      names(ct) = colnames(data)
      if (any(class(data) == "dist")) {
        names(ct) = colnames(as.matrix(data))
      }
      colorList = list()
      colorM = rbind()
      thisPal <- c(
        "#A6CEE3",
        "#1F78B4",
        "#B2DF8A",
        "#33A02C",
        "#FB9A99",
        "#E31A1C",
        "#FDBF6F",
        "#FF7F00",
        "#CAB2D6",
        "#6A3D9A",
        "#FFFF99",
        "#B15928",
        "#bd18ea",
        "#2ef4ca",
        "#f4cced",
        "#f4cc03",
        "#05188a",
        "#e5a25a",
        "#06f106",
        "#85848f",
        "#000000",
        "#076f25",
        "#93cd7f",
        "#4d0776",
        "#ffffff"
      )
      colBreaks = NA
      tmyPal=NULL
      if(is.null(tmyPal) == TRUE) {
        colBreaks = 10
        tmyPal = myPal(colBreaks)
      } else {
        colBreaks = length(tmyPal)
      }
      pc = fm
      pc = pc[hc$order,]
      pc = rbind(pc, 0)
      colorList = setClusterColors(res[[K - 1]][[3]], ct, thisPal, colorList)
      stats::heatmap(pc,
                     Colv = as.dendrogram(hc), 
                     Rowv = NA, 
                     symm = FALSE, 
                     scale = "none",
                     col = tmyPal,
                     na.rm = TRUE, 
                     labRow = F,
                     labCol = F, 
                     mar = c(5,5),
                     main = paste("Consensus matrix k=", K, sep = ""),
                     ColSideColors=colorList[[1]]
      ) 
      legend("topright", legend = unique(ct), fill = unique(colorList[[1]]),
             horiz = FALSE)  
    } else{
      pc<-res$res$pc
      ct<-res$res$consensusClass
      hc<-res$res$consensusTree
      tmyPal<-res$res$tmyPal
      clrs<-res$res$clrs
      sc=res$res$sc
      K<-length(unique(ct))
      stats::heatmap(x=pc,
                     Colv = as.dendrogram(hc),
                     Rowv = NA,
                     symm = FALSE,
                     scale = "none",
                     col = tmyPal,
                     na.rm = TRUE,
                     labRow = F,
                     labCol = F,
                     mar = c(5,5),
                     main = paste("Consensus matrix k=", K, sep = ""),
                     ColSideColors=clrs[[1]]
      )
      legend("topright", legend = unique(ct), fill = unique(clrs[[1]]), horiz = FALSE)
    }
    dev.off()
  }
)

observeEvent(input$cbabt, {
  output$cdfing <-  renderPlot(
     {
       res<-cbarun()
       clustermethod<-isolate({input$clustermethod})
       if(clustermethod =="M3C"){
         require(ggpubr)
         ggarrange(plotlist=res$plots,nrow=2,ncol=2)
       } else{
         CDF(res$ml)
       }
    }
  )
})

observeEvent(input$cbabt, {
  output$cdf <- renderUI({
    plotOutput("cdfing",
               width =   paste0(isolate({input$cdfpwidth}),"%"),
               height = isolate({input$cdfpheight}))
  })})



output$savecdf <- downloadHandler(
  
  filename = function() {
    paste0("CDF-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cdfwidth}),height = isolate({input$cdfheight})) # open the pdf device
    res<-cbarun()
    clustermethod<-isolate({input$clustermethod})
    if(clustermethod =="M3C"){
      library(ggpubr)
      print(ggarrange(plotlist=res$plots,nrow=2,ncol=2))
    } else{
      print(CDF(res$ml))
    }
    dev.off()
  }
)

observeEvent(input$cbabt, {
  output$silhouetteing <-  renderPlot(
    {
      res<-cbarun()
      clustermethod<-isolate({input$clustermethod})
      if(clustermethod =="M3C"){
        require(ggpubr)
        K=length(unique(res$assignments))
        sil=silhouette_SimilarityMatrix(res$realdataresults[[K]]$assignments, res$realdataresults[[K]]$consensus_matrix)
        # sil = silhouette(res$realdataresults[[K]]$assignments, res$realdataresults[[K]]$consensus_matrix)
        plot(sil, col = 2:(K + 1),main="Silhouette plot")
      } else{
        sil = silhouette_SimilarityMatrix(res$res$consensusClass, res$res$consensusMatrix)
        plot(sil, col = 2:(length(unique(res$res$consensusClass)) + 1),main="Silhouette plot")
      }
    }
  )
})

observeEvent(input$cbabt, {
  output$silhouette <- renderUI({
    plotOutput("silhouetteing",
               width =   paste0(isolate({input$silhouettepwidth}),"%"),
               height = isolate({input$silhouetteppheight}))
  })})

output$saveSilhouette <- downloadHandler(
  filename = function() {
    paste0("Silhouette-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$Silhouettewidth}),height = isolate({input$Silhouetteheight})) # open the pdf device
    res<-cbarun()
    clustermethod<-isolate({input$clustermethod})
    if(clustermethod =="M3C"){
      require(ggpubr)
      K=length(unique(res$assignments))
      sil=silhouette_SimilarityMatrix(res$realdataresults[[K]]$assignments, res$realdataresults[[K]]$consensus_matrix)
      plot(sil, col = 2:(K + 1),main="Silhouette plot")
    } else{
      sil = silhouette_SimilarityMatrix(res$res$consensusClass, res$res$consensusMatrix)
      plot(sil, col = 2:(length(unique(res$res$consensusClass)) + 1),main="Silhouette plot")
    }
    dev.off()
  }
)

cbares<-eventReactive(input$cbabt,{
  input$cbabt
  data <- isolate({
    cbaprerun()
  })
  clustermethod<-isolate({input$clustermethod})
  res<-isolate({cbarun()})
  if(clustermethod == "M3C"){
    K<-length(unique(res$assignments))
    assign<-res$realdataresults[[K]]$assignments
    assign<-as.data.frame(assign)
    assign$assign<-paste("Subtype", assign$assign)
  }else{
    assign<-res$res$consensusClass 
    assign<-as.data.frame(assign)
    assign$assign<-paste("Subtype", assign$assign)
  }
  list(assign=assign,data=data)
})

observeEvent(input$cbabt, {
  output$cbaassign <-  DT::renderDT({
    DT::datatable(cbares()$assign,options=list(scrollX=TRUE))
  })
})

output$savecbaassign <- downloadHandler(
  filename = function() {
    paste0("Subtype-assignment-in-the-training-set-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(cbares()$assign, file)
  }
)


validst<-eventReactive(input$cbavalidationbt,{
  input$cbavalidationbt
  assign <- isolate({
    cbares()
  })
  assign<-assign$assign    
  cbavalimethod<-isolate({input$cbavalimethod})     
  
  withProgress(message = "Preprocessing the input data",
               detail = "This may take a while...",
               value = 3,
               {
  preres<-isolate({cbaprerun()})
  # browser()
  expres<-as.data.frame(t(preres$expres))
  clinical<-preres$clinical

  if(identical(row.names(expres),row.names(assign))){
    expres$subtype<-assign$assign
  }else{stop("The rownames of expres and assign are not identical, please check!")}
  
  if(identical(row.names(clinical),row.names(assign))){
    clinical$Subtype<-assign$assign
  }else{
    index<-intersect(row.names(assign),row.names(clinical))
    clinical<-clinical[index,]
    assign<-assign[index,,drop=F]
    clinical$Subtype<-assign$assign
  }
  
  #
  traindata<-list(expres=preres$data$expres,clinical=clinical)
  #
  
  expres<-split(expres,expres$subtype)

  fun<-function(x){
    x[,which(names(x)!="subtype")]
  }

  expres<-lapply(expres,fun)
  expres<-lapply(expres,colMeans)
  
  expres<-do.call(cbind,expres)
  expres<-as.data.frame(expres)
  # expres<-expres[which(row.names(expres)!="subtype"),]
  
  newdata1<-isolate({rawdata2()})
  # newdata<-exprs(newdata1)
  newdata<-newdata1$expres
  
  index<-intersect(row.names(expres),row.names(newdata))
  newdata<-newdata[index,]

  # names(expres)<-paste("Subtype",names(expres))
# print(cbavalimethod)
  validassign<-assignSubtype(ImmuneMW = expres, dataGE = newdata,method=cbavalimethod)
  
  # validclinical<-pData(newdata1)
  validclinical<-newdata1$clinical
  if(identical(row.names(validassign),row.names(validclinical))){
    validclinical$Subtype<-validassign$Subtype
  }else{
    index<-intersect(row.names(validassign),row.names(validclinical))
    validclinical<-validclinical[index,]
    validassign<-validassign[index]
    validclinical$Subtype<-validassign$Subtype
  }
  validata<-list(expres=newdata1$expres,clinical=validclinical)
  list(traindata=traindata,validata=validata,validassign=validassign)
               })
                       })


observeEvent(input$cbavalidationbt, {
  disable("cbavalidationbt")
})

observeEvent(input$cbavalidationbt, {
  updateCollapse(session, "collapsecbavalidate", open = "Subtype Validation", close = "Descriptions and parameters")
})

observeEvent(validst(), {
  shinyjs::show("cbavalid_wrapper")
  # shinyjs::show("cdf_wrapper")
  # shinyjs::show("silhouette_wrapper")
  # shinyjs::show("cbaassign_wrapper")
  enable("cbavalidationbt")
})


observeEvent(input$cbavalidationbt, {
  output$cbavalid <-  DT::renderDT({
    DT::datatable(validst()$validassign,options=list(scrollX=TRUE))
  })
})

output$savecbavalid <- downloadHandler(
  filename = function() {
    paste0("Subtype-assignment-in-the-validation-set-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(validst()$validassign, file)
  }
)



# observe({
#   if(
#     is.null(input$Learner) ||
#     input$Learner == "" ||
#     is.null(input$validat) ||
#     input$validat == "" ||
#     is.null(input$nrtime) ||
#     input$nrtime == "" ||
#     is.null(input$nrstatus) ||
#     input$nrstatus == ""
#   ){
#     disable("nrbt")
#   }
#   else{
#     enable("nrbt")
#   }
# })
# 

observeEvent(input$page_before_cba, {
  newtab <- switch(input$tabs, "meta" = "cba","cba" = "meta")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})
observeEvent(input$page_after_cba, {
  newtab <- switch(input$tabs, "cba" = "cbaclinical","cbaclinical" = "cba")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})





