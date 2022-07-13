observe({
  data <- rawdata()$clinical
  updateSelectizeInput(session,
                       'wgvariable',
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
                       selected = "wgvariable")
})


pres<-eventReactive(input$datacleanbt,{
  input$datacleanbt
  data <- isolate({
    rawdata()
  })

  withProgress(message = "Preprocessing the input data",
               detail = "This may take a while...",
               value = 3,
               {
  clinical <- data$clinical
  expres <- data$expres
  index <- intersect(names(expres), row.names(clinical))
  clinical <- clinical[index, ]
  expres <- expres[ ,index]
  variable<-isolate({input$wgvariable})
  clinical1<-clinical
  clinical<-subset(clinical,select=variable,drop = FALSE)

  if(FALSE %in% lapply(clinical,is.numeric)){
    x<- as.data.frame(which(lapply(clinical,is.numeric)==FALSE))
    if(length(row.names(x) == 1)) {
      createAlert(
        session,
        "wgcnapremess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content =  paste(
          "The variable ",
          paste(row.names(x)),
          "You selected is not numeric variable, please re-select numeric variables or recode is as numeric variable"
        ),
        append = T,
        dismiss=T
      )
    } else{
      createAlert(
        session,
        "wgcnapremess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content =  paste(
          "The variable ",
          paste(row.names(x), collapse = " and "),
          "You selected are not numeric variables, please re-select numeric variables or recode them as numeric variables"
        ),
        append = T,
        dismiss=T,
      )
    }
    return(NULL)
  } else {
  closeAlert(session, "wgcnapremess")
  clinical<-na.omit(clinical)

  expres<-expres[,row.names(clinical)]
  topvar<-isolate({input$topvar})
  selectvar<-isolate({input$setopvar})
  thresholdZ.k<-isolate({input$ZK})
  res<-preproc(expres,clinical,thresholdZ.k,topvar,selectvar)
  }
  index<-intersect(row.names(clinical),row.names(clinical1))
  clinical1<-clinical1[index,]
  res$clinical1<-clinical1
  res
})

})

observe({
  if(is.null(input$wgvariable) || input$wgvariable == ""|| length(input$wgvariable)<2){
    disable("datacleanbt")
  }
  else{
    enable("datacleanbt")
  }
})

# output$wgcnapre <- renderUI(includeHTML("./WGCNA_preprocess.html"))
observeEvent(input$datacleanbt,{
  updateCollapse(session, "collapseWGCNApreprocess", open = "Data Preprocessing", close = "Descriptions and parameters for data preprocessing")

})

observeEvent(input$datacleanbt, {
  output$sampletreeing  <- renderPlot({
    closeAlert(session, "wgcnapremess")
    plotDendroAndColors(
      pres()$sampleTree,
      groupLabels = names(pres()$datColors),
      colors = pres()$datColors,
      main = "Sample dendrogram and trait heatmap"
    )
  })
})

observeEvent(input$datacleanbt, {
  output$sampletree<- renderUI({
    plotOutput("sampletreeing",
               width = paste0(isolate({input$sampwidth}), "%"),
               height = isolate({input$samheight}))
  })})

output$downloadsampletree <- downloadHandler(

  filename = function() {
    paste0("Sample-dendrogram-and-trait-heatmap-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$sampletreewidthdl}),height = isolate({input$sampletreeheightdl})) # open the pdf device
    print(plotDendroAndColors(
      pres()$sampleTree,
      groupLabels = names(pres()$datColors),
      colors = pres()$datColors,
      main = "Sample dendrogram and trait heatmap"
    ))
    dev.off()
  }
)

observeEvent(input$datacleanbt, {
  output$datacleandescrip <- renderText({
    paste(pres()$removed.samples[[2]]," samples/sample were/was detected as outliers and ", pres()$removed.samples[[1]], " samples were retained for subsequent analysis")

  })

})

# observeEvent(input$datacleanbt, {
#
#   disable("setopvar")
#   disable("topvar")
#   disable("ZK")
#   disable("wgvariable")
# })
observeEvent(pres(), {
  shinyjs::show("sampletree_wrapper")
  enable("setopvar")
  enable("topvar")
  enable("ZK")
  enable("wgvariable")
  shinyjs::show("wgcna1_wrapper")
})


WGCNArun<-eventReactive(input$WGCNAbt,{
  input$WGCNAbt
  expres<-isolate({pres()$expres})
  clinical<-isolate({pres()$clinical})
  clinical1<-isolate({pres()$clinical1})

  withProgress(message = "Constructing network and detecting module",
               detail = "This may take a while...",
               value = 3,
               {
  autowgcna(
    expres=expres,
    clinical=clinical,
    clinical1=clinical1,
    RsquaredCut=isolate({input$RsquaredCut}),
    networkType=  isolate({input$networkType}),
    corType= isolate({input$corType}),
    TOMType=isolate({input$TOMType}),
    deepSplit=isolate({input$deepSplit}),
    detectCutHeight=isolate({input$detectCutHeight}),
    minModuleSize=isolate({input$minModuleSize}),
    reassignThreshold=isolate({input$reassignThreshold}),
    mergeCutHeight=isolate({input$mergeCutHeight}),
    numericLabels=isolate({input$numericLabels}),
    pamRespectsDendro=isolate({input$pamRespectsDendro})
  )
               })
})



observeEvent(input$WGCNAbt, {
  output$sftdistributioning  <- renderPlot({
    closeAlert(session, "wgcnapremess")
    powers = 1:20
    sft<-WGCNArun()$sft
    par(mfrow = c(2,2));
    cex1 = 0.9;
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red");
    abline(h=isolate({input$RsquaredCut}),col="red")
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    plot(sft$fitIndices[,1], sft$fitIndices[,6],
         xlab="Soft Threshold (power)",ylab="Median Connectivity", type="n",
         main = paste("Median Connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,6], labels=powers, cex=cex1,col="red")
    plot(sft$fitIndices[,1], sft$fitIndices[,7],
         xlab="Soft Threshold (power)",ylab="Max Connectivity", type="n",
         main = paste("Max Connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,7], labels=powers, cex=cex1,col="red")

  })

})

observeEvent(input$WGCNAbt,{
  updateCollapse(session, "collapseWGCNAnetwork", open = "Network construction and module detection", close = "Descriptions and parameters for network construction and module detection")
  
})

observeEvent(input$WGCNAbt, {

  output$sftdistribution<- renderUI({
    plotOutput("sftdistributioning",
               width = paste0(isolate({input$sftwidth}), "%"),
               height = isolate({input$sftheight}))
  })})



observeEvent(input$WGCNAbt,{
  output$sftdiscrp <- renderText(
    {
      if(!is.na(WGCNArun()$sft$powerEstimate)){
        paste("Soft-thresholding power",WGCNArun()$sft$powerEstimate,"is selected as the soft-thresholding power at Rsquared Cutoff:", isolate({input$RsquaredCut}), sep=" ")
      } else{
        showModal(modalDialog(
          title = "Please note!",
          paste("No soft threshold is identified,  for which the scale free topology fit R^2 exceeds Rsquared Cutoff:",isolate({input$RsquaredCut}),sep=" ,"),
          easyClose = TRUE,
          footer = NULL
        ))
      }
    }
  )
})





output$downloadsftdistribution <- downloadHandler(

  filename = function() {
    paste0("Soft-threshold-selection-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$sftdistributionwidthdl}),height = isolate({input$sftdistributionheightdl})) # open the pdf device
    powers = 1:20
          sft<-WGCNArun()$sft
          par(mfrow = c(2,2));
          cex1 = 0.9;
          plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
               xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
               main = paste("Scale independence"));
          text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
               labels=powers,cex=cex1,col="red");
          abline(h=isolate({input$RsquaredCut}),col="red")
          plot(sft$fitIndices[,1], sft$fitIndices[,5],
               xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
               main = paste("Mean connectivity"))
          text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
          plot(sft$fitIndices[,1], sft$fitIndices[,6],
               xlab="Soft Threshold (power)",ylab="Median Connectivity", type="n",
               main = paste("Median Connectivity"))
          text(sft$fitIndices[,1], sft$fitIndices[,6], labels=powers, cex=cex1,col="red")
          plot(sft$fitIndices[,1], sft$fitIndices[,7],
               xlab="Soft Threshold (power)",ylab="Max Connectivity", type="n",
               main = paste("Max Connectivity"))
          text(sft$fitIndices[,1], sft$fitIndices[,7], labels=powers, cex=cex1,col="red")

    dev.off()
  }
)



observeEvent(input$WGCNAbt,{
  output$modengraming <- renderPlot(

    {
       closeAlert(session, "wgcnapremess")
      if(!is.na(WGCNArun()$sft$powerEstimate)){
        net<-WGCNArun()$net
        mergedColors = labels2colors(net$colors)
        plotDendroAndColors(net$dendrograms[[1]],
                            mergedColors[net$blockGenes[[1]]],
                            "Module colors",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05)
      } else{
        showModal(modalDialog(
          title = "Please note!",
          paste(
            "No soft threshold is identified,  for which the scale free topology fit R^2 exceeds Rsquared Cutoff:",
            isolate({
              input$RsquaredCut
            }),
            sep = " ,"
          ),
          easyClose = TRUE,
          footer = NULL
        ))

      }

    }
  )
})



observeEvent(input$WGCNAbt, {
  output$modengram<- renderUI({
    plotOutput("modengraming",
               width = paste0(isolate({input$MDwidth}), "%"),
               height = isolate({input$MDheight}))
  })})


output$downloadmodengram <- downloadHandler(

  filename = function() {
    paste0("Module-dendrogram-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$modengramwidthdl}),height = isolate({input$modengramheightdl})) # open the pdf device
          net<-WGCNArun()$net
          mergedColors = labels2colors(net$colors)
          print(plotDendroAndColors(net$dendrograms[[1]],
                              mergedColors[net$blockGenes[[1]]],
                              "Module colors",
                              dendroLabels = FALSE, hang = 0.03,
                              addGuide = TRUE, guideHang = 0.05))
    dev.off()
  }
)
observeEvent(WGCNArun(), {
  shinyjs::show("sftdistribution_wrapper")
  shinyjs::show("modengram_wrapper")
  shinyjs::show("wgcna2_wrapper")
})



observe({
  updateSelectizeInput(session,
                       'trait',
                       choices = {
                         input$wgvariable
                       },
                       server = TRUE,
                       selected = NULL)
})

observe({
  MEs<-WGCNArun()$MEs
  updateSelectizeInput(session,
                       'module',
                       choices = {
                         substring(names(MEs), 3)
                       },
                       server = TRUE,
                       selected = NULL)
})

observe({
  MEs<-WGCNArun()$MEs
  updateSelectizeInput(session,
                       'output',
                       choices = {
                         c(substring(names(MEs), 3),"Non-grey modules")
                       },
                       server = TRUE,
                       selected = NULL)
})


netsummary<-eventReactive( input$MTRbt,{
  input$MTRbt
  clinical<-isolate({WGCNArun()$clinical})
  clinical1<-isolate({WGCNArun()$clinical1})
  expres<-isolate({WGCNArun()$expres})
  trait<-isolate({input$trait})
  moduleColors<-isolate(WGCNArun()$moduleColors)
  module = isolate({input$module})
  netres(clinical=clinical,clinical1=clinical1,expres=expres,trait=trait,moduleColors=moduleColors, module=module,output=isolate({input$output}))
})


observeEvent(input$MTRbt,{
  updateCollapse(session, "collapseWGCNAmtr", open = "Module-trait relationships", close = "Descriptions and parameters for module-trait relationships")
  
})


observeEvent(input$MTRbt, {
  output$MTRploting  <- renderPlot({

    MTR(WGCNArun(),isolate({input$BMar}), isolate({input$LMar}), isolate({input$TMar}), isolate({input$RMar}))
  })
})

observeEvent(input$MTRbt, {

  output$MTRplot<- renderUI({
    plotOutput("MTRploting",
               width = paste0(isolate({input$MTRwidth}), "%"),
               height = isolate({input$MTRheight}))
  })})


output$downloadMTR <- downloadHandler(
  filename = function() {
    paste0("Module-trait-relationships-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$MTRwidthdl}),height = isolate({input$MTRheightdl})) # open the pdf device
    print(MTR(WGCNArun(),isolate({input$BMar}), isolate({input$LMar}), isolate({input$TMar}), isolate({input$RMar})))
    dev.off()
  }
)



observeEvent(input$MTRbt, {
  output$GSploting  <- renderPlot({

    verboseScatterplot(abs(netsummary()$geneModuleMembership[netsummary()$moduleGenes, netsummary()$column]),
                       abs(netsummary()$geneTraitSignificance[netsummary()$moduleGenes, 1]),
                       xlab = paste("Module Membership in", isolate({input$module}), "module"),
                       ylab = paste("Gene significance for",isolate({input$trait}),sep=" "),
                       main = paste("Module membership vs gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = isolate({input$module}))
  })
})

observeEvent(input$MTRbt, {

  output$GSplot<- renderUI({
    plotOutput("GSploting",
               width = paste0(isolate({input$GVMwidth}), "%"),
               height = isolate({input$GVMheight}))
  })})

output$downloadGM<- downloadHandler(

  filename = function() {
    paste0("GS-vs-MM-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$GMwidthdl}),height = isolate({input$GMheightdl})) # open the pdf device
    print(verboseScatterplot(abs(netsummary()$geneModuleMembership[netsummary()$moduleGenes, netsummary()$column]),
                              abs(netsummary()$geneTraitSignificance[netsummary()$moduleGenes, 1]),
                              xlab = paste("Module Membership in", isolate({input$module}), "module"),
                              ylab = paste("Gene significance for",isolate({input$trait}),sep=" "),
                              main = paste("Module membership vs gene significance\n"),
                              cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = isolate({input$module})))
    dev.off()
  }
)




observeEvent(input$MTRbt, {
    output$netsum <-  DT::renderDT({
    DT::datatable(netsummary()$Output,options=list(scrollX=TRUE))
  }
  )
})

observeEvent(netsummary(), {
  shinyjs::show("MTR_wrapper")
  shinyjs::show("GM_wrapper")
  shinyjs::show("netsum_wrapper")
})

output$downloadnetsum <- downloadHandler(
  filename = function() {
    paste0("Output-genes-of-WGCNA-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(netsummary()$Output, file)
  }
)



observeEvent(input$page_after_WGCNA, {
  newtab <- switch(input$tabs, "wgcna" = "msurv","msurv" = "wgcna")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_before_WGCNA, {
  newtab <- switch(input$tabs,"dataset" = "wgcna","wgcna" = "dataset")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})



#
