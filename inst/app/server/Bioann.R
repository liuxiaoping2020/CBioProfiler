
observeEvent(input$MTRbt, {
  updateSelectizeInput(
    session,
    "biomakers",
    choices ="Network hub genes",
    server = TRUE,
    selected = "Network hub genes"
  )
})

observeEvent(input$DEGvisbt, {
  updateSelectizeInput(
    session,
    "biomakers",
    choices ="Significant DEGs",
    server = TRUE,
    selected = "Significant DEGs"
  )
})
observeEvent(input$msurvopbt, {
  updateSelectizeInput(
    session,
    "biomakers",
    choices ="Significant SRGs",
    server = TRUE,
    selected = "Significant SRGs"
  )
})
ve1 <- reactive({
  req(input$nrbt) && req(input$msurvopbt)
})

observeEvent(ve1(),{
  updateSelectizeInput(
    session,
    "biomakers",
    choices =c("Significant SRGs","Genes from benchmark experiment"),
    server = TRUE,
    selected = "Significant SRGs"
  )
})

ve2 <- reactive({
  req(input$nrbt) && req(input$DEGvisbt)
})

observeEvent(ve2(),{
  updateSelectizeInput(
    session,
    "biomakers",
    choices =c("Significant DEGs","Genes from benchmark experiment"),
    server = TRUE,
    selected = "Significant DEGs"
  )
})
ve3 <- reactive({
  req(input$nrbt) && req(input$MTRbt)
})

observeEvent(ve3(),{
  updateSelectizeInput(
    session,
    "biomakers",
    choices =c("Network hub genes","Genes from benchmark experiment"),
    server = TRUE,
    selected = "Network hub genes"
  )
})

enrun<-eventReactive(input$enrichbt,{
  input$enrichbt

  library(org.Hs.eg.db)

  if(isolate({input$biomakers})=="Significant DEGs"){
  DEGen<-isolate({DEGop()$Output})
  genean <- as.data.frame(org.Hs.egALIAS2EG)
  names(genean) <- c("geneid", "gene")
  DEGen$gene<-row.names(DEGen)
  DEGen <- merge(genean, DEGen, by = "gene")
  DEGen<-DEGen[order(DEGen$LogFC,decreasing=T),]
  gene<-DEGen$geneid
  geneList = DEGen$LogFC
  }else if(isolate({input$biomakers})=="Significant SRGs"){
    DEGen<-isolate({msurvsig()$Output})
    genean <- as.data.frame(org.Hs.egALIAS2EG)
    names(genean) <- c("geneid", "gene")
    DEGen$gene<-row.names(DEGen)
    DEGen <- merge(genean, DEGen, by = "gene")
    DEGen<-DEGen[order(DEGen$Coefficient,decreasing =T),]
    gene<-DEGen$geneid
    geneList = DEGen$Coefficient
  }else if(isolate({input$biomakers})=="Network hub genes"){
    DEGen<-isolate({netsummary()$Output})
    genean <- as.data.frame(org.Hs.egALIAS2EG)
    names(genean) <- c("geneid", "gene")
    DEGen$gene<-row.names(DEGen)
    DEGen <- merge(genean, DEGen, by = "gene")
    DEGen<-DEGen[order(DEGen[,2],decreasing=T),]
    gene<-DEGen$geneid
    geneList = DEGen[,2]
  } else if(isolate({input$biomakers})=="Genes from benchmark experiment"){
    DEGen<-isolate({nestrun()$listidx})
    genean <- as.data.frame(org.Hs.egALIAS2EG)
    names(genean) <- c("geneid", "gene")
    DEGen$gene<-row.names(DEGen)
    DEGen <- merge(genean, DEGen, by = "gene")
    DEGen<-DEGen[order(DEGen$coeff,decreasing=T),]
    gene<-DEGen$geneid
    geneList = DEGen$coeff
  }

  names(geneList) = as.character(DEGen$geneid)
  geneList = sort(geneList, decreasing = TRUE)
  ana<-isolate({input$FEAana})
  ont<-isolate({input$ont})
  pAdjustMethod<-isolate({input$pAdjustMethod})
  minGSSize<-isolate({input$minGSSize})
  maxGSSize<-isolate(input$maxGSSize)
  pvalueCutoff<-isolate({input$pvalueCutoff})
  qvalueCutoff<-isolate({input$qvalueCutoff})
  readable<-isolate({input$readable})
  gseamethod<-isolate({input$gseamethod})
  nPerm<-isolate({input$nPerm})
  method<-isolate(input$FEAmethod)
  withProgress(message = "Performing functional enrichment analysis",
               detail = "This may take a while...",
               value = 3,
               {
                 FEAres<-FEA(
                   ana=ana,
                   method=method,
                   geneList=geneList,
                   nPerm=nPerm,
                   minGSSize=minGSSize,
                   maxGSSize=maxGSSize,
                   pvalueCutoff=pvalueCutoff,
                   pAdjustMethod=pAdjustMethod,
                   qvalueCutoff=qvalueCutoff,
                   gseamethod=gseamethod,
                   gene=gene,
                   readable=readable,
                   ont= ont
                 )
                 if(nrow(FEAres@result)==0){
                   createAlert(
                     session,
                     "FEAmess",
                     "exampleAlert",
                     title = "Please note!",
                     style =  "danger",
                     content = paste(
                       "No gene sets/functional terms meet the significance cutoff you specified"),
                     append = T,
                     dismiss=T
                   )
                   return(NULL)
                 } else{
                   list(FEAres=FEAres,geneList=geneList)
                 }

               })

})


observe({
  if (is.null(input$biomakers) ||
      input$biomakers == "" ||
      is.null(input$FEAana)||
      input$FEAana == ""
      ) {
    disable("enrichbt")
  }
  else{
    enable("enrichbt")
  }
})

observeEvent(input$enrichbt, {
  disable("enrichbt")
})


observeEvent(input$enrichbt, {
  updateCollapse(session, "collapsebioan", open = "Functional enrichment analysis", close = "Descriptions and parameters")
  
})



observeEvent(is.null(enrun()), {
  enable("enrichbt")
})

observeEvent(input$enrichbt, {
  output$DEGenplotting <- renderPlot({
    if(isolate({
      input$FEAmethod
    }) == "ORA" &
    isolate({
      input$enplottype
    }) %in% c(
      "Heatmap",
      "Ridgeline plot",
      "Geneset enrichment plot1",
      "Geneset enrichment plot2"
    )) {
      showModal(modalDialog(
        title = "Please note!",
        paste(
          "The 'Plot type' you choose does not support the visualization of the result of the 'Functional analysis method",
          isolate({
            input$FEAmethod
          }), "you specified",

          sep = " "
        ),
        easyClose = TRUE,
        footer = NULL
      ))
    } else if (isolate({
      input$FEAmethod
    }) == "GSEA" &
    isolate({
      input$enplottype
    }) == "Bar plot"

    ){
      showModal(modalDialog(
        title = "Please note!",
        paste(
          "The 'Plot type' you choose does not support the visualization of the result of the 'Functional analysis method",
          isolate({
            input$FEAmethod
          }), "you specified",

          sep = " "
        ),
        easyClose = TRUE,
        footer = NULL
      ))
    } else{
      res<-isolate({enrun()$FEAres})
      enplottype<-isolate({input$enplottype})
      geneList<-isolate({enrun()$geneList})
      if (enplottype == "Bar plot") {
        barplot(
          res,
          color = isolate({input$barcolor}),
          showCategory = isolate({input$barshowCategory}),
          font.size = isolate({input$barfont.size}),
          label_format = isolate({input$barlabel_format})
        )
      } else if (enplottype == 'Dot plot') {
        dotplot(
          object = res,
          x = isolate({input$dotxvar}),
          color = isolate({input$dotcolor}),
          showCategory = isolate({input$dotshowCategory}),
          size = NULL,
          split = NULL,
          font.size = isolate({input$dotfont.size}),
          title = "",
          label_format = isolate({input$dotlabel_format}),
          orderBy = "x"
        )
      } else if (enplottype == 'Enrichment Map') {
        res1 <- pairwise_termsim(res)
        emapplot(
          x = res1,
          showCategory = isolate({input$emshowCategory}),
          color = isolate({input$emcolor}),
          layout = isolate({input$emlayout}),
          node_scale = NULL,
          line_scale = NULL,
          min_edge = isolate({input$emmin_edge}),
          node_label_size = NULL,
          cex_label_category = isolate({input$emcex_label_category}),
          cex_category = isolate({input$emcex_category}),
          cex_line = isolate({input$emcex_line})
        )
      } else if (enplottype == 'Gene-concept network') {
        cnetplot(
          x = res,
          showCategory = isolate({input$cneshowCategory}),
          foldChange = NULL,
          layout = "kk",
          colorEdge = isolate({input$cnecolorEdge}),
          circular = isolate({input$cnecircular}),
          node_label = isolate({input$cnenode_label}),
          cex_category = isolate({input$cnecex_category}),
          cex_gene = isolate({input$cnecex_gene}),
          node_label_size = NULL,
          cex_label_category = isolate({input$cnecex_label_category}),
          cex_label_gene = isolate({input$cnecex_label_gene})
        )
      } else if (enplottype == 'Heatmap') {
        heatplot(x = res,
                 showCategory = isolate({input$heatshowCategory}),
                 foldChange = geneList)
      } else if (enplottype == 'Ridgeline plot') {
        ridgeplot(
          x = res,
          showCategory = isolate({input$regshowCategory}),
          fill = isolate({input$regfill}),
          core_enrichment = isolate({input$regcore_enrichment}),
          label_format = isolate({input$reglabel_format})
        )
      } else if (enplottype == 'Geneset enrichment plot1') {
        gseaplot(
          x = res,
          geneSetID = isolate({input$gsgeneSetID1}),
          by = isolate({input$gsby}),
          title = res$Description[isolate({input$gsgeneSetID1})],
          color = isolate({input$gscolor}),
          color.line = isolate({input$gscolor.line}),
          color.vline = isolate({input$gscolor.vline})
        )
      } else if (enplottype == 'Geneset enrichment plot2') {
        gseaplot2(
          x = res,
          geneSetID = isolate({input$gsegeneSetID2}),
          title = res$Description[isolate({input$gsegeneSetID2})],
          color = isolate({input$gsecolor}),
          base_size = isolate({input$gsebase_size}),
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = isolate({input$gsepvalue_table}),
          ES_geom = isolate({input$gseES_geom})
        )
      }
    }

  })
})

observeEvent(input$enrichbt, {
  output$DEGenplot <- renderUI({
    plotOutput("DEGenplotting",
               width = paste0(isolate({input$enriwidth}), "%"),
               height = isolate({input$enriheight}))
  })})

observeEvent(enrun(), {
  shinyjs::show("enrichtable_wrapper")
  shinyjs::show("DEGenplot_wrapper")
  enable("enrichbt")
})

output$downloadDEGenplot <- downloadHandler(
  filename = function() {
    paste0("DEG-enrichment-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$DEGenplotwidthdl}),height = isolate({input$DEGenplotheightdl})) # open the pdf device
    if(isolate({
      input$FEAmethod
    }) == "ORA" &
    isolate({
      input$enplottype
    }) %in% c(
      "Heatmap",
      "Ridgeline plot",
      "Geneset enrichment plot1",
      "Geneset enrichment plot2"
    )) {
      showModal(modalDialog(
        title = "Please note!",
        paste(
          "The 'Plot type' you choose does not support the visualization of the result of the 'Functional analysis method",
          isolate({
            input$FEAmethod
          }), "you specified",

          sep = " "
        ),
        easyClose = TRUE,
        footer = NULL
      ))
    } else if (isolate({
      input$FEAmethod
    }) == "GSEA" &
    isolate({
      input$enplottype
    }) == "Bar plot"

    ){
      showModal(modalDialog(
        title = "Please note!",
        paste(
          "The 'Plot type' you choose does not support the visualization of the result of the 'Functional analysis method",
          isolate({
            input$FEAmethod
          }), "you specified",

          sep = " "
        ),
        easyClose = TRUE,
        footer = NULL
      ))
    } else{
      res<-isolate({enrun()$FEAres})
      enplottype<-isolate({input$enplottype})
      geneList<-isolate({enrun()$geneList})
      if (enplottype == "Bar plot") {
        print(barplot(
          res,
          color = isolate({input$barcolor}),
          showCategory = isolate({input$barshowCategory}),
          font.size = isolate({input$barfont.size}),
          label_format = isolate({input$barlabel_format})
        ))
      } else if (enplottype == 'Dot plot') {
        print(dotplot(
          object = res,
          x = isolate({input$dotxvar}),
          color = isolate({input$dotcolor}),
          showCategory = isolate({input$dotshowCategory}),
          size = NULL,
          split = NULL,
          font.size = isolate({input$dotfont.size}),
          title = "",
          label_format = isolate({input$dotlabel_format}),
          orderBy = "x"
        ))
      } else if (enplottype == 'Enrichment Map') {
        res1 <- pairwise_termsim(res)
        print(emapplot(
          x = res1,
          showCategory = isolate({input$emshowCategory}),
          color = isolate({input$emcolor}),
          layout = isolate({input$emlayout}),
          node_scale = NULL,
          line_scale = NULL,
          min_edge = isolate({input$emmin_edge}),
          node_label_size = NULL,
          cex_label_category = isolate({input$emcex_label_category}),
          cex_category = isolate({input$emcex_category}),
          cex_line = isolate({input$emcex_line})
        ))
      } else if (enplottype == 'Gene-concept network') {
        print(cnetplot(
          x = res,
          showCategory = isolate({input$cneshowCategory}),
          foldChange = NULL,
          layout = "kk",
          colorEdge = isolate({input$cnecolorEdge}),
          circular = isolate({input$cnecircular}),
          node_label = isolate({input$cnenode_label}),
          cex_category = isolate({input$cnecex_category}),
          cex_gene = isolate({input$cnecex_gene}),
          node_label_size = NULL,
          cex_label_category = isolate({input$cnecex_label_category}),
          cex_label_gene = isolate({input$cnecex_label_gene})
        ))
      } else if (enplottype == 'Heatmap') {
        print(heatplot(x = res,
                 showCategory = isolate({input$heatshowCategory}),
                 foldChange = geneList))
      } else if (enplottype == 'Ridgeline plot') {
        print(ridgeplot(
          x = res,
          showCategory = isolate({input$regshowCategory}),
          fill = isolate({input$regfill}),
          core_enrichment = isolate({input$regcore_enrichment}),
          label_format = isolate({input$reglabel_format})
        ))
      } else if (enplottype == 'Geneset enrichment plot1') {

        #geneSetID
          print(gseaplot(
            x = res,
            geneSetID = isolate({input$gsgeneSetID1}),
            by = isolate({input$gsby}),
            title = res$Description[isolate({input$gsgeneSetID1})],
            color = isolate({input$gscolor}),
            color.line = isolate({input$gscolor.line}),
            color.vline = isolate({input$gscolor.vline})
          ))
        # )
        # print(gseaplot(
        #   x = res,
        #   geneSetID = isolate({input$gsgeneSetID}),
        #   by = isolate({input$gsby}),
        #   title = res$Description[gsgeneSetID],
        #   color = isolate({input$gscolor}),
        #   color.line = isolate({input$gscolor.line}),
        #   color.vline = isolate({input$gscolor.vline})
        # ))
        # print(p)
      } else if (enplottype == 'Geneset enrichment plot2') {
        geneSetID = isolate({input$gsegeneSetID})
        print(gseaplot2(
          x = res,
          geneSetID = isolate({input$gsegeneSetID2}),
          title = res$Description[isolate({input$gsegeneSetID2})],
          color = isolate({input$gsecolor}),
          base_size = isolate({input$gsebase_size}),
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = isolate({input$gsepvalue_table}),
          ES_geom = isolate({input$gseES_geom})
        ))
      }
    }
    dev.off()



  }
)

observeEvent(input$enrichbt, {
  output$enrichtable <-  DT::renderDT({
    DT::datatable(enrun()$FEAres@result,options=list(scrollX=TRUE))
  })
})

output$downloadenrichtable <- downloadHandler(
  filename = function() {
    paste0("DEG-enrichment-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(enrun()$FEAres@result, file)
  }
)


observeEvent(input$page_before_bioan, {
  newtab <- switch(input$tabs, "drug" = "bioan","bioan" = "drug")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_bioan, {
  newtab <- switch(input$tabs, "bioan" = "meta","meta" = "bioan")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})


# ######################################################
# observeEvent(input$page_before_KM, {
#   newtab <- switch(input$tabs, "clinical" = "KM","KM" = "clinical")
#   updateTabItems(session, "tabs", newtab)
#   shinyjs::runjs("window.scrollTo(0, 50)")
# })
# 
# observeEvent(input$page_after_KM, {
#   newtab <- switch(input$tabs, "KM" = "CoxPH","CoxPH" = "KM")
#   updateTabItems(session, "tabs", newtab)
#   shinyjs::runjs("window.scrollTo(0, 50)")
# })






