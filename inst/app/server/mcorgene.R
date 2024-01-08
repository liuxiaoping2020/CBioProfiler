observe({
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'mcorgene',
                       choices = {
                         if (!is.null(data)) {
                           if (class(data) ==  class(data.frame()))
                             choices <- as.character(rownames(data))
                           if (class(data) !=  class(data.frame()))
                             choices <-
                               as.character(rownames(as.data.frame(data)))
                           choices
                         }
                       },
                       server = TRUE,
                       selected = "A1BG")
})



mcor<- eventReactive(input$mcorgenebt,{
  input$mcorgenebt
  withProgress(message = "Performing correlation analsis",
               detail = "This may take a while...",
               value =5,
               {
  data <- isolate({
    rawdata()
  })
  data <- lapply(data, as.data.frame)
  expres<-as.data.frame(t(data$expres))
  gene<-isolate({input$mcorgene})
  geneexp<-subset(expres,select=gene)

  expres[,gene]<-NULL
  mcormethod<-isolate({input$mcormethod})
  library(WGCNA)
  res<-corAndPvalue(geneexp, expres,method = mcormethod)
  res<-merge(as.data.frame(t(res$cor)),as.data.frame(t(res$p)),by=0)
  names(res)<-c("Gene","R","Pvalue")

  corsigcut<-isolate({input$corsigcut})

  if(isolate({input$padjust})==TRUE){
    mcorrectmethod<-isolate({input$mcorrectmethod})
    res$padjust<-p.adjust(res$Pvalue,method = mcorrectmethod)
    names(res)<-c("Gene","R","Pvalue","Padjusted")
    
    if(isolate({input$mcorcutmeth})=='Threshold value'){
      tab<-res[res$Padjusted<corsigcut,]
    }else{
      tab<-res[order(res$Padjusted,decreasing = FALSE),]
      tab<-tab[1:isolate({input$mcortp}),,drop=F]
    }
    if(nrow(tab)==0) {
      createAlert(
        session,
        "mcorgenemess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste(
          "No genes meet the significance cutoff (",
          corsigcut,
          ") at the correction method",
          paste("'", mcorrectmethod, "'"),
          "you specified",
          sep = " "
        ),
        append = T,
        dismiss=T
      )
      return(NULL)
    } else{

    closeAlert(session, "mcorgenemess")
    tab$Gene<-as.factor(tab$Gene)
    barplot<-ggplot(data=tab)+geom_bar(aes(x=Gene,y=R, fill=Padjusted), stat='identity')+
      coord_flip() +scale_fill_gradient(low="red", high = "blue") +xlab("")+
      ylab("Correlation coefficient")+
      theme(axis.text.x=element_text(color="black",size=48),axis.text.y=element_text(color="black", size=48))+
      scale_y_continuous(expand=c(0, 0))+
      scale_x_discrete(expand=c(0,0))+theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


    bubbleplot<-ggplot(data=tab,aes(x=R,y=Gene))+geom_point(aes(col=R,size=Padjusted))+
      scale_color_gradient(low="blue",high="red")+
      labs(color="Correlation coefficient", size="P value")+theme_bw()+
      ylab("")+xlab("Correlation coefficient")+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),
            axis.text=element_text(size=12,color="black"),legend.text=element_text(size=10))
    }

  } else{
    if(isolate({input$mcorcutmeth})=='Threshold value'){
      tab<-res[res$Pvalue<corsigcut,]
    }else{
      tab<-res[order(res$Pvalue,decreasing = FALSE),]
      tab<-tab[1:isolate({input$mcortp}),,drop=F]
      
    }

    
    if(nrow(tab)==0) {
      createAlert(
        session,
        "mcorgenemess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste(
          "No genes meet the significance cutoff of P value <",
          corsigcut,
          "you specified",
          sep = " "
        ),
        append = T,
        dismiss=T
      )
      return(NULL)
    } else {

    closeAlert(session, "mcorgenemess")
    tab$Gene<-as.factor(tab$Gene)
    barplot<-ggplot(data=tab)+geom_bar(aes(x=Gene,y=R, fill=Pvalue), stat='identity')+
      coord_flip() +scale_fill_gradient(low="red", high = "blue") +xlab("")+
      ylab("Correlation coefficient")+
      theme(axis.text.x=element_text(color="black",size=48),axis.text.y=element_text(color="black", size=48))+
      scale_y_continuous(expand=c(0, 0))+
      scale_x_discrete(expand=c(0,0))+theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

    bubbleplot<-ggplot(data=tab,aes(x=R,y=Gene))+geom_point(aes(col=R,size=Pvalue))+
      scale_color_gradient(low="blue",high="red")+
      labs(color="Correlation coefficient", size="P value")+theme_bw()+
      ylab("")+xlab("Correlation coefficient")+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),
            axis.text=element_text(size=12,color="black"),legend.text=element_text(size=10))
    }
  }
  res<-list(tab,barplot,bubbleplot)
  names(res)<-c("tab","barplot","bubbleplot")
  res
 })


})

observe({
  if(is.null(input$mcorgene) || input$mcorgene == ""){
    disable("mcorgenebt")
  }
  else{
    enable("mcorgenebt")
  }
})
observeEvent(input$mcorgenebt, {
  disable("mcorgenebt")
})

observeEvent(input$mcorgenebt, {
  updateCollapse(session, "collapsemcorgene", open = "Most correlated genes", close = "Descriptions and parameters")
})

observeEvent(input$mcorgenebt, {
  output$mcortable <-  DT::renderDT({
    mcor()$tab
  }
  ,rownames = FALSE)
})



observeEvent(input$mcorgenebt, {
  output$mcorbarploting  <- renderPlot({
     closeAlert(session, "mcorgenemess")
    mcor()$barplot
  })
})

observeEvent(input$mcorgenebt, {
  output$mcorbarplot<- renderUI({
    plotOutput("mcorbarploting",
               width = paste0(isolate({input$mcorbarwidth}), "%"),
               height = isolate({input$mcorbarheight}))
  })})



observeEvent(input$mcorgenebt, {
  output$mcorbubploting  <- renderPlot({
    closeAlert(session, "mcorgenemess")
    mcor()$bubbleplot
  })
})

observeEvent(input$mcorgenebt, {
  output$mcorbubplot<- renderUI({
    plotOutput("mcorbubploting",
               width = paste0(isolate({input$mcorbubwidth}), "%"),
               height = isolate({input$mcorbubheight}))
  })})

observeEvent(input$mcorgenebt, {
  shinyjs::show("mcorgenebarplot_wrapper")
  shinyjs::show("mcorgenebubbleplot_wrapper")
  shinyjs::show("mcorgenetable_wrapper")
  enable("mcorgenebt")
})


output$downmcorgenetable <- downloadHandler(
  filename = function() {
    paste0("most-correlated-gene-table-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(mcor()$tab, file)
  }
)

output$downmcorgenebarplot <- downloadHandler(
  filename = function() {
    paste0("Gene-correlation-barplot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$mcorgenebarwidthdl}),height = isolate({input$mcorgenebarheightdl})) # open the pdf device
    print(mcor()$barplot)
    dev.off()
  }
)

output$downmcorgenebubbleplot <- downloadHandler(

  filename = function() {
    paste0("Gene-correlation-bubbleplot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$mcorgenebubblewidthdl}),height = isolate({input$mcorgenebubbleheightdl})) # open the pdf device
    print(mcor()$bubbleplot)
    dev.off()
  }
)


observeEvent(input$page_before_mcorgene, {
  newtab <- switch(input$tabs, "mcorgene" = "SurvROC","SurvROC" = "mcorgene")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_mcorgene, {
  newtab <- switch(input$tabs, "gene" = "mcorgene","mcorgene" = "gene")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})


