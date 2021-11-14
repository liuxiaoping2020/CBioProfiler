observe({
  data <- rawdata()$expres
  updateSelectizeInput(session, 'gene1', choices = {
    if (!is.null(data)) {
      if (class(data) ==  class(data.frame())) choices <- as.character(rownames(data))
      if (class(data) !=  class(data.frame())) choices <- as.character(rownames(as.data.frame(data)))
      choices
    }
   },server = TRUE,select="A1BG",
  )
 }
)

observe({
  data <- rawdata()$expres
  updateSelectizeInput(session, 'gene2', choices = {
    if (!is.null(data)) {
      if (class(data) ==  class(data.frame())) choices <- as.character(rownames(data))
      if (class(data) !=  class(data.frame())) choices <- as.character(rownames(as.data.frame(data)))
      choices
    }
  },server = TRUE,select="TP53",
  )
}
)

ggcor<-eventReactive(input$genecorbt,{
  input$genecorbt
  data <- isolate({
    rawdata()
  })
  data <- lapply(data, as.data.frame)
  gene1<- isolate({input$gene1})
  gene2<- isolate({input$gene2})
  method<-isolate({input$cormethods})
  name<-c(gene1,gene2)
  expres<-data$expres
  cordata<-as.data.frame(t(expres[name,]))
  scatter<-ggscatter(data=cordata, y = name[1],x = name[2],
                      size = 1,
                      add = "reg.line", conf.int = TRUE,
                      add.params = list(color = "blue", fill = "gray"),
                      cor.coef = TRUE, cor.method = method,cor.coef.size=4.5,
                      cor.coeff.args = list(label.y=max(cordata[,1])+describe(cordata[,1])$se,label.sep = ","),
                      xlab =paste("Normalized expression of",name[2],sep=" ") , ylab = paste("Normalized expression of",name[1],sep=" ") )
 cort<-cor.test(cordata[,1],cordata[,2],method=method)
 res<-list(scatter,cort)
})
observe({
  if(is.null(input$gene1) || input$gene1 == ""||is.null(input$gene2) || input$gene2 == ""){
    disable("genecorbt")
  }
  else{
    enable("genecorbt")
  }
})


observeEvent(input$genecorbt, {
  output$genecorploting  <- renderPlot({
    ggcor()[[1]]
  })
})

observeEvent(input$genecorbt, {
  output$genecorplot<- renderUI({
    plotOutput("genecorploting",
               width = paste0(isolate({input$genecorwidth}), "%"),
               height = isolate({input$genecorheight}))
  })})

observeEvent(input$genecorbt, {
  output$genecorsummary <-  renderPrint(
    {ggcor()[[2]]}
  )
})

observeEvent(input$genecorbt, {
  shinyjs::show("genecor_wrapper")
})

output$downgenecor <- downloadHandler(

  filename = function() {
    paste0("Gene-gene-correlation",Sys.Date(),'.pdf')
  },
  content = function(file) {

    pdf(file,width = isolate({input$genecorwidthdl}),height = isolate({input$genecorheightdl}))
    print(ggcor()[[1]])
    dev.off()
  }
)

observeEvent(input$page_before_gene, {
  newtab <- switch(input$tabs, "mcorgene" = "gene","gene" = "mcorgene")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_gene, {
  newtab <- switch(input$tabs, "gene" = "genediff","genediff" = "gene")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})







