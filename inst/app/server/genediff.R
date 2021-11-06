observeEvent(rawdata(),{
  data <- rawdata()$expres
  updateSelectizeInput(session,
                       'genediff',
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


observeEvent(rawdata(),{
  data <- rawdata()$clinical
  updateSelectizeInput(session, "genediffgroup", choices = {
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
   selected = "Grade"
  )
})


observeEvent(rawdata(),{
  data <- rawdata()$clinical
  updateSelectizeInput(session, "genediffsub", choices = {
    if (!is.null(data)) {
      if (class(data) ==  class(data.frame()))
        choices <- as.character(colnames(data))
      if (class(data) !=  class(data.frame()))
        choices <-
          as.character(colnames(as.data.frame(data)))
      choices
    }
  },server = TRUE,
  selected = NULL
  )
})

glevel<-reactive({
  req(input$genediffgroup)
  data <- rawdata()$clinical
  data<-as.data.frame(data)

  group<-paste(input$genediffgroup)
  g<-as.factor(data[,group])
  unique(g)

})

observeEvent(glevel(),{
  updateSelectInput(session, "comparefer", choices = glevel())
})

genedifffun<-reactive({
  input$genediffbt
  data <- isolate({
    rawdata()
  })
  data <- lapply(data, as.data.frame)
  expres<-data$expres
  clinical<-data$clinical
  clinical[clinical == ""]<- NA
  group1 <-isolate({input$genediffgroup})
  subgroup<-isolate({input$genediffsub})

  clinical<-clinical[complete.cases(clinical[,group1]),]
  index<-intersect(colnames(expres),rownames(clinical))
  expres<-expres[,index]
  clinical<-clinical[index,]
  gene1 <-isolate({input$genediff})

  gene<-as.numeric(expres[gene1,])
  group<-clinical[,group1]
  da<-as.data.frame(cbind(gene,group))
  da$gene<-as.numeric(da$gene)
  da$group<-as.factor(da$group)
  subgroup<-clinical[,subgroup]
  da$subgroup<-as.factor(subgroup)

  if(length(levels(factor(da$group)))==1){
    createAlert(
      session,
      "genediffmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The group variable you selected has only one level (only one value), which is not suitable for T/wilcox test. Please check!",
      append = T,
      dismiss = T
    )
    return(NULL)
  }else if(length(levels(factor(da$subgroup)))==1){
    createAlert(
      session,
      "genediffmess",
      "exampleAlert",
      title = "Please note!",
      style =  "danger",
      content = "The subgroup variable you selected has only one level (only one value), which is not suitable for T/wilcox test. Please check!",
      append = T,
      dismiss = T
    )
    return(NULL)
  }else{

  genediff(
    da=da,
    x="group",
    y="gene",
    subgroup=isolate({input$gdsubgroup}),
    plottype=isolate({input$genediffplottype}),
    subgroupname="subgroup",
    xlab=group1,
    ylab=paste("Relative expression of", gene1,sep=" "),
    genediffp=isolate({input$genediffp}),
    compar = isolate({input$comparisons}),
    genediffmethod=isolate({input$genediffmethod}),
    comparefer=isolate({input$comparefer})
  )
  }
}
)

observe({
  if(is.null(input$genediff) || input$genediff == ""||is.null(input$genediffgroup) || input$genediffgroup == "" ){
    disable("genediffbt")
  }
  else{
    enable("genediffbt")
  }
})



observeEvent(input$genediffbt, {
  shinyjs::show("genediff_wrapper")
})

observeEvent(input$genediffbt, {
  output$genediffplot <-  renderPlot(
    width = isolate({input$genediffbarwidth}),
    height = isolate({input$genediffbarheight}),
    res = 96,
    {genedifffun()[[1]]}
  )
})

observeEvent(input$genediffbt, {
  output$genediffploting  <- renderPlot({
    genedifffun()
  })
})

observeEvent(input$genediffbt, {
  output$genediffplot<- renderUI({
    plotOutput("genediffploting",
               width = paste0(isolate({input$genediffbarwidth}), "%"),
               height = isolate({input$genediffbarheight}))
  })})

output$downloadgenediff<- downloadHandler(

  filename = function() {
    paste0("Gene-expression-difference-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$genediffwidthdl}),height = isolate({input$genediffheightdl})) # open the pdf device
    print(genedifffun())
    dev.off()
  }
)

observeEvent(input$page_before_genediff, {
  newtab <- switch(input$tabs, "genediff" = "gene","gene" = "genediff")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_genediff, {
  newtab <- switch(input$tabs, "immune" = "genediff","genediff" = "immune")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})




