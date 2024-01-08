output$dataset = DT::renderDataTable(dataset, options=list(scrollX=TRUE),server = FALSE)

selected.rows <- reactiveValues(index = NULL)
observe({
  selected.rows$index <- input$dataset_rows_selected
})

metarun<-eventReactive(input$metabt,
                       withProgress(
                         message = "Inputting the data",
                         detail = "It may take a while, please be patient",
                         value = 5,{
                           library(meta)
                           library(Biobase)
                           library(CuratedCancerPrognosisData)
                           # if(input$inputType == "Public dataset"){
                           s=selected.rows$index
                           s<-as.numeric(s)
                           # print(s)
                           # print(mode(s))
                           dataset<-dataset$Dataset[s]
                           # dataset<-DatasetID()
                           dataset<-as.vector(dataset)
                           # print(mode(dataset))
                           # print(dataset)
                           gene<-isolate({input$metamarker})
                           time<-isolate({input$metatime})
                           status<-isolate({input$metastatus})
                           # print(dataset)
                           # dat<-list()
                           # for(i in 1:length(dataset)){
                           #   data(list=dataset[i])
                           #   print(dataset[i])
                           #   dat[[i]]<-get(dataset[i])
                           # }
                           # names(dat)<-dataset
                           dat <- lapply(dataset, function(x) {
                             data(list=x)
                             get(x)
                           })
                           
                           names(dat) <- dataset
                           
                           note1<-c()
                           for(i in 1:length(dat)){
                             if(gene %in% row.names(exprs(dat[[i]]))==F){
                               note1[i]<-dataset[i]
                             }}
                           note2<-c()
                           for(i in 1:length(dat)){
                             if(!all(c(time,status) %in% names(pData(dat[[i]])))){
                               note2[i]<-dataset[i]
                             }
                           }

                           if(length(note1)!=0 & length(note2)!=0){
                             createAlert(
                               session,
                               "metamess",
                               "exampleAlert",
                               title = "Please note!",
                               style =  "info",
                               content = paste("Please note! Either the", time ,"or the", status,"was not found in the clinical data of the dataset:", paste(na.omit(note2),collapse = ", "),", and",gene, "is not found in dataset",paste(na.omit(note1),collapse = ", ")),
                               append = T,
                               dismiss=T
                             )
                           }
                           if(length(note1)!=0 & length(note2)==0){
                             createAlert(
                               session,
                               "metamess",
                               "exampleAlert",
                               title = "Please note!",
                               style =  "info",
                               content = paste("Please note!",gene, "is not found in dataset",paste(na.omit(note1),collapse = ", ")),
                               append = T,
                               dismiss=T
                             )
                           }
                           if(length(note1)==0 & length(note2)!=0){
                             createAlert(
                               session,
                               "metamess",
                               "exampleAlert",
                               title = "Please note!",
                               style =  "info",
                               content = paste("Either the", time ,"or the", status,"was not found in the clinical data of the dataset:", paste(na.omit(note2), collapse= ", "))
,
                               append = T,
                               dismiss=T
                             )
                           }
                           metagene(dat=dat,gene=gene,time=time,status=status)
                           })
)


observe({
  if(
    is.null(selected.rows$index)||
    # length(DatasetID())==0||
    is.null(input$metamarker) ||
    input$metamarker == ""||
    is.null(input$metatime)||
    input$metatime == ""  ||
    is.null(input$metastatus)||
    input$metastatus == ""  
    
  ){
    disable("metabt")
  }
  else{
    enable("metabt")
  }
})
observeEvent(input$metabt, {
  disable("metabt")
})

observeEvent(input$metabt, {
  updateCollapse(session, "collapsedatainputmeta", open = "Meta-analysis of biomarker", close = "Select dataset for meta-analysis")
})

observeEvent(metarun(), {
  shinyjs::show("meta_wrapper")
  # shinyjs::show("metatable_wrapper")
  enable("metabt")
})

# observeEvent(input$metabt, {
#   output$metatable <-  DT::renderDT({
#     DT::datatable(metarun()$table,options=list(scrollX=TRUE))
#   })
# })

# output$savemetatable <- downloadHandler(
#   filename = function() {
#     paste0("Cytotoxic-activity-table-",Sys.Date(), ".csv")
#   },
#   content = function(file) {
#     write.csv(metarun()$table, file)
#   }
# )


observeEvent(input$metabt,{
  output$metaploting  <- renderPlot({
    closeAlert(session, "metamess")
    forest(metarun(),fixed=F)
  })
})

observeEvent(input$metabt, {
  output$metaplot<- renderUI({
    plotOutput("metaploting",
               width = paste0(isolate({input$metawidth}), "%"),
               height = isolate({input$metaheight}))
  })})

output$downloadmeta <- downloadHandler(
  filename = function() {
    paste0("Forest-plot-of-meta-analysis-result-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$metawidthdl}),height = isolate({input$metaheightdl}))
    print(forest(metarun(),fixed=F))
    dev.off()
  }
)

observeEvent(input$page_before_meta, {
  newtab <- switch(input$tabs, "meta" = "bioan","bioan" = "meta")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_meta, {
  newtab <- switch(input$tabs, "cba" = "meta","meta" = "cba")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









