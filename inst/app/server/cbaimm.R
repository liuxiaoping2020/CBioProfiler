
value <- reactiveValues(data = NULL)

observeEvent(input$cbavalidationbt, {
  value$data <- isolate({validst()})
})

# observe({
#   if(!is.null(isolate({validst()}))){
#     enable("cbaimmunebt")
#   }else{
#     disable("cbaimmunebt")
#   }
# })

# observe({
#   if(is.null(value$data)){
#     shinyjs::disable("cbaimmunebt")
#   }else{
#     shinyjs::enable("cbaimmunebt")
#   }
# })

cbaimmrun<-eventReactive(input$cbaimmunebt,{
  input$cbaimmunebt
  res<-isolate({validst()})
  cohort<-isolate({input$cbaimmcor})
  if(cohort=="Training set"){
    data<-res$traindata
  }else{
    data<-res$validata
  }
  
  expres<-data$expres
  clinical<-data$clinical
  withProgress(message = "Calculating immune cell composition",
               detail = "This may take a while...",
               value = 3,
               {
                 geneset<-isolate({input$cbageneset})
                 require(ConsensusTME)
                 if(geneset=="Bindea"){
                   idx<-intersect(row.names(expres),as.vector(unlist(methodSignatures$Bindea)))
                 }else if(geneset=="Danaher"){
                   idx<-intersect(row.names(expres),as.vector(unlist(methodSignatures$Danaher)))
                 }else if(geneset=="Davoli"){
                   idx<-intersect(row.names(expres),as.vector(unlist(methodSignatures$Davoli)))
                 }else if(geneset=="MCP.Counter"){
                   idx<-intersect(row.names(expres),as.vector(unlist(methodSignatures$MCP.Counter)))
                 }else if(geneset=="xCell"){
                   idx<-intersect(row.names(expres),as.vector(unlist(methodSignatures$xCell)))
                 }
                 if(length(idx)==0){
                   createAlert(
                     session,
                     "cbaimmunemess",
                     "exampleAlert",
                     title = "Please note!",
                     style =  "danger",
                     content = paste(
                       "No genes in the reference immune geneset could be matched to the identifiers in the expression data you selected, please Check!!"
                     ),
                     append = T,
                     dismiss=T
                   )
                   return(NULL)
                 } else {
                   score<-ims(exp=as.matrix(expres),geneset=geneset)
                   score<-as.data.frame(t(score))

                   index<-intersect(row.names(clinical),row.names(score))
                   name<-gsub("_"," ",names(score))

                   clinical<-clinical[index,]
                   score<-score[index,]
                   score$Subtype<-clinical$Subtype
                   compare<-as.data.frame(combn(unique(as.vector(clinical$Subtype)),2))
                   comp<-list()
                   for(i in 1:length(compare)){
                     comp[[i]]<-compare[[i]]
                   }
                   barp(data = score, 
                        paired = isolate({input$cbaimmpair}), 
                        plottype = isolate({input$cbaimmplot}),
                        method =isolate({input$cbaimmeth}),
                        comp=comp,
                        name=name
                        )
                   
                   #ggarrange(plotlist=pp,common.legend = T,align = "hv",labels="AUTO")
                 }
          })
})



observe({
  if(
    # is.null(input$cbavalidationbt)||
    # input$cbavalidationbt||
    is.null(input$cbageneset)||
    input$cbageneset == ""   ||
    is.null(input$cbaimmeth) ||
    input$cbaimmeth == ""    ||
    is.null(input$cbaimmpair)||
    input$cbaimmpair == ""   ||
    is.null(input$cbaimmplot)||
    input$cbaimmplot==""
  ){
    disable("cbaimmunebt")
  }
  else{
    enable("cbaimmunebt")
  }
})

observeEvent(input$cbaimmunebt, {
  disable("cbaimmunebt")
})

observeEvent(input$cbaimmunebt, {
  updateCollapse(session, "cbacollapseimmune",
                 open = "Immune cell infiltration distribution accross different subtypes", 
                 close = "Descriptions and parameters")
})

observe({
  if(is.null(value$data)){
    disable("cbaimmunebt")
  }else{
    enable("cbaimmunebt")
  }
})


observeEvent(cbaimmrun(), {
  shinyjs::show("cbaimmuneplot_wrapper")
  enable("cbaimmunebt")
})

observeEvent(is.null(cbaimmrun()), {
  enable("cbaimmunebt")
})

observeEvent(input$cbaimmunebt,{
  output$cbaimmuneploting  <- renderPlot({
    closeAlert(session, "cbaimmunemess")
    #cbaimmrun()
    ggarrange(plotlist=cbaimmrun(),common.legend = T,align = "hv",labels="AUTO")
  })
})

observeEvent(input$cbaimmunebt, {
  output$cbaimmuneplot<- renderUI({
    plotOutput("cbaimmuneploting",
               width = paste0(isolate({input$cbaimmwidth}), "%"),
               height = isolate({input$cbaimmheight}))
  })})

output$cbasaveimune <- downloadHandler(
  filename = function() {
    paste0("Immune-correlation-plot-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$cbaimmpltwidth}),height = isolate({input$cbaimmpltheight}))
    #print(cbaimmrun())
    print(ggarrange(plotlist=cbaimmrun(),common.legend = T,align = "hv",labels="AUTO"))
    dev.off()
  }
)

observeEvent(input$page_before_cbaimmune, {
  newtab <- switch(input$tabs, "cbaDEG" = "cbaimmune","cbaimmune" = "cbaDEG")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_cbaimmune, {
  newtab <- switch(input$tabs, "cbaimmune" = "cbastemness","cbastemness" = "cbaimmune")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})








