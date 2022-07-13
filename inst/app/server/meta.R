
metacancer <- reactive({
  dataset[dataset$Cancer==input$metacancer,]
})

observeEvent(metacancer(), {
  choices <- metacancer()$Dataset
  updateSelectInput(session, "metaaccession", choices = choices,selected=NULL)
})

metadata<-eventReactive(input$metadatabt,
                       withProgress(
                         message = "Inputting the data",
                         detail = "It may take a while, please be patient",
                         value = 5,{
                           # if(input$inputType == "Public dataset"){
                           dataset<-isolate({input$metaaccession})
                             # data(list=accession)
                             # dat<-get(accession)
                             # clinical<-pData(dat)
                             # clinical[which(lapply(as.data.frame(clinical),smna)==nrow(clinical))]<-NULL
                             # expres<-as.data.frame(exprs(dat))
                             # row.names(expres)<-make.names(row.names(expres))
                             # dat<-list(expres=expres,clinical=clinical)
                             # dat$clinical[dat$clinical == ""]<- NA
                             # 
                             # dat
                             dat<-list()
                             for(i in 1:length(dataset)){
                               data(list=dataset[i])
                               da<-get(dataset[i])
                               clinical<-pData(da)
                               clinical[which(lapply(as.data.frame(clinical),smna)==nrow(clinical))]<-NULL
                               expres<-as.data.frame(exprs(da))
                               row.names(expres)<-make.names(row.names(expres))
                               dat[[i]]<-list(expres=expres,clinical=clinical)
                               dat[[i]]$clinical[dat[[i]]$clinical == ""]<- NA
                               
                             }
                             names(dat)<-dataset
                             dat
                           # } 
                           # else{
                           #   req(input$expfile)
                           #   expres <- read.csv(input$expfile$datapath,
                           #                      header = T,row.names=1)
                           #   row.names(expres)<-make.names(row.names(expres))
                           #   req(input$clin)
                           #   clinical<- read.csv(input$clin$datapath,
                           #                       header = T,row.names=1)
                           #   clinical[which(lapply(as.data.frame(clinical),smna)==nrow(clinical))]<-NULL
                           #   dat<-list(expres,clinical)
                           #   names(dat)<-c("expres","clinical")
                           #   dat$clinical[dat$clinical == ""]<- NA
                           #   dat
                           # }
                           
                           })
)


observe({
  data <- metadata()[[1]]$expres
  print(dim(data))
  updateSelectizeInput(session,
                       'metamarker',
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
                       selected = "TP53")
})

metarun<-eventReactive(input$metabt,{
  input$metabt
  dataset<-isolate({input$metaaccession})
  dat<-isolate({metadata()})
  print(dataset)
  gene<-isolate({input$metamarker})
  time<-isolate({input$metatime})
  status<-isolate({input$metastatus})
  
  check<-list()
  for(i in 1:length(dataset)){
    name<-row.names(exprs(dat[[i]]))
    check[[i]]<-gene %in% name
  }
  check<-unlist(check)

  
  withProgress(message = "Performing meta-analysis of the biomarker you selected",
               detail = "This may take a while...",
               value = 3,
               {
                 if(length(check==T)==1){
                     createAlert(
                       session,
                       "metamess",
                       "exampleAlert",
                       title = "Please note!",
                       style =  "danger",
                       content = paste(gene,"only exists in one gene expression dataset, which does not meet the requirements of meta-analysis, please check!!!"),
                       append = T,
                       dismiss=T
                     )
                   return(NULL)
                   # stop(paste(gene,"only exists in one gene expression dataset, which does not meet the requirements of meta-analysis, please check!!!"))
                 }else if(any(F %in% check )){
                   createAlert(
                     session,
                     "metamess",
                     "exampleAlert",
                     title = "Please note!",
                     style =  "info",
                     content = paste("Please note!",gene,"is not found in the gene expression profile of", paste(names(dat)[which(check==F)])),
                     append = T,
                     dismiss=T
                   )
                   # warning(paste("Please note!",gene,"is not found in the gene expression profile of", paste(names(dat)[which(check==F)])))
                   func<-function(dt,gene){
                     if(gene %in% row.names(exprs(dt))==F){
                       stop(paste(gene, "is not found in",paste0(dataset[i],","), " please check!!!"))
                     }else{
                       expres<-exprs(dt)
                       expres<-as.data.frame(expres)
                       clinical<-pData(dt)
                       clinical<-clinical[complete.cases(clinical[,time]) & clinical[,time]>0,]
                       index<-intersect(row.names(clinical),names(expres))
                       clinical<-clinical[index,]
                       expres<-expres[,index]
                       clinical$gene<-as.numeric(expres[gene,])
                       ph<-coxph(Surv(clinical[[time]], clinical[[status]]) ~ gene, clinical) 
                       ph<-summary(ph)
                       coeffs <- c(ph$coefficients[1,1:2], ph$conf.int[1,3:4], ph$coefficients[1,5])
                       names(coeffs)<-c('coef','HR','lower95','upper95','pValue')
                     }
                     return(coeffs)
                   }
                   res<-lapply(dat,func,gene=gene)
                   res<-do.call(rbind,res)
                   res<-as.data.frame(res)
                   gen<-function(es,cil,ciu,studlab){
                     require(meta)
                     TE<-log(es)
                     seTE<-(log(ciu)-log(cil))/1.96/2
                     meta<-metagen(TE,seTE,sm="HR",studlab=studlab)
                     return(meta)
                   }  
                   meta<-gen(es=res$HR,cil=res$lower95,ciu=res$upper95,studlab=row.names(res)) 
                 }else {
                   func<-function(dt,gene){
                     if(gene %in% row.names(exprs(dt))==F){
                       stop(paste(gene, "is not found in",paste0(dataset[i],","), " please check!!!"))
                     }else{
                       expres<-exprs(dt)
                       expres<-as.data.frame(expres)
                       clinical<-pData(dt)
                       clinical<-clinical[complete.cases(clinical[,time]) & clinical[,time]>0,]
                       index<-intersect(row.names(clinical),names(expres))
                       clinical<-clinical[index,]
                       expres<-expres[,index]
                       clinical$gene<-as.numeric(expres[gene,])
                       ph<-coxph(Surv(clinical[[time]], clinical[[status]]) ~ gene, clinical) 
                       ph<-summary(ph)
                       coeffs <- c(ph$coefficients[1,1:2], ph$conf.int[1,3:4], ph$coefficients[1,5])
                       names(coeffs)<-c('coef','HR','lower95','upper95','pValue')
                     }
                     return(coeffs)
                   }
                   res<-lapply(dat,func,gene=gene)
                   res<-do.call(rbind,res)
                   res<-as.data.frame(res)
                   gen<-function(es,cil,ciu,studlab){
                     require(meta)
                     TE<-log(es)
                     seTE<-(log(ciu)-log(cil))/1.96/2
                     meta<-metagen(TE,seTE,sm="HR",studlab=studlab)
                     return(meta)
                   }  
                   meta<-gen(es=res$HR,cil=res$lower95,ciu=res$upper95,studlab=row.names(res)) 
                 }
                 meta
                 # gmgene<-c("GZMA","PRF1")
                 # if(length(setdiff(gmgene,row.names(expres)))!=0){
                 #   createAlert(
                 #     session,
                 #     "metamess",
                 #     "exampleAlert",
                 #     title = "Please note!",
                 #     style =  "danger",
                 #     content = paste(setdiff(gmgene),names(expres),"is not found in the gene expression profile"),
                 #     append = T,
                 #     dismiss=T
                 #   )
                 #   return(NULL)
                 # }else{
                 #   cormeta(expres=expres,gene=gene, method=method)
                 # }
                 
               })
})


observe({
  if(
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
  updateCollapse(session, "collapsedatainputmeta", open = "Meta-analysis of biomarker", close = "Descriptions and parameters")
})



observeEvent(metarun(), {
  shinyjs::show("meta_wrapper")
  # shinyjs::show("metatable_wrapper")
})

# observeEvent(input$metabt, {
#   output$metatable <-  DT::renderDT({
#     DT::datatable(metarun()$table,options=list(scrollX=TRUE))
#   })
# })
# 
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
    forest(metarun())
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
    paste0("Forest-plot-of-meta-analysis-",Sys.Date(),'.pdf')
  },
  content = function(file) {
    pdf(file,width = isolate({input$metawidthdl}),height = isolate({input$metaheightdl}))
    print(forest(metarun()))
    dev.off()
  }
)

observeEvent(input$page_before_metaness, {
  newtab <- switch(input$tabs, "meta" = "stemness","stemness" = "meta")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_metaness, {
  newtab <- switch(input$tabs, "bioan" = "cbameta","cbameta" = "bioan")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})









