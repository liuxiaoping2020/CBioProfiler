options(shiny.maxRequestSize=10000*1024^2)
cancer <- reactive({
  filter(dataset,Cancer == input$cancer)
})
observeEvent(cancer(), {
  choices <- cancer()$Dataset
  updateSelectInput(session, "accession", choices = choices,selected=NULL)
})


Rawdata<-eventReactive(input$data,
                       withProgress(
                         message = "Inputting the data",
                         detail = "It may take a while, please be patient",
                         value = 5,{
                       if(input$inputType == "Public dataset"){
                         accession<-isolate({input$accession})
                         data(list=accession)
                         dat<-get(accession)
                         clinical<-pData(dat)
                         clinical[which(lapply(as.data.frame(clinical),smna)==nrow(clinical))]<-NULL
                         expres<-as.data.frame(exprs(dat))
                         row.names(expres)<-make.names(row.names(expres))
                         dat<-list(expres=expres,clinical=clinical)
                         dat$clinical[dat$clinical == ""]<- NA

                         dat
                       } else{
                         req(input$expfile)
                         expres <- read.csv(input$expfile$datapath,
                                            header = T,row.names=1)
                         row.names(expres)<-make.names(row.names(expres))
                         req(input$clin)
                         clinical<- read.csv(input$clin$datapath,
                                             header = T,row.names=1)
                         clinical[which(lapply(as.data.frame(clinical),smna)==nrow(clinical))]<-NULL
                         dat<-list(expres,clinical)
                         names(dat)<-c("expres","clinical")
                         dat$clinical[dat$clinical == ""]<- NA
                         dat
                       }})
)


observeEvent(input$data, {
  disable("data")
})

observeEvent(Rawdata(), {
  shinyjs::show("clinicalinput.table")
  shinyjs::show("express.table")
  enable("data")
})

observeEvent(input$data,{
  output$viewclinical <-  DT::renderDT({
    DT::datatable(Rawdata()$clinical,options=list(scrollX=TRUE))
  })

})

output$saveclinicalinput <- downloadHandler(
  filename = function() {
    paste0(isolate({input$accession}),"-clinical-data-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(Rawdata()$clinical, file)
  }
)

observeEvent(input$data,{
  output$viewexpres <-  DT::renderDT({
    DT::datatable(Rawdata()$expres,options=list(scrollX=TRUE))
  })
})

output$saveexpresinput <- downloadHandler(
  filename = function() {
    paste0(isolate({input$accession}),"-gene-expression-data-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(Rawdata()$expres, file)
  }
)


###############################################################################################################################
cancer2 <- reactive({
  filter(dataset,Cancer == input$cancer2)
})
observeEvent(cancer2(), {
  choices2 <- cancer2()$Dataset
  updateSelectInput(session, "accession2", choices = choices2,selected=NULL)
})

Rawdata2<-eventReactive(input$data2,
                        withProgress(
                          message = "Inputting the data",
                          detail = "It may take a while, please be patient",
                          value = 5,{
                        if(input$inputType2 == "Public dataset"){
                          accession2<-isolate({input$accession2})
                          data(list=accession2)
                          dat2<-get(accession2)
                          clinical2<-pData(dat2)
                          clinical2[which(lapply(as.data.frame(clinical2),smna)==nrow(clinical2))]<-NULL

                          expres2<-as.data.frame(exprs(dat2))
                          row.names(expres2)<-make.names(row.names(expres2))
                          dat2<-list(expres=expres2,clinical=clinical2)
                          dat2$clinical[dat2$clinical == ""]<- NA
                          dat2
                        } else{
                          req(input$expfile2)
                          expres2 <- read.csv(input$expfile2$datapath,
                                              header = T,row.names=1)
                          row.names(expres2)<-make.names(row.names(expres2))
                          req(input$clin2)
                          clinical2<- read.csv(input$clin2$datapath,
                                               header = T,row.names=1)
                          clinical2[which(lapply(as.data.frame(clinical2),smna)==nrow(clinical2))]<-NULL
                          dat2<-list(expres2,clinical2)
                          names(dat2)<-c("expres","clinical")
                          dat2$clinical[dat2$clinical == ""]<- NA
                          dat2
                        }})
)

observeEvent(input$data2, {
  disable("data2")
})

observeEvent(Rawdata2(), {
  shinyjs::show("viewclinical2_wrapper")
  shinyjs::show("viewexpress2_wrapper")
  enable("data2")
})

observeEvent(input$data2,{
  output$viewclinical2 <-  DT::renderDT({
    DT::datatable(Rawdata2()$clinical,options=list(scrollX=TRUE))
  })
})

output$saveviewclinical2 <- downloadHandler(
  filename = function() {
    paste0("Clinical-data-of-validation-set-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(Rawdata2()$clinical, file)
  }
)

observeEvent(input$data2,{
  output$viewexpress2 <-DT::renderDT({
    DT::datatable(Rawdata2()$expres,options=list(scrollX=TRUE))
  })
})

output$saveviewexpress2 <- downloadHandler(
  filename = function() {
    paste0("Gene-expression-data-of-validation-set-",Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(Rawdata2()$expres, file)
  }
)

inpt <- reactive({
  req(input$data) && req(input$data2)
})

raw<-eventReactive(inpt(),{
  req(isolate({Rawdata()}))
  req(isolate({Rawdata2()}))

  index<-intersect(row.names(Rawdata()$expres),row.names(Rawdata2()$expres))
if(length(index)==0){
  createAlert(
    session,
    "inputmess",
    "exampleAlert",
    title = "Please note!",
    style =  "danger",
    content = paste(
      "No common genes found between discovery set and external validation set regarding the gene expression data, please check !!!"
    ),
    append = T,
    dismiss=T
  )
  return(NULL)
}else{
  Rawexpres<-Rawdata()$expres
  Raw2expres<-Rawdata2()$expres
  Rawclin<-Rawdata()$clinical
  Raw2clin<-Rawdata2()$clinical
  Rawexpres<-Rawexpres[index,]
  Raw2expres<-Raw2expres[index,]
  rawdata<-list(expres=Rawexpres,clinical= Rawclin)
  rawdata2<-list(expres=Raw2expres,clinical= Raw2clin)

  list(rawdata=rawdata,rawdata2=rawdata2)
}
  })


rawdata<-eventReactive(inpt(),
  {raw()$rawdata}
)

rawdata2<-eventReactive(inpt(),
 { raw()$rawdata2}
)
###################################################################################################################################

observeEvent(input$page_after_introduction, {
  newtab <- switch(input$tabs, "welcome1" = "dataset","dataset" = "welcome1")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_before_datainput, {
  newtab <- switch(input$tabs, "dataset" = "welcome1","welcome1" = "dataset")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_datainput, {
  newtab <- switch(input$tabs, "dataset" = "wgcna","wgcna" = "dataset")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})
  














