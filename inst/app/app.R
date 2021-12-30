suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyBS)
  library(table1)
  library(survminer)
  library(survival)
  library(DT)
  library(psych)
  library(survivalROC)
  library(tidyverse)
  library(WGCNA)
  library(ggpubr)
  library(limma)
  library(colourpicker)
  library(randomForestSRC)
  library(mlr)
  library(shinyjs)
  library(enrichplot)
  library(rms)
  library(TCGAbiolinks)
  library(ReactomePA)
  library(clusterProfiler)
  library(msigdbr)
  library(caret)
  library(Biobase)
  library(CuratedCancerPrognosisData)
  library("ggplotify")
  library(ConsensusTME)
})
sumz<-function(x){
  p<-sum(na.omit(x)==0)/length(x)
  return(p)
}
rmz<-function(x,per){
  idx<-apply(x,1,sumz)
  
  index<-which(idx>per)
  if(length(index)>0){
    x<-x[-index,]
  }else{
    x
  }
  
  return(x)
}

stem<-function(exp,gene,method){
  library(TCGAbiolinks)
  exp<-as.data.frame(exp)
  stemness <- TCGAanalyze_Stemness(stemSig = PCBC_stemSig,
                                   dataGE = exp)
  stemness$gene<-as.numeric(exp[gene,])
  stemness<-subset(stemness,select=c(StemnessScore,gene  ))
  scatter<-ggscatter(data=stemness, y = "StemnessScore",x = "gene",
                     size = 1,
                     add = "reg.line", conf.int = TRUE,
                     add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                     cor.coef = TRUE, cor.method = method,cor.coef.size=4.5,
                     cor.coeff.args = list(label.y=max(stemness[,1])+describe(stemness[,1])$se,label.sep = ","),
                     xlab =paste("Normalized expression of",gene,sep=" ") , ylab = "Stemness score" )
  cort<-cor.test(stemness$StemnessScore,stemness$gene,method=method)
  names(stemness)[2]<-gene
  res<-list(scatter=scatter,cort=cort,stemnesstable=stemness)
  return(res)
}

ims<-function(exp,geneset){
  require(ConsensusTME)
  if(geneset=="Bindea"){
    score<-geneSetEnrichment(exp,methodSignatures$Bindea)
  }else if(geneset=="Danaher"){
    score<-geneSetEnrichment(exp,methodSignatures$Danaher)
  }else if(geneset=="Davoli"){
    score<-geneSetEnrichment(exp,methodSignatures$Davoli)
  }else if(geneset=="MCP.Counter"){
    score<-geneSetEnrichment(exp,methodSignatures$MCP.Counter)
  }else if(geneset=="xCell"){
    score<-geneSetEnrichment(exp,methodSignatures$xCell)
  }
  return(as.data.frame(score))
}


plotcor<-function(score,exp,gene,select,normmethod,fdrcutoff, type,plottype){
  CorMatrix <- function(cor,p) {
    ut <- upper.tri(cor)
    data.frame(row = rownames(cor)[row(cor)[ut]] ,
               column = rownames(cor)[col(cor)[ut]],
               cor =(cor)[ut],
               p = p[ut] ) }
  score<-as.data.frame(t(score))
  data<-cbind(as.numeric(exp[gene,]),score)
  names(data)[1]<-"idx"
  require(Hmisc)
  res <- rcorr(as.matrix(data))
  res<-CorMatrix (res$r, res$P)
  res<-res[which(res$row=="idx"),]
  res<-na.omit(res)
  res$Adjusted.P<-p.adjust(res$p,method=normmethod)
  if(select==T){
    res<-res[res$Adjusted.P<fdrcutoff,]
  }
  
  if(type=="pearson"){
    ylab<-"Pearson's correlation coefficient"
  }else{
    ylab<-"Spearman's correlation coefficient"
  }
  require(ggplot2)
  if(plottype=="bar plot"){
    correplot<-ggplot(data=res)+geom_bar(aes(x=column,y=cor, fill=Adjusted.P), stat='identity')+
      coord_flip() +scale_fill_gradient(low="red", high = "blue") +xlab("")+
      ylab(ylab)+
      theme(axis.text.x=element_text(color="black",size=48),axis.text.y=element_text(color="black", size=48))+
      scale_y_continuous(expand=c(0, 0))+
      scale_x_discrete(expand=c(0,0))+theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }else{
    correplot<-ggplot(data=res,aes(x=cor,y=column))+geom_point(aes(col=cor,size=Adjusted.P))+
      scale_color_gradient(low="blue",high="red")+
      labs(color="Correlations", size="Adjusted.P")+theme_bw()+
      ylab("")+xlab(ylab)+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),
            axis.text=element_text(size=12,color="black"),legend.text=element_text(size=10))
  }
  res<-subset(res,select=c( column ,       cor   ,         p  , Adjusted.P))
  names(res)[1]<-"Cell type"
  list(cordata=res,corplot=correplot,score=score)
}

 
nomo<-function(train,time,status,variable,varlabels,yearpoint){
  train<-train[complete.cases(train[,time]) & train[,time]>0,]
  if(max(train[,time])>600){
    by<- 365
    time.inc=1095
  } else{
    time.inc=36
    by<-12
  }
  years<-seq(from=0,to=max(train[,time]),by=by)
  years<-years[-1]
  
  train$time<-train[,time]
  train$status<-train[,status]
  train<-train[,c("time","status",variable)]
  varlabels<-c("Time","Status",varlabels)
  for(i in 1:length(varlabels)){
    label(train[[i]]) <-varlabels[i]
  }
  dd<<-datadist(train)
  options(datadist="dd")
  Srv = Surv(train$time, train$status)
  var<-paste(variable,collapse ="+")
  fomu<-as.formula(paste("Srv~",var,seq=""))
  f2 <- cph(fomu, surv=T,data = train,x=T, y=T,time.inc=time.inc)
  med <- Quantile(f2)
  surv <- Survival(f2)
  yearid<-as.numeric(gsub("-year","",yearpoint))
  times<-years[yearid]
  if(length(yearid)==1) {
    fun <- list(function(x)
      surv(times[1], x))
    funlabel <- paste(yearpoint[1], "Survival Probability")
  } else if (length(yearid) == 2) {
    fun <- list(function(x)
      surv(times[1], x),
      function(x)
        surv(times[2], x))
    funlabel <-
      c(
        paste(yearpoint[1], "Survival Probability"),
        paste(yearpoint[2], "Survival Probability")
      )
  }else if (length(yearid) == 3) {
    fun <- list(function(x)
      surv(times[1], x),
      function(x)
        surv(times[2], x),
      function(x)
        surv(times[3], x))
    funlabel <-
      c(
        paste(yearpoint[1], "Survival Probability"),
        paste(yearpoint[2], "Survival Probability"),
        paste(yearpoint[3], "Survival Probability")
      )
  }
  nom <<- nomogram(f2, fun=fun,funlabel=funlabel)
  return(nom)
}

geom_effect <- function(mapping = NULL,
                        data = NULL,
                        stat = "identity",
                        position = ggstance::position_dodgev(height = 0.5),
                        ...,
                        fatten = 2,
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomEffect,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      fatten = fatten,
      na.rm = na.rm,
      ...
    )
  )
}

GeomEffect <- ggproto("GeomEffect", Geom,
                      default_aes = aes(
                        colour = "black", size = 0.6, linetype = 1, shape = 21,
                        fill = "black", alpha = NA, stroke = 1, filled = TRUE
                      ),
                      
                      
                      draw_key = function(data, params, size) {
                        if (is.character(data$shape)) {
                          data$shape <- translate_shape_string(data$shape)
                        }
                        grid::pointsGrob(
                          0.5, 0.5,
                          pch = data$shape,
                          gp = grid::gpar(
                            col = scales::alpha(data$colour, data$alpha),
                            # fill = scales::alpha(data$fill, data$alpha),
                            fill = scales::alpha(data$colour, data$alpha),
                            fontsize = data$size * .pt + data$stroke * .stroke / 2,
                            lwd = data$stroke * .stroke / 2
                          )
                        )
                      },
                      
                      # TODO: guide for significance, when specified
                      
                      required_aes = c("x", "y", "xmin", "xmax"),
                      
                      # TODO: check for the number of shapes?
                      
                      draw_panel = function(data,
                                            panel_params,
                                            coord,
                                            fatten = 2) {
                        ggstance::GeomPointrangeh$draw_panel(
                          # @Ilari Indeed, the transform() is needed here for fatten to go through
                          transform(data, fatten = fatten) %>%
                            dplyr::mutate(
                              fill = dplyr::case_when(
                                is.na(.data$filled) ~ "#00000000",
                                !.data$filled ~ "white",
                                TRUE ~ .data$colour
                              )
                            ),
                          panel_params,
                          coord
                        )
                      }
)

geom_stripes <- function(mapping = NULL,
                         data = NULL,
                         stat = "identity",
                         position = "identity",
                         ...,
                         show.legend = NA,
                         inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomStripes,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(...)
  )
}

GeomStripes <- ggplot2::ggproto("GeomStripes", ggplot2::Geom,
                                required_aes = c("y"),
                                
                                default_aes = ggplot2::aes(
                                  xmin = -Inf, xmax = Inf,
                                  odd = "#22222222", even = "#00000000",
                                  alpha = NA, colour = "black", linetype = "solid", size = NA
                                ),
                                
                                
                                draw_key = ggplot2::draw_key_rect,
                                
                                draw_panel = function(data, panel_params, coord) {
                                  ggplot2::GeomRect$draw_panel(
                                    data %>%
                                      dplyr::mutate(
                                        y = round(.data$y),
                                        ymin = .data$y - 0.5,
                                        ymax = .data$y + 0.5
                                      ) %>%
                                      dplyr::select(
                                        .data$xmin, .data$xmax,
                                        .data$ymin, .data$ymax,
                                        .data$odd, .data$even,
                                        .data$alpha, .data$colour, .data$linetype, .data$size
                                      ) %>%
                                      unique() %>%
                                      dplyr::arrange(.data$ymin) %>%
                                      dplyr::mutate(
                                        .n = dplyr::row_number(),
                                        fill = dplyr::if_else(
                                          .data$.n %% 2L == 1L,
                                          true = .data$odd,
                                          false = .data$even
                                        )
                                      ) %>%
                                      dplyr::select(-.data$.n, -.data$odd, -.data$even),
                                    panel_params,
                                    coord
                                  )
                                }
)



internalvc<-function(split,time,status,variable,ratio,bootstrap,yearpoint,res ){
  if(split==T){
    data<-res$Trainsurv
    data<-data[complete.cases(data[,time]) & data[,time]>0,]
    if(max(data[,time])>600){
      by<- 365
      time.inc=1095
    } else{
      time.inc=36
      by<-12
    }
    years<-seq(from=0,to=max(data[,time]),by=by)
    years<-years[-1]
    years<-years[-length(years)]
    yearid<-as.numeric(gsub("-year","",yearpoint))
    pred.time<-years[yearid]
    resvalid<-valid(data=res$Trainsurv,
                    time=time,
                    status=status,
                    variable=variable,
                    ratio=ratio,
                    bootstrap=bootstrap,
                    time.inc=time.inc,
                    pred.time=pred.time
    )
    resvalidpred<-validpred(
      traindata = res$Trainsurv,
      testdata = res$Testsurv,
      time = time,
      status=status,
      variable=variable,
      ratio=ratio,
      bootstrap=bootstrap,
      time.inc=time.inc,
      pred.time=pred.time
    )
    train<-as.data.frame(resvalid$cidx)
    train$dataset<-rep("Training set",nrow(train))
    test<-as.data.frame(resvalidpred$cidx)
    test$dataset<-rep("Test set",nrow(test))
    da<-rbind(train,test)
    names(da)<-c("C-index","Dataset")
    p <-
      ggboxplot(
        da,
        x = "Dataset",
        y = "C-index",
        fill = "Dataset",
        xlab = NULL,
        ylab = "C-index",
        font.y = 12,
        font.xtickslab = 12,
        font.ytickslab = 12,
        add = "mean_sd"
      ) + theme(
        axis.title.x = element_blank(),
        legend.position = "none"
      )
    idex<-lapply(list(resvalid$cidx,resvalidpred$cidx),sumcidx)
    validation<-list(validtrain=resvalid,validtest=resvalidpred)
    result<-list(p=p,idex=idex,validation=validation)
    
  } else if(split==F){
    
    resvalid<-valid(data=res$Trainsurv,
                    time=time,
                    status=status,
                    variable=variable,
                    ratio=ratio,
                    bootstrap=bootstrap,
                    time.inc=time.inc,
                    pred.time=pred.time
    )
    train<-as.data.frame(resvalid$cidx)
    train$dataset<-rep("Training set",nrow(train))
    da<-train
    names(da)<-c("C-index","Dataset")
    p <-
      ggboxplot(
        da,
        x = "Dataset",
        y = "C-index",
        fill = "Dataset",
        xlab = NULL,
        ylab = "C-index",
        font.y = 12,
        font.xtickslab = 12,
        font.ytickslab = 12,
        add = "mean_sd"
      ) +theme(
        axis.title.x = element_blank(),
        legend.position = "none"
      )
    idex<-lapply(list(resvalid$cidx),sumcidx)
    validation<-list(validtrain=resvalid)
    result<-list(p=p,idex=idex,validation=validation)
  }
  return(result)
}
smna<-function(x){
  sum(is.na(x))
}

valid<-function(data,time,status,variable,ratio,bootstrap,time.inc,pred.time){
  require(rms)
  require(survival)
  if(any(variable %in% names(data))==F){
    stop("At least one variable you provided is not found in data")
  }
  data$time<-data[,time]
  data$status<-data[,status]
  data<-data[complete.cases(data$time)&data$time>0,]
  fomu<-as.formula(paste(paste("Surv(time,status)~"),paste(variable,collapse = "+")))
  cidx = cindex_boostrapping(data,fomu,ratio,bootstrap)
  dd<<-datadist(data)
  options(datadist="dd")
  f3 <- cph(fomu, surv=T,data = data,x=T, y=T,time.inc=time.inc)
  calibration<-list()
  for (i in 1:length(pred.time)){
    calibration[[i]] <- calibrate(f3, cmethod='KM', method="boot", u=pred.time[i], m=floor(nrow(data)/3), B=bootstrap)
  }
  res<-list(cidx=cidx,calibration=calibration)
  return(res)
}

validpred<-function(traindata,testdata,time,status,variable,ratio,bootstrap,time.inc,pred.time){
  traindata$time<-traindata[,time]
  traindata$status<-traindata[,status]
  traindata<-traindata[complete.cases(traindata$time)&traindata$time>0,]
  fomu<-as.formula(paste(paste("Surv(time,status)~"),paste(variable,collapse = "+")))
  model<-coxph(fomu,traindata)
  testdata$time<-testdata[,time]
  testdata$status<-testdata[,status]
  testdata<-testdata[complete.cases(testdata$time)&testdata$time>0,]
  pred<-predict(model,newdata=testdata,type="lp")
  testdata$lp<-pred
  fomu<-as.formula(paste(paste("Surv(time,status)~"),"lp"))
  cidx = cindex_boostrapping(testdata,fomu,ratio,bootstrap)
  dd<<-datadist(testdata)
  options(datadist="dd")
  f3 <- cph(fomu, surv=T,data = testdata,x=T, y=T,time.inc=time.inc)
  calibration<-list()
  for (i in 1:length(pred.time)){
    calibration[[i]] <- calibrate(f3, cmethod='KM', method="boot", u=pred.time[i], m=floor(nrow(testdata)/3), B=bootstrap)
  }
  res<-list(cidx=cidx,calibration=calibration)
  return(res)
}

cindex_boostrapping = function(inputdata, formula, sampling_ratio, iteration){
  
  sample_size = dim(inputdata)[1]
  sample_idx = c(1:sample_size)
  subset_sample_size = round(sample_size* sampling_ratio)
  cidx = matrix(NA,iteration,1)
  
  for( iter in c(1:iteration)){
    subset_sample_idx = sample(sample_idx, subset_sample_size, replace = FALSE, prob = NULL)
    
    subset_sample_data = inputdata[subset_sample_idx,]
    
    cresult = coxph(formula,data=subset_sample_data)
    
    cidx[iter,1] = cresult$concordance['concordance']
  }
  return(cidx)
}


sumcidx<-function(x){
  m<-mean(x)
  s<-sd(x)
  w<-paste0("C-index: ",round(m,3),", SD: ",round(s,3))
  return(w)
}


extvalid<-function(extdata,variable,res,ratio,boostrap,time,status,yearpoint){
  require(mlr)
  if(any(res$`Selected variables` %in% row.names(extdata$expres)) ==F ){
    stop("Please check! Not all genes derived from previous benchmark experiment are found in the gene expression data you provided!")
  }
  
  modeltrain<-res$`Fitted model`
  testset<-extdata$expres
  testsurv<-extdata$clinical
  testsurv$time<-testsurv[[time]]
  testsurv$status<-testsurv[[status]]
  testsurv<-testsurv[complete.cases(testsurv$time)&testsurv$time>0,]
  
  data<-testsurv
  data<-data[complete.cases(data[,time]) & data[,time]>0,]
  if(max(data[,time])>600){
    by<- 365
    time.inc=1095
  } else{
    time.inc=36
    by<-12
  }
  years<-seq(from=0,to=max(data[,time]),by=by)
  years<-years[-1]
  years<-years[-length(years)]
  yearid<-as.numeric(gsub("-year","",yearpoint))
  pred.time<-years[yearid]
  #
  index<-intersect(names(testset),row.names(testsurv))
  if(length(index)==0){
    stop("The sample names of the clinical data and gene expression data are not identical")
  }
  testset<-testset[,index]
  testsurv<-testsurv[index,]
  testset1<-testset
  testset<-testset[modeltrain$features,]
  testset<-as.data.frame(t(testset))
  testset1<-as.data.frame(t(testset1))
  testset$time<-testsurv$time
  testset$status<-testsurv$status
  names(testset)<-gsub("-","_",names(testset))
  names(testset)<-gsub("/","_",names(testset))
  names(testset)<-gsub(" ","_",names(testset))
  
  test.task <- makeSurvTask(data = testset, target = c("time","status"))
  
  if(res$lrnid=="RandomForestSRC"){
    model.predicted<-predict(modeltrain,testset1)
    test.chf <- model.predicted$chf
    testsurv$Risk<-rowSums(test.chf)
  } else if (res$lrnid=="Glmboost"){
    modelprediction.test <- predict(modeltrain, test.task)
    testsurv$Risk<-modelprediction.test$data$response
  } else if(res$lrnid=="ElasticNet"){
    modelprediction.test <- predict(modeltrain, test.task)
    testsurv$Risk<-modelprediction.test$data$response
    
  } else if(res$lrnid=="Ridge"){
    
    modelprediction.test <- predict(modeltrain, test.task)
    testsurv$Risk<-modelprediction.test$data$response
  }else if (res$lrnid=="Lasso"){
    modelprediction.test <- predict(modeltrain, test.task)
    testsurv$Risk<-modelprediction.test$data$response
  } else if(res$lrnid == "Coxboost"){
    modelprediction.test <- predict(modeltrain, test.task)
    testsurv$Risk<-modelprediction.test$data$response
  }
  resvalid<-validpred(
    traindata = res$Trainsurv,
    testdata = testsurv,
    time = time,
    status=status,
    variable=variable,
    ratio=ratio,
    bootstrap=boostrap,
    # time.inc=time.inc,
    pred.time=pred.time
  )
  external<-as.data.frame(resvalid$cidx)
  external$dataset<-rep("External set",nrow(external))
  da<-external
  names(da)<-c("C-index","Dataset")
  p <-
    ggboxplot(
      da,
      x = "Dataset",
      y = "C-index",
      fill = "Dataset",
      xlab = NULL,
      ylab = "C-index",
      font.y = 12,
      font.xtickslab = 12,
      font.ytickslab = 12,
      add = "mean_sd"
    ) +theme(
      axis.title.x = element_blank(),
      legend.position = "none"
    )
  idex<-lapply(list(resvalid$cidx),sumcidx)
  validation<-list(validexternal=resvalid)
  result<-list(p=p,idex=idex,validation=validation)
  
  return(result)
}
unicox<-function(time,status,data){
  require(survival)
  Surv<-Surv(time,status)
  coxph <- c()
  for (i in seq_len(nrow(data))) {
    gene<- unlist(data[i,])
    coxtest <- coxph(Surv ~ gene,data = data)
    summcph <- summary(coxtest)
    coeffs <- c(summcph$coefficients[1,1:2], summcph$conf.int[1,3:4],
                summcph$coefficients[1,5])
    coxph <- rbind(coxph, coeffs)
  }
  colnames(coxph) <- c('Coefficient','HR','Lower95','Upper95','PValue')
  rownames(coxph) <- rownames(data)
  return(coxph)
}

fp<-function(tab,xticks){
  forestplot(
    tab[, c(6,5)],
    hrzl_lines = gpar(col = "#444444"),
    boxsize = .35,line.margin =0.1,
    mean = tab[, "HR"],
    lower = tab[, "Lower95"],
    upper = tab[, "Upper95"],
    zero = 1, ci.vertices = T,
    ci.vertices.height = 0.2,
    xticks = c(0, 1, xticks),lineheight = unit(0.5, "cm"),
    
    col = fpColors(box = "darkblue",lines = "darkblue",
                   zero = "gray",summary = "black"),
    xlab = "Hazard ratio",clip=c(0,max(xticks)),
    graphwidth = unit(50, "mm"),
    is.summary = c(TRUE, rep(FALSE, (nrow(tab) - 1)),T),
    txt_gp = fpTxtGp(cex = 1, label = gpar(fontfamily = "serif")),
    colgap = unit(3.5, "mm")
  )
}

unicoxforestp<-function(data,xlim){
  #data should be a result object of unicox data
  data<-round(data,3)
  
  data$p<- format.pval(data$PValue, digits = 3, eps = 0.001)
  data$label<-row.names(data)
  
  data$method<-rep("Univariate",nrow(data))
  data$text<-paste(data$HR,paste("(",data$Lower95,"-",paste(paste(data$Upper95,",",sep=""),data$p,sep=" "),")",sep=""),sep=" ")
  
  f<-rep(NA,ncol(data))
  data<-rbind(f,data)
  data$method[1]<-"Univariate"
  data$UCI2<-ifelse(data$Upper95>xlim,xlim,NA)
  LfLabels<-data.frame(x=1,
                       y=length(unique(data$label)),
                       lab="HR (95% CI, P value)")
  LfLabel2<-data.frame(x=0.5,
                       y=length(unique(data$label)),
                       lab="Variable")
  p1<-ggplot(data = data, aes(y = label)) +
    geom_stripes()+
    ggplot2::theme(panel.background = element_rect(colour = NA,fill = NA))+
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    )+
    theme(plot.margin = margin(5.5, -5.5, 5.5, 5.5, "pt"))+xlab(" ")+
    geom_text(aes(x = 0.5, label = label),vjust=0.5,  hjust=0,fontface="bold",size=4.8,) +
    geom_text(aes(x = 1, label = text),vjust=0.5,  hjust=1, size=4.8
    ) +
    geom_text(data=LfLabels,aes(x,y,label=lab, fontface="bold"),hjust=1,size=5)+
    geom_text(data=LfLabel2,aes(x,y,label=lab, fontface="bold"),hjust=0,size=5)+
    geom_hline(aes(yintercept=length(unique(data$label))-0.5),linetype=2, size=1)+
    geom_hline(aes(yintercept=length(unique(data$label))+0.5),linetype=1, size=1.5)+
    theme(axis.line.x= element_line("black", size= 1.5))
  
  p2 <-  ggplot2::ggplot(data, aes(x = as.numeric(HR), y = label))+
    ggplot2::theme(
      panel.background = element_rect(colour = NA,fill = NA))+theme(
        axis.title.x = element_text(face="bold",size=13),
        axis.text.x = element_text(face="bold",size=13),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
    xlab("Hazard ratio")+
    scale_x_continuous(breaks  = c(0,1,xlim)
    )+coord_cartesian(xlim=c(0,xlim))+
    geom_stripes() +
    theme(plot.margin = margin(5.5, 5.5, 5.5, -4.5, "pt"))+
    geom_line(data = data.frame(x=c(1,1),y = c(0.5,nrow(data)-0.5)),aes(x=x,y=y))+
    geom_effect(ggplot2::aes(xmin = as.numeric(Lower95), xmax =as.numeric(Upper95)),size=1.5
    )+
    geom_segment(aes(x = as.numeric(Lower95),yend=label,xend = as.numeric(UCI2)+xlim/18),lineend = "round", # See available arrow types in example above
                 linejoin = "round",
                 arrow = arrow(length = unit(0.5, "cm"),type="closed"))+
    geom_hline(aes(yintercept=length(unique(data$label))-0.5),linetype=2, size=1)+
    geom_hline(aes(yintercept=length(unique(data$label))+0.5),linetype=1, size=1.5)+
    theme(axis.line.x= element_line("black", size= 1.5))
  p<-ggarrange(p1,p2,nrow = 1,align="h")
  return(p)
}

coxforestp<-function(unicox,multicox,legend.pos,xlim,varname){
  if(!identical(row.names(unicox),row.names(multicox))){
    stop("The rownames of unicox result and multicox result are not identical")
  }
  if(!is.null(varname)){
    if(length(varname)!=length(row.names(unicox))){
      stop("The length of varname is not identical to the variables in the cox model")
    }
  }
  unicox<-round(unicox,3)
  unicox<-as.data.frame(unicox)
  unicox$Pvalue<- format.pval(unicox$`P value`, digits = 3, eps = 0.001)
  multicox<-round(multicox,3)
  multicox<-as.data.frame(multicox)
  multicox$Pvalue<- format.pval(multicox$`P value`, digits = 3, eps = 0.001)
  
  # unicox$label<-row.names(unicox)
  # multicox$label<-row.names(multicox)
  unicox$method<-rep("Univariate",nrow(unicox))
  multicox$method<-rep("Multivariate",nrow(multicox))
  if(varname==""||is.null(varname)){
    unicox$label<-row.names(unicox)
    multicox$label<-row.names(multicox)
    data<-rbind(unicox,multicox)
  } else{
    row.names(unicox)<-varname
    row.names(multicox)<-varname
    unicox$label<-row.names(unicox)
    multicox$label<-row.names(multicox)
    data<-rbind(unicox,multicox)
  }
  
  data$p<-data$Pvalue
  data$text<-paste(data$HR,paste("(",data$LCI,"-",paste(paste(data$UCI,",",sep=""),data$p,sep=" "),")",sep=""),sep=" ")
  data$UCI2<-ifelse(data$UCI>xlim,xlim,NA)
  f<-rep(NA,ncol(data))
  data<-rbind(f,data)
  data$method[1]<-"Univariate"
  
  LfLabels<-data.frame(x=1,
                       y=length(unique(data$label)),
                       lab="HR (95% CI, P value)")
  LfLabel2<-data.frame(x=0.5,
                       y=length(unique(data$label)),
                       lab="Variable")
  p1<-ggplot(data = data, aes(y = label)) +
    geom_stripes()+
    ggplot2::theme(panel.background = element_rect(colour = NA,fill = NA))+
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    )+
    theme(plot.margin = margin(15.5, -4, 5.5, 5.5, "pt"))+xlab(" ")+
    geom_text(aes(x = 0.5, label = label),vjust=0.5,  hjust=0,fontface="bold",size=5,) +
    geom_text(aes(x = 1, label = text,group=method),vjust=0.5,  hjust=1, size=5,
              position = ggstance::position_dodgev(height = 0.5)) +
    geom_text(data=LfLabels,aes(x,y,label=lab, fontface="bold"),hjust=1,size=5)+
    geom_text(data=LfLabel2,aes(x,y,label=lab, fontface="bold"),hjust=0,size=5)+
    geom_hline(aes(yintercept=length(unique(data$label))-0.5),linetype=2, size=1)+
    geom_hline(aes(yintercept=length(unique(data$label))+0.5),linetype=1, size=1.5)+
    theme(axis.line.x= element_line("black", size= 1.5))
  
  p2 <-  ggplot2::ggplot(data, aes(x = as.numeric(HR), y = label))+
    ggplot2::theme(
      panel.background = element_rect(colour = NA,fill = NA))+theme(
        axis.title.x = element_text(face="bold",size=13),
        axis.text.x = element_text(face="bold",size=13),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
    xlab("Hazard ratio")+
    scale_x_continuous(breaks  = c(0,1,xlim)
    )+coord_cartesian(xlim=c(0,xlim))+
    geom_stripes() +theme(legend.position = legend.pos,
                          legend.title = element_blank(),
                          # legend.key.size = unit(0.15, "inches"),
                          legend.margin = margin(-4, 0, 0, 0),
                          legend.text = element_text( face = "bold"),
                          legend.background = element_blank(),
                          legend.box.background = element_blank(),
                          legend.key = element_blank()
    )+
    theme(plot.margin = margin(15.5, 5.5, 5.5, -5.5, "pt"))+
    geom_line(data = data.frame(x=c(1,1),y = c(0.5,length(unique(data$label))-0.5)),aes(x=x,y=y))+
    geom_effect(ggplot2::aes(xmin = as.numeric(LCI), xmax =as.numeric(UCI),
                             colour = method, shape =method,fill=method),size=2,
                position = ggstance::position_dodgev(height = 0.5))+
    geom_point(aes(x = as.numeric(UCI2)+xlim/18,group=method,fill=method,colour=method), shape = 23,
               position = ggstance::position_dodgev(height = 0.5), size = 5.5, show.legend = F) +
    geom_hline(aes(yintercept=length(unique(data$label))-0.5),linetype=2, size=1)+
    geom_hline(aes(yintercept=length(unique(data$label))+0.5),linetype=1, size=1.5)+
    theme(axis.line.x= element_line("black", size= 1.5))
  p<-ggarrange(p1,p2,nrow = 1,align="h")
  return(p)
}


genediff<-function(da,x,y,subgroup,plottype,subgroupname,xlab,ylab,genediffp,genediffmethod,compar,comparefer){
  require(ggpubr)
  if(subgroup ==F ){
    if(plottype=="Bar plot"){
      p<-ggbarplot(
        da,
        x = x,
        y = y,
        fill = x,
        add = "mean_sd",
        ylab = ylab,
        xlab = xlab
      ) + theme(legend.position = 'none')
    } else if(plottype=="Box plot"){
      p<-ggboxplot(
        da,
        x = x,
        y = y,
        fill = x,
        add = "mean_sd",
        ylab = ylab,
        xlab = xlab
      ) + theme(legend.position = 'none')
    } else{
      p<-ggviolin(
        da,
        x = x,
        y = y,
        fill = x,
        add = "boxplot", add.params = list(fill = "white"),
        color=x,
        ylab = ylab,
        xlab = xlab
      ) + theme(legend.position = 'none')
    }
  } else{
    if(plottype=="Bar plot"){
      p<-ggbarplot(
        da,
        x = x,
        y = y,
        fill = x,
        add = "mean_sd",
        facet.by = subgroupname,
        ylab = ylab,
        xlab = xlab
      ) + theme(legend.position = 'none')
    }else if(plottype=="Box plot"){
      
      p<-ggboxplot(
        da,
        x = x,
        y = y,
        fill = x,
        add = "mean_sd",
        facet.by = subgroupname,
        ylab = ylab,
        xlab = xlab
      ) + theme(legend.position = 'none')
    } else{
      p<-ggviolin(
        da,
        x = x,
        y = y,
        fill = x,
        add = "boxplot", add.params = list(fill = "white"),
        color=x,
        facet.by =subgroupname,
        ylab = ylab,
        xlab = xlab
      ) + theme(legend.position = 'none')
    }
  }
  if(genediffp== TRUE){
    if(genediffmethod== "T"){
      if(compar == "Pairwise comparison"){
        my_comparisons = combn(levels(da$group),2,simplify = F)
        p<-p+stat_compare_means(comparisons = my_comparisons,method = "t.test")
      }else{
        p<-p+stat_compare_means(ref.group = comparefer,method = "t.test")
      }
    }else{
      if(compar == "Pairwise comparison"){
        my_comparisons = combn(levels(da$group),2,simplify = F)
        p<-p+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
      }else{
        p<-p+stat_compare_means(ref.group = comparefer,method = "wilcox.test")
      }
      
    }
  }else{
    p
  }
  return(p)
}

htm <-function(
  DEG,
  
  sigDEG,
  group,
  color,
  DEGscale,
  expres,
  coltitsize,
  ClusterR,
  cluster_row_slices,
  cludistanrow,
  clumethodrow,
  Rdend_side,
  showRname,
  Rnameside,
  showFDR,
  showFC,
  ClusterC,
  heatname,
  cluster_column_slices,
  cludistancol,
  clumethodcol,
  Cdend_side,
  showCname,
  Cnameside,
  heatColors){
  
  normalize <- function(x) {
    return((x - min(x)) / (max(x) - min(x)))
  }
  
  
  heat <- expres[row.names(sigDEG), ]
  
  if (DEGscale == 'Scale') {
    heat <- t(scale(t(heat), center = F, scale = T))
  } else if (DEGscale == 'Center') {
    heat <- t(scale(t(heat), center = T, scale = F))
  } else if (DEGscale == 'Log') {
    heat <- t(log(t(heat)))
  } else if (DEGscale == "Z-score") {
    heat <- t(scale(t(heat), center = T, scale = T))
  } else if (DEGscale == "0-1 normalization") {
    heat <- t(normalize(t(heat)))
  }
  
  FDR1 <- log10(as.numeric(sigDEG$AdjustedP))
  
  FC <- sigDEG$LogFC
  
  require(ComplexHeatmap)
  if (showFDR == TRUE & showFC == FALSE) {
    row_an <-
      rowAnnotation(Log10adjP = anno_barplot(FDR1),
                    annotation_name_rot = 90)
  } else if (showFC == TRUE & showFDR == FALSE) {
    row_an <-
      rowAnnotation(LogFC = anno_barplot(FC), annotation_name_rot = 90)
  } else if (showFC == TRUE & showFDR == TRUE) {
    row_an <-
      rowAnnotation(
        Log10adjP = anno_barplot(FDR1),
        LogFC = anno_barplot(FC),
        annotation_name_rot = 90
      )
  } else{
    row_an <- NULL
  }
  
  Heatmap(
    heat,
    name = heatname,
    col = heatColors,
    border = F,
    column_title = levels(group),
    column_title_gp = gpar(
      fill = color,
      col = "black",
      fontsize = 15
    ),
    
    cluster_rows = ClusterR,
    cluster_row_slices = cluster_row_slices,
    clustering_distance_rows = cludistanrow,
    clustering_method_rows = clumethodrow,
    row_dend_side = Rdend_side,
    show_row_names = showRname,
    row_names_side = Rnameside,
    right_annotation = row_an,
    
    cluster_columns = ClusterC,
    cluster_column_slices = cluster_column_slices,
    clustering_distance_columns = cludistancol,
    clustering_method_columns	= clumethodcol,
    column_dend_side = Cdend_side,
    show_column_names = showCname,
    column_names_side = Cnameside,
    column_split = group
  )
}




voca<-function(DEG, vocaXlab,vocaYlab,vocacutp,vocacutfc,legendPosition){
  require(EnhancedVolcano)
  if (vocaXlab == "Default") {
    xlab <- bquote( ~ Log[2] ~ "fold change")
  } else{
    xlab <- vocaXlab
  }
  
  if (vocaYlab == "Default") {
    ylab <- bquote( ~ -Log[10] ~ italic(adjustP))
  } else{
    ylab <- vocaYlab
    
  }
  
  EnhancedVolcano(
    toptable = DEG,
    lab = row.names(DEG),
    x = "LogFC",
    y = "AdjustedP",
    selectLab = NULL,
    xlim = c(min(DEG$LogFC, na.rm = TRUE) - 1,
             max(DEG$LogFC, na.rm = TRUE) + 1),
    ylim = c(0, max(-log10(DEG$AdjustedP), na.rm = TRUE) + 5),
    
    xlab = xlab,
    ylab = ylab,
    axisLabSize = 18,
    title = NULL,
    subtitle = NULL,
    caption = paste0('Total = ', nrow(DEG), ' variables'),
    titleLabSize = 18,
    subtitleLabSize = 14,
    captionLabSize = 14,
    pCutoff = vocacutp,
    FCcutoff = vocacutfc,
    
    cutoffLineType = 'longdash',
    cutoffLineCol = 'black',
    cutoffLineWidth = 0.4,
    pointSize = 2.0,
    labSize = 3.0,
    labCol = 'black',
    labFace = 'plain',
    labhjust = 0.5,
    labvjust = 1.5,
    boxedLabels = FALSE,
    shape = 19,
    shapeCustom = NULL,
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    colCustom = NULL,
    colAlpha = 1 / 2,
    colGradient = NULL,
    colGradientBreaks = c(pCutoff, 1.0),
    colGradientLabels = c('0', '1.0'),
    colGradientLimits = c(0, 1.0),
    .legend = c('NS', 'Log2 FC', 'Adjust P', 'Adjust P & Log2 FC'),
    legendLabels = c(
      'NS',
      expression(Log[2] ~ FC),
      'Adjust P',
      expression(Adjust ~ P ~ and ~ log[2] ~ FC)
    ),
    legendPosition = legendPosition,
    legendLabSize = 14,
    legendIconSize = 4.0,
    shade = NULL,
    shadeLabel = NULL,
    shadeAlpha = 1 / 2,
    shadeFill = 'grey',
    shadeSize = 0.01,
    shadeBins = 2,
    drawConnectors = FALSE,
    widthConnectors = 0.5,
    typeConnectors = 'closed',
    endsConnectors = 'first',
    lengthConnectors = unit(0.01, 'npc'),
    colConnectors = 'grey10',
    hline = NULL,
    hlineType = 'longdash',
    hlineCol = 'black',
    hlineWidth = 0.4,
    vline = NULL,
    vlineType = 'longdash',
    vlineCol = 'black',
    vlineWidth = 0.4,
    gridlines.major = F,
    gridlines.minor = F,
    border = 'partial',
    borderWidth = 0.8,
    borderColour = 'black'
  )
}

MAplt <-
  function(DEG,
           DEGmethod,
           MAstopmeth,
           MAcutp,
           MAfc,
           Topgene,
           MAgenesym,
           MAXlab,
           MAYlab,
           MAlegendPosition) {
    
    
    if (DEGmethod == "edgeR-LRT" | DEGmethod == "edgeR-QLF") {
      DEG <- subset(DEG, select = c(logCPM, LogFC, AdjustedP))
      names(DEG) <- c("baseMeanLog2" , "log2FoldChange" , "padj")
    } else if (DEGmethod == "DESeq2") {
      DEG <- subset(DEG, select = c(BaseMean, LogFC, AdjustedP))
      names(DEG) <- c("baseMean" , "log2FoldChange" , "padj")
    } else{
      DEG <- subset(DEG, select = c(AveExpr, LogFC, AdjustedP))
      names(DEG) <- c("baseMean" , "log2FoldChange" , "padj")
    }
    if (MAstopmeth == "Adjusted P") {
      MAstopmeth1 <- "padj"
    } else{
      MAstopmeth1 <- "fc"
    }
    
    ggmaplot(
      DEG,
      fdr = MAcutp,
      fc =MAfc,
      genenames = row.names(DEG),
      detection_call = NULL,
      size = 0.4,
      alpha = 1,
      font.label = c(12, "plain", "black"),
      label.rectangle = FALSE,
      palette = c("#B31B21", "#1465AC", "darkgray"),
      top = Topgene,
      select.top.method = MAstopmeth1,
      label.select = MAgenesym,
      main = NULL,
      xlab = MAXlab,
      ylab = MAYlab,
      ggtheme = theme_classic(),
      legend = MAlegendPosition
    )
  }




padjplt<-function(DEG,col){
  as.ggplot(function() hist(
    DEG$AdjustedP,
    col = col,
    border = "white",
    xlab = "Adjusted P value",
    ylab = "Number of genes",
    main = "Adjusted P value distribution",
    plot=T
  ))
  
}


GO<-function(method,
             gene,
             geneList,
             ont,
             pAdjustMethod,
             minGSSize,
             maxGSSize,
             pvalueCutoff,
             qvalueCutoff,
             readable,
             gseamethod#,
             # nPerm
){
  
  if(method=="ORA"){
    GOres <- enrichGO(
      gene          = gene,
      OrgDb         = org.Hs.eg.db,
      ont           = ont,
      pAdjustMethod = pAdjustMethod,
      minGSSize     = minGSSize,
      maxGSSize     = maxGSSize,
      pvalueCutoff  = pvalueCutoff,
      qvalueCutoff  = qvalueCutoff,
      readable      = readable
    )
  } else {
    GOres <- gseGO(
      geneList     = geneList ,
      OrgDb        = org.Hs.eg.db,
      keyType      = "ENTREZID",
      ont          = ont,
      # nPerm        = nPerm,
      minGSSize    = minGSSize,
      maxGSSize    = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      verbose      = FALSE,
      by           = gseamethod )
  }
  return(GOres)
}

KEGG<-function(method,
               gene,
               geneList,
               # ont,
               pAdjustMethod,
               minGSSize,
               maxGSSize,
               pvalueCutoff,
               qvalueCutoff,
               readable,
               gseamethod){
  if(method=="ORA"){
    KEGGres<-enrichKEGG(
      gene          = gene,
      organism      = "hsa",
      keyType       = "kegg",
      pvalueCutoff  = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      # universe,
      minGSSize     = minGSSize,
      maxGSSize     = maxGSSize,
      qvalueCutoff  = qvalueCutoff,
      use_internal_data = FALSE
    )
  } else if(method=="GSEA"){
    KEGGres <- gseKEGG(
      geneList     = geneList,
      organism     = 'hsa',
      exponent     = 1,
      
      minGSSize    = minGSSize,
      maxGSSize    = maxGSSize,
      pAdjustMethod= pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      verbose      = FALSE,
      by           = gseamethod )
  } else if(method=="Module ORA"){
    KEGGres<-enrichMKEGG(
      gene        =gene,
      organism    = "hsa",
      keyType     = "kegg",
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      
      minGSSize   = minGSSize,
      maxGSSize   = maxGSSize,
      qvalueCutoff = qvalueCutoff
    )
  } else{
    KEGGres<-gseMKEGG(
      geneList=geneList,
      organism = "hsa",
      keyType = "kegg",
      exponent = 1,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      eps = 1e-10,
      
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      verbose = TRUE,
      seed = FALSE,
      by = gseamethod
    )
  }
  return(KEGGres)
}
msig <-
  function(gene,
           method,
           pvalueCutoff,
           pAdjustMethod,
           minGSSize,
           maxGSSize,
           qvalueCutoff,
           TERM2GENE,
           geneList,
           
           gseamethod) {
    require(msigdbr)
    gene_sets = msigdbr(species = "Homo sapiens", category = "H")
    msigdbr_t2g = gene_sets %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
    if (method == "ORA") {
      msigres <- enricher(
        gene          = gene,
        pvalueCutoff  = pvalueCutoff,
        pAdjustMethod = pAdjustMethod,
        
        minGSSize     = minGSSize ,
        maxGSSize     = maxGSSize,
        qvalueCutoff  = qvalueCutoff,
        TERM2GENE     = msigdbr_t2g,
        TERM2NAME     = NA
      )
    } else {
      msigres <- GSEA(
        geneList      = geneList,
        exponent      = 1,
        
        minGSSize     = minGSSize,
        maxGSSize     = maxGSSize,
        eps           = 1e-10,
        pvalueCutoff  = pvalueCutoff,
        pAdjustMethod = pAdjustMethod,
        TERM2GENE     = msigdbr_t2g,
        TERM2NAME     = NA,
        verbose       = F,
        seed          = FALSE,
        by            = gseamethod
      )
    }
    return(msigres)
  }

Path <-
  function(gene,
           method,
           pvalueCutoff,
           pAdjustMethod,
           qvalueCutoff,
           minGSSize,
           maxGSSize,
           readable,
           geneList,
           
           gseamethod) {
    if (method == "ORA") {
      pathres <-
        enrichPathway(
          gene = gene,
          organism = "human",
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          qvalueCutoff = qvalueCutoff,
          
          minGSSize = minGSSize,
          maxGSSize = maxGSSize,
          readable = readable
        )
    } else{
      pathres <- gsePathway(
        geneList  = geneList,
        organism  = "human",
        exponent = 1,
        
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = pAdjustMethod,
        verbose = TRUE,
        seed = FALSE,
        by  = gseamethod
      )
    }
    return(pathres)
    
  }

FEA <-
  function(ana,
           method,
           geneList,
           nPerm,
           minGSSize,
           maxGSSize,
           pvalueCutoff,
           pAdjustMethod,
           qvalueCutoff,
           gseamethod,
           gene,
           readable,
           ont) {
    if(ana=="GO"){
      res<-GO(
        method=method,
        gene=gene,
        geneList=geneList,
        ont=ont,
        pAdjustMethod=pAdjustMethod,
        minGSSize=minGSSize,
        maxGSSize=maxGSSize,
        pvalueCutoff=pvalueCutoff,
        qvalueCutoff=qvalueCutoff,
        readable=readable,
        gseamethod=gseamethod#,
        
      )
    } else if(ana=="KEGG"){
      res<-KEGG(
        method=method,
        gene=gene,
        geneList=geneList,
        pAdjustMethod=pAdjustMethod,
        minGSSize=minGSSize,
        maxGSSize=maxGSSize,
        pvalueCutoff=pvalueCutoff,
        qvalueCutoff=qvalueCutoff,
        readable=readable,
        gseamethod=gseamethod#,
        
      )
    }else if (ana=="MSigDb"){
      res<-msig(
        gene=gene,
        method=method,
        pvalueCutoff=pvalueCutoff,
        pAdjustMethod=pAdjustMethod,
        minGSSize=minGSSize,
        maxGSSize=maxGSSize,
        qvalueCutoff=qvalueCutoff,
        TERM2GENE=TERM2GENE,
        geneList=geneList,
        
        gseamethod=gseamethod)
    }else if (ana=="Reactome Pathway"){
      res<-Path(
        gene=gene,
        method=method,
        pvalueCutoff=pvalueCutoff,
        pAdjustMethod=pAdjustMethod,
        qvalueCutoff=qvalueCutoff,
        minGSSize=minGSSize,
        maxGSSize=maxGSSize,
        readable=readable,
        geneList=geneList,
        
        gseamethod=gseamethod
      )
    }
  }



FS<-function (Data, value)
{
  vars = apply(Data, 1, var)
  index = sort(vars, decreasing = TRUE, index.return = TRUE)
  if (value > nrow(Data)) {
    value = nrow(Data)
    cat("Warning: the feature selection number is beyond the original feature numnber")
  }
  cutoff = index$x[value]
  index = index$ix[1:value]
  selectData = Data[index, ]
  selectData
}
preproc<-function(expres,clinical,thresholdZ.k, topvar,selectvar){
  if(selectvar==TRUE){
    expres = FS(expres, value = topvar)
  }
  gsg = goodSamplesGenes(t(expres), verbose = 3);
  if (!gsg$allOK)
  {
    expres = expres[gsg$goodGenes,gsg$goodSamples]
  }
  
  require(flashClust)
  A=adjacency(expres,type="distance")
  k=as.numeric(apply(A,2,sum))-1
  Z.k=scale(k)
  thresholdZ.k=thresholdZ.k
  outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
  sampleTree = flashClust(as.dist(1-A), method = "average")
  traitColors=data.frame(numbers2colors(clinical,signed=FALSE))
  dimnames(traitColors)[[2]]=paste(names(clinical),"C",sep="")
  datColors=data.frame(outlierC=outlierColor,traitColors)
  
  removed.samples= Z.k<thresholdZ.k | is.na(Z.k)
  index<-which(removed.samples==T)
  expres<-expres[,-index]
  clinical<-clinical[-index,drop=F,]
  removed.samples<-table(removed.samples)
  res<-list(expres=expres,clinical=clinical,sampleTree=sampleTree,datColors=datColors,removed.samples=removed.samples)
  return(res)
}


autowgcna <-
  function(expres,
           clinical,
           networkType,
           RsquaredCut,
           corType,
           TOMType,
           deepSplit,
           detectCutHeight,
           minModuleSize,
           reassignThreshold,
           mergeCutHeight,
           numericLabels,
           pamRespectsDendro,
           clinical1
  ){
    powers = 1:20
    sft = pickSoftThreshold(
      t(expres),
      powerVector = powers,
      verbose = 5,
      RsquaredCut=RsquaredCut,
      networkType = networkType
    )
    cor <- WGCNA::cor
    if(!is.na(sft$powerEstimate)){
      net = blockwiseModules(
        datExpr=t(as.data.frame(expres)),############
        power = sft$powerEstimate,
        corType=corType,
        TOMType = TOMType,
        networkType= networkType,
        deepSplit=deepSplit,
        detectCutHeight=detectCutHeight,
        minModuleSize = minModuleSize,
        reassignThreshold = reassignThreshold,
        mergeCutHeight = 0.25,
        numericLabels = TRUE,
        pamRespectsDendro = FALSE,
        saveTOMs = F,
        verbose = 3
      )
      cor<-stats::cor
      moduleLabels = net$colors
      moduleColors = labels2colors(net$colors)
      
      MEs0 = moduleEigengenes(t(expres), moduleColors)$eigengenes
      MEs = orderMEs(MEs0)
      
      geneTree = net$dendrograms[[1]]
      res <-
        list(
          sft = sft,
          net = net,
          expres = expres,
          clinical = clinical,
          clinical1=clinical1,
          moduleLabels = moduleLabels,
          moduleColors = moduleColors,
          MEs = MEs,
          geneTree = geneTree
        )
    } else{
      res<-sft
    }
    return(res)
  }

MTR<-function(xx,BMar, LMar, TMar, RMar){
  moduleLabels = xx$moduleLabels
  moduleColors = xx$moduleColors
  
  geneTree = xx$geneTree
  nGenes = nrow(xx$expres);
  nSamples = ncol(xx$expres);
  
  MEs = xx$MEs
  moduleTraitCor = cor(MEs, xx$clinical, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(BMar, LMar, TMar, RMar));
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(xx$clinical),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
}



netres <- function(clinical,
                   clinical1,
                   expres,
                   trait,
                   moduleColors,
                   module,
                   output) {
  MEs0 = moduleEigengenes(t(expres), moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  nSamples <- ncol(expres)
  traitdata = as.data.frame(clinical[, trait])
  names(traitdata) = trait
  
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(t(expres), MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep = "")
  names(MMPvalue) = paste("p.MM", modNames, sep = "")
  geneTraitSignificance = as.data.frame(cor(t(expres), traitdata, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", names(traitdata), sep = "")
  names(GSPvalue) = paste("p.GS.", names(traitdata), sep = "")
  
  geneInfo0 = data.frame(moduleColor = moduleColors,
                         geneTraitSignificance,
                         GSPvalue)
  modOrder = order(-abs(cor(MEs, traitdata, use = "p")))
  
  for (mod in 1:ncol(geneModuleMembership)) {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]])
    
    names(geneInfo0) = c(oldNames,
                         paste("MM.", modNames[modOrder[mod]], sep = ""),
                         paste("p.MM.", modNames[modOrder[mod]], sep = ""))
  }
  geneOrder = order(geneInfo0$moduleColor,-abs(geneInfo0[, names(geneTraitSignificance)]))
  
  geneInfo = geneInfo0[geneOrder,]
  if(output!="Non-grey modules"){
    Output<-geneInfo[geneInfo$moduleColor==output,]
  }else{
    Output<-geneInfo[geneInfo$moduleColor!=output,]
  }
  column = match(module, modNames)
  moduleGenes = moduleColors == module
  
  res <- list(geneInfo = geneInfo,
              column = column,
              moduleGenes = moduleGenes,
              moduleColors=moduleColors,
              geneModuleMembership=geneModuleMembership,
              geneTraitSignificance=geneTraitSignificance,
              clinical=clinical1,
              expres=expres,
              # clinical1=clinical1,
              Output=Output
  )
  return(res)
}

#######################################################################################################################
#######################################################################################################################
nestml<-function(expres,clinical,split,endpoint,ratio,ctrl,inner,outer,lnid,seed){
  require(caret)
  require(mlr)
  require(randomForestSRC)
  expres$time<-clinical[,endpoint[1]]
  expres$status<-clinical[,endpoint[2]]
  if(split==T){
    set.seed(seed)
    index <- createDataPartition(expres$status, p = ratio, list = FALSE)
    trainset<-expres[index,]
    testset<-expres[-index,]
    trainsurv<-clinical[index,]
    testsurv<-clinical[-index,]
    train.task <- makeSurvTask(data = trainset, target = c("time","status"))
    test.task <- makeSurvTask(data = testset, target = c("time","status"))
  } else {
    train.task <- makeSurvTask(data = expres, target = c("time","status"))
  }
  
  #learn rf
  rf<-makeLearner(cl = "surv.randomForestSRC", predict.type = "response",importance = TRUE,id="RandomForestSRC")
  rf_par <- makeParamSet(makeIntegerParam("mtry", 1,15),#p/3 for regression, other famile is sqrt(p)
                         makeDiscreteParam("nodesize", c(3, 5, 8, 10,15,18,20)),
                         makeDiscreteParam("ntree", c(500,1000,1500,2000)),
                         makeIntegerParam("nodedepth", lower = 5, upper = 20))
  lrn_rf = makeTuneWrapper(rf, resampling = inner,
                           par.set = rf_par, control = ctrl,
                           show.info = FALSE)
  #learn glmboost
  glmboost<-makeLearner(cl = "surv.glmboost", predict.type = "response",id="Glmboost")
  glmboost_par <- makeParamSet(makeIntegerParam("mstop", 1e2, 1e3),
                               makeDiscreteParam("nu", c(0.05, 0.1, 0.3, 0.5, 0.8,1)))
  lrn_glmboost = makeTuneWrapper(glmboost, resampling = inner,
                                 par.set = glmboost_par, control = ctrl,
                                 show.info = FALSE)
  
  # learn coxboost
  coxboost<-makeLearner(cl = "surv.CoxBoost", predict.type = "response",id="Coxboost")
  coxboost_par <- makeParamSet(
    #makeNumericLearnerParam(id = "penalty", lower = 9*sum(train$status==1)),
    makeIntegerLearnerParam(id = "stepno", lower = 50, upper = 200)
  )
  lrn_coxboost = makeTuneWrapper(coxboost, resampling = inner,
                                 par.set = coxboost_par, control = ctrl,
                                 show.info = F)
  #learn en
  en<-makeLearner(cl = "surv.glmnet",predict.type = "response",id="ElasticNet")
  en_par <- makeParamSet(
    makeNumericLearnerParam(id="alpha", lower = 0, upper = 1),
    makeNumericLearnerParam(id="s",lower=0.001,upper=30)
    #makeNumericLearnerParam(id = "s", lower = 0)
  )
  lrn_en = makeTuneWrapper(en, resampling = inner,
                           par.set = en_par, control = ctrl,
                           show.info = FALSE)
  #learn ridge
  ridge<-makeLearner(cl = "surv.glmnet",predict.type = "response",alpha = 0,id="Ridge")
  ridge_par <- makeParamSet(makeNumericParam("s", lower = 0, upper = 20))
  lrn_ridge = makeTuneWrapper(ridge, resampling = inner,
                              par.set = ridge_par, control = ctrl,
                              show.info = FALSE)
  #learn lasso
  lasso<-makeLearner(cl = "surv.glmnet",predict.type = "response",alpha = 1,id="Lasso")
  lasso_par <- makeParamSet(makeNumericParam("s", lower = 0, upper = 20))
  lrn_lasso = makeTuneWrapper(lasso, resampling = inner,
                              par.set = lasso_par, control = ctrl,
                              show.info = FALSE)
  
  lrns<-list(lrn_rf=lrn_rf,lrn_glmboost=lrn_glmboost,lrn_coxboost=lrn_coxboost,lrn_en=lrn_en,lrn_ridge=lrn_ridge,lrn_lasso=lrn_lasso)
  
  idx<-match(lnid,c("RandomForestSRC","Glmboost","Coxboost","ElasticNet","Ridge","Lasso"))
  lrns<-lrns[idx]
  set.seed(seed)
  res = benchmark(lrns, train.task, outer, keep.extract = T)
  
  #fit the best-performing model, get predictions and conduct feature selection
  outmeasure<-print(res)
  if(all(outmeasure$cindex.test.mean %in%  NaN)){
    return(NULL)
  }else{
    
    bestlearner<-outmeasure[which.max(outmeasure$cindex.test.mean),]
    
    tune.res = getBMRTuneResults(res, task.ids = bestlearner$task.id,
                                 learner.ids = bestlearner$learner.id)
    fun<-function(x){x$y}
    if(split==T) {
      bestinner<-which.max(lapply(tune.res$trainset[[1]],fun))
    }else {
      bestinner<-which.max(lapply(tune.res$expres[[1]],fun))
    }
    lrnid<-gsub(".tuned","",bestlearner$learner.id)
    print(paste(lrnid, "performs best among the learns you selected"))
    if(split==T) {
      if(lrnid == "RandomForestSRC") {
        finalmodel <-
          setHyperPars(rf, par.vals = tune.res$trainset[[1]][[bestinner]]$x)
        modeltrain <- mlr::train(learner=finalmodel, task=train.task)
        
        modeltrain <- getLearnerModel(modeltrain)
        
        selvar<-max.subtree(modeltrain)$topvars
        listidx<-as.data.frame(max.subtree(modeltrain)$order)
        listidx<-listidx[selvar,1,drop=F]
        names(listidx)<-"coeff"
        
        
        model.predicted<-predict(modeltrain,testset)
        train.chf<-modeltrain$chf
        test.chf <- model.predicted$chf
        testsurv$Risk<-rowSums(test.chf)
        trainsurv$Risk<-rowSums(train.chf)
        
        
      } else if (lrnid == "Glmboost") {
        finalmodel <-
          setHyperPars(glmboost, par.vals = tune.res$trainset[[1]][[bestinner]]$x)
        modeltrain <- mlr::train(learner=finalmodel, task=train.task)
        modelprediction.test <- predict(modeltrain, test.task)
        modelprediction.train <- predict(modeltrain, train.task)
        
        testsurv$Risk<-modelprediction.test$data$response
        trainsurv$Risk<-modelprediction.train$data$response
        modeldata <- getLearnerModel(modeltrain)
        selvar<-as.data.frame(coef(modeldata))
        selvar<-selvar[selvar[,1]!=0,,drop=F]
        listidx<-as.data.frame(selvar)[-1,]
        names(listidx)<-"coeff"
        selvar<-row.names(selvar)[-1]
        
      } else if (lrnid == "ElasticNet") {
        finalmodel <-
          setHyperPars(en, par.vals = tune.res$trainset[[1]][[bestinner]]$x)
        modeltrain <- mlr::train(learner=finalmodel, task=train.task)
        modelprediction.test <- predict(modeltrain, test.task)
        modelprediction.train <- predict(modeltrain, train.task)
        
        testsurv$Risk<-modelprediction.test$data$response
        trainsurv$Risk<-modelprediction.train$data$response
        modeldata <- getLearnerModel(modeltrain)
        selvar <-
          as.data.frame(as.matrix(coef(modeldata,
                                       s = tune.res$trainset[[1]][[bestinner]]$x$s,
                                       alpha = tune.res$trainset[[1]][[bestinner]]$x$alpha)))
        selvar<-selvar[selvar[,1]!=0,,drop=F]
        listidx<-as.data.frame(selvar)
        names(listidx)<-"coeff"
        selvar<-row.names(selvar)
      } else if (lrnid == "Ridge") {
        finalmodel <-
          setHyperPars(ridge, par.vals = tune.res$trainset[[1]][[bestinner]]$x)
        modeltrain <- mlr::train(learner=finalmodel, task=train.task)
        modelprediction.test <- predict(modeltrain, test.task)
        modelprediction.train <- predict(modeltrain, train.task)
        
        testsurv$Risk<-modelprediction.test$data$response
        trainsurv$Risk<-modelprediction.train$data$response
        modeldata <- getLearnerModel(modeltrain)
        selvar <-
          as.data.frame(as.matrix(coef(modeldata,
                                       s = tune.res$trainset[[1]][[bestinner]]$x$s)))
        selvar<-selvar[selvar[,1]!=0,,drop=F]
        listidx<-as.data.frame(selvar)
        names(listidx)<-"coeff"
        selvar<-row.names(selvar)
      } else if (lrnid == "Lasso") {
        finalmodel <-
          setHyperPars(lasso, par.vals = tune.res$trainset[[1]][[bestinner]]$x)
        modeltrain <- mlr::train(learner=finalmodel, task=train.task)
        modelprediction.test <- predict(modeltrain, test.task)
        modelprediction.train <- predict(modeltrain, train.task)
        
        testsurv$Risk<-modelprediction.test$data$response
        trainsurv$Risk<-modelprediction.train$data$response
        modeldata <- getLearnerModel(modeltrain)
        selvar <-
          as.data.frame(as.matrix(coef(modeldata,
                                       s = tune.res$trainset[[1]][[bestinner]]$x$s)))
        selvar<-selvar[selvar[,1]!=0,,drop=F]
        listidx<-as.data.frame(selvar)
        names(listidx)<-"coeff"
        selvar<-row.names(selvar)
        
      } else if (lrnid == "Coxboost") {
        finalmodel <-
          setHyperPars(coxboost, par.vals = tune.res$trainset[[1]][[bestinner]]$x)
        modeltrain <- mlr::train(learner=finalmodel, task=train.task)
        modelprediction.test <- predict(modeltrain, test.task)
        modelprediction.train <- predict(modeltrain, train.task)
        
        testsurv$Risk<-modelprediction.test$data$response
        trainsurv$Risk<-modelprediction.train$data$response
        modeldata <- getLearnerModel(modeltrain)
        selvar <-
          as.data.frame(as.matrix(coef(modeldata,
                                       s = tune.res$trainset[[1]][[bestinner]]$x$s)))
        
        selvar<-selvar[selvar$V1!=0,,drop=F]
        listidx<-as.data.frame(selvar)
        names(listidx)<-"coeff"
        selvar<-row.names(selvar)
      }
      
      res<-list(res,selvar, modeltrain, trainset, testset,lrnid,expres,clinical,testsurv,trainsurv,listidx)
      names(res)<-c("Benchmark result","Selected variables","Fitted model","Training set", "Test set","lrnid","expres","clinical","Testsurv","Trainsurv","listidx")
    } else {
      if(lrnid == "RandomForestSRC") {
        finalmodel <-
          setHyperPars(rf, par.vals = tune.res$expres[[1]][[bestinner]]$x)
        modeltrain <- mlr::train(learner=finalmodel, task=train.task)
        
        modeltrain <- getLearnerModel(modeltrain)
        
        selvar<-max.subtree(modeltrain)$topvars
        listidx<-as.data.frame(max.subtree(modeltrain)$order)
        listidx<-listidx[selvar,1,drop=F]
        names(listidx)<-"coeff"
        
        train.chf<-modeltrain$chf
        
        clinical$Risk<-rowSums(train.chf)
        
        
      } else if (lrnid == "Glmboost") {
        finalmodel <-
          setHyperPars(glmboost, par.vals = tune.res$expres[[1]][[bestinner]]$x)
        modeltrain <- mlr::train(learner=finalmodel, task=train.task)
        
        modelprediction.train <- predict(modeltrain, train.task)
        
        clinical$Risk<-modelprediction.train$data$response
        modeldata <- getLearnerModel(modeltrain)
        selvar<-as.data.frame(coef(modeldata))
        selvar<-selvar[selvar[,1]!=0,,drop=F]
        listidx<-as.data.frame(selvar)[-1,]
        names(listidx)<-"coeff"
        selvar<-row.names(selvar)[-1]
        
      } else if (lrnid == "ElasticNet") {
        
        finalmodel <-
          setHyperPars(en, par.vals = tune.res$expres[[1]][[bestinner]]$x)
        modeltrain <- mlr::train(learner=finalmodel, task=train.task)
        
        modelprediction.train <- predict(modeltrain, train.task)
        
        clinical$Risk   <- modelprediction.train$data$response
        modeldata <- getLearnerModel(modeltrain)
        selvar <-
          as.data.frame(as.matrix(coef(modeldata,
                                       s = tune.res$expres[[1]][[bestinner]]$x$s,
                                       alpha = tune.res$expres[[1]][[bestinner]]$x$alpha)))
        selvar<-selvar[selvar[,1]!=0,,drop=F]
        listidx<-as.data.frame(selvar)
        names(listidx)<-"coeff"
        selvar<-row.names(selvar)
      } else if (lrnid == "Ridge") {
        finalmodel <-
          setHyperPars(ridge, par.vals = tune.res$expres[[1]][[bestinner]]$x)
        modeltrain <- mlr::train(learner=finalmodel, task=train.task)
        
        modelprediction.train <- predict(modeltrain, train.task)
        
        clinical$Risk<-modelprediction.train$data$response
        modeldata <- getLearnerModel(modeltrain)
        selvar <-
          as.data.frame(as.matrix(coef(modeldata,
                                       s = tune.res$expres[[1]][[bestinner]]$x$s)))
        selvar<-selvar[selvar[,1]!=0,,drop=F]
        listidx<-as.data.frame(selvar)
        names(listidx)<-"coeff"
        selvar<-row.names(selvar)
      } else if (lrnid == "Lasso") {
        finalmodel <-
          setHyperPars(lasso, par.vals = tune.res$expres[[1]][[bestinner]]$x)
        modeltrain <- mlr::train(learner=finalmodel, task=train.task)
        
        modelprediction.train <- predict(modeltrain, train.task)
        
        clinical$Risk<-modelprediction.train$data$response
        modeldata <- getLearnerModel(modeltrain)
        selvar <-
          as.data.frame(as.matrix(coef(modeldata,
                                       s = tune.res$expres[[1]][[bestinner]]$x$s)))
        selvar<-selvar[selvar[,1]!=0,,drop=F]
        listidx<-as.data.frame(selvar)
        names(listidx)<-"coeff"
        selvar<-row.names(selvar)
        
      } else if (lrnid == "Coxboost") {
        finalmodel <-
          setHyperPars(coxboost, par.vals = tune.res$expres[[1]][[bestinner]]$x)
        modeltrain <- mlr::train(learner=finalmodel, task=train.task)
        
        modelprediction.train <- predict(modeltrain, train.task)
        
        clinical$Risk<-modelprediction.train$data$response
        
        modeldata <- getLearnerModel(modeltrain)
        selvar <-
          as.data.frame(as.matrix(coef(modeldata,
                                       s = tune.res$expres[[1]][[bestinner]]$x$s)))
        selvar<-selvar[selvar$V1!=0,,drop=F]
        listidx<-as.data.frame(selvar)
        names(listidx)<-"coeff"
        selvar<-row.names(selvar)
      }
      
      res<-list(res,selvar,modeltrain,lrnid,expres,clinical,listidx)
      names(res)<-c("Benchmark result","Selected variables", "Fitted model","lrnid","expres","clinical","listidx")
    }
    
    return(res)
  }
  
}


cv<-function(expres,clinical,split,endpoint,ratio,lnid,inner,ctrl,rdesc,seed){
  require(caret)
  require(mlr)
  require(randomForestSRC)
  
  rf<-makeLearner(cl = "surv.randomForestSRC", predict.type = "response",importance = TRUE,id="RandomForestSRC")
  rf_par <- makeParamSet(makeIntegerParam("mtry", 1,15),#p/3 for regression, other famile is sqrt(p)
                         makeDiscreteParam("nodesize", c(3, 5, 8, 10,15,18,20)),
                         makeDiscreteParam("ntree", c(500,1000,1500,2000)),
                         makeIntegerParam("nodedepth", lower = 5, upper = 20))
  rfm<-list(rf,rf_par)
  
  glmboost<-makeLearner(cl = "surv.glmboost", predict.type = "response",id="Glmboost")
  glmboost_par <- makeParamSet(makeIntegerParam("mstop", 1e2, 1e3),
                               makeDiscreteParam("nu", c(0.05, 0.1, 0.3, 0.5, 0.8,1)))
  glmboostm<-list(glmboost,glmboost_par)
  
  coxboost<-makeLearner(cl = "surv.CoxBoost", predict.type = "response",id="Coxboost")
  coxboost_par <- makeParamSet(
    
    makeIntegerLearnerParam(id = "stepno", lower = 50, upper = 200)
  )
  coxboostm<-list(coxboost,coxboost_par)
  
  en<-makeLearner(cl = "surv.glmnet",predict.type = "response",id="ElasticNet")
  en_par <- makeParamSet(
    makeNumericLearnerParam(id="alpha", lower = 0, upper = 1),
    makeNumericLearnerParam(id="s",lower=0.001,upper=30)
    
  )
  enm<-list(en,en_par)
  
  ridge<-makeLearner(cl = "surv.glmnet",predict.type = "response",alpha = 0,id="Ridge")
  ridge_par <- makeParamSet(makeNumericParam("s", lower = 0, upper = 20))
  ridgem<-list(ridge,ridge_par)
  
  lasso<-makeLearner(cl = "surv.glmnet",predict.type = "response",alpha = 1,id="Lasso")
  lasso_par <- makeParamSet(makeNumericParam("s", lower = 0, upper = 20))
  lassom<-list(lasso,lasso_par)
  
  mset<-list(RandomForestSRC=rfm,
             Glmboost=glmboostm,
             Coxboost=coxboostm,
             ElasticNet=enm,
             Ridge=ridgem,
             Lasso=lassom
  )
  
  
  if(identical(row.names(clinical),row.names(expres))){
    expres$time<-clinical[,endpoint[1]]
    expres$status<-clinical[,endpoint[2]]
    
  } else (stop("The rownames of clinical dataframe and gene expression dataframe are not identical."))
  
  if(split==T){
    set.seed(seed)
    index <- createDataPartition(expres$status, p = ratio, list = FALSE)
    
    trainset<-expres[index,]
    testset<-expres[-index,]
    
    trainsurv<-clinical[index,]
    testsurv<-clinical[-index,]
    train.task <- makeSurvTask(data = trainset, target = c("time","status"))
    test.task <- makeSurvTask(data = testset, target = c("time","status"))
    
  } else {
    train.task <- makeSurvTask(data = expres, target = c("time","status"))
  }
  
  idx<-match(lnid,c("RandomForestSRC","Glmboost","Coxboost","ElasticNet","Ridge","Lasso"))
  mset<-mset[idx]
  
  modeltune<-function(x,task,inner, ctrl){
    tune<-tuneParams(learner=x[[1]], task = task, resampling = inner, par.set = x[[2]], control = ctrl,show.info=T)
    model <- setHyperPars(x[[1]], par.vals = tune$x)
    return(model)
  }
  
  modelcv<-lapply(mset,FUN=modeltune,task=train.task,inner=inner, ctrl=ctrl)
  if(split==F){
    set.seed(seed)
    res = benchmark(modelcv, train.task, rdesc, keep.extract = T)
    outmeasure<-print(res)
    #fit the best-performing model, get predictions and conduct feature selection
    
    if(all(outmeasure$cindex.test.mean %in%  NaN)){
      res<-NULL
    }else{
      
      finalmodel <-modelcv[[which.max(outmeasure$cindex.test.mean)]]
      modeltrain <- mlr::train(learner=finalmodel, task=train.task)
      if(finalmodel$id=="RandomForestSRC"){
        modeltrain <- getLearnerModel(modeltrain)
        selvar<-max.subtree(modeltrain)$topvars
        listidx<-as.data.frame(max.subtree(modeltrain)$order)
        listidx<-listidx[selvar,1,drop=F]
        names(listidx)<-"coeff"
        train.chf<-modeltrain$chf
        
        clinical$Risk<-rowSums(train.chf)
      }else if(finalmodel$id=="Glmboost"){
        modelprediction.train <- predict(modeltrain, train.task)
        clinical$Risk<-modelprediction.train$data$response
        
        modeldata <- getLearnerModel(modeltrain)
        selvar<-as.data.frame(coef(modeldata))
        selvar<-selvar[selvar[,1]!=0,,drop=F]
        listidx<-as.data.frame(selvar)[-1,]
        names(listidx)<-"coeff"
        selvar<-row.names(selvar)[-1]
      } else if(finalmodel$id=="ElasticNnet"){
        modelprediction.train <- predict(modeltrain, train.task)
        
        clinical$Risk<-modelprediction.train$data$response
        modeldata <- getLearnerModel(modeltrain)
        selvar <-
          as.data.frame(as.matrix(coef(modeldata,
                                       s = finalmodel$par.vals$s,
                                       alpha =finalmodel$par.vals$alpha)))
        selvar<-selvar[selvar[,1]!=0,,drop=F]
        listidx<-as.data.frame(selvar)
        names(listidx)<-"coeff"
        selvar<-row.names(selvar)
      } else if(finalmodel$id== "Ridge"){
        modelprediction.train <- predict(modeltrain, train.task)
        
        clinical$Risk<-modelprediction.train$data$response
        modeldata <- getLearnerModel(modeltrain)
        selvar <-
          as.data.frame(as.matrix(coef(modeldata,
                                       s = finalmodel$par.vals$s)))
        selvar<-selvar[selvar[,1]!=0,,drop=F]
        listidx<-as.data.frame(selvar)
        names(listidx)<-"coeff"
        
        selvar<-row.names(selvar)
      }else if(finalmodel$id=="Lasso"){
        modelprediction.train <- predict(modeltrain, train.task)
        
        clinical$Risk<-modelprediction.train$data$response
        modeldata <- getLearnerModel(modeltrain)
        selvar <-
          as.data.frame(as.matrix(coef(modeldata,
                                       s = finalmodel$par.vals$s)))
        selvar<-selvar[selvar[,1]!=0,,drop=F]
        listidx<-as.data.frame(selvar)
        names(listidx)<-"coeff"
        
        selvar<-row.names(selvar)
      }else if(finalmodel$id=="Coxboost"){
        modelprediction.train <- predict(modeltrain, train.task)
        
        clinical$Risk<-modelprediction.train$data$response
        modeldata <- getLearnerModel(modeltrain)
        selvar <-
          as.data.frame(as.matrix(coef(modeldata,
                                       s = finalmodel$par.vals$s)))
        selvar<-selvar[selvar$V1!=0,,drop=F]
        listidx<-as.data.frame(selvar)
        names(listidx)<-"coeff"
        selvar<-row.names(selvar)
      }
      
      res<-list(res,selvar, modeltrain,finalmodel$id,expres,clinical,listidx)
      names(res)<-c("Benchmark result","Selected variables", "Fitted model", "lrnid","expres","clinical","listidx")
      
    }}else{
      set.seed(seed)
      res = benchmark(modelcv, test.task, rdesc, keep.extract = T)
      outmeasure<-print(res)
      if(all(outmeasure$cindex.test.mean %in%  NaN)){
        res<-NULL
      }else{
        
        finalmodel <-modelcv[[which.max(outmeasure$cindex.test.mean)]]
        
        modeltrain <- mlr::train(learner=finalmodel, task=train.task)
        if(finalmodel$id=="RandomForestSRC"){
          modeltrain <- getLearnerModel(modeltrain)
          selvar<-max.subtree(modeltrain)$topvars
          listidx<-as.data.frame(max.subtree(modeltrain)$order)
          listidx<-listidx[selvar,1,drop=F]
          names(listidx)<-"coeff"
          model.predicted<-predict(modeltrain,testset)
          train.chf<-modeltrain$chf
          test.chf <- model.predicted$chf
          testsurv$Risk<-rowSums(test.chf)
          trainsurv$Risk<-rowSums(train.chf)
        } else if (finalmodel$id=="Glmboost"){
          modelprediction.test <- predict(modeltrain, test.task)
          modelprediction.train <- predict(modeltrain, train.task)
          
          testsurv$Risk<-modelprediction.test$data$response
          trainsurv$Risk<-modelprediction.train$data$response
          modeldata <- getLearnerModel(modeltrain)
          selvar<-as.data.frame(coef(modeldata))
          selvar<-selvar[selvar[,1]!=0,,drop=F]
          listidx<-as.data.frame(selvar)[-1,]
          names(listidx)<-"coeff"
          selvar<-row.names(selvar)[-1]
        } else if(finalmodel$id=="ElasticNet"){
          modelprediction.test <- predict(modeltrain, test.task)
          modelprediction.train <- predict(modeltrain, train.task)
          
          testsurv$Risk<-modelprediction.test$data$response
          trainsurv$Risk<-modelprediction.train$data$response
          modeldata <- getLearnerModel(modeltrain)
          selvar <-
            as.data.frame(as.matrix(coef(modeldata,
                                         s = finalmodel$par.vals$s,
                                         alpha = finalmodel$par.vals$alpha)))
          selvar<-selvar[selvar[,1]!=0,,drop=F]
          listidx<-as.data.frame(selvar)
          names(listidx)<-"coeff"
          selvar<-row.names(selvar)
        } else if(finalmodel$id=="Ridge"){
          modelprediction.test <- predict(modeltrain, test.task)
          modelprediction.train <- predict(modeltrain, train.task)
          
          testsurv$Risk<-modelprediction.test$data$response
          trainsurv$Risk<-modelprediction.train$data$response
          modeldata <- getLearnerModel(modeltrain)
          selvar <-
            as.data.frame(as.matrix(coef(modeldata,
                                         s = finalmodel$par.vals$s)))
          selvar<-selvar[selvar[,1]!=0,,drop=F]
          listidx<-as.data.frame(selvar)
          names(listidx)<-"coeff"
          selvar<-row.names(selvar)
        }else if (finalmodel$id=="Lasso"){
          modelprediction.test <- predict(modeltrain, test.task)
          modelprediction.train <- predict(modeltrain, train.task)
          
          testsurv$Risk<-modelprediction.test$data$response
          trainsurv$Risk<-modelprediction.train$data$response
          modeldata <- getLearnerModel(modeltrain)
          selvar <-
            as.data.frame(as.matrix(coef(modeldata,
                                         s = finalmodel$par.vals$s)))
          selvar<-selvar[selvar[,1]!=0,,drop=F]
          listidx<-as.data.frame(selvar)
          names(listidx)<-"coeff"
          
          selvar<-row.names(selvar)
        } else if(finalmodel$id == "Coxboost"){
          
          modelprediction.test <- predict(modeltrain, test.task)
          modelprediction.train <- predict(modeltrain, train.task)
          
          testsurv$Risk<-modelprediction.test$data$response
          trainsurv$Risk<-modelprediction.train$data$response
          modeldata <- getLearnerModel(modeltrain)
          selvar <-
            as.data.frame(as.matrix(coef(modeldata,
                                         s = finalmodel$par.vals$s)))
          selvar<-selvar[selvar$V1!=0,,drop=F]
          listidx<-as.data.frame(selvar)
          names(listidx)<-"coeff"
          selvar<-row.names(selvar)
          
        }
        
        res<-list(res,selvar, modeltrain, trainset, testset, finalmodel$id,expres,clinical,testsurv,trainsurv,listidx)
        names(res)<-c("Benchmark result","Selected variables", "Fitted model","Training set", "Test set","lrnid","expres","clinical","Testsurv","Trainsurv","listidx")
      }
    }
  return(res)
  
}

survROC<-function(data,bmtime,bmstatus,marker,method,predyear,cutpoint){
  require(survivalROC)
  
  data[complete.cases(data[,bmtime]) & data[,bmtime]>0,]
  if(max(data[,bmtime])>600){
    by<- 365
  } else{
    by<-12
  }
  y<-strsplit(predyear,"-",fixed = T)
  years<-c()
  yearid<-c()
  for(i in 1:length(y)){
    yearid[i]<-as.numeric(y[[i]][1])
    years[i]<-yearid[i]*by
  }
  
  sumROC<-list()
  for (i in 1:length(years)){
    sumROC[[i]] <- survivalROC(Stime = data[,bmtime],status = data[,bmstatus],marker = data[[marker]],
                               predict.time =years[i],method = method,span = 0.25*nrow(data)^(-0.20))
  }
  sumAUC<-list()
  for (i in 1:length(sumROC)){
    sumAUC[[i]]<-sumROC[[i]]$AUC
  }
  ROCdata<-c()
  for(i in 1:length(sumROC)){
    predict.time<-sumROC[[i]]$predict.time
    TP<-sumROC[[i]]$TP
    FP<-sumROC[[i]]$FP
    auc<-sumROC[[i]]$AUC
    tmp<-c(predict.time,TP,FP,auc)
    ROCdata<-rbind(ROCdata,tmp)
  }
  survivalROC_helper <- function(t) {
    survivalROC(Stime = data[,bmtime],status = data[,bmstatus],marker = data[[marker]],
                predict.time =t,method = method,span = 0.25*nrow(data)^(-0.20))
  }
  
  time<-years
  survivalROC_data <- data_frame(t = time) %>%
    mutate(survivalROC = map(t, survivalROC_helper),
           
           auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
           
           df_survivalROC = map(survivalROC, function(obj) {
             as_data_frame(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    unnest() %>%
    arrange(t, FP, TP)
  survivalROC_data1<-mutate(survivalROC_data,auc = sprintf("%.3f",auc))
  survivalROC_data1$years<-survivalROC_data1$t/by
  
  survivalROC_data1<-unite(survivalROC_data1,year, years,auc,sep = " year AUC: " )
  AUC =factor(survivalROC_data1$year)
  
  ROC.1<-sumROC[[which.max(sumAUC)]]
  
  dot <- data.frame(TP = ROC.1$TP[with(ROC.1, which.min((1-TP)^2+ FP^2))],
                    FP = ROC.1$FP[with(ROC.1, which.min((1-TP)^2+ FP^2))])
  dot <- rbind(c(1,0),dot)
  if(cutpoint==T){
    cutoff<- ROC.1$cut.values[with(ROC.1, which.min((1-TP)^2+ FP^2))]
    ROC.plot<-ggplot(survivalROC_data1,mapping = aes(x = FP, y = TP)) +
      geom_path(aes(color= AUC))+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
      theme_bw() +scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+xlab("False positive rate")+
      ylab("True positive rate") +
      theme(legend.position = c(0.7,0.2))+
      geom_path(mapping = aes(x = FP,y = TP),data = dot)+
      annotate("text",x = dot$FP[2] + 0.05,y = dot$TP[2],label = paste0("Cutoff: ",round(cutoff,3)))
  } else{
    ROC.plot<-ggplot(survivalROC_data1,mapping = aes(x = FP, y = TP)) +
      geom_path(aes(color= AUC))+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
      theme_bw() +scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+xlab("False positive rate")+
      ylab("True positive rate") +
      theme(legend.position = c(0.7,0.2))
  }
  return(ROC.plot)
}

km <- function(data,
               time,
               status,
               marker,
               groupby,
               ratio,
               value,
               high,
               low,
               survxlab,
               survP,
               survRT,
               survCI,
               color1,
               color2
){
  require(survminer)
  require(survival)
  data[complete.cases(data[,time]) & data[,time]>0,]
  data<-data[complete.cases(data[,marker]),]
  if(groupby=="Percentage"){
    data$KMgroup <-ifelse(
      data[[marker]] <= quantile(data[[marker]], ratio),
      paste(marker, low, sep = " "),
      paste(marker, high, sep = " ")
    )
  } else{
    data$KMgroup <-ifelse(
      data[[marker]] <= value,
      paste(marker, low, sep = " "),
      paste(marker, high, sep = " ")
    )
  }
  
  data$time <- data[, time]
  data$status <- data[, status]
  fit <- survfit(Surv(time, status) ~ KMgroup, data = data)
  level <- levels(factor(data$KMgroup))
  c <-
    ifelse(
      level[1] == paste(marker, low, sep = " "),
      color1,
      color2
    )
  
  d <-
    ifelse(
      level[2] == paste(marker, high, sep = " "),
      color2,
      color1
    )
  
  e <-
    ifelse(
      level[1] == paste(marker, low, sep = " "),
      low, 
      high 
    )
  
  f <-
    ifelse(
      level[2] == paste(marker, low, sep = " "),
      low,
      high 
    )
  
  ggsurvplot(
    fit,
    data = data,
    risk.table = survRT,
    risk.table.height = 0.3,
    risk.table.y.text = FALSE,
    risk.table.title = "",
    main =NULL,
    palette = c(c, d),
    pval = survP,
    pval.method = T,
    conf.int = survCI,
    risk.table.y.text.col = T,
    legend = c(0.8, 0.90),
    legend.title = "",
    xlab = survxlab,
    legend.labs = c(e, f)
  )
}

cox2<-function(data,time,status,feature,maxtick,varname){
  require(survival)
  require(forestplot)
  data$time<-data[,time]
  data$status<-data[,status]
  Surv<-Surv(data$time, data$status)
  feature<-paste(feature,collapse ="+")
  fomu<-as.formula(paste("Surv","~",feature,sep=''))
  fit<-coxph(fomu,data=data)
  sumfit<-summary(fit)
  a<-cbind(sumfit$conf.int[,1:4],sumfit$coefficients[,5])
  colnames(a)<-c("HR","exp(-coef)","LCI","UCI","Pvalue")
  multicox <- as.data.frame(a)
  multicox<-round(multicox,3)
  multicox$Pvalue<- format.pval(multicox$Pvalue, digits = 3, eps = 0.001)
  unicox <- c()
  for (i in 1:length(feature)) {
    fomu <- as.formula(paste("Surv~", feature[i], sep = ""))
    cox <- coxph(fomu, data = data)
    cox <- summary(cox)
    conf.int <- cox$conf.int
    coef1 <- cox$coef
    a <- cbind(conf.int, coef1[, 5])
    colnames(a) <- c("HR", "exp(-coef)", "LCI", "UCI", "Pvalue")
    unicox <- rbind(unicox, a)
  }
  unicox<-round(unicox,3)
  unicox<-as.data.frame(unicox)
  unicox$Pvalue<- format.pval(unicox$Pvalue, digits = 3, eps = 0.001)
  
  index<-intersect(rownames(multicox),row.names(unicox))
  
  
  if(is.null(varname)){
    unicox<-unicox[index,]
    multicox<-multicox[index,]
  }else{
    unicox<-unicox[index,]
    multicox<-multicox[index,]
    
    if(length(varname)==length(index)){
      row.names(unicox)<-varname
      row.names(multicox)<-varname}else{
        stop("Error: The length of variable names does not equal to the clinical variables included in the CoxPH model")
      }
  }
  
  fun<-function(x){
    x$summary<-paste(x$HR,"(",x$LCI,"-",x$UCI,", ",x[,5],")",sep="")
    return(x)
  }
  unicox<-fun(unicox)
  multicox<-fun(multicox)
  
  tab<-cbind(unicox,multicox)
  tab$summar<-paste(tab[,6],"\n",tab[,12],sep="")
  
  names(tab)<-c("HR.x","exp(-coef).x","LCI.x","UCI.x","P value.x","summary.x","HR.y","exp(-coef).y","LCI.y","UCI.y","P value.y","summary.y" ,"summary")
  tab$Row.names<-row.names(tab)
  
  tablecox<-subset(tab,select=c(summary.x,summary.y))
  tab<-subset(tab,select=c(Row.names,summary, HR.x,LCI.x,UCI.x,HR.y,LCI.y,UCI.y))
  
  names(tablecox)<-c("Univariate HR(95% CI,Pvalue)","Multivariable HR(95% CI,Pvalue)")
  
  cn<-c("Variable", "HR(95% CI, P value)",NA,NA,NA,NA,NA,NA)
  tab<-rbind(cn,tab)
  fun1<-function(x){
    x$HR.x<-as.numeric(x$HR.x)
    x$LCI.x<-as.numeric(x$LCI.x)
    x$UCI.x<-as.numeric(x$UCI.x)
    x$HR.y<-as.numeric(x$HR.y)
    x$LCI.y<-as.numeric(x$LCI.y)
    x$UCI.y<-as.numeric(x$UCI.y)
    return(x)
  }
  tab<-fun1(tab)
  
  plotcox<-grid.grabExpr(print(forestplot(
    tab[, 1:2],
    hrzl_lines = gpar(col = "#444444"),
    legend = c("Univariate", "Multivariable"),
    legend_args = fpLegend(pos = list(x = .4, y = 0.9)),
    fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
    boxsize = .15,line.margin = .26,
    mean = cbind(tab[, "HR.x"], tab[, "HR.y"]),
    lower = cbind(tab[, "LCI.x"], tab[, "LCI.y"]),
    upper = cbind(tab[, "UCI.x"], tab[, "UCI.y"]),
    zero = 1, ci.vertices = T,
    ci.vertices.height = 0.08,
    xticks = c(0, 1, maxtick),clip = c(0,maxtick),
    lineheight = unit(1.8, "cm"),
    col = fpColors(box = c("darkblue", "darkred"), zero = "gray"),
    xlab = "Hazard ratio",
    graphwidth = unit(50, "mm"),
    is.summary = c(TRUE, rep(FALSE, (nrow(tab) - 1))),
    txt_gp = fpTxtGp(cex = 1, label = gpar(fontfamily = "serif")),
    colgap = unit(3.5, "mm"),mar = unit(rep(-2, times = 4), "mm")
  )))
  
  sumcox<-sumfit
  res<-list(tablecox,plotcox,sumcox)#,return_data = T)
  names(res)<-c("tablecox","plotcox","sumcox")
  
  return(res)
}



cox<-function(data,time,status,feature,maxtick,varname,legend.pos){
  require(survival)
  require(forestplot)
  data$time<-data[,time]
  data$status<-data[,status]
  Surv<-Surv(data$time, data$status)
  feature1<-paste(feature,collapse ="+")
  fomu<-as.formula(paste("Surv","~",feature1,sep=''))
  fit<-coxph(fomu,data=data)
  sumfit<-summary(fit)
  a<-cbind(sumfit$conf.int[,1:4],sumfit$coefficients[,5])
  colnames(a)<-c("HR","exp(-coef)","LCI","UCI","P value")
  multicox <- as.data.frame(a)
  
  unicox <- c()
  for (i in 1:length(feature)) {
    fomu <- as.formula(paste("Surv~", feature[i], sep = ""))
    cox <- coxph(fomu, data = data)
    cox <- summary(cox)
    conf.int <- cox$conf.int
    coef1 <- cox$coef
    a <- cbind(conf.int, coef1[, 5])
    colnames(a) <- c("HR", "exp(-coef)", "LCI", "UCI", "P value")
    unicox <- rbind(unicox, a)
  }
  
  unicox<-as.data.frame(unicox)
  
  index<-intersect(rownames(multicox),row.names(unicox))
  unicox<-unicox[index,]
  multicox<-multicox[index,]
  plotcox<-coxforestp(unicox=unicox, multicox=multicox, legend.pos=legend.pos, xlim=maxtick, varname=varname)
  
  
  sumcox<-sumfit
  unicox<-round(unicox,3)
  unicox<-as.data.frame(unicox)
  unicox$Pvalue<- format.pval(unicox$`P value`, digits = 3, eps = 0.001)
  
  
  multicox<-round(multicox,3)
  
  multicox$Pvalue<- format.pval(multicox$`P value`, digits = 3, eps = 0.001)
  unicox$summaryx<-paste(unicox$HR,paste("(",unicox$LCI,"-",paste(paste(unicox$UCI,",",sep=""),unicox$Pvalue,sep=" "),")",sep=""),sep=" ")
  multicox$summaryy<-paste(multicox$HR,paste("(",multicox$LCI,"-",paste(paste(multicox$UCI,",",sep=""),multicox$Pvalue,sep=" "),")",sep=""),sep=" ")
  
  tab<-merge(unicox,multicox,by=0)
  
  rownames(tab)<-tab$Row.names
  
  tablecox<-subset(tab,select=c(summaryx,summaryy))
  names(tablecox)<-c("Univariate HR(95% CI,Pvalue)","Multivariable HR(95% CI,Pvalue)")
  res<-list(tablecox,plotcox,sumcox)#,return_data = T)
  names(res)<-c("tablecox","plotcox","sumcox")
  return(res)
}
cox1<-function(data,time,status,feature,maxtick,split,marker,varname,legend.pos){
  if (varname == ""){
    varname <- NULL
  } else{
    
    varname <- strsplit(varname, "|", fixed = T)[[1]]
  }
  feature <- c(feature,marker)
  if (split == F) {
    res<-cox(
      data = data$clinical,
      time = time,
      status = status,
      feature = feature,
      maxtick = maxtick,
      varname = varname,
      legend.pos=legend.pos
    )
  } else{
    res<-lapply(
      list(data$Trainsurv, data$Testsurv),
      FUN = cox,
      time = time,
      status = status,
      feature = feature,
      maxtick = maxtick,
      varname = varname,
      legend.pos=legend.pos
    )
  }
  return(res)
}


predm <- function(expclin,
                  expval,
                  time,
                  status,
                  lnid,
                  modeltrain) {
  
  expclin <- expclin[complete.cases(expclin[[time]]) & expclin[[time]] > 0, ]
  index <- intersect(row.names(expclin), names(expval))
  
  expval <- as.data.frame(t(expval))
  
  expval <- expval[index, ]
  expclin <- expclin[index, ]
  selvar<-modeltrain$features
  expval1<-expval
  if (all(selvar %in% names(expval))) {
    expval <- expval[, selvar]
    
  } else{
    stop(paste(
      setdiff(selvar, names(expval)),
      "is not in the expression profile you provided"
    ))
    
  }
  expval$time <- expclin[[time]]
  expval$status <- expclin[[status]]
  expval1$time <- expclin[[time]]
  expval1$status <- expclin[[status]]
  
  exp.task <- makeSurvTask(data = expval, target = c("time", "status"))
  
  if (lnid == "RandomForestSRC") {
    
    model.predicted <- predict(modeltrain, expval1)
    test.chf <- model.predicted$chf
    expclin$Risk <- rowSums(test.chf)
    
  } else{
    
    modelprediction.test <- predict(modeltrain, exp.task)
    expclin$Risk <- modelprediction.test$data$response
    
  }
  return(expclin)
}

#ui.R

dataset<-read.csv("dataset.csv",header=T,row.names=1)
header <- dashboardHeader(
  title = "CBioExplorer",
  titleWidth = 250
)

sidebar <-  dashboardSidebar(
  width = 250,
  sidebarMenu(
    id = "tabs",
    
    tags$div("Introduction",
             style= "font-size: 1.5em;
                 margin-top: 6px;
                 padding:0 0.25em;
                 text-align: center;
                 background: rgba(255, 255, 255, 0);
                 color: orange"),
    menuItem("Introduction", tabName = "welcome1", icon = icon("home"), selected = T),
    
    tags$div("Data",
             style= "font-size: 1.5em;
                 margin-top: 6px;
                 padding:0 0.25em;
                 text-align: center;
                 background: rgba(255, 255, 255, 0);
                 color: #FF5151"),
    menuSubItem("Data input", tabName = "dataset",icon=icon("database")),
    tags$div("Analysis",
             style = "margin-top: 6px;
                          font-size: 1.5em;
                          padding: 0 1.25em;
                          text-align: center;
                          background: rgba(255, 255, 255, 0);
                          color: #0AC71B"),
    
    menuItem("Dimensionality reduction", icon=icon("tree"),
             menuSubItem("WGCNA",tabName = "wgcna",icon=NULL),
             menuSubItem("Survival related genes",tabName = "msurv",icon=NULL),
             menuSubItem("Differentially expressed genes",tabName = "DEG",icon=NULL)
    ),
    menuItem("Benchmark experiment", icon=icon("laptop"),
             
             menuSubItem("Benchmark experiment",tabName = "nestr",icon=NULL)),
    
    menuItem("Prediction model", icon=icon("notes-medical"),
             menuSubItem("Construct model",tabName = "bmpm",icon=NULL),
             menuSubItem("Validate model",tabName = "valmo",icon=NULL),
             menuSubItem("Nomogram",tabName = "nomo",icon=NULL)
    ),
    menuItem("Clinical annotation", icon=icon("layer-group"),
             menuSubItem("Correlation with clinical features",tabName = "clinical",icon=NULL),
             menuSubItem("Kaplan-Meier curve",tabName = "KM",icon=NULL),
             menuSubItem("CoxPH model",tabName = "CoxPH",icon=NULL),
             menuSubItem("Time-dependent ROC",tabName = "SurvROC",icon=NULL),
             menuSubItem("Most correlated genes",tabName = "mcorgene",icon=NULL),
             menuSubItem("Correlation with specific gene",tabName = "gene",icon=NULL),
             menuSubItem("Gene expression in different groups",tabName = "genediff",icon=NULL),
             menuSubItem("Correlation with Immune infiltration",tabName = "immune",icon=NULL),
             menuSubItem("Correlation with stemness score",tabName = "stemness",icon=NULL)
    )
    
    ,
    menuItem("Biological annotation", icon=icon("dna"),
             menuSubItem("Biological annotation",tabName = "bioan",icon=NULL)
    ) ,
    
    tags$div("Source",
             style= "margin-top: 6px;
                         font-size: 1.5em;
                         padding: 0 1.25em;
                         text-align: center;
                         background: rgba(255, 255, 255, 0);
                        color: #909CFF"),
    
    menuSubItem("Dataset source", href = "https://liuxiaoping2020.github.io/CBioExplorerDatasource/",icon=icon("external-link")),
    menuSubItem("CBioExplorer standalone app", href = "https://gitee.com/liuxiaoping2020/CBioExplorer",icon=icon("external-link")),
    menuSubItem("R package reference", href = "https://liuxiaoping2020.github.io/CBioExplorer_reference/",icon=icon("external-link")),
    
    tags$div("Tutorial",
             style= "margin-top: 6px;
                         font-size: 1.5em;
                         padding: 0 1.25em;
                         text-align: center;
                         background: rgba(255, 255, 255, 0);
                        color: deepskyblue"),
    
    menuSubItem("Tutorial" , href = "https://github.com/liuxiaoping2020/CBioExplorer_tutorial/blob/main/CBioExplorer_tutorial.pdf", icon = icon("book")),
    
  tags$div("Contact",
               style= "margin-top: 6px;
                           font-size: 1.5em;
                           padding: 0 1.25em;
                           text-align: center;
                           background: rgba(255, 255, 255, 0);
                          color: lightgreen"),
      tags$div("Xing-Huan Wang",
               style= "margin-top: 6px;
                           font-size: 1em;
                           padding: 0 1.25em;
                           text-align: left;
                           background: rgba(255, 255, 255, 0);
                          color: white"),

      tags$div("Email: wangxinghuan@whu.edu.cn",
               style= "margin-top: 6px;
                           font-size: 1em;
                           padding: 0 1.25em;
                           text-align: left;
                           background: rgba(255, 255, 255, 0);
                          color: white"),
      br(),
 tags$div("Sheng Li",
               style= "margin-top: 6px;
                           font-size: 1em;
                           padding: 0 1.25em;
                           text-align: left;
                           background: rgba(255, 255, 255, 0);
                          color: white"),

      tags$div("Email: lisheng-znyy@whu.edu.cn",
               style= "margin-top: 6px;
                           font-size: 1em;
                           padding: 0 1.25em;
                           text-align: left;
                           background: rgba(255, 255, 255, 0);
                          color: white"),
br(),
      tags$div("Xiao-Ping Liu",
               style= "margin-top: 6px;
                           font-size: 1em;
                           padding: 0 1.25em;
                           text-align: left;
                           background: rgba(255, 255, 255, 0);
                          color: white"),

      tags$div("Email: liuxiaoping@whu.edu.cn",
               style= "margin-top: 6px;
                           font-size: 1em;
                           padding: 0 1.25em;
                           text-align: left;
                           background: rgba(255, 255, 255, 0);
                          color: white")
    
    
  )
)

body <-  dashboardBody(
  
  tabItems(
#intro.R
    
    
    tabItem(
      tabName = "welcome1",
      fluidRow(column(
        12,
        HTML('<img style="width: 25%; display: block; margin-left: auto; margin-right: auto;" src="logo.png"/>'),
        # br(),
        HTML('<div align="justify">
                          <h3>Introduction </h3>
                          CBioExplorer (Cancer Biomarker Explorer) was developed to facilitate researchers and clinicians to screen, characterize, annotate and translate cancer biomarkers from molecular level to clinical settings more comfortably with graphical user interfaces (GUI). The whole pipeline of CBioExplorer includes data collection, data curation, dimensionality reduction using three methods of WGCNA, univariate Cox proportional hazards regression model, differentially expressed gene analysis, benchmark experiment with 6 machine learning learners (Lasso, Ridge, Elastic net, Glmboost, Coxboost, Randomforest) using cross validation and nested cross validation based on R package mlr, prediction model construction using Cox proportional hazards regression model and nomogram, clinical annotation using a variety of clinical approaches, and biological annotation using over-representation analysis (ORA) and gene set enrichment analysis (GSEA). The overview of CBioExplorer is summarized below:
                          </div>'),
        br(),
        HTML('<img style="width: 85%; display: block; margin-left: auto; margin-right: auto;" src="overview.png"/>'),
        HTML('<div align="justify">
                          <h3>Notes</h3>
                          <ul>
                          <li>Thanks for considering CBioExplorer for your study. Due to the limited computing power of the server, when multiple users use CBioExplorer at the same time, the response of the program may become slower, please be patient and only click the button once and wait until one step done.</li>
                          <li>If you plan to use CBioExplorer for high-iterative nested cross validation calculations, we strongly recommend that you download the CBioExplorer source code to your R software for corresponding calculations. This is caused by the limited computing power of the server. We apologize to you for this.</li>
                          <li>The App will be disconnected from our server after an hour if there is no mouse action on the web browser.</li>
                          </ul>
                          </div>'),
        
        
        h3("Citation"),
        p("To be added"),
        h3("Licence"),
        p("Open source under GPLV3.0. Both CBioExplorer software and curated cancer gene expression data are free for academic and non-commercial use."),
        
        h3("Issue report & feedback"),
        p("The user guide for CBioExplorer can be found "
          , a("here.", href="https://github.com/liuxiaoping2020/CBioExplorer_tutorial/blob/main/CBioExplorer_tutorial.pdf")
          , style="padding-left: 0em"),
        p("If you have questions to raise or are experiencing difficulties using the CBioExplorer, please report it at "
          , a("CBioExplorer Github issues.", href="https://github.com/liuxiaoping2020/CBioExplorer/issues")
          , style="padding-left: 0em"),
        
        h3("Contact"),
        p("This web app is developed and maintained by Xiao-Ping Liu at Wang's Lab, Department of Urology, Zhongnan hospital of WUhan University"),
        p("If you have any questions, comments, or suggestions, please feel free to contact the developer at liuxiaoping@whu.edu.cn"),
        br(),
        
        p("All rights reserved.",style = "display: flex; flex-direction: column; align-items: center; width: 100%")
        
        
        
        
        
      )
      ),
      fluidRow(
        div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
            HTML('<i>Data input</i>'),
            actionButton(inputId = 'page_after_introduction',label = '',icon = icon('arrow-right'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
        )
      )
    ),
#datainput.R
    tabItem(tabName = "dataset",
            
            
            fluidRow(
              column(9,
                     box(
                       title="Overview of discovery data",width = NULL,
                       solidHeader=T,collapsible=T,collapsed=F,status="success",
                       bsAlert("inputmess"),
                       tabBox(width = NULL,
                              tabPanel("Clinical data", DT::dataTableOutput('viewclinical'),
                                       useShinyjs(),
                                       fluidRow(column(4,
                                                       hidden(div(id = "clinicalinput.table",
                                                                  downloadButton('saveclinicalinput', 'Download clinical data', class = "butt2")
                                                       ))
                                       ))
                              ),
                              tabPanel("Gene expression data",DT::dataTableOutput('viewexpres'),
                                       useShinyjs(),
                                       fluidRow(column(4,
                                                       hidden(div(id = "express.table",
                                                                  downloadButton('saveexpresinput', 'Download gene expression data', class = "butt2")
                                                       ))
                                       ))
                                       
                                       
                              )))
                     # ,
                     # box(
                     #   title="Overview of gene expression data",width = NULL, solidHeader=T,collapsible=T,collapsed=T,status="success",
                     # # bsCollapsePanel("Overview of gene expression data",style = "info",
                     #   DT::dataTableOutput('viewexpres'),
                     #   useShinyjs(),
                     #   fluidRow(column(4,
                     #                   hidden(div(id = "express.table",
                     #                              downloadButton('saveexpresinput', 'Download gene expression data', class = "butt2")
                     #                   ))
                     #   ))
                     # )
                     # )
              ),
              column(3,
                     box(width = NULL,status = "danger", solidHeader=T,title = "Input discovery set",
                         selectInput("inputType", "Dataset type",
                                     choices = c("Public dataset", "Customized dataset"),
                                     selected = "Public dataset"),
                         conditionalPanel(
                           condition = "input.inputType == 'Customized dataset'",
                           # box(title = "Upload a gene expression matrix", width = NULL,
                           # h3('Upload a gene expression matrix'),
                           # solidHeader = TRUE, collapsible = TRUE,
                           # fileInput("expfile" ,"Gene expression matrix (.csv) with gene in row and sample in column.", accept = ".csv")
                           # )
                           # ,
                           # bsTooltip("expfile", "A matrix (.csv) with gene in row and sample in column.","left"),
                           div(
                             div(
                               # edit1
                               style="width:85%; display:inline-block; vertical-align: middle;",
                               fileInput("expfile", label = h4("Gene expression matrix (.csv)"),
                                         accept = ".csv")
                             ),
                             div(
                               # edit2
                               style="display:inline-block; vertical-align: middle;",
                               bsButton("q1", label = "", icon = icon("question"),
                                        style = "default"
                               ),
                               bsPopover(id = "q1", title = NULL,
                                         content = paste0("A matrix (.csv) with gene in row and sample in column."),
                                         placement = "left",
                                         trigger = "hover",
                                         options = list(container = "body")
                               )
                             )
                           ),
                           div(
                             div(
                               # edit1
                               style="width:85%; display:inline-block; vertical-align: middle;",
                               fileInput("clin", label = h4("Clinical data matrix (.csv)"),
                                         accept = ".csv")
                             ),
                             div(
                               # edit2
                               style="display:inline-block; vertical-align: middle;",
                               bsButton("q2", label = "", icon = icon("question"),
                                        style = "default"
                               ),
                               bsPopover(id = "q2", title = NULL,
                                         content = paste0("A matrix (.csv) with sample in row and clinical feature in column."),
                                         placement = "left",
                                         trigger = "hover",
                                         options = list(container = "body")
                               )
                             )
                           )
                           
                           # box(title = "Upload clinical data",width = NULL,
                           # h4('Upload clinical data'),
                           # solidHeader = TRUE, collapsible = TRUE,
                           # fileInput("clin", "Clinical data matrix", accept = ".csv"),
                           # bsTooltip("clin", "A matrix (.csv) with sample in row and clinical feature in columns.","left")
                           # )
                         ),
                         
                         
                         conditionalPanel(
                           condition = "input.inputType == 'Public dataset'",
                           selectInput("cancer", "Cancer", choices = unique(dataset$Cancer),selected=NULL),
                           selectInput("accession", "Accession", choices = NULL,selected=NULL)
                         ),
                         actionButton("data",
                                      "Submit dataset",
                                      style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                      icon = icon("picture-o"))
                     ))
            ),
            
            ###################################################################################################################
            
            fluidRow(
              column(9,
                     box(
                       title="Overview of external validation data",width = NULL,
                       solidHeader=T,collapsible=T,collapsed=F,status="success",
                       
                       tabBox(width = NULL,
                              # tabPanel("Clinical data", DT::dataTableOutput('viewclinical'),
                              #          useShinyjs(),
                              #          fluidRow(column(4,
                              #                          hidden(div(id = "clinicalinput.table",
                              #                                     downloadButton('saveclinicalinput', 'Download clinical data', class = "butt2")
                              #                          ))
                              #          ))
                              # ),
                              
                              tabPanel("Clinical data of external set", DT::dataTableOutput('viewclinical2'),
                                       fluidRow(column(4,align="left",
                                                       hidden(div(id = "viewclinical2_wrapper",
                                                                  downloadButton('saveviewclinical2', 'Download table', class = "butt2")
                                                       ))
                                       ))
                              ),
                              
                              # tabPanel("Gene expression data",DT::dataTableOutput('viewexpres'),
                              #          useShinyjs(),
                              #          fluidRow(column(4,
                              #                          hidden(div(id = "express.table",
                              #                                     downloadButton('saveexpresinput', 'Download gene expression data', class = "butt2")
                              #                          ))
                              #          ))
                              #
                              #
                              # )
                              tabPanel("Gene expression data of external set", DT::dataTableOutput('viewexpress2'),
                                       fluidRow(column(4,align="left",
                                                       hidden(div(id = "viewexpress2_wrapper",
                                                                  downloadButton('saveviewexpress2', 'Download table', class = "butt2")
                                                       ))
                                       ))
                              )
                       ))
                     # ,
                     # box(
                     #   title="Overview of gene expression data",width = NULL, solidHeader=T,collapsible=T,collapsed=T,status="success",
                     # # bsCollapsePanel("Overview of gene expression data",style = "info",
                     #   DT::dataTableOutput('viewexpres'),
                     #   useShinyjs(),
                     #   fluidRow(column(4,
                     #                   hidden(div(id = "express.table",
                     #                              downloadButton('saveexpresinput', 'Download gene expression data', class = "butt2")
                     #                   ))
                     #   ))
                     # )
                     # )
              ),
              column(3,
                     box(
                       title = "Input validation set",
                       width = NULL,
                       status = "danger",
                       solidHeader = T,
                       collapsible = T,
                       collapsed = F,
                       
                       selectInput("inputType2", "Dataset type",
                                   choices = c("Public dataset", "Customized dataset"),
                                   selected = "Public dataset"),
                       conditionalPanel(
                         condition = "input.inputType2 == 'Customized dataset'",
                         
                         div(
                           div(
                             # edit1
                             style="width:85%; display:inline-block; vertical-align: middle;",
                             fileInput("expfile2", label = h4("Gene expression matrix (.csv)"),
                                       accept = ".csv")
                           ),
                           div(
                             # edit2
                             style="display:inline-block; vertical-align: middle;",
                             bsButton("q5", label = "", icon = icon("question"),
                                      style = "default"
                             ),
                             bsPopover(id = "q5", title = NULL,
                                       content = paste0("A matrix (.csv) with gene in row and sample in column."),
                                       placement = "left",
                                       trigger = "hover",
                                       options = list(container = "body")
                             )
                           )
                         ),
                         div(
                           div(
                             style="width:85%; display:inline-block; vertical-align: middle;",
                             fileInput("clin2", label = h4("Clinical data matrix (.csv)"),
                                       accept = ".csv")
                           ),
                           div(
                             style="display:inline-block; vertical-align: middle;",
                             bsButton("q6", label = "", icon = icon("question"),
                                      style = "default"
                             ),
                             bsPopover(id = "q6", title = NULL,
                                       content = paste0("A matrix (.csv) with sample in row and clinical feature in column."),
                                       placement = "left",
                                       trigger = "hover",
                                       options = list(container = "body")
                             )
                           )
                         )
                         
                         
                       ),
                       conditionalPanel(
                         condition = "input.inputType2 == 'Public dataset'",
                         selectInput("cancer2", "Cancer", choices = unique(dataset$Cancer),selected=NULL),
                         selectInput("accession2", "Accession", choices = NULL,selected=NULL)
                       ),
                       actionButton("data2",
                                    "Submit dataset",
                                    style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                    icon = icon("picture-o"))
                     )
              )
            ),
            
            fluidRow(
              div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
                  actionButton(inputId = 'page_before_datainput',label = '',icon = icon('arrow-left'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
                  HTML('<i>Introduction</i>')
              ),
              div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
                  HTML('<i>Dimensionality reduction: WGCNA</i>'),
                  actionButton(inputId = 'page_after_datainput',label = '',icon = icon('arrow-right'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
              )
            )
            
    ),
    
#ClinicalCorrelation.R
    tabItem(tabName = "clinical",
            fluidRow(column(
              9,
              box(
                title = "Table for the correlations between gene expression and clinical features ",
                width = NULL,
                solidHeader = T,
                collapsible = F,
                status = "success",
                align="center",
                htmlOutput("table1")
              )
            ),
            column(
              3,
              box(
                title = "Correlations between gene expression and clinical features",
                width = NULL,
                status = "danger",
                solidHeader = T,
                collapsible = F,
                
                selectizeInput(
                  "table1bygene",
                  label = "Official gene symbol",
                  choices = NULL,
                  multiple = F,
                  selected = "TP53"
                ),
                bsTooltip("table1bygene", "Input a gene with official gene symbol","left"),
                
                selectizeInput(
                  "tbgroupby",
                  label = "Group by",
                  choices = c("Percentage","Value"),
                  multiple = F,
                  selected = "Percentage"
                ),
                conditionalPanel(
                  condition = "input.tbgroupby == 'Percentage'",
                  sliderInput("grouppercent", "Cutoff percentage", min = 0, max = 1, value = 0.5),
                  bsTooltip("grouppercent", "Input a cutoff percentage (range:0-1) to categorize the samples into low and high expression group regarding to your interested gene. Example if 0.25: 0%-25% = Low, 25%-100% high","left")
                  
                ),
                conditionalPanel(
                  condition = "input.tbgroupby == 'Value'",
                  numericInput('tabgpvalue', "Cutoff value", value=1,step = 1),
                  bsTooltip("tabgpvalue", "Input a specific cutoff value to categorize the samples into low and high expression group regarding to your interested gene.","left"),
                  
                  
                ),
                selectizeInput(
                  "feature",
                  label = "Select clinical features",
                  choices = NULL,
                  multiple = T
                ),
                bsTooltip("feature", "Select clinical features you want to include in the table","left"),
                useShinyjs(),
                actionButton(
                  "table1bt",
                  "Draw table1",
                  style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                  icon = icon("picture-o")
                ))
              
              
            )),
            fluidRow(
              div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
                  actionButton(inputId = 'page_before_clinical',label = '',icon = icon('arrow-left'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
                  HTML('<i>Prediction model: Nomogram</i>')
              ),
              div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
                  HTML('<i>Kaplan-Meier curve</i>'),
                  actionButton(inputId = 'page_after_clinical',label = '',icon = icon('arrow-right'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
              )
            )
    ),

#KM.R
    tabItem(tabName = "KM",
            fluidRow(column(9,
                            box(
                              title = "Kaplan-Meier plot",
                              width = NULL,
                              solidHeader = T,
                              collapsible = T,
                              collapsed = F,
                              bsAlert("kmmess"),
                              status = "success",
                              align = "center",
                              uiOutput('KMplot')#,
                              
                              ,
                              fluidRow(column(4,align="left",
                                              useShinyjs(),
                                              hidden(div(id = "km_wrapper",
                                                         
                                                         splitLayout(
                                                           
                                                           
                                                           numericInput("kmwidth","Figure width",value = 10),
                                                           numericInput("kmheight","Figure height",value = 10)),
                                                         
                                                         downloadButton('downloadKM', 'Download figure', class = "butt2")
                                              ))
                              ))
                              
                            )
            ),
            column(3,box(
              title = "Kaplan-Meier plot",
              width = NULL,
              status = "danger",
              solidHeader = T,
              collapsible = T,
              
              selectizeInput(
                "KMgene",
                label = "Official gene symbol",
                choices = NULL,
                multiple = F,
                selected = "TP53"
              ),
              bsTooltip("KMgene", "Input a gene with official gene symbol","left"),
              selectizeInput(
                "kmgroupby",
                label = "Group by",
                choices = c("Percentage","Value"),
                multiple = F,
                selected = "Percentage"
              ),
              
              conditionalPanel(
                condition = "input.kmgroupby == 'Percentage'",
                sliderInput("genecut", "Cutoff percentage", min = 0, max = 1, value = 0.5),
                bsTooltip("genecut", "Input a cutoff (range:0-1) to categorize the samples into low and high expression group regarding to your interested gene. Example if 0.25: 0%-25% = Low, 25%-100% high","left"),
                
              ),
              conditionalPanel(
                condition = "input.kmgroupby == 'Value'",
                numericInput('kmgpvalue', "Cutoff value", value=1,step = 1),
                bsTooltip("kmgpvalue", "Input a specific cutoff value to categorize the samples into low and high expression group regarding to your interested gene.","left")
                
                
              ),
              selectizeInput(
                "survivaltime",
                label = "Survival time",
                choices = NULL,
                multiple = F,
                selected = "OS.time"
              ),
              bsTooltip("survivaltime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),
              
              selectizeInput(
                "survivalstatus",
                label = "Survival status",
                choices = NULL,
                multiple = F,
                selected = "OS"
              ),
              bsTooltip("survivalstatus", "Select survival time column. Example: OS, RFS, PFS","left"),
              box(title="Customized setting",width=NULL,solidHeader=F,collapsible=T,status=NULL,collapsed =T,
                  textInput("survxlab","Label of X axis",placeholder="Overall survival (months)",value="Overall survival (months)"),
                  bsTooltip("survxlab", "Define the label of X axis","left"),
                  checkboxInput("survP", "Show p-value?", value = TRUE, width = NULL),
                  checkboxInput("survRT", "Show risk table?", value = TRUE, width = NULL),
                  checkboxInput("survCI", "Show confidence interval?", value = F, width = NULL),
                  colourpicker::colourInput("kmcolor1", "Color of group 1", value = "#00AFBB", showColour = "background",closeOnClick = TRUE),
                  colourpicker::colourInput("kmcolor2", "Color of group 2", value = "#FC4E07", showColour = "background",closeOnClick = TRUE),
              ),
              box(title="Size control",width=NULL,solidHeader=F,collapsible=T,status=NULL,collapsed =T,
                  sliderInput("survwidth", "Plot Width (%)", min = 0, max = 100, value =50),
                  sliderInput("survheight", "Plot Height (px)", min = 0, max = 1200, value = 500)),
              useShinyjs(),
              actionButton(
                "KMplotbt",
                "Draw KM curve",
                style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                icon = icon("picture-o")
              )
              
            ))
            ),
            
            fluidRow(
              div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
                  actionButton(inputId = 'page_before_KM',label = '',icon = icon('arrow-left'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
                  HTML('<i>Correlation with clinical features</i>')
              ),
              div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
                  HTML('<i>CoxPH model</i>'),
                  actionButton(inputId = 'page_after_KM',label = '',icon = icon('arrow-right'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
              )
            )
    ),
    
#CoxPH.R
    tabItem(tabName = "CoxPH",
            fluidRow(
              column(9,
                     box(title="Cox proportional hazard regression model",width = NULL,
                         solidHeader=T,collapsible=T,status="success",#align="center",
                         bsAlert("CoxPHmess"),
                         tabBox(width = NULL,
                                tabPanel("Forestplot", align="center",
                                         
                                         uiOutput('Coxforestplot'),
                                         useShinyjs(),
                                         fluidRow(column(4,align="left",
                                                         hidden(div(id = "mybox_wrapper",
                                                                    splitLayout(
                                                                      numericInput("FPwidth","Figure width",value = 10),
                                                                      numericInput("FPheight","Figure height",value = 10)),
                                                                    downloadButton('saveforest', 'Download figure', class = "butt2")
                                                         ))
                                         ))
                                ),
                                
                                tabPanel("Table", dataTableOutput('Coxtable'),
                                         useShinyjs(),
                                         fluidRow(column(4,
                                                         hidden(div(id = "mybox_wrapper.table",
                                                                    downloadButton('saveforesttable', 'Download table', class = "butt2")
                                                         ))
                                         ))
                                ),
                                tabPanel("Summary", verbatimTextOutput("Coxsummary"),align="left",
                                         
                                         
                                         tags$head(
                                           tags$style(
                                             "#Coxsummary{ font-size:12px; font-style:Arial;width: 1000px; max-width: 215%;background: ghostwhite;}"
                                           )
                                         )
                                )
                         )
                     )
                     
              ),
              column(3,box(
                title = "Cox proportional hazards regression model",
                width = NULL,
                status = "danger",
                solidHeader = T,
                collapsible = T,
                
                selectizeInput(
                  "coxclinvar",
                  label = "Select clinical features",
                  choices = NULL,
                  multiple = T
                ),
                bsTooltip("coxclinvar", "Select clinical variables to included to CoxPH model","left"),
                selectizeInput(
                  "coxgene",
                  label = "Official gene symbol",
                  choices = NULL,
                  multiple = T
                ),
                
                bsTooltip("coxgene", "Input one or more genes with official gene symbol that you want to included in the CoxPH model.","left"),
                
                selectizeInput("CoxPHtime",label = "Survival time",choices = NULL,multiple = F,selected = "OS.time"),
                bsTooltip("CoxPHtime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),
                
                selectizeInput("CoxPHstatus",label = "Survival status",choices = NULL, multiple = F,selected = "OS"),
                bsTooltip("CoxPHstatus", "Select survival status column. Example: OS, RFS, PFS","left"),
                
                box(title = "Size control",width = NULL,
                    solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                    sliderInput("Coxwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
                    sliderInput("Coxheight", "Plot Height (px)", min = 0, max = 1000, value = 430),
                    sliderInput("coxfiglg", "Figure legend position", value = c(0.5, 0.9), min = 0, max = 1),
                    bsTooltip("coxfiglg", "Set the coordinates of the legend box. Their values should be between 0 and 1. c(0,0) corresponds to the 'bottom left' and c(1,1) corresponds to the 'top right' position","left"),
                    textInput("coxPHvarname", "Variable names", placeholder="Age,Gender,Stage,Grade..."),
                    bsTooltip("coxPHvarname", "Define the variable names you interested,the number of variable names should be identical with the number of variables you included in the CoxPH model and use '|' to separate multiple variable names","left"),
                    numericInput("maxtick", "Maximum of xticks", 5, min = 1, max = 20,step=1),
                    bsTooltip("maxtick", "Define the maximum of xticks and clip to adjust the size of the forestplot","left")
                    
                ),
                
                useShinyjs(),
                actionButton("CoxPHbt",
                             "Perform CoxPH model",
                             style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                             icon = icon("picture-o"))
              )
              
              )
            ),
            fluidRow(
              div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
                  actionButton(inputId = 'page_before_CoxPH',label = '',icon = icon('arrow-left'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
                  HTML('<i>Kaplan-Meier curve</i>')
              ),
              div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
                  HTML('<i>Time-dependent ROC</i>'),
                  actionButton(inputId = 'page_after_CoxPH',label = '',icon = icon('arrow-right'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
              )
            )
    ),
    
    
#SurvROC.R
    tabItem (tabName = "SurvROC",
             fluidRow(column(
               9,
               box(
                 title = "Time-dependent ROC analysis",
                 width = NULL,
                 height = NULL,
                 solidHeader = T,
                 collapsible = F,
                 status = "success",
                 align="center",
                 uiOutput("SurvROCplot"),
                 fluidRow(column(4,align="left",
                                 useShinyjs(),
                                 hidden(div(id = "survROC_wrapper",
                                            splitLayout(
                                              numericInput("survROCwidth.dl","Figure width",value = 10),
                                              numericInput("survROCheight.dl","Figure height",value = 10)),
                                            downloadButton('downloadsurvROC', 'Download figure', class = "butt2")
                                 ))
                 ))
               )
             ),
             column(
               3,
               box(
                 title = "Survival ROC plot",
                 width = NULL,
                 status = "danger",
                 solidHeader = T,
                 collapsible = T,
                 selectizeInput(
                   "SurvROCgene",
                   label = "Official gene symbol",
                   choices = NULL,
                   multiple = F,
                   selected= "TP53"
                 ),
                 bsTooltip("SurvROCgene", "Input a gene with official gene symbol", "left"),
                 selectizeInput(
                   "SurvROCtime",
                   label = "Survival time",
                   choices = NULL,
                   multiple = F,
                   selected= "OS.time"
                 ),
                 bsTooltip(
                   "SurvROCtime",
                   "Select survival time column. Example: OS.time, RFS.time, PFS.time",
                   "left"
                 ),
                 
                 selectizeInput(
                   "SurvROCstatus",
                   label = "Survival status",
                   choices = NULL,
                   multiple = F,
                   selected= "OS"
                 ),
                 bsTooltip(
                   "SurvROCstatus",
                   "Select survival status column. Example: OS, RFS, PFS",
                   "left"
                 ),
                 selectizeInput(
                   "predictyear",
                   label = "Prediction years",
                   choices = NULL,
                   multiple = T
                 ),
                 bsTooltip(
                   "predictyear",
                   "Define the time points in years you want to predict based on time dependent ROC analysis. the longest time point should not the max of survival (relapse) duration",
                   "left"
                 ),
                 selectInput('SurvROCmethod', 'Method for survival ROC', choices =c("KM","NNE"), selected = "NNE"),
                 checkboxInput("cutoff", "Show optimal cutoff ?", value = F, width = NULL),
                 
                 box(
                   title = "Size control",
                   width = NULL,
                   solidHeader = TRUE,
                   collapsible = TRUE,
                   collapsed = T,
                   sliderInput("SurvROCwidth","Plot Width (%)",min = 0,max = 100,value = 50),
                   sliderInput("SurvROCheight",
                               "Plot Height (px)",  min = 0,max = 1000,value = 410)
                 ),
                 useShinyjs(),
                 actionButton(
                   "SurvROCbt",
                   "Submit",
                   style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                   icon = icon("picture-o")
                   
                 )
               )
             )
             ),
             fluidRow(
               div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
                   actionButton(inputId = 'page_before_SurvROC',label = '',icon = icon('arrow-left'),
                                style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
                   HTML('<i>CoxPH model</i>')
               ),
               div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
                   HTML('<i>Most correlated genes</i>'),
                   actionButton(inputId = 'page_after_SurvROC',label = '',icon = icon('arrow-right'),
                                style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
               )
             )
             
    ),
    
#genecor.R
    tabItem(tabName = "gene",
            fluidRow(column(9,
                            box(
                              title = "The correlation between two genes",
                              width = NULL,
                              solidHeader = T,
                              collapsible = F,
                              status = "success",
                              align="center",
                              tabBox(width = NULL,
                                     
                                     tabPanel(
                                       "Scatter plot",
                                       uiOutput('genecorplot'),
                                       useShinyjs(),
                                       fluidRow(column(4,align="left",
                                                       hidden(div(id = "genecor_wrapper",
                                                                  splitLayout(
                                                                    numericInput("genecorwidthdl","Figure width",value = 10),
                                                                    numericInput("genecorheightdl","Figure height",value = 10)),
                                                                  downloadButton('downgenecor', 'Download figure', class = "butt2")
                                                       ))
                                       ))
                                     ),
                                     tabPanel(
                                       "Summary",
                                       verbatimTextOutput("genecorsummary"),align="left",
                                       tags$head(
                                         tags$style(
                                           "#genecorsummary{ font-size:12px; font-style:Arial;width: 1000px; max-width: 215%;background: ghostwhite;}"
                                         )
                                       )
                                     )
                                     
                              )
                            )
            ),
            column(
              3,
              box(
                title = "Gene-gene correlation analysis",width = NULL,status = "danger",solidHeader=T,
                collapsible = T,
                selectizeInput("gene1",label = "Gene 1",choices = NULL,multiple = F),
                selectizeInput("gene2",label = "Gene 2",choices = NULL,multiple = F),
                radioButtons(
                  "cormethods","Correlation methods",
                  c("Pearson's correlation" = "pearson",
                    "Spearman's correlation" = "spearman",
                    "Kendall's correlation" = "kendall"
                  )
                ),
                box(title = "Size control",width = NULL,
                    solidHeader = F, collapsible = TRUE, collapsed = TRUE,
                    sliderInput("genecorwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
                    sliderInput("genecorheight", "Plot Height (px)", min = 0, max = 1000, value = 400)),
                actionButton(
                  "genecorbt","Submit",
                  style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                  icon = icon("picture-o")
                )
              )
            )),
            fluidRow(
              div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
                  actionButton(inputId = 'page_before_gene',label = '',icon = icon('arrow-left'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
                  HTML('<i>Most correlated genes</i>')
              ),
              div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
                  HTML('<i>Gene expression in different groups</i>'),
                  actionButton(inputId = 'page_after_gene',label = '',icon = icon('arrow-right'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
              )
            )
    ),
    
    
#mcorgene.R
    tabItem(
      tabName = "mcorgene",
      fluidRow(column(9,box(
        title = "Most correlated genes",
        width = NULL,
        solidHeader = T,
        collapsible = T,
        status = "success",
        bsAlert("mcorgenemess"),
        align="center",
        tabBox(width = NULL,
               tabPanel("Bar plot", uiOutput('mcorbarplot'),
                        useShinyjs(),
                        fluidRow(column(4,align="left",
                                        hidden(div(id = "mcorgenebarplot_wrapper",
                                                   
                                                   splitLayout(
                                                     
                                                     numericInput("mcorgenebarwidthdl","Figure width",value = 10),
                                                     numericInput("mcorgenebarheightdl","Figure width",value = 10)),
                                                   downloadButton('downmcorgenebarplot', 'Download figure', class = "butt2")
                                        ))
                        ))
               ),
               tabPanel("Bubble plot", uiOutput('mcorbubplot') ,
                        useShinyjs(),
                        fluidRow(column(4,align="left",
                                        hidden(div(id = "mcorgenebubbleplot_wrapper",
                                                   splitLayout(
                                                     numericInput("mcorgenebubblewidthdl","Figure width",value = 10),
                                                     numericInput("mcorgenebubbleheightdl","Figure height",value = 10)),
                                                   downloadButton('downmcorgenebubbleplot', 'Download figure', class = "butt2")
                                        ))
                        ))
                        
                        
               ),
               tabPanel("Table", dataTableOutput('mcortable'),
                        useShinyjs(),
                        fluidRow(column(4,align="left",
                                        hidden(div(id = "mcorgenetable_wrapper",
                                                   downloadButton('downmcorgenetable', 'Download table', class = "butt2")
                                        ))
                        ))
                        
                        
               )
        )
        
      )),
      column(3,(box(
        title = "Most correlated gene analysis",
        width = NULL,
        status = "danger",
        solidHeader = T,
        collapsible = T,
        selectInput('mcorgene', 'Official gene symbol', choices = NULL, selected = NULL),
        bsTooltip("mcorgene", "Select an official gene symbol you interested in ","left"),
        selectizeInput('mcormethod', "Correlation method", choices =c("pearson","spearman","kendall") , selected = "pearson"),
        bsTooltip("mcormethod", "Specify the correction method","left"),
        checkboxInput("padjust", "Adjust P value ?", value = TRUE, width = NULL),
        conditionalPanel(
          condition = "input.padjust == true ",
          selectInput('mcorrectmethod', "Correction method", choices =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr") , selected = "fdr"),
          bsTooltip("mcormethod", "Specify the correction method for P values. For more details, please refer to https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust","left")
        ),
        numericInput("corsigcut","Significance criterion", value=0.05,max =0.1, step =0.001),
        bsTooltip("corsigcut", "Specify the criterion for significant correlation ","left"),
        box(title = "Siz control",width = NULL,#status = "info",
            solidHeader = F,collapsible = TRUE,collapsed=TRUE,
            box(title = "Bar plot",width = NULL,collapsible = TRUE,collapsed=TRUE,
                sliderInput("mcorbarwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
                sliderInput("mcorbarheight", "Plot Height (px)", min = 0, max = 1000, value = 430)),
            box(title = "Bubble plot",width = NULL,collapsible = TRUE,collapsed=TRUE,
                sliderInput("mcorbubwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
                sliderInput("mcorbubheight", "Plot Height (px)", min = 0, max = 1000, value = 430))),
        useShinyjs(),
        
        actionButton("mcorgenebt",
                     "Submit",
                     style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                     icon = icon("picture-o"))#,
        
      )
      )
      
      )
      ),
      
      fluidRow(
        div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
            actionButton(inputId = 'page_before_mcorgene',label = '',icon = icon('arrow-left'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
            HTML('<i>Time-dependent ROC</i>')
        ),
        div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
            HTML('<i>Correlation with specific gene</i>'),
            actionButton(inputId = 'page_after_mcorgene',label = '',icon = icon('arrow-right'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
        )
      )
    ),
    
#genediff.R
    tabItem(
      tabName="genediff",
      fluidRow(column(9,box(
        title="Gene expression in different groups",
        width = NULL,
        solidHeader = T,
        status = "success",
        align="center",
        bsAlert("genediffmess"),
        uiOutput('genediffplot'),
        useShinyjs(),
        fluidRow(column(4,align="left",
                        hidden(div(id = "genediff_wrapper",
                                   splitLayout(
                                     numericInput("genediffwidthdl","Figure width",value = 10),
                                     numericInput("genediffheightdl","Figure width",value = 10)),
                                   downloadButton('downloadgenediff', 'Download figure', class = "butt2")
                        ))
        ))
      )
      ),
      column(3,box(
        title="Gene expression in different groups",
        width = NULL,
        status = "danger",
        solidHeader = T,
        collapsible = T,
        selectizeInput('genediff', 'Official gene symbol', choices = NULL, selected = "A1BG"),
        bsTooltip("genediff", "Select an official gene symbol you are interested in ","left"),
        selectizeInput('genediffgroup', "Group",choices = NULL, selected = "Grade"),
        bsTooltip("genediffgroup", "Select a clinical variable to divide the gene expression into different groups","left"),
        checkboxInput("gdsubgroup", "Perform subgroup analysis ?", value = F, width = NULL),
        conditionalPanel(
          condition = "input.gdsubgroup == true ",
          selectizeInput('genediffsub', "Subgroup", choices =NULL , selected = "Sex"),
          bsTooltip("genediffsub", "Select a clinical variable as a subgroup indicator to compare the expression of the interested gene in subgroups","left")
        ),
        
        checkboxInput("genediffp", "Show P value ?", value = TRUE, width = NULL),
        conditionalPanel(
          condition = "input.genediffp == true",
          selectizeInput('genediffmethod', "Statistical test", choices=c("T", "Wilcoxon"), selected = NULL),
          bsTooltip("genediffmethod", "Select a statistical test for the comparisons, 'T' means 'T test', and 'wilcoxon' means 'wilcoxon test' ","left"),
          selectizeInput('comparisons', "Comparisons", choices=c("Pairwise comparison", "Comparisons to a reference"), selected = NULL),
          conditionalPanel(
            condition = "input.comparisons == 'Comparisons to a reference'",
            selectizeInput('comparefer', "Comparison reference", choices=NULL, selected = NULL),
          )
        ),
        selectizeInput('genediffplottype', "Plot type", choices=c("Box plot", "Bar plot", "Violin plot"), selected = "Box plot",multiple=F),
        
        box(title = "Size control",width = NULL,status = NULL,
            solidHeader = T,collapsible = TRUE,collapsed=TRUE,
            sliderInput("genediffbarwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
            sliderInput("genediffbarheight", "Plot Height (px)", min = 0, max = 1000, value = 430)
        ),
        
        actionButton("genediffbt",
                     "Submit",
                     style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                     icon = icon("picture-o"))
        
        
      )
      )
      
      ),
      fluidRow(
        div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
            actionButton(inputId = 'page_before_genediff',label = '',icon = icon('arrow-left'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
            HTML('<i>Correlation with specific gene</i>')
        ),
        div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
            HTML('<i>Correlation with Immune infiltration</i>'),
            actionButton(inputId = 'page_after_genediff',label = '',icon = icon('arrow-right'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
        )
      )
      
    ),
    #
    
#unicox.R
    tabItem(
      tabName = "msurv",fluidRow(column(9,box(
        title="Univariate CoxPH table",
        width = NULL,
        solidHeader = T,
        status = "success",
        align="center",
        bsAlert("msurvmess"),
        DT::dataTableOutput('unicoxtable'),
        fluidRow(column(4,align="left",
                        hidden(div(id = "unicoxtable_wrapper",
                                   downloadButton('downloadunicoxtable', 'Download table', class = "butt2")
                        ))
        ))
      )
      ),
      column(3,box(title="Univariate CoxPH",
                   width = NULL,
                   status = "danger",
                   solidHeader = T,
                   collapsible = T,
                   selectizeInput('msurvtime', "Survival time",choices = NULL, selected = "OS.time"),
                   bsTooltip("msurvtime", "Select the survival time column for univariate CoxPH analysis","left"),
                   selectizeInput('msurvstatus', "Survival status",choices = NULL, selected = "OS"),
                   bsTooltip("msurvstatus", "Select the survival Status column for univariate CoxPH analysis","left"),
                   # selectizeInput('msurvcor', "Adjusting method", choices =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr") , selected = "fdr"),
                   # bsTooltip("msurvcor", "Specify the correction method","left"),
                   actionButton("msurvbt",
                                "Submit",
                                style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                icon = icon("picture-o"))
      ))),
      hidden(div(id = "unicoxanalysis_wrapper",
                 fluidRow(column(9,
                                 box(
                                   title="Significant univariate CoxPH result",
                                   width = NULL,
                                   solidHeader = T,
                                   status = "success",
                                   align="center",
                                   tabBox(width=NULL,
                                          tabPanel(title="Forestplot",value="sigfor",uiOutput('sigfor'),
                                                   fluidRow(column(4,align="left",
                                                                   hidden(div(id = "unicoxForestplot_wrapper",
                                                                              splitLayout(
                                                                                numericInput("unicoxFPwidthdl","Figure width",value = 10),
                                                                                numericInput("unicoxFPheightdl","Figure height",value = 10)),
                                                                              downloadButton('downloadunicoxForestplot', 'Download figure', class = "butt2")
                                                                   ))
                                                   ))
                                          ),
                                          tabPanel(title="Significant CoxPH table",value="sigunt",DT::dataTableOutput('sigcoxtable'),
                                                   fluidRow(column(4,align="left",
                                                                   hidden(div(id = "sigcoxtable_wrapper",
                                                                              downloadButton('downloadsigcoxtable', 'Download table', class = "butt2")
                                                                   ))
                                                   ))
                                          )
                                   )
                                 )
                                 
                                 
                 ),
                 column(3,box(
                   title="Significant univariate CoxPH",
                   width = NULL,
                   status = "danger",
                   solidHeader = T,
                   collapsible = T,
                   selectizeInput('msurvcor', "Adjusting method", choices =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr") , selected = "fdr"),
                   bsTooltip("msurvcor", "Specify the correction method","left"),
                   selectizeInput('msurvsel', "Significance cutoff method", choices =c("P value", "Adjusted P value") , selected = "Adjusted P value"),
                   bsTooltip("msurvsel", "Specify the method to select to most signifcant survival related genes.","left"),
                   conditionalPanel(
                     condition = "input.msurvsel == 'Adjusted P value'",
                     numericInput("msurvap","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005)),
                   conditionalPanel(
                     condition = "input.msurvsel == 'P value'",
                     numericInput("msurvp","P value cutoff",0.05,min = 0,max = 1,step =0.005)
                   ),
                   box(title="Size control",width=NULL,solidHeader=F,collapsible=T,collapsed =T,
                       sliderInput("msurvwidth", "Width (px)", min = 0, max = 100, value = 50),
                       sliderInput("msurvheight", "Height (px)", min = 0, max = 3500, value = 700),
                       sliderInput("msurvxtick", "xticks", min = 0, max = 20, value = 5),
                       bsTooltip("msurvxtick", "Specified x-axis largest tick mark","left")
                   ),
                   
                   actionButton("msurvopbt",
                                "Submit",
                                style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                icon = icon("picture-o"))
                   
                 )
                 )
                 )
                 
      )
      ),
      
      fluidRow(
        div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
            actionButton(inputId = 'page_before_unicox',label = '',icon = icon('arrow-left'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
            HTML('<i>Dimensionality reduction: WGCNA</i>')
        ),
        div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
            HTML('<i>Dimensionality reduction: Differentially expressed genes</i>'),
            actionButton(inputId = 'page_after_unicox',label = '',icon = icon('arrow-right'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
            
        )
      )
    ),
    
#DEG.R
    tabItem(
      tabName = "DEG",
      fluidRow(
        column(9,
               box(
                 title = "Differentially expressed gene table",width = NULL,align="center",
                 bsAlert("DEGmsg"),
                 solidHeader=T,collapsible=T,status="success", DT::dataTableOutput('DEGtable'),
                 fluidRow(column(4,align="left",
                                 hidden(div(id = "DEG_wrapper",
                                            
                                            downloadButton('downloadDEGtable', 'Download table', class = "butt2")
                                 ))
                 ))
               )),
        
        column(3,box(title="Differentially expressed gene",width=NULL,
                     solidHeader=T,
                     collapsible=T,status="danger",collapsed = F,
                     # checkboxInput("filterDEG", "Filter low expression ?", value = F, width = NULL),
                     # conditionalPanel(condition = "input.filterDEG == true",
                     #                  numericInput(inputId = "filnumb",label = "Minimun CPM reads",value=1,min=0,max=50,step = 1),
                     #                  bsTooltip("filnumb", "Set the minimun CPM reads of genes to be retained for subsequent analysis, CPM reads row sums <n","left"),
                     #                  numericInput(inputId = "leasamp",label="Least number of samples",value=5,min=1,step=1),
                     #                  bsTooltip("leasamp", "Set the least samples with minimun CPM counts","left")
                     # ),
                     #
                     selectizeInput('DEGmethod', "Method", choices ="limma", selected = "limma"),
                     bsTooltip("DEGmethod", "'limma' means conducting moderated contrast t-test for each gene in limma","left"),
                     # conditionalPanel(
                     #   condition = "input.DEGmethod == 'edgeR-LRT' ",
                     #   selectizeInput('DEGnorm', "Normalization method", choices =c("TMM","TMMwsp","RLE","upperquartile","none"), selected = "TMM")
                     #
                     # ),
                     # conditionalPanel(
                     #   condition = "input.DEGmethod == 'edgeR-QLF'",
                     #   selectizeInput('DEGnorm', "Normalization method", choices =c("TMM","TMMwsp","RLE","upperquartile","none"), selected = "TMM")
                     #
                     # ),
                     selectizeInput('DEGfactor', "Factor/Condition", choices =NULL , selected = "Sex"),
                     bsTooltip("DEGfactor", "A categorical variable to define the group for DEG analysis","left"),
                     
                     selectizeInput('DEGpadj', "P adjust methods", choices =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none") , selected = "fdr"),
                     bsTooltip("DEGpadj", "Specify the correction method for P values. For more details, please refer to https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust","left"),
                     # box(
                     #  selectizeInput(
                     #    'mlDEGsel',
                     #    "DEG cutoff method",
                     #    choices = c("Adjusted P", "LogFC", "Adjusted P & LogFC"),
                     #    selected = "LogFC"
                     #  ),
                     #  conditionalPanel(
                     #    condition = "input.mlDEGsel == 'Adjusted P'",
                     #    numericInput("mlgeneselp","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005)),
                     #  conditionalPanel(
                     #    condition = "input.mlDEGsel == 'LogFC'",
                     #    numericInput("mlgeneselogFC","LogFC cutoff",2,min = 0,max = 1,step =0.005)
                     #  ),
                     #  conditionalPanel(
                     #    condition = "input.mlDEGsel == 'Adjusted P & LogFC'",
                     #    numericInput("mlgeneselp1","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005),
                     #    numericInput("mlgeneselogFC1","LogFC cutoff",0.05,min = 0,max = 1,step =0.005)
                     #  # )
                     # ),
                     actionButton("DEGbt",
                                  "Submit",
                                  style = "background-color: #000080;
                                          color: #FFFFFF;
                                          margin-left: auto;
                                          margin-right: auto;
                                          width: 100%",
                                  icon = icon("picture-o"))
        )
        )
      ),
      
      hidden(div(id = "DEG1_wrapper",
                 fluidRow(column(9,
                                 box(
                                   title = "DEG output",
                                   width = NULL,
                                   solidHeader = T,
                                   collapsible = T,
                                   collapsed = F,
                                   status = "success",
                                   bsAlert("DEGvismess"),
                                   align="center",
                                   tabBox(
                                     width=NULL,
                                     tabPanel(title="Output genes",value="DEGog",DT::dataTableOutput('DEGog'),
                                              fluidRow(column(4,align="left",
                                                              hidden(div(id = "DEGog_wrapper",
                                                                         
                                                                         downloadButton('downloadopgene', 'Download table', class = "butt2")
                                                              ))
                                              ))
                                              
                                              
                                     ),
                                     tabPanel(title="DEG Visualization",value="DEGvis",uiOutput('DEGvis'),
                                              fluidRow(column(4,align="left",
                                                              hidden(div(id = "DEGvis_wrapper",
                                                                         splitLayout(
                                                                           numericInput("DEGviswidthdl","Figure width",value = 10),
                                                                           numericInput("DEGvisheightdl","Figure height",value = 10)),
                                                                         downloadButton('downloadDEGvis', 'Download figure', class = "butt2")
                                                              ))
                                              ))
                                              
                                     )
                                     
                                   )
                                   
                                   
                                 )
                 ),column(3,
                          box(title="DEG output",solidHeader=T,collapsible=T,width=NULL,collapsed = F,status = "danger",
                              selectizeInput(
                                'mlDEGsel',
                                "DEG cutoff method",
                                choices = c("Adjusted P", "LogFC", "Adjusted P & LogFC"),
                                selected = "LogFC"
                              ),
                              conditionalPanel(
                                condition = "input.mlDEGsel == 'Adjusted P'",
                                numericInput("mlgeneselp","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005)),
                              conditionalPanel(
                                condition = "input.mlDEGsel == 'LogFC'",
                                numericInput("mlgeneselogFC","LogFC cutoff",2,min = 0,max = 100,step =0.5)
                              ),
                              conditionalPanel(
                                condition = "input.mlDEGsel == 'Adjusted P & LogFC'",
                                numericInput("mlgeneselp1","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005),
                                numericInput("mlgeneselogFC1","LogFC cutoff",2,min = 0,max = 100,step =0.5)
                              ),
                              selectizeInput(
                                'DEGvismeth',
                                "Visualization method",
                                choices = c("Heatmap", "Volcano plot", "MA plot","Adjusted P plot"),
                                selected = "Heatmap"
                              ),
                              bsTooltip("DEGvismeth", "Choose visualization method for DEGs, which includes Heatmap, Volcano plot, MA plot,Adjusted P plot.","left"),
                              
                              conditionalPanel(
                                condition= "input.DEGvismeth == 'Heatmap'",
                                
                                bsCollapse(id = "collheat",open=NULL,
                                           bsCollapsePanel("Basic setting",style = "info",
                                                           textInput("heatname", "Heatmap name", value = "Expression", width = NULL, placeholder = NULL),
                                                           colourpicker::colourInput("heatcol1", "Min colour", "blue",closeOnClick = TRUE,allowTransparent = TRUE),
                                                           colourpicker::colourInput("heatcol2", "Median colour", "white",closeOnClick = TRUE,allowTransparent = TRUE),
                                                           colourpicker::colourInput("heatcol3", "Max colour", "red",closeOnClick = TRUE,allowTransparent = TRUE),
                                                           bsTooltip("heatname", "Name of the heatmap. By default the heatmap name is used as the title of the heatmap legend.","left"),
                                                           checkboxInput("heatnorm", "Normalize the heatmap ?", value = TRUE, width = NULL),
                                                           conditionalPanel(
                                                             condition = "input.heatnorm == true ",
                                                             selectizeInput('DEGscale', "Normalization method", choices =c("Scale","Center","Log","Z-score","0-1 normalization") , selected = "Scale"),
                                                             bsTooltip("DEGscale", "Select a normalzation method for the heatmap ","left")
                                                           )),
                                           bsCollapsePanel("Row setting",style = "info",
                                                           checkboxInput("ClusterR", "Cluster on rows ?", value = TRUE, width = NULL),
                                                           conditionalPanel(
                                                             condition = "input.ClusterR == true",
                                                             checkboxInput("cluster_row_slices", "Cluster on row slice ?", value = TRUE, width = NULL),
                                                             bsTooltip("cluster_row_slices", "If rows are split into slices, whether perform clustering on the slice means?","left"),
                                                             selectizeInput('cludistanrow', "Cluster distance on rows", choices =c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall") , selected ="euclidean"),
                                                             selectizeInput('clumethodrow', "Cluster method on rows", choices =c("ward.D", "ward.D2", "single", "complete", "average","mcquitty","median","centroid") , selected ="complete"),
                                                             selectizeInput('Rdend_side', "Row dendrogram side", choices =c("right","left") , selected ="left"),
                                                           ),
                                                           checkboxInput("showRname", "Show row names ?", value = TRUE, width = NULL),
                                                           conditionalPanel(
                                                             condition = "input.showRname == true",
                                                             selectizeInput('Rnameside', "Row name side", choices =c("right","left") , selected ="right"),
                                                           ),
                                                           checkboxInput("showFDR", "Show adjusted P value ?", value = F, width = NULL),
                                                           checkboxInput("showFC", "Show logFC ?", value = F, width = NULL)
                                           ),
                                           bsCollapsePanel("column title setting",style = "info",
                                                           # uiOutput('manySliders'),
                                                           numericInput("coltitsize","Column title font size",20,min = 1,max = 50,step = 5)
                                           ),
                                           bsCollapsePanel("Column  setting",style = "info",
                                                           checkboxInput("ClusterC", "Cluster on columns ?", value = TRUE, width = NULL),
                                                           conditionalPanel(
                                                             condition = "input.ClusterC == true",
                                                             checkboxInput("cluster_column_slices", "Cluster on column slices?", value = TRUE, width = NULL),
                                                             bsTooltip("cluster_column_slices", "If columns are split into slices, whether perform clustering on the slice means?","left"),
                                                             selectizeInput('cludistancol', "Cluster distance on rows", choices =c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall") , selected ="euclidean"),
                                                             selectizeInput('clumethodcol', "Cluster method on rows", choices =c("ward.D", "ward.D2", "single", "complete", "average","mcquitty","median","centroid") , selected ="complete"),
                                                             selectizeInput('Cdend_side', "Row dendrogram side", choices =c("top","bottom") , selected ="top")
                                                           ),
                                                           checkboxInput("showCname", "Show column names ?", value = F, width = NULL),
                                                           conditionalPanel(
                                                             condition = "input.showCname == true",
                                                             selectizeInput('Cnameside', "Column name side", choices =c("top","bottom") , selected ="bottom"),
                                                           )
                                           )
                                           # ,
                                           #                             bsCollapsePanel("Size control",style = "info",
                                           #                                             sliderInput("heatwidth", "Heatmap Width (%)", min = 0, max = 1000, value = 840),
                                           #                                             sliderInput("heatheight", "Heatmap Height (px)", min = 0, max = 1000, value = 430)
                                )
                                # )
                              ),
                              conditionalPanel(
                                condition =  "input.DEGvismeth == 'Volcano plot'",
                                bsCollapse(id = "collvolca",open=NULL,
                                           bsCollapsePanel("Cutoff for significance",style = "info",
                                                           numericInput("vocacutp","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005),
                                                           numericInput("vocacutfc","LogFC cutoff",2,min = 0,max = 1,step =0.005)
                                           ),
                                           bsCollapsePanel("Axis label",style = "info",
                                                           textInput("vocaXlab", label = "X-axis label", value =  "Default"),
                                                           textInput("vocaYlab", label = "Y-axis label", value = "Default"),
                                                           selectInput("legendPosition","Legend position", choices =c("top","bottom","left", "right") , selected ="top")
                                           )
                                           
                                           # bsCollapsePanel("Size control",style = "info",
                                           #                 sliderInput("Volcanowidth", "Heatmap Width", min = 0, max = 1000, value = 450),
                                           #                 sliderInput("Volcanoheight", "Heatmap Height (px)", min = 0, max = 1000, value = 430)
                                )
                                
                              ),
                              conditionalPanel(
                                condition= "input.DEGvismeth == 'MA plot'",
                                bsCollapse(id = "collMA",open=NULL,
                                           bsCollapsePanel("Cutoff for significance",style = "info",
                                                           numericInput("MAcutp","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005),
                                                           numericInput("MAfc","LogFC cutoff",2,min = 0,max = 1,step =0.005)
                                           ),
                                           bsCollapsePanel("Gene selection",style = "info",
                                                           selectInput("MASelmeth","Selection methods", choices =c("Top genes","Cutoff","Gene symbol"), selected ="Top genes"),
                                                           numericInput("Topgene","Top gene",20,min = 0,max = Inf,step =1),
                                                           selectInput("MAstopmeth","Select top method", choices =c("Adjusted P","LogFC"), selected ="Adjusted P"),
                                                           selectizeInput("MAgenesym","Gene symbol", choices =NULL , selected =NULL,multiple = T)
                                           ),
                                           bsCollapsePanel("Axis label",style = "info",
                                                           textInput("MAXlab", label = "X-axis label", value =  "Log2 mean expression"),
                                                           textInput("MAYlab", label = "Y-axis label", value = "Log2 fold change"),
                                                           selectInput("MAlegendPosition","Legend position", choices =c("top","bottom","left", "right") , selected ="top")
                                           )
                                           # ,
                                           # bsCollapsePanel("Size control",style = "info",
                                           #                 sliderInput("MAwidth", "MA-plot Width", min = 0, max = 1000, value = 450),
                                           #                 sliderInput("MAheight", "MA-plot Height (px)", min = 0, max = 1000, value = 430)
                                )
                                
                                
                              ),
                              conditionalPanel(
                                condition= "input.DEGvismeth == 'Adjusted P plot'",
                                colourpicker::colourInput("padjcol", "Color of the histogram", value = "gray")
                              ),
                              box(title="Size control",width=NULL,solidHeader=F,collapsible=T,collapsed =T,
                                  sliderInput("padjwidth", "Width (px)", min = 0, max = 100, value =50),
                                  sliderInput("padjheight", "Height (px)", min = 0, max = 1000, value = 430))
                              ,
                              actionButton("DEGvisbt",
                                           "Submit",
                                           style = "background-color: #000080;
              color: #FFFFFF;
              margin-left: auto;
              margin-right: auto;
              width: 100%",
                                           icon = icon("picture-o"))
                          )
                 )
                 ))),
      
      fluidRow(
        div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
            actionButton(inputId = 'page_before_DEG',label = '',icon = icon('arrow-left'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
            HTML('<i>Dimensionality reduction: Survival related genes</i>')
        ),
        div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
            HTML('<i>Benchmark experiment</i>'),
            actionButton(inputId = 'page_after_DEG',label = '',icon = icon('arrow-right'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
            
        )
      )
    ),
    
    
    
#WGCNA.R
    tabItem(
      tabName = "wgcna",
      fluidPage(
        column(9,
               box(
                 title = "Data Preprocessing",
                 width = NULL,
                 height=NULL,
                 solidHeader = T,
                 collapsible = T,
                 collapsed = F,
                 status = "success",
                 align="center",
                 bsAlert("wgcnapremess"),
                 
                 uiOutput('sampletree'),
                 textOutput(outputId = "datacleandescrip"),
                 useShinyjs(),
                 fluidRow(column(4,align="left",
                                 hidden(div(id = "sampletree_wrapper",
                                            splitLayout(
                                              numericInput("sampletreewidthdl","Figure width",value = 10),
                                              numericInput("sampletreeheightdl","Figure height",value = 10)),
                                            downloadButton('downloadsampletree', 'Download figure', class = "butt2")
                                 ))
                 ))
               )
        ),
        column(3,
               box(title = "Data preprocess",
                   solidHeader = T,
                   collapsible = T, width = NULL, collapsed = F, status = "danger",
                   checkboxInput("setopvar", "Select most variable genes ?", value = T, width = NULL),
                   conditionalPanel(
                     condition = "input.setopvar == true ",
                     numericInput("topvar","Most variable genes to include",2000,min = 1,max = 15000,step =100),
                     bsTooltip("topvar", "Selection top n variable genes to be included for network constuction based on the variances of genes cross all samples","left")
                   ),
                   numericInput("ZK","Z.K",-2.5,min = -10,max = 10,step =1),
                   bsTooltip("ZK", "Cutoff in sample network for outlier detection, please refer to Horvath S (2011) Weighted Network Analysis. Applications in Genomics and Systems Biology. Springer Book. ISBN: 978-1-4419-8818-8 for more details","left"),
                   selectizeInput(
                     'wgvariable',
                     "Numeric clinical variable",
                     choices = NULL,
                     selected = NULL,
                     multiple=T
                   ),
                   bsTooltip("wgvariable", "Select 2 or more numeric clinical variables to be included for correlation analysis","left"),
                   box(title="Size control",solidHeader=F,collapsible=T,width=NULL,collapsed = T,
                       sliderInput("sampwidth", "Heatmap Width (%)", min = 0, max = 100, value = 100),
                       sliderInput("samheight", "Heatmap Height (px)", min = 0, max = 2000, value = 400)
                   ),
                   
                   actionButton("datacleanbt",
                                "Submit",
                                style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                icon = icon("picture-o"))
               )
        )
      ),
      hidden(div(id = "wgcna1_wrapper",
                 fluidPage(column(9,
                                  box(
                                    title = "Network construction and module detection",
                                    width = NULL,
                                    solidHeader = T,
                                    collapsible = T,
                                    collapsed = F,
                                    align="center",
                                    status = "success",
                                    tabBox(
                                      width=NULL,
                                      tabPanel(
                                        title = "Soft threshold selection",
                                        value = "sft",
                                        uiOutput('sftdistribution'),
                                        textOutput(outputId = "sftdiscrp"),
                                        useShinyjs(),
                                        fluidRow(column(4,align="left",
                                                        hidden(div(id = "sftdistribution_wrapper",
                                                                   splitLayout(
                                                                     numericInput("sftdistributionwidthdl","Figure width",value = 10),
                                                                     numericInput("sftdistributionheightdl","Figure height",value = 10)),
                                                                   downloadButton('downloadsftdistribution', 'Download figure', class = "butt2")
                                                        ))
                                        ))
                                      ),
                                      tabPanel(title="Module dendrogram",value="modeng",uiOutput('modengram'),
                                               fluidRow(column(4,align="left",
                                                               hidden(div(id = "modengram_wrapper",
                                                                          splitLayout(
                                                                            numericInput("modengramwidthdl","Figure width",value = 10),
                                                                            numericInput("modengramheightdl","Figure height",value = 10)),
                                                                          downloadButton('downloadmodengram', 'Download figure', class = "butt2")
                                                               ))
                                               ))
                                      )
                                    )
                                  )
                 ),
                 column(3,box(title="Network & module detection",
                              solidHeader=T,
                              collapsible=T,width=NULL,collapsed = F,status = "danger",
                              bsCollapse(id = "Networkmoduledetection",open="Network construction",
                                         bsCollapsePanel("Network construction",style = "info",
                                                         
                                                         
                                                         selectizeInput(
                                                           'networkType',
                                                           "NetworkType",
                                                           choices = c("unsigned", "signed", "signed hybrid"),
                                                           selected = "unsigned",
                                                           multiple=F),
                                                         bsTooltip("networkType", "network type. Allowed values are (unique abbreviations of) 'signed','unsigned','signed hybrid'. More detail can be found at https://cran.r-project.org/web/packages/WGCNA/WGCNA.pdf","left"),
                                                         numericInput("RsquaredCut","Rsquared Cutoff",0.9,min = 0,max = 1,step =0.1),
                                                         bsTooltip("RsquaredCut", "Desired minimum scale free topology fitting index R^2.","left"),
                                                         selectizeInput(
                                                           'corType',
                                                           "Correlation Type",
                                                           choices = c("pearson", "bicor"),
                                                           selected = "pearson",
                                                           multiple=F
                                                         )
                                                         ,
                                                         bsTooltip("corType", "Specifying the correlation type, 'pearson' means 'Pearson's Correlation','bicor' means 'Bidweight midcorrelation' ","left"),
                                                         selectizeInput(
                                                           'TOMType',
                                                           "TOM Type",
                                                           choices = c("none", "unsigned", "signed", "signed Nowick", "unsigned 2", "signed 2","signed Nowick 2"),
                                                           selected ="unsigned",
                                                           multiple=F
                                                         ),
                                                         bsTooltip("TOMType", "Specifying the  Topology Overlap Matrix (TOM) type' ","left")
                                                         
                                         ),
                                         
                                         
                                         bsCollapsePanel("Module detection",style = "info",
                                                         numericInput("deepSplit","Deep split",2,min = 0,max = 4,step =1),
                                                         bsTooltip("deepSplit", "Provides a simplified control over how sensitive module detection should be to module splitting, with 0 least and 4 most sensitive' ","left"),
                                                         
                                                         numericInput("detectCutHeight","Maximum joining heights",0.995,min = 0,max = 1,step =0.001),
                                                         bsTooltip("detectCutHeight", "For method=='tree' it defaults to 0.99. For method=='hybrid' it defaults to 99% of the range between the 5th percentile and the maximum of the joining heights on the dendrogram' ","left"),
                                                         numericInput("minModuleSize","Minimum Module size",30,min = 1,max = 1000,step =5),
                                                         bsTooltip("minModuleSize", "Minimum module size for module detection' ","left"),
                                                         numericInput("reassignThreshold","P value  threshold", 0,min = 0,max = 1,step =0.0001),
                                                         bsTooltip("reassignThreshold", "p-value ratio threshold for reassigning genes between modules' ","left"),
                                                         numericInput("mergeCutHeight","Threshold to Merge Modules", 0.25,min = 0,max = 1,step =0.0001),
                                                         bsTooltip("mergeCutHeight", "Dendrogram cut height for module merging.","left"),
                                                         checkboxInput("numericLabels", "Should modules be labeled numbers ?", value = T, width = NULL),
                                                         bsTooltip("numericLabels", "Should the returned modules be labeled by colors (FALSE), or by numbers (TRUE)?","left"),
                                                         checkboxInput("pamStage", "PAM-like stage ?", value = F, width = NULL),
                                                         bsTooltip("pamStage", "Only used for method 'hybrid'. If TRUE, the second (PAM-like) stage will be performed.","left"),
                                                         checkboxInput("pamRespectsDendro", "PAM stage will respect the dendrogram ?", value = F, width = NULL),
                                                         bsTooltip("pamRespectsDendro", "Only used for method 'hybrid'. If TRUE, the PAM stage will respect the dendrogram in the sense that objects and small clusters will only be assigned to clusters that belong to the same branch that the objects or small clusters being assigned belong to.","left")
                                         ),
                                         
                                         bsCollapsePanel("Size control",style = "info",
                                                         bsCollapse(id = "WGCNASizecontrol",open=NULL,
                                                                    bsCollapsePanel("Soft threshold selection",style = "info",
                                                                                    sliderInput("sftwidth", "Heatmap Width (%)", min = 0, max = 100, value = 50),
                                                                                    sliderInput("sftheight", "Heatmap Height (px)", min = 0, max = 1000, value = 600)),
                                                                    bsCollapsePanel("Module dendrogram",style = "info",
                                                                                    sliderInput("MDwidth", "Heatmap Width (%)", min = 0, max = 100, value = 50),
                                                                                    sliderInput("MDheight", "Heatmap Height (px)", min = 0, max = 1000, value = 450))
                                                                    
                                                                    
                                                         )))
                              ,
                              
                              actionButton("WGCNAbt",
                                           "Submit",
                                           style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                           icon = icon("picture-o"))
                              
                 )
                 )
                 
                 )
      )),
      hidden(div(id = "wgcna2_wrapper",
                 fluidPage(column(9,
                                  box(
                                    title = "Module-trait relationships",
                                    width = NULL,
                                    solidHeader = T,
                                    collapsible = T,
                                    collapsed = F,
                                    status = "success",
                                    align="center",
                                    tabBox(
                                      width=NULL,
                                      tabPanel(title="Module - trait relationships",value="MTR",uiOutput('MTRplot'),
                                               fluidRow(column(4,align="left",
                                                               hidden(div(id = "MTR_wrapper",
                                                                          splitLayout(
                                                                            numericInput("MTRwidthdl","Figure width",value = 10),
                                                                            numericInput("MTRheightdl","Figure height",value = 10)),
                                                                          downloadButton('downloadMTR', 'Download figure', class = "butt2")
                                                               ))
                                               ))
                                      ),
                                      tabPanel(title="GS vs MM",value="GM",uiOutput('GSplot'),
                                               fluidRow(column(4,align="left",
                                                               hidden(div(id = "GM_wrapper",
                                                                          splitLayout(
                                                                            numericInput("GMwidthdl","Figure width",value = 10),
                                                                            numericInput("GMheightdl","Figure height",value = 10)),
                                                                          downloadButton('downloadGM', 'Download figure', class = "butt2")
                                                               ))
                                               ))
                                      ),
                                      tabPanel(title="Output genes",value="netsum",DT::dataTableOutput('netsum'),
                                               fluidRow(column(4,align="left",
                                                               hidden(div(id = "netsum_wrapper",
                                                                          downloadButton('downloadnetsum', 'Download table', class = "butt2")
                                                               ))
                                               ))
                                      )
                                    )
                                  )
                 ),
                 column(3,
                        box(title="Module-trait relationships",
                            solidHeader=T,
                            collapsible=T,width=NULL,collapsed =F,status = "danger",
                            selectizeInput(
                              'trait',
                              "Clinical trait",
                              choices = NULL,
                              selected = NULL,
                              multiple=F
                            ),
                            bsTooltip("trait", "Select the clinical trait you are interested and to identify correlated modules and genes","left"),
                            selectizeInput(
                              'module',
                              "Module",
                              choices = NULL,
                              selected = NULL,
                              multiple=F
                            ),
                            bsTooltip("module", "Select the module you are interested","left"),
                            selectizeInput(
                              'output',
                              "Output ",
                              choices = NULL,
                              selected = NULL,
                              multiple=F
                            ),
                            bsTooltip("output", "Output genes in modules for downstream analysis","left"),
                            box(title="Size control",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                                box(title="Module trait relationships",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                                    sliderInput("MTRwidth", "Heatmap Width (%)", min = 0, max = 100, value = 50),
                                    sliderInput("MTRheight", "Heatmap Height (px)", min = 0, max = 1000, value = 750),
                                    sliderInput("BMar", "Bottom margin", min = 0, max = 20, value = 5),
                                    sliderInput("LMar", "Left margin", min = 0, max = 20, value = 5),
                                    sliderInput("TMar", "Top margin", min = 0, max = 20, value = 2),
                                    sliderInput("RMar", "Right margin", min = 0, max = 20, value = 2)
                                ),
                                box(title="GS vs MM",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                                    sliderInput("GVMwidth", "Heatmap Width (%)", min = 0, max = 100, value = 50),
                                    sliderInput("GVMheight", "Heatmap Height (px)", min = 0, max = 1000, value = 600)
                                )
                            ),
                            actionButton("MTRbt",
                                         "Submit",
                                         style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                         icon = icon("picture-o"))
                        )
                 )
                 )
      )
      ),
      
      fluidRow(
        div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
            actionButton(inputId = 'page_before_WGCNA',label = '',icon = icon('arrow-left'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
            HTML('<i>Data input</i>')
        ),
        div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
            HTML('<i>Dimensionality redction: Survival related genes</i>'),
            actionButton(inputId = 'page_after_WGCNA',label = '',icon = icon('arrow-right'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
        )
      )
      
    ),
    
# nestr.R
    tabItem(
      tabName = "nestr",
      fluidPage(
        column(9,
               box(
                 title = "Benchmark experiment output",
                 width = NULL,
                 height=NULL,
                 solidHeader = T,
                 collapsible = T,
                 collapsed = F,
                 align="center",
                 status = "success",
                 bsAlert("nestrmess"),
                 tabBox(width = NULL,
                        tabPanel("Cindex comparison", align="center",
                                 uiOutput('modelcomp'),
                                 textOutput(outputId = "modeldescrip"),
                                 fluidRow(column(4,align="left",
                                                 hidden(div(id = "modelcomp_wrapper",
                                                            splitLayout(
                                                              numericInput("modelcompwidthdl","Figure width",value = 10),
                                                              numericInput("modelcompheightdl","Figure height",value = 10)),
                                                            downloadButton('downloadmodelcomp', 'Download figure', class = "butt2")
                                                 ))
                                 ))
                        ),
                        tabPanel("Biomarker output", dataTableOutput('biomarkout'),
                                 useShinyjs(),
                                 fluidRow(column(4,align="left",
                                                 hidden(div(id = "biomarkout_wrapper.table",
                                                            downloadButton('savebiomarkout', 'Download table', class = "butt2")
                                                 ))
                                 ))
                        )
                        
                 )
               )
        ),
        column(3,
               box(title="Benchmark experiment setting",
                   solidHeader=T,
                   collapsible=T,
                   width=NULL,
                   collapsed = F,
                   status = "danger",
                   selectizeInput(
                     'Learner',
                     "Survival learner",
                     choices = c("Lasso","ElasticNet","Ridge","Glmboost","Coxboost","RandomForestSRC"),
                     selected ="Lasso",
                     multiple=T
                   ),
                   bsTooltip("Learner", "Select the survivals you want to train you data","left"),
                   selectizeInput(
                     'validat',
                     "Method",
                     choices = c("Cross-validation","Nested cross-validation"),
                     selected ="Cross-validation",
                     multiple=F
                   ),
                   bsTooltip("validat", "Select a method for model and feature selection","left"),
                   textInput("seed","Random Seed",value="123",placeholder = "1234"),
                   bsTooltip("seed", "Set the random seed for benchmark experiment","left"),
                   selectizeInput("nrtime",label = "Survival time",choices = NULL,multiple = F,selected = "OS.time"),
                   bsTooltip("nrtime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),
                   selectizeInput("nrstatus",label = "Survival status",choices = NULL,multiple = F,selected = "OS"),
                   bsTooltip("nrstatus", "Select survival time column. Example: OS, RFS, PFS","left"),
                   selectizeInput(
                     'optimization',
                     "Optimization algorithm",
                     choices = c("Grid search","Random search"),
                     selected ="Grid search",
                     multiple=F
                   ),
                   bsTooltip("optimization", "Choose an optimization algorithm to search an appropriate set of parameters from a given search space for given learners.","left"),
                   conditionalPanel(
                     condition = "input.optimization == 'Grid search'",
                     numericInput("resolution","Resolution",10,min = 5,max = 50,step =1),
                     bsTooltip("resolution", "Resolution of the grid for each numeric/integer parameter in par.set. For vector parameters, it is the resolution per dimension. Either pass one resolution for all parameters, or a named vector. See ParamHelpers::generateGridDesign. Default is 10","left")
                   ),
                   conditionalPanel(
                     condition = "input.optimization == 'Random search'",
                     numericInput("maxit","Maxit",100,min = 10,max = 300,step =5),
                     bsTooltip("maxit", "Number of iterations for random search. Default is 100.","left")
                   ),
                   numericInput("ifold","Inner fold",5,min = 3,max = 20,step =1),
                   bsTooltip("ifold", "Set the fold number for inner K-fold cross-validation","left"),
                   conditionalPanel(
                     condition = "input.validat == 'Nested cross-validation'",
                     numericInput("ofold","Outer fold",5,min = 3,max = 20,step =1),
                     bsTooltip("ofold", "Set the fold number for outer K-fold cross-validation","left")
                   ),
                   conditionalPanel(
                     condition = "input.validat == 'Cross-validation'",
                     numericInput("biter","Bootstrap iterations",100,min = 50,max = 1000,step =100),
                     bsTooltip("biter", "Set the iteration for bootstrap validation","left")
                   ),
                   
                   checkboxInput("split", "Split the dataset ?", value = T, width = NULL),
                   bsTooltip("split", "Whether split the whole dataset into internal training and test set based on stratified sampling.","left"),
                   conditionalPanel(
                     condition = "input.split == true",
                     numericInput("sratio","Split ratio",0.6,min = 0,max = 1,step =0.1),
                     bsTooltip("sratio", "0-1. The percentage of samples that goes to training","left")
                   ),
                   box(title="Size control",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                       sliderInput("nrwidth", "Width (%)", min = 0, max = 100, value = 50),
                       sliderInput("nrheight", "Height (px)", min = 0, max = 1000, value = 450)
                   ),
                   actionButton("nrbt",
                                "Submit",
                                style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                icon = icon("picture-o"))
                   
               )
        )
      ),
      
      fluidRow(
        div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
            actionButton(inputId = 'page_before_nestr',label = '',icon = icon('arrow-left'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
            HTML('<i>Dimensionality reduction: Differentially expressed genes</i>')
        ),
        div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
            HTML('<i>Prediction model: Construct model</i>'),
            actionButton(inputId = 'page_after_nestr',label = '',icon = icon('arrow-right'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
            
        )
      )
    ),
#predmo.R
    tabItem(tabName = "bmpm",
            fluidRow(column(9,
                            box(
                              title = "Prediction model",
                              width = NULL,
                              solidHeader = T,
                              collapsible = F,
                              status = "success",
                              bsAlert("bmmess"),
                              align="center",
                              tabBox(width = NULL,
                                     tabPanel("Survival ROC curve", uiOutput('bmsurvROC'),
                                              useShinyjs(),
                                              fluidRow(column(4,align="left",
                                                              hidden(div(id = "bmsurvROC_wrapper",
                                                                         splitLayout(
                                                                           numericInput("bmsurvROCwidthdl","Figure width",value = 10),
                                                                           numericInput("bmsurvROCheightdl","Figure height",value = 10)),
                                                                         downloadButton('savebmsurvROC', 'Download figure', class = "butt2")
                                                              ))
                                              ))
                                     ),
                                     tabPanel("KM plot", uiOutput('bmKMplot'),
                                              fluidRow(column(4,align="left",
                                                              hidden(div(id = "bmKMplot_wrapper",
                                                                         splitLayout(
                                                                           numericInput("bmKMplotwidthdl","Figure width",value = 10),
                                                                           numericInput("bmKMplotheightdl","Figure height",value = 10)),
                                                                         downloadButton('savebmKMplot', 'Download figure', class = "butt2")
                                                              ))
                                              ))
                                     ),
                                     tabPanel("Forestplot",uiOutput('bmCoxforest'),
                                              fluidRow(column(4,align="left",
                                                              hidden(div(id = "bmCoxforest_wrapper",
                                                                         splitLayout(
                                                                           numericInput("bmCoxforestwidthdl","Figure width",value = 10),
                                                                           numericInput("bmCoxforestheightdl","Figure height",value = 10)),
                                                                         downloadButton('savebmCoxforest', 'Download figure', class = "butt2")
                                                              ))
                                              ))
                                     ),
                                     tabPanel("CoxPH table", dataTableOutput('bmCoxtable'),
                                              fluidRow(column(4,align="left",
                                                              hidden(div(id = "bmCoxtable_wrapper",
                                                                         downloadButton('savebmCoxtable', 'Download table', class = "butt2")
                                                              ))
                                              ))
                                     )#,
                                     
                                     
                              )
                            )
            ),
            column(3,box(
              title = "Prediction model",
              width = NULL,
              status = "danger",
              solidHeader = T,
              collapsible = T,
              selectizeInput(
                "bmtime",
                label = "Survival time",
                choices = NULL,
                multiple = F,
                selected = "OS.time"
              ),
              bsTooltip("bmtime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),
              
              selectizeInput(
                "bmstatus",
                label = "Survival status",
                choices = NULL,
                multiple = F,
                selected = "OS"
              ),
              bsTooltip("bmstatus", "Select survival status column. Example: OS, RFS, PFS","left"),
              
              box(title="Survival ROC setting",solidHeader = F,collapsible = T,width = NULL,collapsed = F,
                  selectizeInput(
                    "bmpredictyear",
                    label = "Prediction years",
                    choices = NULL,
                    multiple = T,
                    selected="1-year"
                  ),
                  bsTooltip(
                    "bmpredictyear",
                    "Define the time points in years you want to predict based on time dependent ROC analysis. The longest time point should not exceed the max of survival (relapse) duration",
                    "left"
                  ),
                  selectInput('bmSurvROCmethod', 'Method for survival ROC', choices =c("KM","NNE"), selected = "NNE"),
                  checkboxInput("bmcutoff", "Show optimal cutoff ?", value = F, width = NULL),
                  box(title = "Size control",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                      sliderInput("bmSurvROCwidth","Plot Width (%)",min = 0,max = 100, value = 50),
                      sliderInput("bmSurvROCheight","Plot Height (px)",min = 0,max = 1000,value = 600))
              ),
              
              box(title = "KM plot setting",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                  selectizeInput(
                    "bmkmgroupby",
                    label = "Group by",
                    choices = c("Percentage","Value"),
                    multiple = F,
                    selected = "Percentage"
                  ),
                  conditionalPanel(
                    condition = "input.bmkmgroupby == 'Percentage'",
                    sliderInput("riskcut", "Cutoff percentage", min = 0, max = 1, value = 0.5),
                    bsTooltip("riskcut", "Input a cutoff (range:0-1) to categorize the samples into low and high risk group. Example if 0.25: 0%-25% = Low, 25%-100% high","left"),
                    
                  ),
                  conditionalPanel(
                    condition = "input.bmkmgroupby == 'Value'",
                    numericInput('bmkmgpvalue', "Cutoff value", value=1,step = 1),
                    bsTooltip("bmkmgpvalue", "Input a specific cutoff value to categorize the samples into low and high risk group.","left")
                  ),
                  
                  box(title="Customized setting",width=NULL,solidHeader=F,collapsible=T,status=NULL,collapsed =T,
                      textInput("bmsurvxlab","Label of X axis",placeholder="Overall survival (months)",value="Overall survival (months)"),
                      bsTooltip("bmsurvxlab", "Define the label of X axis","left"),
                      checkboxInput("bmsurvP", "Show p-value?", value = TRUE, width = NULL),
                      checkboxInput("bmsurvRT", "Show risk table?", value = TRUE, width = NULL),
                      checkboxInput("bmsurvCI", "Show confidence interval?", value = F, width = NULL),
                      colourpicker::colourInput("bmkmcolor1", "Color of group 1", value = "#00AFBB", showColour = "background",closeOnClick = TRUE),
                      colourpicker::colourInput("bmkmcolor2", "Color of group 2", value = "#FC4E07", showColour = "background",closeOnClick = TRUE)),
                  box(title="Size control",width=NULL,solidHeader=F,collapsible=T,status=NULL,collapsed =T,
                      sliderInput("bmsurvwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
                      sliderInput("bmsurvheight", "Plot Height (px)", min = 0, max = 1200, value = 650))
              ),
              box(title = "CoxPH model setting",solidHeader = F,collapsible = T,width = NULL,collapsed = F,
                  selectizeInput("bmcoxclinvar",label = "Clinical covariates",choices = NULL,multiple = T),
                  bsTooltip("bmcoxclinvar", "Select clinical variables to included to CoxPH model","left"),
                  
                  box(title = "Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                      sliderInput("bmCoxwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
                      sliderInput("bmCoxheight", "Plot Height (px)", min = 0, max = 2500, value = 550),
                      textInput("bmcoxfeat", "Variable names", placeholder="Age|Gender|Stage|Grade"),
                      bsTooltip("bmcoxfeat", "Define the variable names you interested,the number of variable names should be identical with the number of variables you included in the CoxPH model and use '|' to separate multiple variable names","left"),
                      numericInput("bmmaxtick", "Maximum of xticks", 5, min = 1, max = 20,step=1),
                      bsTooltip("bmmaxtick", "Define the maximum of xticks and clip to adjust the size of the forestplot","left"),
                      sliderInput("bmlegend.pos", "Figure legend position", value = c(0.5, 0.9), min = 0, max = 1),
                      bsTooltip("bmlegend.pos", "Set the coordinates of the legend box. Their values should be between 0 and 1. c(0,0) corresponds to the 'bottom left' and c(1,1) corresponds to the 'top right' position","left")
                      
                  )
                  
              ),
              actionButton(
                "pmbt",
                "Submit",
                style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                icon = icon("picture-o")
              )
              
            ))
            ),
            
            fluidRow(
              div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
                  actionButton(inputId = 'page_before_bmpm',label = '',icon = icon('arrow-left'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
                  HTML('<i>Benchmark experiment</i>')
              ),
              div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
                  HTML('<i>Prediction model: Validate model</i>'),
                  actionButton(inputId = 'page_after_bmpm',label = '',icon = icon('arrow-right'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
                  
              )
            )
    ),
    
    
# valmo.R
    tabItem(
      tabName = "valmo",
      
      fluidRow(column(9,
                      box(
                        title = "Validation model",
                        width = NULL,
                        solidHeader = T,
                        collapsible = T,
                        status = "success",
                        bsAlert("valmess"),
                        align="center",
                        tabBox(width = NULL,
                               tabPanel("Survival ROC curve", uiOutput('valsurvROC'),
                                        fluidRow(column(4,align="left",
                                                        hidden(div(id = "valsurvROC_wrapper",
                                                                   splitLayout(
                                                                     numericInput("valsurvROCwidthdl","Figure width",value = 10),
                                                                     numericInput("valsurvROCheightdl","Figure height",value = 10)),
                                                                   downloadButton('savevalsurvROC', 'Download figure', class = "butt2")
                                                        ))
                                        ))
                               ),
                               tabPanel("KM plot", uiOutput('valKMplot'),
                                        fluidRow(column(4,align="left",
                                                        hidden(div(id = "valKMplot_wrapper",
                                                                   splitLayout(
                                                                     numericInput("valKMplotwidthdl","Figure width",value = 10),
                                                                     numericInput("valKMplotheightdl","Figure height",value = 10)),
                                                                   downloadButton('savevalvalKMplot', 'Download figure', class = "butt2")
                                                        ))
                                        ))
                               ),
                               tabPanel("Forestplot", uiOutput('valCoxforest'),
                                        fluidRow(column(4,align="left",
                                                        hidden(div(id = "valCoxforest_wrapper",
                                                                   splitLayout(
                                                                     numericInput("valCoxforestwidthdl","Figure width",value = 10),
                                                                     numericInput("valCoxforestheightdl","Figure height",value = 10)),
                                                                   downloadButton('savevalvalCoxforest', 'Download figure', class = "butt2")
                                                        ))
                                        ))
                               ),
                               tabPanel("CoxPH table", dataTableOutput('valCoxtable'),
                                        fluidRow(column(4,align="left",
                                                        hidden(div(id = "valCoxtable_wrapper",
                                                                   
                                                                   downloadButton('savevalCoxtable', 'Download table', class = "butt2")
                                                        ))
                                        ))
                               )
                        )
                      )),
               column(3,
                      box(
                        title = "Valiation model",
                        width = NULL,
                        status = "danger",
                        solidHeader = T,
                        collapsible = T,
                        selectizeInput(
                          "Valtime",
                          label = "Survival time",
                          choices = NULL,
                          multiple = F,
                          selected = "OS.time"
                        ),
                        bsTooltip("Valtime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),
                        
                        selectizeInput(
                          "Valstatus",
                          label = "Survival status",
                          choices = NULL,
                          multiple = F,
                          selected = "OS"
                        ),
                        bsTooltip("Valstatus", "Select survival status column. Example: OS, RFS, PFS","left"),
                        
                        box(title="Survival ROC setting",solidHeader = F,collapsible = T,width = NULL,collapsed = F,
                            selectizeInput(
                              "Valpredictyear",
                              label = "Prediction years",
                              choices = NULL,
                              multiple = T,
                              selected="1-year"
                            ),
                            bsTooltip(
                              "Valpredictyear",
                              "Define the time points in years you want to predict based on time dependent ROC analysis. the longest time point should not the max of survival (relapse) duration",
                              "left"
                            ),
                            selectInput('ValSurvROCmethod', 'Method for survival ROC', choices =c("KM","NNE"), selected = "NNE"),
                            checkboxInput("Valcutoff", "Show optimal cutoff ?", value = F, width = NULL),
                            box(title="Size control",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                                sliderInput("ValSurvROCwidth","Plot Width (%)",min = 0,max = 100, value = 50),
                                sliderInput("ValSurvROCheight","Plot Height (px)",min = 0,max = 1000,value = 600))
                        ),
                        
                        box(title = "KM plot setting",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                            
                            selectizeInput(
                              "Valkmgroupby",
                              label = "Group by",
                              choices = c("Percentage","Value"),
                              multiple = F,
                              selected = "Percentage"
                            ),
                            conditionalPanel(
                              condition = "input.Valkmgroupby == 'Percentage'",
                              sliderInput("Valriskcut", "Cutoff", min = 0, max = 1, value = 0.5),
                              bsTooltip("Valriskcut", "Input a cutoff (range:0-1) to categorize the samples into low and high risk groups.","left"),
                            ),
                            conditionalPanel(
                              condition = "input.Valkmgroupby == 'Value'",
                              numericInput('Valkmgpvalue', "Cutoff value", value=1,step = 1),
                              bsTooltip("Valkmgpvalue", "Input a specific cutoff value to categorize the samples into low and high risk group.","left")
                            ),
                            
                            box(title = "Customized setting",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                                textInput("Valsurvxlab","Label of X axis",placeholder="Overall survival (months)",value="Overall survival (months)"),
                                bsTooltip("Valsurvxlab", "Define the label of X axis","left"),
                                checkboxInput("ValsurvP", "Show p-value?", value = TRUE, width = NULL),
                                checkboxInput("ValsurvRT", "Show risk table?", value = TRUE, width = NULL),
                                checkboxInput("ValsurvCI", "Show confidence interval?", value = F, width = NULL),
                                colourpicker::colourInput("Valkmcolor1", "Color of group 1", value = "#00AFBB", showColour = "background",closeOnClick = TRUE),
                                colourpicker::colourInput("Valkmcolor2", "Color of group 2", value = "#FC4E07", showColour = "background",closeOnClick = TRUE)),
                            
                            box(title = "Size control",solidHeader = F,collapsible = T,width = NULL,collapsed = T,
                                sliderInput("Valsurvwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
                                sliderInput("Valsurvheight", "Plot Height (px)", min = 0, max = 1200, value = 650))
                        ),
                        box(title = "CoxPH model",solidHeader = F,collapsible = T,width = NULL,collapsed = F,
                            selectizeInput("Valcoxclinvar",label = "Clinical Covariates" ,choices = NULL,multiple = T),
                            bsTooltip("Valcoxclinvar", "Select clinical variables to included to CoxPH model","left"),
                            
                            box(title = "Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                sliderInput("ValCoxwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
                                sliderInput("ValCoxheight", "Plot Height (px)", min = 0, max = 2500, value = 430),
                                textInput("Valcoxfeat", "Variable names", placeholder="Age|Gender|Stage|Grade"),
                                bsTooltip("Valcoxfeat", "Define the variable names you interested,the number of variable names should be identical with the number of variables you included in the CoxPH model and use '|' to separate multiple variable names","left"),
                                numericInput("Valmaxtick", "maximum of xticks", 5, min = 1, max = 20,step=1),
                                bsTooltip("Valmaxtick", "Define the maximum of xticks and clip to adjust the size of the forestplot","left"),
                                sliderInput("valegend.pos", "Figure legend position", value = c(0.5, 0.9), min = 0, max = 1),
                                bsTooltip("valegend.pos", "Set the coordinates of the legend box. Their values should be between 0 and 1. c(0,0) corresponds to the 'bottom left' and c(1,1) corresponds to the 'top right' position","left")
                            )
                        ),
                        actionButton(
                          "Valbt",
                          "Submit",
                          style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                          icon = icon("picture-o")
                        )
                      )
               )
               # )
               # )
      ),
      fluidRow(
        div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
            actionButton(inputId = 'page_before_valmo',label = '',icon = icon('arrow-left'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
            HTML('<i>Prediction model: Construct model</i>')
        ),
        div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
            HTML('<i>Prediction model: Nomogram</i>'),
            actionButton(inputId = 'page_after_valmo',label = '',icon = icon('arrow-right'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
            
        )
      )
    ),
    
    
# nomo.R
    tabItem(
      tabName = "nomo",
      fluidRow(
        column(9,
               box(
                 title = "Nomogram",
                 width = NULL,height=NULL,
                 solidHeader = T,
                 collapsible = T,
                 bsAlert("nomomess"),
                 status = "success",
                 align="center",
                 uiOutput('nomoplot'),
                 fluidRow(column(4,align="left",
                                 hidden(div(id = "nomogram_wrapper",
                                            splitLayout(
                                              numericInput("nomogramwidthdl","Figure width",value = 10),
                                              numericInput("nomogramheightdl","Figure heigth",value = 6)),
                                            downloadButton('savenomogram', 'Download figure', class = "butt2")
                                 ))
                 ))
               )),
        column(3,
               box(
                 title = "Nomogram",
                 width = NULL,
                 status = "danger",
                 solidHeader = T,
                 collapsible = T,
                 selectizeInput(
                   "nomotime",
                   label = "Survival time",
                   choices = NULL,
                   multiple = F,
                   selected = "OS.time"
                 ),
                 bsTooltip("nomotime", "Select survival time column. Example: OS.time, RFS.time, PFS.time","left"),
                 
                 selectizeInput(
                   "nomostatus",
                   label = "Survival status",
                   choices = NULL,
                   multiple = F,
                   selected = "OS"
                 ),
                 bsTooltip("nomostatus", "Select survival status column. Example: OS, RFS, PFS","left"),
                 
                 selectizeInput(
                   "nomovar",
                   label = "Clinical variable",
                   choices = NULL,
                   multiple = T,
                   selected = "Risk"
                 ),
                 bsTooltip("nomovar", "Select survival variable you want to included in the nomogram","left"),
                 textInput("nomovarlab", "Variable names", placeholder="Age|Gender|Stage|Grade"),
                 bsTooltip("nomovarlab", "Define the variable names you interested,the number of variable names should be identical with the number of variables you included in the CoxPH model to contrust the nomogram and use '|' to separate multiple variable names","left"),
                 selectizeInput(
                   "nomoypoint",
                   label = "Prediction years",
                   choices = NULL,
                   multiple = T,
                   options = list(maxItems = 3),
                   selected="1-year"
                 ),
                 bsTooltip(
                   "nomoypoint",
                   "Define the time points in years you want to predict based on the nomogram. The longest time point should not exceed the max of survival (relapse) duration",
                   "left"
                 ),
                 box(title = "Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                     sliderInput("nomowidth", "Plot Width (%)", min = 0, max =100, value = 50),
                     sliderInput("nomoheight", "Plot Height (px)", min = 0, max = 2500, value = 430)
                 ),
                 actionButton(
                   "nomogrambt",
                   "Submit",
                   style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                   icon = icon("picture-o"))
               )
        )
      ),
      hidden(div(id = "invalcal_wrapper",
                 fluidRow(
                   column(9,
                          box(
                            title = "Internal validation and calibration",
                            width = NULL,height=NULL,
                            solidHeader = T,
                            collapsible = T,
                            bsAlert("innomomess"),
                            status = "success",
                            align="center",
                            tabBox(width = NULL,
                                   tabPanel("Internal validation", uiOutput('invalid'),
                                            textOutput(outputId = "invaliddescrip"),
                                            fluidRow(column(4,align="left",
                                                            hidden(div(id = "invalid_wrapper",
                                                                       splitLayout(
                                                                         numericInput("invalidwidthdl","Figure width",value = 10),
                                                                         numericInput("invalidheightdl","Figure height",value = 10)),
                                                                       downloadButton('saveinvalid', 'Download figure', class = "butt2")
                                                            ))
                                            ))),
                                   tabPanel("Internal calibration", uiOutput('incalib'),
                                            fluidRow(column(4,align="left",
                                                            hidden(div(id = "incalib_wrapper",
                                                                       splitLayout(
                                                                         numericInput("incalibwidthdl","Figure width",value = 10),
                                                                         numericInput("incalibheightdl","Figure heigth",value = 10)),
                                                                       downloadButton('saveincalib', 'Download figure', class = "butt2")
                                                            ))
                                            )))
                            )
                            
                          )
                   ),column(3,
                            box(
                              title = "Internal validation and calibration",
                              width = NULL,
                              status = "danger",
                              solidHeader = T,
                              collapsible = T,
                              selectizeInput(
                                "invalidypoint",
                                label = "Prediction years",
                                choices = NULL,
                                multiple = T,
                                options = list(maxItems = 3),
                                selected="1-year"
                              ),
                              bsTooltip(
                                "invalidypoint",
                                "Define the time points in years you want to predict based on the nomogram. The longest time point should not exceed the max of survival (relapse) duration",
                                "left"
                              ),
                              numericInput("inreps","Bootstrap iterations",1000,min = 50,max=3000,step = 100),
                              bsTooltip("inreps", "Number of bootstrap resamplings","left"),
                              numericInput("inratio","Ratio",0.8,min = 0.1,max=1,step = 0.1),
                              bsTooltip("inratio", "Ratio of resamplings","left"),
                              box(title = "C-index Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                  sliderInput("invalidwidth", "Plot Width (%)", min = 0, max =100, value = 50),
                                  sliderInput("invalidheight", "Plot Height (px)", min = 0, max = 2500, value = 430)
                              ),
                              box(title = "Calibration plot Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                  sliderInput("incalibwidth", "Plot Width (%)", min = 0, max =100, value = 50),
                                  sliderInput("incalibheight", "Plot Height (px)", min = 0, max = 2500, value = 430),
                                  checkboxInput("incalibsubtitle","Show subtitle ?",value = F)
                              ),
                              actionButton(
                                "invalidbt",
                                "Submit",
                                style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                icon = icon("picture-o"))
                            )
                   )
                   
                   
                   
                 )
      )),
      hidden(div(id = "exvalcal_wrapper",
                 fluidRow(
                   column(9,
                          box(
                            title = "External validation and calibration",
                            width = NULL,height=NULL,
                            solidHeader = T,
                            collapsible = T,
                            bsAlert("exnomomess"),
                            status = "success",
                            align="center",
                            tabBox(width = NULL,
                                   tabPanel("External validation", uiOutput('exvalid'),
                                            textOutput(outputId = "exvaliddescrip"),
                                            fluidRow(column(4,align="left",
                                                            hidden(div(id = "exvalid_wrapper",
                                                                       splitLayout(
                                                                         numericInput("exvalidwidthdl","Figure width",value = 10),
                                                                         numericInput("exvalidheightdl","Figure height",value = 10)),
                                                                       downloadButton('saveexvalid', 'Download figure', class = "butt2")
                                                            ))
                                            ))),
                                   tabPanel("External calibration", uiOutput('extcalib'),
                                            fluidRow(column(4,align="left",
                                                            hidden(div(id = "excalib_wrapper",
                                                                       splitLayout(
                                                                         numericInput("excalibwidthdl","Figure width",value = 10),
                                                                         numericInput("excalibheightdl","Figure height",value = 10)),
                                                                       downloadButton('saveexcalib', 'Download figure', class = "butt2")
                                                            ))
                                            )))
                                   
                                   
                            )
                          )
                   ),column(3,
                            box(
                              title = "External validation and calibration",
                              width = NULL,
                              status = "danger",
                              solidHeader = T,
                              collapsible = T,
                              
                              
                              splitLayout(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                          cellWidths = c("0%","50%", "50%"),
                                          selectizeInput(
                                            "extvaltime",
                                            label = "Survival time",
                                            choices = NULL,
                                            multiple = F,
                                            selected = "OS.time"
                                          ),
                                          selectizeInput(
                                            "extvalstatus",
                                            label = "Survival status",
                                            choices = NULL,
                                            multiple = F,
                                            selected = "OS"
                                          )),
                              bsTooltip("extvaltime", "Select survival time column. Example: OS.time, RFS.time, PFS.time, etc.","left"),
                              bsTooltip("extvalstatus", "Select survival status column. Example: OS, RFS, PFS, etc.","left"),
                              
                              splitLayout(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                          cellWidths = c("0%","50%", "50%"),
                                          selectizeInput(
                                            "exvalidypoint",
                                            label = "Prediction years",
                                            choices = NULL,
                                            multiple = T,
                                            options = list(maxItems = 3),
                                            selected="1-year"
                                          ),
                                          selectizeInput(
                                            "extvalvar",
                                            label = "Clinical variable",
                                            choices = NULL,
                                            multiple = T,
                                            selected = "Risk"
                                          )
                              ),
                              bsTooltip("extvalvar", "Select survival variable you want to included","left"),
                              bsTooltip(
                                "exvalidypoint",
                                "Define the time points in years you want to predict based on the nomogram. The longest time point should not exceed the max of survival (relapse) duration",
                                "left"
                              ),
                              splitLayout(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                          cellWidths = c("0%","50%", "50%"),
                                          numericInput("exreps","Boostrap iterations",1000,min = 50,max=3000,step = 100),
                                          numericInput("exratio","Ratio",0.8,min = 0.1,max=1,step = 0.1)),
                              bsTooltip("exreps", "Number of resamplings","left"),
                              bsTooltip("exratio", "Ratio of resamplings","left"),
                              
                              box(title = "C-index Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                  sliderInput("exvalidwidth", "Plot Width (%)", min = 0, max =100, value = 50),
                                  sliderInput("exvalidheight", "Plot Height (px)", min = 0, max = 2500, value = 430)
                              ),
                              box(title = "Calibration plot Size control",width = NULL,solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                  sliderInput("excalibwidth", "Plot Width (%)", min = 0, max =100, value = 50),
                                  sliderInput("excalibheight", "Plot Height (px)", min = 0, max = 2500, value = 430),
                                  checkboxInput("excalibsubtitle","Show subtitle ?",value = F)
                              ),
                              actionButton(
                                "exvalidbt",
                                "Submit",
                                style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                icon = icon("picture-o"))
                            )
                   )
                   
                   
                   
                 )
      )),
      
      fluidRow(
        div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
            actionButton(inputId = 'page_before_nomo',label = '',icon = icon('arrow-left'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
            HTML('<i>Prediction model: Validate model</i>')
        ),
        div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
            HTML('<i>Correlation with clinical features</i>'),
            actionButton(inputId = 'page_after_nomo',label = '',icon = icon('arrow-right'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
        )
      )
    ),
    
    
# Bioann.R

    tabItem(
      tabName = "bioan",
      fluidRow(column(9,
                      box(
                        title = "Functional Enrichment analysis",
                        width = NULL,
                        solidHeader = T,
                        collapsible = T,
                        collapsed = F,
                        align="center",
                        status = "success",
                        bsAlert("FEAmess"),
                        tabBox(width = NULL,
                               tabPanel("Enrichment table", dataTableOutput('enrichtable'),
                                        fluidRow(column(4,align="left",
                                                        hidden(div(id = "enrichtable_wrapper",
                                                                   downloadButton('downloadenrichtable', 'Download table', class = "butt2")
                                                        ))
                                        ))
                               ),
                               tabPanel(title="Enrichment plot",value="DEGen",uiOutput('DEGenplot'),
                                        fluidRow(column(4,align="left",
                                                        hidden(div(id = "DEGenplot_wrapper",
                                                                   splitLayout(
                                                                     numericInput("DEGenplotwidthdl","Figure width",value = 10),
                                                                     numericInput("DEGenplotheightdl","Figure height",value = 10)),
                                                                   downloadButton('downloadDEGenplot', 'Download figure', class = "butt2")
                                                        ))
                                        ))
                               )
                        )
                      )),
               
               column(3,
                      box(title="Enrichment analysis",
                          solidHeader=T,
                          collapsible=T,
                          width=NULL,
                          collapsed = F,
                          status = "danger",
                          bsCollapse(id = "colEnrichment",open="Enrichment method",
                                     bsCollapsePanel("Enrichment method",style = "info",
                                                     
                                                     selectizeInput(
                                                       'biomakers',
                                                       "Biomarkers",
                                                       choices = NULL,
                                                       multiple = F
                                                     ),
                                                     bsTooltip("biomakers", "'Significant DEGs' means significantly differential experssion genes at the cutoff you specified in the 'Differentially expressed genes' module; 'Significant SRGs' means significant survival-related genes at the cutoff you specified in the 'Survival related genes' module; 'Network hub genes' means genes in one or all non-grey modules you selected in the 'WGCNA' module; 'Genes from benchmark experiment' means genes derived from benchmark experiment based on cross-validation or nested cross-validation.","left"),
                                                     selectizeInput(
                                                       'FEAana',
                                                       "Functional analysis",
                                                       choices = c("GO", "KEGG", "MSigDb","Reactome Pathway"),
                                                       selected = "GO"
                                                     ),
                                                     bsTooltip("FEAana", "Define the functional enrichment analysis methods,GO, gene ontology analysis; KEGG, Kyoto Encyclopedia of Genes and Genomes analysis; MsigDb, Molecular Signatures Database analysis; Reactome Pathway, Reactome Pathway analysis","left"),
                                                     
                                                     conditionalPanel(
                                                       condition = "input.FEAana == 'GO'",
                                                       selectizeInput(
                                                         'ont',
                                                         "GO categories",
                                                         choices = c("BP", "CC", "MF","ALL"),
                                                         selected = "BP"
                                                       ),
                                                       bsTooltip("ont", "Define the subcategory of gene ontology.BP, biological process; CC, cellular component;MF, molecular function; 'All' means performing functional GO analysis based on the all three subcategories","left"),
                                                       checkboxInput("readable", "Map gene ID to gene Name ?", value = F, width = NULL)
                                                     ),
                                                     
                                                     selectizeInput(
                                                       'FEAmethod',
                                                       "Functional analysis method",
                                                       choices = c("ORA", "GSEA"),
                                                       selected = "ORA"
                                                     ),
                                                     
                                                     conditionalPanel(
                                                       condition = "input.FEAmethod == 'GSEA'",
                                                       selectizeInput(
                                                         'gseamethod',
                                                         "GSEA algorithm",
                                                         choices = c("DOSE", "fgsea"),
                                                         selected = "fgsea"
                                                       )
                                                     ),
                                                     
                                                     numericInput("minGSSize","Minimal size of genes",10,min = 1,max = 50,step =1),
                                                     bsTooltip("minGSSize", "Minimal size of genes annotated by functional term for testing","left"),
                                                     numericInput("maxGSSize","Maximal size of genes",500,min = 10,max = 100,step =1),
                                                     bsTooltip("maxGSSize", "Maximal size of genes annotated by functional term for testing","left"),
                                                     numericInput("pvalueCutoff","P value Cutoff",0.05,min = 0,max = 1,step =0.05),
                                                     bsTooltip("pvalueCutoff", "Adjusted pvalue cutoff on enrichment tests to report","left"),
                                                     selectizeInput(
                                                       'pAdjustMethod',
                                                       "P Adjust Method",
                                                       choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                                       selected = "BH"
                                                     ),
                                                     numericInput("qvalueCutoff","q value Cutoff",0.2,min = 0,max = 1,step =0.05)
                                                     
                                     ),
                                     bsCollapsePanel("Visualization",style = "info",
                                                     selectizeInput(
                                                       'enplottype',
                                                       "Plot type",
                                                       choices = c("Bar plot", "Dot plot", "Gene-concept network", "Heatmap", "Enrichment Map", "Ridgeline plot", "Geneset enrichment plot1", "Geneset enrichment plot2"),
                                                       selected = "Dot plot"
                                                     ),
                                                     bsTooltip("enplottype", "Please note! Only 'Dot plot', 'Gene-concept network',  'Enrichment Map' are suitable for visualizing the results of ORA analysis, while only 'Dot plot', 'Gene-concept network', 'Heatmap', 'Enrichment Map', 'Ridgeline plot', 'Geneset enrichment plot1', 'Geneset enrichment plot2 are suitable for visualizing the result of GSEA","left"),
                                                     conditionalPanel(
                                                       
                                                       condition = "input.enplottype == 'Bar plot'",
                                                       
                                                       selectizeInput(
                                                         'barcolor',
                                                         "Color",
                                                         choices = c('pvalue', 'p.adjust', 'qvalue'),
                                                         selected = "p.adjust"
                                                       ),
                                                       bsTooltip("barcolor", "The value that the color of the bar map to","left"),
                                                       numericInput("barshowCategory","Show Category",5,min = 1,max = Inf,step =1),
                                                       bsTooltip("barshowCategory", "Number of categories to show","left"),
                                                       numericInput("barfont.size","Font size",12,min = 1,max = Inf,step =1),
                                                       bsTooltip("barfont.size", "Font size of the text","left"),
                                                       
                                                       numericInput("barlabel_format","Label format",30,min = 1,max = Inf,step =1)
                                                       
                                                     ),
                                                     conditionalPanel(
                                                       
                                                       condition = "input.enplottype == 'Dot plot'",
                                                       selectizeInput(
                                                         'dotcolor',
                                                         "Color",
                                                         choices = c('pvalue', 'p.adjust', 'qvalue'),
                                                         selected = "p.adjust"
                                                       ),
                                                       selectizeInput(
                                                         'dotxvar',
                                                         "variable for x-axis",
                                                         choices = c("GeneRatio", 'Count'),
                                                         selected = "GeneRatio"
                                                       ),
                                                       bsTooltip("dotcolor", "The value that the color of the dot map to","left"),
                                                       numericInput("dotshowCategory","Show Category",5,min = 1,max = Inf,step =1),
                                                       bsTooltip("dotshowCategory", "Number of categories to show","left"),
                                                       
                                                       numericInput("dotfont.size","Font size",12,min = 1,max = Inf,step =1),
                                                       bsTooltip("dotfont.size", "Font size of the text","left"),
                                                       
                                                       numericInput("dotlabel_format","Label format",30,min = 1,max = Inf,step =1)
                                                       
                                                       
                                                     ),
                                                     
                                                     conditionalPanel(
                                                       
                                                       condition = "input.enplottype == 'Enrichment Map'",
                                                       
                                                       numericInput("emshowCategory","Show Category",30,min = 1,max = Inf,step =1),
                                                       bsTooltip("emshowCategory", "Number of categories to show","left"),
                                                       selectizeInput(
                                                         'emcolor',
                                                         "Color",
                                                         choices = c("pvalue", "p.adjust" , "qvalue"),
                                                         selected = "p.adjust"
                                                       ),
                                                       bsTooltip("emcolor", "Variable that used to color enriched terms", "left"),
                                                       selectizeInput(
                                                         'emlayout',
                                                         "Layout",
                                                         choices = c('nicely','star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' , 'lgl'),
                                                         selected = "nicely"
                                                       ),
                                                       bsTooltip("emlayout", "Set the layout of the plot", "left"),
                                                       numericInput("emmin_edge","Min edge",0.2,min = 0,max = 1,step =0.1),
                                                       bsTooltip("emmin_edge", "Minimum percentage of overlap genes to display the edge, should between 0 and 1, default value is 0.2", "left"),
                                                       numericInput("emcex_label_category","Cex label category",1,min = 1,max = Inf,step =1),
                                                       bsTooltip("emcex_label_category", "Scale of category node label size", "left"),
                                                       numericInput("emcex_category","Cex category",1,min = 1,max = Inf,step =1),
                                                       bsTooltip("emcex_category", "Number indicating the amount by which plotting category nodes should be scaled relative to the default", "left"),
                                                       numericInput("emcex_line","cex line",1,min = 1,max = Inf,step =1),
                                                       bsTooltip("emcex_line", "scale of line width", "left")
                                                       
                                                     ),
                                                     
                                                     conditionalPanel(
                                                       
                                                       condition = "input.enplottype == 'Gene-concept network'",
                                                       
                                                       numericInput("cneshowCategory","Show Category",5,min = 1,max = Inf,step =1),
                                                       bsTooltip("cneshowCategory", "Number of categories to show","left"),
                                                       checkboxInput("cnecolorEdge", "Coloring Edge ?", value = T, width = NULL),
                                                       bsTooltip("cnecolorEdge", "whether coloring edge by enriched terms","left"),
                                                       checkboxInput("cnecircular", "Circular layout ?", value = F, width = NULL),
                                                       bsTooltip("cnecircular", "whether using circular layout ?","left"),
                                                       selectizeInput(
                                                         'cnenode_label',
                                                         "Node label",
                                                         choices = c('category', 'gene', 'all' , 'none'),
                                                         selected = "all"
                                                       ),
                                                       bsTooltip("cnenode_label", "select which labels to be displayed", "left"),
                                                       numericInput("cnecex_category","Cex category",1,min = 1,max = Inf,step =1),
                                                       bsTooltip("cnecex_category", "Number indicating the amount by which plotting category nodes should be scaled relative to the default.", "left"),
                                                       numericInput("cnecex_gene","Cex gene",1,min = 1,max = Inf,step =1),
                                                       bsTooltip("cnecex_gene", "Number indicating the amount by which plotting gene nodes should be scaled relative to the default", "left"),
                                                       numericInput("cnecex_label_category","Cex label category",1,min = 1,max = Inf,step =1),
                                                       bsTooltip("cnecex_label_category", "Scale of category node label size", "left"),
                                                       numericInput("cnecex_label_gene","Cex label gene",1,min = 1,max = Inf,step =1),
                                                       bsTooltip("cnecex_label_gene", "scale of gene node label size", "left")
                                                       
                                                     )
                                                     ,
                                                     conditionalPanel(
                                                       
                                                       condition = "input.enplottype == 'Heatmap'",
                                                       
                                                       numericInput("heatshowCategory","Show Category",5,min = 1,max = Inf,step =1),
                                                       bsTooltip("heatshowCategory", "Number of categories to show","left")
                                                       
                                                     ),
                                                     conditionalPanel(
                                                       
                                                       condition = "input.enplottype == 'Ridgeline plot'",
                                                       
                                                       numericInput("regshowCategory","Show Category",30,min = 1,max = Inf,step =1),
                                                       bsTooltip("regshowCategory", "Number of categories to show","left"),
                                                       selectizeInput(
                                                         'regfill',
                                                         "Fill",
                                                         choices = c("pvalue", "p.adjust" , "qvalue"),
                                                         selected = "p.adjust"
                                                       ),
                                                       bsTooltip("regfill", "Variable that used to color enriched terms", "left"),
                                                       checkboxInput("regcore_enrichment", "Coloring Edge ?", value = T, width = NULL),
                                                       bsTooltip("regcore_enrichment", "Whether only using core_enriched genes", "left"),
                                                       numericInput("reglabel_format","Label format",30,min = 1,max = Inf,step =1),
                                                       bsTooltip("reglabel_format", "A numeric value sets wrap length, alternatively a custom function to format axis labels", "left")
                                                       
                                                     ),
                                                     conditionalPanel(
                                                       
                                                       condition = "input.enplottype == 'Geneset enrichment plot1'",
                                                       
                                                       numericInput("gsgeneSetID1","GeneSet ID",1,min = 1,max = 500,step =1),
                                                       bsTooltip("gsgeneSetID1", "The numeric geneSet ID", "left"),
                                                       selectizeInput(
                                                         'gsby',
                                                         "By",
                                                         choices = c("runningScore", "preranked", "all"),
                                                         selected ="all"
                                                       ),
                                                       bsTooltip("gsby", "plot the geneset enrichment plot by runningScore, preranked or both", "left"),
                                                       
                                                       colourpicker::colourInput("gscolor", "Color", value = "black", showColour = "background",closeOnClick = TRUE),
                                                       bsTooltip("gscolor", "Color of line segments", "left"),
                                                       colourpicker::colourInput("gscolor.line", "Score Line color", value = "green", showColour = "background",closeOnClick = TRUE),
                                                       bsTooltip("gscolor.line", "Color of running enrichment score line", "left"),
                                                       colourpicker::colourInput("gscolor.vline", "Vertical line color", value = "#FA5860", showColour = "background",closeOnClick = TRUE),
                                                       bsTooltip("gscolor.vline", "color of vertical line which indicating the maximum/minimal running enrichment score", "left")
                                                       
                                                     ),
                                                     conditionalPanel(
                                                       
                                                       condition= "input.enplottype == 'Geneset enrichment plot2'",
                                                       
                                                       numericInput("gsegeneSetID2","GeneSet ID",1,min = 1,max = Inf,step =1),
                                                       bsTooltip("gsegeneSetID2", "The numeric geneSet ID", "left"),
                                                       colourpicker::colourInput("gsecolor", "Color", value = "black", showColour = "background",closeOnClick = TRUE),
                                                       bsTooltip("gsecolor", "Color of line segments", "left"),
                                                       numericInput("gsebase_size","Font size",11,min = 1,max = Inf,step =1),
                                                       bsTooltip("gsebase_size", "Base font size", "left"),
                                                       checkboxInput("gsepvalue_table", "Add pvalue table ?", value = F, width = NULL),
                                                       selectizeInput(
                                                         'gseES_geom',
                                                         "Enrichment score geom",
                                                         choices = c("line", "dot"),
                                                         selected ="line"
                                                       )
                                                       
                                                     ),
                                                     box(title="Size control",solidHeader=F,collapsible=T,width=NULL,collapsed = T,
                                                         sliderInput("enriwidth", "Heatmap Width (%)", min = 0, max = 100, value = 50),
                                                         sliderInput("enriheight", "Heatmap Height (px)", min = 0, max = 1000, value = 430)
                                                         
                                                     )
                                                     
                                     )
                          ),
                          actionButton("enrichbt",
                                       "Submit",
                                       style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                       icon = icon("picture-o"))
                      ))),
      fluidRow(
        div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
            actionButton(inputId = 'page_before_bioan',label = '',icon = icon('arrow-left'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
            HTML('<i>Correlation with stemness score</i>')
        )
        
        
        
        
      )
    ),
    
    
# immune.R
    tabItem(
      tabName = "immune",
      fluidRow(
        column(9,box(
          title="Correlation with Immune cell infiltration",width = NULL,
          solidHeader=T,collapsible=T,status="success",collapsed = F,
          bsAlert("immunemess"),
          tabBox(width = NULL,
                 tabPanel("Correlation plot", align="center",
                          uiOutput('immuneplot'
                          )
                          ,
                          useShinyjs(),
                          fluidRow(column(4,align="left",
                                          hidden(div(id = "immuneplot_wrapper",
                                                     splitLayout(
                                                       numericInput("immpltwidth","Figure width",value = 10),
                                                       numericInput("immpltheight","Figure height",value = 10)),
                                                     downloadButton('saveimune', 'Download figure', class = "butt2")
                                          ))
                          ))
                 ),
                 
                 tabPanel("Immune cell composition matrix", dataTableOutput('immcelltable'),
                          useShinyjs(),
                          fluidRow(column(4,
                                          hidden(div(id = "immunecell_wrapper",
                                                     downloadButton('saveimmcelltable', 'Download table', class = "butt2")
                                          ))
                          ))
                 ),
                 tabPanel("Correlation table", dataTableOutput('immcortable'),
                          useShinyjs(),
                          fluidRow(column(4,
                                          hidden(div(id = "immcor_wrapper",
                                                     downloadButton('saveimmcortable', 'Download table', class = "butt2")
                                          ))
                          ))
                 )
          )
        )
        ),
        column(3,box(
          title = "Correlation with Immune cell infiltration",
          width = NULL,
          status = "danger",
          solidHeader = T,
          collapsible = T,
          collapsed = F,
          selectizeInput(
            "geneset",
            label = "Reference immune geneset",
            choices = c("Bindea","Danaher","Davoli","MCP.Counter","xCell"),
            multiple = F
          ),
          bsTooltip("geneset", "Immune gene signatures from Bindea and colleagues, Danaher and colleagues, Davoli and colleagues, MCP-Counter, and xCell, respectively","left"),
          selectizeInput(
            "immgene",
            label = "Official gene symbol",
            choices = NULL,
            multiple = F
          ),
          bsTooltip("immgene", "Input one gene (biomarker candidate) with official gene symbol that you want to analyze the correlation between the gene and tummor microenvironment.","left"),
          selectizeInput(
            "immmethod",
            label = "Correlation method",
            choices = c("pearson", "spearman"),
            multiple = F,
            selected="spearman"
          ),
          selectizeInput(
            "immnorm",
            label = "Adjust method for P",
            choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),
            multiple = F,
            selected="bonferroni"
          ),
          checkboxInput("corselect", "Select most correlated cells ?", value = T, width = NULL),
          bsTooltip("corselect", "Select the most correlated immune cells based on the cutoff of adjusted p value? ","left"),
          conditionalPanel(
            condition = "input.corselect == true ",
            numericInput("immcutoff","Adjusted P cutoff",0.05,min = 0,max = 1,step =0.005)
          ),
          
          selectizeInput(
            "immpltype",
            label = "Plot type",
            choices = c("bar plot","bubble plot"),
            multiple = F,
            selected="bar plot"
          ),
          box(title = "Size control",width = NULL,
              solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
              sliderInput("immwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
              sliderInput("immheight", "Plot Height (px)", min = 0, max = 1000, value = 430)
          ),
          actionButton("immunebt",
                       "Submit",
                       style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                       icon = icon("picture-o"))
        )
        )
      ),
      fluidRow(
        div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
            actionButton(inputId = 'page_before_immune',label = '',icon = icon('arrow-left'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
            HTML('<i>Gene expression in different groups</i>')
        ),
        div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
            HTML('<i>Correlation with stemness score</i>'),
            actionButton(inputId = 'page_after_immune',label = '',icon = icon('arrow-right'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
        )
      )
    ),
    
# stemness.R
    tabItem(
      tabName = "stemness",
      fluidRow(
        column(9,box(
          title="Correlation with stemness score",width = NULL,
          solidHeader=T,collapsible=T,status="success",collapsed = F,
          bsAlert("stemnessmess"),
          tabBox(width = NULL,
                 tabPanel("Correlation plot", align="center",
                          uiOutput('stemnessplot'
                          )
                          ,
                          useShinyjs(),
                          fluidRow(column(4,align="left",
                                          hidden(div(id = "stemnessplot_wrapper",
                                                     splitLayout(
                                                       numericInput("stempltwidth","Figure width",value = 10),
                                                       numericInput("stempltheight","Figure height",value = 10)),
                                                     downloadButton('savestemness', 'Download figure', class = "butt2")
                                          ))
                          ))
                 ),
                 
                 tabPanel("Stemness table", dataTableOutput('stemtable'),
                          useShinyjs(),
                          fluidRow(column(4,
                                          hidden(div(id = "stemtable_wrapper",
                                                     downloadButton('savestemtable', 'Download table', class = "butt2")
                                          ))
                          ))
                 )
                 
          )
        )
        ),
        column(3,box(
          title = "Correlation with Stemness score",
          width = NULL,
          status = "danger",
          solidHeader = T,
          collapsible = T,
          collapsed = F,
          selectizeInput(
            "stemgene",
            label = "Official gene symbol",
            choices = NULL,
            multiple = F
          ),
          bsTooltip("stemgene", "Input one gene (biomarker candidate) with official gene symbol that you want to analyze the correlation between the gene and tummor microenvironment ","left"),
          selectizeInput(
            "stemmethod",
            label = "Correlation method",
            choices = c("pearson", "spearman"),
            multiple = F,
            selected="spearman"
          ),
          box(title = "Size control",width = NULL,
              solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
              sliderInput("stemwidth", "Plot Width (%)", min = 0, max = 100, value = 50),
              sliderInput("stemheight", "Plot Height (px)", min = 0, max = 1000, value = 430)
          ),
          actionButton("stembt",
                       "Submit",
                       style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                       icon = icon("picture-o"))
        )
        )
      ),
      fluidRow(
        div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
            actionButton(inputId = 'page_before_stemness',label = '',icon = icon('arrow-left'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
            HTML('<i>Correlation with Immune infiltration</i>')
        ),
        div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
            HTML('<i>Biological annotationf</i>'),
            actionButton(inputId = 'page_after_stemness',label = '',icon = icon('arrow-right'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
        )
      )
    )
)
)


ui<-shinyUI(
  bootstrapPage(
    
    dashboardPage(title="CBioExplorer",
                  skin = "green",
                  header,
                  sidebar,
                  body)
    
  )
)

server <- function(input, output, session) {
  
  session$onSessionEnded(stopApp)
  server.path <- ifelse(system.file("app", package = "CBioExplorer") == "",
                        "server",
                        file.path(system.file("app", package = "CBioExplorer"),"server"))
  eset.path <- ifelse(system.file("app", package = "CBioExplorer") == "",
                      "eset",
                      file.path(system.file("app", package = "CBioExplorer"),"eset"))
  
  #datainput.R
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
  
  #ClinicalCorrelation.R
  observe({
    data <- rawdata()$clinical
    updateSelectizeInput(session, 'feature', choices = {
      if (!is.null(data)) {
        if (class(data) ==  class(data.frame())) choices <- as.character(colnames(data))
        if (class(data) !=  class(data.frame())) choices <- as.character(colnames(as.data.frame(data)))
        choices
      }
    }, server = TRUE)
  }
  )
  
  observe({
    data <- rawdata()$expres
    updateSelectizeInput(session,
                         'table1bygene',
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
  
  
  percent<-function(x,per){
    n<-length(x)
    half <- (n + 1L) %/% (1/per)
    if(n %% (1/per) == 0) mean(sort(x, partial = half + 0L:1L)[half + 0L:1L])
    else   sort(x, partial = half)[half]
  }
  
  
  drawtable1<-eventReactive(input$table1bt,{
    
    input$table1bt
    data<-isolate({rawdata()})
    data<-lapply(data,as.data.frame)
    index<-intersect(colnames(data$expres),row.names(data$clinical))
    expres<-data$expres[,index]
    clinical<-data$clinical[index,]
    clinical[clinical == ""]<- NA
    gene<-isolate({input$table1bygene})
    clinical$gene<-as.numeric(expres[gene,])
    grouppercent<-isolate({input$grouppercent})
    if(isolate({input$tbgroupby=="Percentage"})){
      clinical$Group<-ifelse(clinical$gene<=quantile(clinical$gene,grouppercent),paste(gene,"low expression group",sep=" "),paste(gene,"high expression group", sep=" "))
    }else{
      clinical$Group<-ifelse(clinical$gene<=isolate({input$tabgpvalue}),paste(gene,"low expression group",sep=" "),paste(gene,"high expression group", sep=" "))
    }
    if("Age" %in% names(clinical)){
      units(clinical$Age) <- "years"
    }
    if("OS.time" %in% names(clinical)){
      if(max(na.omit(clinical$OS.time))>600){
        units(clinical$OS.time)="Day"
      } else{
        units(clinical$OS.time)="Month"
        
      }
    }
    if("RFS.time" %in% names(clinical)){
      if(max(na.omit(clinical$RFS.time))>600){
        units(clinical$RFS.time)="Day"
      } else{
        units(clinical$RFS.time)="Month"
      }
    }
    
    if("PFS.time" %in% names(clinical)){
      if(max(na.omit(clinical$PFS.time))>600){
        units(clinical$PFS.time)="Day"
      } else{
        units(clinical$PFS.time)="Month"
      }
    }
    
    if("DFS.time" %in% names(clinical)){
      if(max(na.omit(clinical$DFS.time))>600){
        units(clinical$DFS.time)="Day"
      } else{
        units(clinical$DFS.time)="Month"
      }
    }
    if("DMFS.time" %in% names(clinical)){
      if(max(na.omit(clinical$DMFS.time))>600){
        units(clinical$DMFS.time)="Day"
      } else{
        units(clinical$DMFS.time)="Month"
      }
    }
    if("DRFS.time" %in% names(clinical)){
      if(max(na.omit(clinical$DRFS.time))>600){
        units(clinical$DRFS.time)="Day"
      } else{
        units(clinical$DRFS.time)="Month"
      }
    }
    
    clinical$Group<-as.factor(clinical$Group)
    clinical$Group<-factor(clinical$Group,levels=c(levels(clinical$Group),"P value"))
    
    
    rndr <- function(x, name, ...) {
      if (length(x) == 0) {
        y <- clinical[[name]]
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
          p <- t.test(y ~ clinical$Group)$p.value
        } else {
          p <- chisq.test(table(y, droplevels(clinical$Group)))$p.value
        }
        s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
        s
      } else {
        render.default(x=x, name=name, ...)
      }
    }
    rndr.strat <- function(label, n, ...) {
      ifelse(n==0, label, render.strat.default(label, n, ...))
    }
    
    feature<-isolate({input$feature})
    feature<-paste(feature,collapse ="+")
    fomu<-as.formula(paste("~",feature,"|Group",sep=''))
    table1(fomu,data=clinical,droplevels=F, render=rndr, render.strat=rndr.strat,overall = F)
    
  }
  )
  
  observe({
    if(is.null(input$table1bygene) || input$table1bygene == "" ||is.null(input$feature) ||input$feature == "" ){
      disable("table1bt")
    }
    else{
      enable("table1bt")
    }
  })
  
  
  
  observeEvent(input$table1bt , {
    
    output$table1 <- renderText({
      
      drawtable1()
    })
  })
  observeEvent(input$table1bt , {
    shinyjs::show("clincor_wrapper")
  })
  
  output$downloadtable1 <- downloadHandler(
    
    filename = function() {
      paste0("Table1-clinical-correlation-",Sys.Date(),'.html')
    },
    content = function(file) {
      saveWidget(print(drawtable1()), file,selfcontained = T)
    }
  )
  
  
  observeEvent(input$page_before_clinical, {
    newtab <- switch(input$tabs, "clinical" = "nomo","nomo" = "clinical")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  observeEvent(input$page_after_clinical, {
    newtab <- switch(input$tabs, "KM" = "clinical","clinical" = "KM")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  #KM.R
  observe({
    data <- rawdata()$expres
    updateSelectizeInput(session,
                         'KMgene',
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
  
  
  observe({
    data <- rawdata()$clinical
    updateSelectizeInput(session,
                         'survivaltime',
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
                         'survivalstatus',
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
  
  plotKM <- eventReactive(input$KMplotbt,{
    input$KMplotbt
    data <- isolate({
      rawdata()
    })
    
    data <- lapply(data, as.data.frame)
    time <- isolate({ input$survivaltime})
    status <- isolate({input$survivalstatus})
    data$clinical <-data$clinical[complete.cases(data$clinical[, time]) & data$clinical[, time] > 0,]
    index <-intersect(colnames(data$expres), row.names(data$clinical))
    clinical <- data$clinical[index,]
    clinical[clinical == ""] <- NA
    expres <- data$expres[, index]
    KMgene <- isolate({
      input$KMgene
    })
    
    clinical$KMgene <- as.numeric(expres[KMgene,])
    clinical<-clinical[complete.cases(clinical$KMgene),]
    
    
    OS.time<-clinical[,time]
    OS<-clinical[,status]
    if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
      createAlert(
        session,
        "kmmess",
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
        "kmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival status column is not numerically in 1, and 0, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    } else if(isolate({input$kmgroupby})=="Value"&& isolate({input$kmgpvalue})>=max(clinical$KMgene) ){
      
      createAlert(
        session,
        "kmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste("The cutoff value you set exceeds the range: (", paste(range(clinical$KMgene),collapse = ", "), ") of the expression level of", paste(KMgene), "you specified, please check!"),
        append = T,
        dismiss=T
      )
      return(NULL)
      
    }else if(isolate({input$kmgroupby})=="Value" && isolate({input$kmgpvalue})<= min(clinical$KMgene)){
      createAlert(
        session,
        "kmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste("The cutoff value you set exceeds the range: (", paste(range(clinical$KMgene),collapse = ", "), ") of the expression level of", paste(KMgene), "you specified, please check!"),
        append = T,
        dismiss=T
      )
      return(NULL)
    }else {
      genecut <- isolate({input$genecut})
      if(isolate({input$kmgroupby})=="Percentage"){
        clinical$KMgroup <-ifelse(
          clinical$KMgene <= quantile(clinical$KMgene, genecut),
          paste(KMgene, "low expression group", sep = " "),
          paste(KMgene, "high expression group", sep = " ")
        )
      }else{
        
        clinical$KMgroup<-ifelse(clinical$KMgene<=isolate({input$kmgpvalue}),paste(KMgene,"low expression group",sep=" "),paste(KMgene,"high expression group", sep=" "))
      }
      
      survxlab <- isolate({input$survxlab})
      survP <- isolate({input$survP})
      survRT <- isolate({input$survRT})
      survCI <- isolate({input$survCI})
      KMcurvename <- isolate({input$KMcurvename})
      clinical$time <- clinical[, time]
      clinical$status <- clinical[, status]
      fit <- survfit(Surv(time, status) ~ KMgroup, data = clinical)
      level <- levels(factor(clinical$KMgroup))
      c <-
        ifelse(
          level[1] == paste(KMgene, "low expression group", sep = " "),
          
          paste(isolate({input$kmcolor2})),
          
          paste(isolate({input$kmcolor1}))
        )
      
      d <-
        ifelse(
          level[2] == paste(KMgene, "high expression group", sep = " "),
          paste(isolate({input$kmcolor1})),
          paste(isolate({input$kmcolor2}))
        )
      
      e <-
        ifelse(
          level[1] == paste(KMgene, "low expression group", sep = " "),
          paste(KMgene, "low expression group", sep = " "),
          paste(KMgene, "high expression group", sep = " ")
        )
      
      f <-
        ifelse(
          level[2] == paste(KMgene, "low expression group", sep = " "),
          paste(KMgene, "low expression group", sep = " "),
          paste(KMgene, "high expression group", sep = " ")
        )
      
      res <- ggsurvplot(
        fit,
        data = clinical,
        risk.table = survRT,
        risk.table.height = 0.3,
        risk.table.y.text = FALSE,
        risk.table.title = "",
        main = "Survival curve",
        palette = c(c, d),
        pval = survP,
        pval.method = T,
        conf.int = survCI,
        risk.table.y.text.col = T,
        legend = c(0.8, 0.90),
        legend.title = "",
        xlab = survxlab,
        legend.labs = c(e, f)
      )
      res
    }
  })
  
  observe({
    if(is.null(input$KMgene) || input$KMgene== "" ||is.null(input$survivaltime) ||input$survivaltime == "" ||is.null(input$survivalstatus) ||input$survivalstatus == ""){
      disable("KMplotbt")
    }
    else{
      enable("KMplotbt")
    }
  })
  
  
  
  observeEvent(input$KMplotbt, {
    output$KMploting <-  renderPlot(
      
      {print(plotKM())}
    )
  })
  
  observeEvent(input$KMplotbt, {
    
    output$KMplot <- renderUI({
      plotOutput("KMploting",
                 
                 width =   paste0(isolate({input$survwidth}),"%"),
                 height = isolate({input$survheight}))
    })})
  observeEvent(input$KMplotbt, {
    shinyjs::show("km_wrapper")
  })
  
  
  
  
  output$downloadKM <- downloadHandler(
    
    filename = function() {
      paste0("K-M-curve-",Sys.Date(),'.pdf')
    },
    
    content = function(file) {
      
      pdf(file,width = isolate({input$kmwidth}), height =isolate({input$kmheight})) # open the pdf device
      print(plotKM())
      dev.off()
    }
  )
  
  observeEvent(input$page_before_KM, {
    newtab <- switch(input$tabs, "clinical" = "KM","KM" = "clinical")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  observeEvent(input$page_after_KM, {
    newtab <- switch(input$tabs, "KM" = "CoxPH","CoxPH" = "KM")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  #CoxPH.R
  observe({
    data <- rawdata()$clinical
    updateSelectizeInput(session,
                         'CoxPHtime',
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
                         'CoxPHstatus',
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
  
  observe({
    data <- rawdata()$clinical
    
    updateSelectizeInput(session, 'coxclinvar', choices = {
      if (!is.null(data)) {
        if (class(data) ==  class(data.frame())) choices <- as.character(colnames(data))
        if (class(data) !=  class(data.frame())) choices <- as.character(colnames(as.data.frame(data)))
        choices
      }
    }, server = TRUE)
  }
  )
  
  observe({
    data <- rawdata()$expres
    updateSelectizeInput(session,
                         'coxgene',
                         choices = {
                           if (!is.null(data)) {
                             if (class(data) ==  class(data.frame()))
                               
                               choices <- as.character(make.names(row.names(data)))
                             if (class(data) !=  class(data.frame()))
                               
                               choices <- as.character(make.names(row.names(as.data.frame(data))))
                             choices
                           }
                         },
                         server = TRUE,
                         selected = "TP53")
  })
  
  
  coxPHmodel <- eventReactive(input$CoxPHbt,{
    input$CoxPHbt
    data <- isolate({
      rawdata()
    })
    
    data <- lapply(data, as.data.frame)
    time <- isolate({input$CoxPHtime})
    status <- isolate({input$CoxPHstatus})
    data$clinical <-data$clinical[complete.cases(data$clinical[, time]) & data$clinical[, time] > 0, ]
    index <- intersect(colnames(data$expres), row.names(data$clinical))
    clinical <- data$clinical[index, ]
    clinical[clinical == ""] <- NA
    expres <- data$expres[, index]
    row.names(expres)<-gsub("-","_",row.names(expres))
    coxclinvar <- isolate({input$coxclinvar})
    
    coxgene <- isolate({input$coxgene})
    
    maxtick<-isolate(input$maxtick)
    
    OS.time<-clinical[,time]
    OS<-clinical[,status]
    
    clinfeat<-subset(clinical,select=coxclinvar)
    name<-names(clinfeat)
    fc<-function(x){
      x<-factor(x)
      x<-length(levels(x))
      return(x)
    }
    clinfeature<-lapply(clinfeat,fc)
    
    if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
      createAlert(
        session,
        "CoxPHmess",
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
        "CoxPHmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival status column is not numerically in 1, and 0, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    } else if(1 %in% clinfeature){
      createAlert(
        session,
        "CoxPHmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste("Variable",paste(name[which(clinfeat==1)]),"you selected has only one level (only one value), which is not accepted by CoxPH model. Please check!"),
        append = T,
        dismiss = T
      )
      return(NULL)
    }else{
      clinical$time<-clinical[,time]
      
      clinical$status<-clinical[,status]
      if(!is.null(coxgene)){
        
        feature<-c(coxclinvar,coxgene)
        coxgene<-subset(as.data.frame(t(expres)),select=coxgene)
        clinical<-merge(clinical,coxgene,by=0)
        
      }else{
        feature<-coxclinvar
        clinical<-clinical
      }
      
      Surv<-Surv(clinical$time, clinical$status)
      feature1<-paste(feature,collapse ="+")
      fomu<-as.formula(paste("Surv","~",feature1,sep=''))
      fit<-coxph(fomu,data=clinical)
      sumfit<-summary(fit)
      
      a<-cbind(sumfit$conf.int[,1:4],sumfit$coefficients[,5])
      colnames(a)<-c("HR","exp(-coef)","LCI","UCI","P value")
      multicox <- as.data.frame(a)
      
      unicox <- c()
      for (i in 1:length(feature)) {
        fomu <- as.formula(paste("Surv~", feature[i], sep = ""))
        cox <- coxph(fomu, data = clinical)
        cox <- summary(cox)
        conf.int <- cox$conf.int
        coef1 <- cox$coef
        a <- cbind(conf.int, coef1[, 5])
        colnames(a) <- c("HR", "exp(-coef)", "LCI", "UCI", "P value")
        unicox <- rbind(unicox, a)
      }
      
      unicox<-as.data.frame(unicox)
      index<-intersect(rownames(multicox),row.names(unicox))
      
      varname<-isolate({input$coxPHvarname})
      var1 <- strsplit(varname, "|", fixed = T)[[1]]
      
      
      if(length(var1)==length(index)){
        unicox<-unicox[index,]
        multicox<-multicox[index,]
        row.names(unicox)<-var1
        row.names(multicox)<-var1
      }else if(varname==""){
        unicox<-unicox[index,]
        multicox<-multicox[index,]
      }else{
        createAlert(
          session,
          "CoxPHmess",
          "exampleAlert",
          title = "Please note!",
          style =  "danger",
          content = "The length of variable names does not equal to the clinical variables included in the CoxPH model",
          append = T,
          dismiss = T
        )
        return(NULL)
      }
      
      plotcox<-coxforestp(unicox=unicox,multicox=multicox,legend.pos=isolate({input$coxfiglg}) ,xlim=isolate({input$maxtick}),varname=NULL)
      sumcox<-sumfit
      unicox<-round(unicox,3)
      unicox<-as.data.frame(unicox)
      unicox$Pvalue<- format.pval(unicox$`P value`, digits = 3, eps = 0.001)
      
      multicox<-round(multicox,3)
      
      multicox$Pvalue<- format.pval(multicox$`P value`, digits = 3, eps = 0.001)
      unicox$summaryx<-paste(unicox$HR,paste("(",unicox$LCI,"-",paste(paste(unicox$UCI,",",sep=""),unicox$Pvalue,sep=" "),")",sep=""),sep=" ")
      multicox$summaryy<-paste(multicox$HR,paste("(",multicox$LCI,"-",paste(paste(multicox$UCI,",",sep=""),multicox$Pvalue,sep=" "),")",sep=""),sep=" ")
      
      tab<-merge(unicox,multicox,by=0)
      rownames(tab)<-tab$Row.names
      
      tablecox<-subset(tab,select=c(summaryx,summaryy))
      names(tablecox)<-c("Univariate HR(95% CI,Pvalue)","Multivariable HR(95% CI,Pvalue)")
      res<-list(tablecox,plotcox,sumcox)
      names(res)<-c("tablecox","plotcox","sumcox")
      
      res
    }
  })
  
  
  
  observe({
    if(is.null(input$coxclinvar) || input$coxclinvar == ""||is.null(input$CoxPHtime) || input$CoxPHtime == "" || is.null(input$CoxPHstatus) || input$CoxPHstatus == ""){
      disable("CoxPHbt")
    }
    else{
      enable("CoxPHbt")
    }
  })
  
  observeEvent(input$CoxPHbt, {
    output$Coxtable <-  DT::renderDT({
      coxPHmodel()$tablecox
    })
  })
  
  
  observeEvent(input$CoxPHbt, {
    shinyjs::show("mybox_wrapper")
  })
  
  
  
  observeEvent(input$CoxPHbt, {
    output$Coxforestploting  <- renderPlot({
      coxPHmodel()$plotcox
    })
  })
  
  observeEvent(input$CoxPHbt, {
    output$Coxforestplot<- renderUI({
      plotOutput("Coxforestploting",
                 width = paste0(isolate({input$Coxwidth}), "%"),
                 height = isolate({input$Coxheight}))
    })})
  
  
  observeEvent(input$CoxPHbt, {
    shinyjs::show("mybox_wrapper.table")
  })
  
  
  observeEvent(input$CoxPHbt, {
    output$Coxsummary <-   renderPrint({
      coxPHmodel()$sumcox
    })
  })
  
  output$saveforest <- downloadHandler(
    
    filename = function() {
      paste0("CoxPH-forestplot-",Sys.Date(),'.pdf')
    },
    content = function(file) {
      
      pdf(file,width = isolate({input$FPwidth}),height = isolate({input$FPheight})) # open the pdf device
      print(coxPHmodel()$plotcox)
      dev.off()
      
    }
  )
  
  output$saveforesttable <- downloadHandler(
    filename = function() {
      paste0("CoxPH-table-",Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(coxPHmodel()$tablecox, file)
    }
  )
  
  observeEvent(input$page_before_CoxPH, {
    newtab <- switch(input$tabs, "CoxPH" = "KM","KM" = "CoxPH")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  observeEvent(input$page_after_CoxPH, {
    newtab <- switch(input$tabs, "SurvROC" = "CoxPH","CoxPH" = "SurvROC")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  
  #SurvROC.R
  
  observeEvent(rawdata(), {
    data <- rawdata()$clinical
    updateSelectInput(session, "SurvROCtime", choices = {
      if (!is.null(data)) {
        if (class(data) ==  class(data.frame()))
          choices <- as.character(colnames(data))
        if (class(data) !=  class(data.frame()))
          choices <-
            as.character(colnames(as.data.frame(data)))
        choices
      }
    },
    selected = "OS.time"
    )
  })
  
  observeEvent(rawdata(), {
    data <- rawdata()$clinical
    updateSelectInput(session, "SurvROCstatus", choices = {
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
    selected = "OS"
    )
  })
  
  observe({
    data <- rawdata()$expres
    updateSelectizeInput(session,
                         'SurvROCgene',
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
  
  yearlabelx <- reactive({
    req(input$SurvROCtime)
    SurvROCtime<-isolate({input$SurvROCtime})
    data <- rawdata()$clinical
    if(SurvROCtime %in% names(data)){
      time<-data[,SurvROCtime]
      
      if(is.numeric(time)==F){
        createAlert(
          session,
          "exnomomess",
          "exampleAlert",
          title = "Please note!",
          style =  "danger",
          content = "The survival time you selected is not numeric data, please check!",
          append = T,
          dismiss=T
        )
        return(NULL)
      } else {
        
        if(max(na.omit(data[,SurvROCtime]))>600){
          by<- 365
        } else{
          by<-12
        }
        
        years<-seq(from=0,to=quantile(na.omit(data[,SurvROCtime]),0.95),by=by)
        years<-years[-1]
        yearlabel<-c()
        for(i in 1:length(years)){
          yearlabel[i]<- paste(i,"year",sep="-")
        }
        yearlabel
      }
    }
  })
  
  observeEvent(yearlabelx(), {
    choices <- yearlabelx()
    updateSelectInput(session, "predictyear", choices = choices)
  })
  
  
  SurvROCplot<-reactive({
    input$SurvROCbt
    data <- isolate({
      rawdata()
    })
    SurvROCgene<-isolate({input$SurvROCgene})
    SurvROCtime<-isolate({input$SurvROCtime})
    
    SurvROCstatus<-isolate({input$SurvROCstatus})
    method<-isolate({input$SurvROCmethod})
    clinical<- data$clinical
    
    OS.time<-clinical[,SurvROCtime]
    OS<-clinical[,SurvROCstatus]
    
    if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
      createAlert(
        session,
        "exnomomess",
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
        "exnomomess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival status column is not numerically in 1, and 0, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    }else{
      
      expres<- data$expres
      clinical<-clinical[complete.cases(clinical[,SurvROCtime])& clinical[,SurvROCtime]>0,]
      index<-intersect(row.names(clinical),colnames(expres))
      clinical<-clinical[index,]
      expres<-expres[,index]
      clinical$SurvROCgene<-as.numeric(expres[SurvROCgene,])
      withProgress(
        message = "Constructing the prediction model",
        detail = "It may take a while, please be patient",
        value = 5,{
          survROC(
            data =clinical,
            bmtime = SurvROCtime,
            bmstatus = SurvROCstatus,
            marker = "SurvROCgene",
            method = method,
            predyear = isolate({input$predictyear}),
            cutpoint=isolate({input$cutoff})
          )
        })
      
      # if(max(clinical[,SurvROCtime])>600){
      #   by<- 365
      # } else{
      #   by<-12
      # }
      # years<-seq(from=0,to=quantile(na.omit(clinical[,SurvROCtime]),0.95),by=by)
      # years<-years[-1]
      # yearlabel<-c()
      # for(i in 1:length(years)){
      #   yearlabel[i]<- paste(i,"year",sep="-")
      # }
      # 
      # sumROC<-list()
      # for (i in 1:length(years)){
      #   sumROC[[i]] <- survivalROC(Stime = clinical[,SurvROCtime],status = clinical[,SurvROCstatus],marker = clinical$SurvROCgene,
      #                              predict.time =years[i],method = method,span = 0.25*nrow(clinical)^(-0.20))
      # }
      # sumAUC<-list()
      # for (i in 1:length(sumROC)){
      #   sumAUC[[i]]<-sumROC[[i]]$AUC
      # }
      # 
      # ROCdata<-c()
      # for(i in 1:length(sumROC)){
      #   predict.time<-sumROC[[i]]$predict.time
      #   TP<-sumROC[[1]]$TP
      #   FP<-sumROC[[i]]$FP
      #   auc<-sumROC[[i]]$AUC
      #   tmp<-c(predict.time,TP,FP,auc)
      #   ROCdata<-rbind(ROCdata,tmp)
      # }
      # 
      # survivalROC_helper <- function(t) {
      #   survivalROC(Stime = clinical[,SurvROCtime],status = clinical[,SurvROCstatus],marker = clinical$SurvROCgene,
      #               predict.time =t,method = method,span = 0.25*nrow(clinical)^(-0.20))
      # }
      # 
      # timeponit<-isolate({input$predictyear})
      # yearid<-as.numeric(substr(timeponit,1,1))
      # time<-years[yearid]
      # 
      # survivalROC_data <- data_frame(t = time) %>%
      #   mutate(survivalROC = map(t, survivalROC_helper),
      #          auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
      #          df_survivalROC = map(survivalROC, function(obj) {
      #            as_data_frame(obj[c("cut.values","TP","FP")])
      #          })) %>%
      #   dplyr::select(-survivalROC) %>%
      #   unnest() %>%
      #   arrange(t, FP, TP)
      # 
      # survivalROC_data1<-mutate(survivalROC_data,auc =sprintf("%.3f",auc))
      # survivalROC_data1$years<-survivalROC_data1$t/by
      # survivalROC_data1<-unite(survivalROC_data1,year, years,auc,sep = " year AUC: " )
      # AUC =factor(survivalROC_data1$year)
      # sumAUC1<-list()
      # for(i in 1:length(yearid)){
      #   sumAUC1[[i]]<-sumAUC[[yearid[i]]]
      # }
      # 
      # sumROC1<-list()
      # for(i in 1:length(yearid)){
      #   sumROC1[[i]]<-sumROC[[yearid[i]]]
      # }
      # 
      # 
      # ROC.1<-sumROC1[[which.max(sumAUC1)]]
      # 
      # dot <- data.frame(TP = ROC.1$TP[with(ROC.1, which.min((1-TP)^2+ FP^2))],
      #                   FP = ROC.1$FP[with(ROC.1, which.min((1-TP)^2+ FP^2))])
      # dot <- rbind(c(1,0),dot)
      # 
      # if(isolate({input$cutoff})==T){
      #   cutoff.imm <- ROC.1$cut.values[with(ROC.1, which.min((1-TP)^2+ FP^2))]
      #   ROC.plot<-ggplot(survivalROC_data1,mapping = aes(x = FP, y = TP)) +
      #     geom_path(aes(color= AUC))+
      #     geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
      #     theme_bw() +scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+xlab("False positive rate")+
      #     ylab("True positive rate") +
      #     theme(legend.position = c(0.7,0.2))+
      #     geom_path(mapping = aes(x = FP,y = TP),data = dot)+
      #     annotate("text",x = dot$FP[2] + 0.05,y = dot$TP[2],label = paste0("Cutoff: ",round(cutoff.imm,3)))
      # } else{
      #   ROC.plot<-ggplot(survivalROC_data1,mapping = aes(x = FP, y = TP)) +
      #     geom_path(aes(color= AUC))+
      #     geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
      #     theme_bw() +scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+xlab("False positive rate")+
      #     ylab("True positive rate") +
      #     theme(legend.position = c(0.7,0.2))
      # }
      # ROC.plot
    }}
    
  )
  
  observe({
    if (is.null(input$SurvROCgene) ||
        input$SurvROCgene == "" ||
        is.null(input$SurvROCtime) ||
        input$SurvROCtime == "" ||
        is.null(input$SurvROCstatus) ||
        input$SurvROCstatus == "" ||
        is.null(input$predictyear) || input$predictyear == "") {
      disable("SurvROCbt")
    }
    else{
      enable("SurvROCbt")
    }
  })
  
  observeEvent(input$SurvROCbt, {
    output$SurvROCplotting <- renderPlot({
      SurvROCplot()
    })
  })
  
  observeEvent(input$SurvROCbt, {
    output$SurvROCplot <- renderUI({
      plotOutput("SurvROCplotting",
                 width = paste0(isolate({input$SurvROCwidth}), "%"),
                 height = isolate({input$SurvROCheight}))
    })})
  
  observeEvent(input$SurvROCbt, {
    shinyjs::show("survROC_wrapper")
  })
  
  output$downloadsurvROC <- downloadHandler(
    
    filename = function() {
      paste0("SurvivalROC-",Sys.Date(),'.pdf')
    },
    content = function(file) {
      pdf(file,width = isolate({input$survROCwidth.dl}),height = isolate({input$survROCheight.dl})) # open the pdf device
      print(SurvROCplot())
      dev.off()
    }
  )
  
  observeEvent(input$page_before_SurvROC, {
    newtab <- switch(input$tabs, "CoxPH" = "SurvROC","SurvROC" = "CoxPH")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  observeEvent(input$page_after_SurvROC, {
    newtab <- switch(input$tabs, "SurvROC" = "mcorgene","mcorgene" = "SurvROC")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  #genecor.R
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
  
  #mcorgene.R
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
                   
                   res<-corAndPvalue(geneexp, expres,method = mcormethod)
                   res<-merge(as.data.frame(t(res$cor)),as.data.frame(t(res$p)),by=0)
                   names(res)<-c("Gene","R","Pvalue")
                   
                   corsigcut<-isolate({input$corsigcut})
                   
                   if(isolate({input$padjust})==TRUE){
                     mcorrectmethod<-isolate({input$mcorrectmethod})
                     res$padjust<-p.adjust(res$Pvalue,method = mcorrectmethod)
                     names(res)<-c("Gene","R","Pvalue","Padjusted")
                     tab<-res[res$Padjusted<corsigcut,]
                     
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
                     tab<-res[res$Pvalue<corsigcut,]
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
  })
  
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
  
  
  #genediff.R
  
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
  
  genedifffun<-eventReactive(input$genediffbt,{
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
  
  #unicox.R
  observeEvent(rawdata(),{
    data <- rawdata()$clinical
    updateSelectizeInput(session, "msurvtime", choices = {
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
    selected = "OS.time"
    )
  }
  )
  
  observeEvent(rawdata(),{
    data <- rawdata()$clinical
    updateSelectizeInput(session, "msurvstatus", choices = {
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
    selected = "OS"
    )
  }
  )
  
  msurvfun<-eventReactive(input$msurvbt,{
    input$msurvbt
    withProgress(message = "Performing Univariate Cox proportional hazards analysis",
                 detail = "This may take a while...",
                 value =5,
                 {
                   data <- isolate({
                     rawdata()})
                   data <- lapply(data, as.data.frame)
                   expres<-data$expres
                   clinical<-data$clinical
                   OS.time<-isolate({input$msurvtime})
                   OS<-isolate({input$msurvstatus})
                   
                   # method<-isolate({input$msurvcor})
                   clinical<-clinical[complete.cases(clinical[,OS.time])&clinical[,OS.time]>0,]
                   index<-intersect(row.names(clinical),names(expres))
                   clinical<-clinical[index,]
                   expres<-expres[,index]
                   OS.time<-clinical[,OS.time]
                   OS<-clinical[,OS]
                   if(mode(OS.time)!="numeric" || mode(OS)!="numeric" ){
                     createAlert(
                       session,
                       "msurvmess",
                       "exampleAlert",
                       title = "Please note!",
                       style =  "danger",
                       content = "Eihter the survival time or survival status column is not numeric data, please check!",
                       append = T,
                       dismiss=T
                     )
                     return(NULL)
                     
                   }else if(any(levels(factor(OS))== c("0", "1"))!=T){
                     createAlert(
                       session,
                       "msurvmess",
                       "exampleAlert",
                       title = "Please note!",
                       style =  "danger",
                       content = "The survival status column is not numerically in 1, and 0, please check!",
                       append = T,
                       dismiss=T
                     )
                     return(NULL)
                   }else{
                     res<-unicox(OS.time,OS,expres)
                     res<-as.data.frame(res)
                     res<-na.omit(res)
                     # res$P.adjusted<-p.adjust(res$PValue, method = method)
                     list(res=res,clinical=clinical,expres=expres)
                   }
                 })
  })
  
  observe({
    if(is.null(input$msurvtime) || input$msurvtime == ""||is.null(input$msurvstatus) || input$msurvstatus == "" ){
      disable("msurvbt")
    }
    else{
      enable("msurvbt")
    }
  })
  
  observeEvent(input$msurvbt, {
    output$unicoxtable <-  DT::renderDT(
      {
        DT::datatable(msurvfun()$res
                      ,options=list(scrollX=TRUE))
      }
    )
  })
  
  observeEvent(input$msurvbt, {
    disable("msurvtime")
    disable("msurvstatus")
    disable("msurvcor")
    disable("msurvbt")
    disable("msurvopbt")
  })
  observeEvent(is.null(msurvfun()), {
    # shinyjs::show("unicoxanalysis_wrapper")
    enable("msurvtime")
    enable("msurvstatus")
    enable("msurvcor")
    enable("msurvbt")
    enable("msurvopbt")
    # shinyjs::show("unicoxtable_wrapper")
  })
  observeEvent(msurvfun(), {
    shinyjs::show("unicoxanalysis_wrapper")
    enable("msurvtime")
    enable("msurvstatus")
    enable("msurvcor")
    enable("msurvbt")
    enable("msurvopbt")
    shinyjs::show("unicoxtable_wrapper")
  })
  
  
  
  output$downloadunicoxtable <- downloadHandler(
    filename = function() {
      paste0("Univariate-CoxPH-analysis-",Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(msurvfun()$res, file)
    }
  )
  #########
  msurvsig<-reactive({
    input$msurvopbt
    res<-isolate({msurvfun()$res})
    method<-isolate({input$msurvcor})
    res$P.adjusted<-p.adjust(res$PValue, method = method)
    if(isolate({input$msurvsel})=="P value"){
      sigres<-res[res$PValue<isolate({input$msurvp}),]
    } else{
      sigres<-res[res$P.adjusted<isolate({input$msurvap}),]
    }
    
    #sigres$P.adjusted<-NULL
    
    if(nrow(sigres)==0){
      createAlert(
        session,
        "msurvmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste(
          "No gene meet the significance cutoff you specified"),
        append = T,
        dismiss=T
      )
      list(res = res,
           clinical = isolate({
             msurvfun()$clinical
           }),
           expres = isolate({
             msurvfun()$expres
           }))
    } else{
      
      list(
        res = res,
        clinical = isolate({
          msurvfun()$clinical
        }),
        expres = isolate({
          msurvfun()$expres
        }),
        Output = sigres
      )
    }
  })
  
  
  observeEvent(input$msurvopbt, {
    output$sigcoxtable <-  DT::renderDT(
      {
        DT::datatable(msurvsig()$Output,
                      options=list(scrollX=TRUE))
      }
      
    )
  })
  
  output$downloadsigcoxtable <- downloadHandler(
    filename = function() {
      paste0("significant-Univariate-CoxPH-analysis-",Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(msurvsig()$Output, file)
    }
  )
  
  forestp<-reactive({
    input$msurvopbt
    withProgress(message = "Drawing forestplot",
                 detail = "This may take a while...",
                 value =5,
                 {
                   res<-isolate({msurvsig()$Output})
                   unicoxforestp(res,xlim = isolate(input$msurvxtick))
                 })
  })
  
  observeEvent(input$msurvopbt, {
    disable("msurvopbt")
    disable("msurvbt")
  })
  
  observeEvent(msurvsig(), {
    shinyjs::show("unicoxForestplot_wrapper")
    enable("msurvopbt")
    enable("msurvbt")
    shinyjs::show("sigcoxtable_wrapper")
  })
  
  
  
  
  observeEvent(input$msurvopbt, {
    output$sigforing  <- renderPlot({
      closeAlert(session, "msurvmess")
      forestp()
    })
  })
  
  observeEvent(input$msurvopbt, {
    output$sigfor <- renderUI({
      plotOutput("sigforing",
                 width = paste0(isolate({input$msurvwidth}), "%"),
                 height = isolate({input$msurvheight}))
    })})
  
  
  output$downloadunicoxForestplot <- downloadHandler(
    
    filename = function() {
      paste0("significant-univariate-forestplot-",Sys.Date(),'.pdf')
    },
    content = function(file) {
      pdf(file,width = isolate({input$unicoxFPwidthdl}),height = isolate({input$unicoxFPheightdl})) # open the pdf device
      print(forestp())
      dev.off()
    }
  )
  
  
  observeEvent(input$page_before_unicox, {
    newtab <- switch(input$tabs, "wgcna" = "msurv","msurv" = "wgcna")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  observeEvent(input$page_after_unicox, {
    newtab <- switch(input$tabs, "DEG" = "msurv","msurv" = "DEG")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  #DEG.R
  observeEvent(rawdata(), {
    data <- rawdata()$clinical
    updateSelectizeInput(
      session,
      "DEGfactor",
      choices = {
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
      selected = "Sex"
    )
  })
  
  DEGana <- eventReactive(input$DEGbt,{
    input$DEGbt
    data <- isolate({ rawdata()})
    DEGmethod <- isolate({input$DEGmethod})
    
    DEGFDR <- isolate({input$DEGFDR})
    DEGFC <- isolate({input$DEGFC})
    DEGfactor <- isolate({input$DEGfactor})
    
    DEGpadj <- isolate({input$DEGpadj})
    mlDEGsel<-isolate({input$mlDEGsel})
    mlgeneselp<-isolate({input$mlgeneselp})
    mlgeneselogFC<-isolate({input$mlgeneselogFC})
    mlgeneselp1<-isolate({input$mlgeneselp1})
    mlgeneselogFC1<-isolate(input$mlgeneselogFC1)
    
    clinical <- data$clinical
    expres <- data$expres
    
    if("" %in% row.names(expres)){
      expres<-expres[-which(row.names(expres)==""),]
    }
    
    clinical[clinical == ""] <- NA
    clinical<-clinical[complete.cases(clinical[,DEGfactor]),]
    
    if(is.numeric(clinical[,DEGfactor])==F){
      clinical[,DEGfactor]<-gsub(' ','',clinical[,DEGfactor])
      clinical[,DEGfactor]<-gsub('-','_',clinical[,DEGfactor])
    }
    
    index <- intersect(names(expres), row.names(clinical))
    clinical <- clinical[index, ]
    expres <- expres[, index]
    
    expres<-as.data.frame(rmz(expres,per=0.9))
    
    group <- factor(clinical[, DEGfactor])
    
    if(is.numeric(clinical[,DEGfactor])){
      createAlert(
        session,
        "DEGmsg",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The group varibale you selected is not categorical variable, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    }else if(length(levels(group))==1 || length(levels(group))>10){
      
      createAlert(
        session,
        "DEGmsg",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The levels group varibale you selected equals to 1 or greater than 10, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    }else{
      withProgress(message = 'Performing differential expression analysis using Limma',
                   detail = 'This may take a while...',
                   value = 3,
                   {
                     expres <- as.matrix(expres)
                     quant <- as.numeric(quantile(expres, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
                     
                     
                     
                     val <- (quant[5] > 100) ||
                       (quant[6] - quant[1] > 50 &&
                          quant[2] > 0)
                     
                     if (val) {
                       expres[which(expres <= 0)] <- NaN
                       expres <- log2(expres)
                     }
                     
                     expres <- as.data.frame(expres)
                     design <- model.matrix( ~ group + 0, expres)
                     colnames(design) <- levels(group)
                     fit <- lmFit(expres, design)
                     Grp <- levels(group)
                     
                     if (length(Grp) == 2) {
                       comp <- paste(Grp[1], Grp[2], sep = "-")
                     } else if (length(Grp) > 2) {
                       comp <- paste(Grp, c(tail(Grp,-1), head(Grp, 1)), sep = "-")
                     }
                     
                     cont.matrix <-makeContrasts(contrasts = comp, levels = design)
                     fit2 <- contrasts.fit(fit, cont.matrix)
                     fit2 <- eBayes(fit2)
                     
                     DEG <- topTable(fit2, adjust.method = DEGpadj,number =(nrow(expres)))
                     names(DEG) <- c("LogFC", "AveExpr", "t", "PValue", "AdjustedP", "B")
                     res<-list(DEG=DEG, group=group, expres=expres,clinical=clinical)
                   })
      
      return(res)
    }
    
  })
  
  observeEvent(input$DEGbt, {
    output$DEGtable <-  DT::renderDT({
      signif(DEGana()$DEG, 3)
    })
  })
  
  observe({
    if(is.null(input$DEGfactor) || input$DEGfactor == ""){
      disable("DEGbt")
    }
    else{
      enable("DEGbt")
    }
  })
  
  observeEvent(input$DEGbt, {
    disable("DEGbt")
  })
  
  observeEvent(DEGana(), {
    shinyjs::show("DEG_wrapper")
    enable("DEGbt")
    shinyjs::show("DEG1_wrapper")
  })
  
  output$downloadDEGtable <- downloadHandler(
    filename = function() {
      paste0("DEG-analysis-result-table-",Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(DEGana()$DEG, file)
    }
  )
  
  
  DEGop <- eventReactive(input$DEGvisbt,{
    input$DEGvisbt
    DEG<-isolate({DEGana()$DEG})
    
    mlDEGsel<-isolate({input$mlDEGsel})
    mlgeneselp<-isolate({input$mlgeneselp})
    mlgeneselogFC<-isolate({input$mlgeneselogFC})
    mlgeneselp1<-isolate({input$mlgeneselp1})
    mlgeneselogFC1<-isolate(input$mlgeneselogFC1)
    if(mlDEGsel=="Adjusted P"){
      sigDEG<-DEG[DEG$AdjustedP<mlgeneselp,]
    } else if(mlDEGsel=="LogFC"){
      sigDEG<-DEG[abs(DEG$LogFC)>mlgeneselogFC,]
    } else{
      sigDEG<-DEG[DEG$AdjustedP<mlgeneselp1 & abs(DEG$LogFC)>mlgeneselogFC1,]
    }
    sigDEG<-sigDEG[complete.cases(sigDEG$AveExpr),]
    if(nrow(sigDEG)==0){
      createAlert(
        session,
        "DEGvismess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste(
          "No genes meet the DEG significance cutoff you specified"),
        append = T,
        dismiss=T
      )
      return(NULL)
    } else{
      list(Output=sigDEG, DEG=DEG, group=isolate({DEGana()$group}), expres=isolate({DEGana()$expres}),clinical=isolate({DEGana()$clinical}))
      
    }
  })
  
  
  observeEvent(input$DEGvisbt, {
    output$DEGog <-  DT::renderDT({
      signif(DEGop()$Output, 3)
    })
  })
  
  observeEvent(DEGop(), {
    shinyjs::show("DEGog_wrapper")
    
    shinyjs::show("DEGvis_wrapper")
    shinyjs::show("DEG2_wrapper")
  })
  output$downloadopgene <- downloadHandler(
    filename = function() {
      paste0("Significantlly-different-expression-gene-table-",Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(DEGop()$Output, file)
    }
  )
  
  
  heatmap<- eventReactive(input$DEGvisbt,{
    input$DEGvisbt
    withProgress(message = "Drawing heatmap",
                 detail = "This may take a while...",
                 value = 3,
                 {
                   DEG<-isolate({DEGop()$DEG})
                   group<-isolate({DEGop()$group})
                   expres<-isolate({DEGop()$expres})
                   sigDEG<-isolate({DEGop()$Output})
                   
                   color <- c()
                   level <- levels(group)
                   htm(
                     DEG=DEG,
                     sigDEG=sigDEG,
                     group=group,
                     color=color,
                     DEGscale=isolate({input$DEGscale}),
                     expres=expres,
                     coltitsize=isolate({input$coltitsize}),
                     ClusterR=isolate({input$ClusterR}),
                     cluster_row_slices=isolate({input$cluster_row_slices}),
                     cludistanrow=isolate({input$cludistanrow}),
                     clumethodrow=isolate({input$clumethodrow}),
                     Rdend_side=isolate({input$Rdend_side}),
                     showRname=isolate({input$showRname}),
                     Rnameside=isolate({input$Rnameside}),
                     showFDR=isolate({input$showFDR}),
                     showFC=isolate({input$showFC}),
                     ClusterC=isolate({input$ClusterC}),
                     heatname=isolate({input$heatname}),
                     cluster_column_slices=isolate({input$cluster_column_slices}),
                     cludistancol=isolate({input$cludistancol}),
                     clumethodcol=isolate({input$clumethodcol}),
                     Cdend_side=isolate({input$Cdend_side}),
                     showCname=isolate({input$showCname}),
                     Cnameside=isolate({input$Cnameside}),
                     heatColors=isolate({input$heatColors})
                   )
                 })
  })
  
  
  vocanaplt<-eventReactive(input$DEGvisbt,{
    input$DEGvisbt
    DEG<-isolate({DEGana()$DEG})
    
    withProgress(message = "Drawing volcano plot",
                 detail = "This may take a while...",
                 value = 3,
                 {
                   
                   voca(
                     DEG = DEG,
                     vocaXlab = isolate({input$vocaXlab}),
                     vocaYlab = isolate({input$vocaYlab}),
                     vocacutp = isolate({input$vocacutp}),
                     vocacutfc = isolate({input$vocacutfc}),
                     legendPosition = isolate({input$legendPosition})
                   )
                 })
    
  })
  
  
  MAplot<-eventReactive(input$DEGvisbt,{
    input$DEGvisbt
    
    DEG<-isolate({DEGana()$DEG})
    withProgress(message = "Drawing MA plot",
                 detail = "This may take a while...",
                 value = 3,
                 {
                   MAplt(
                     DEG=DEG,
                     DEGmethod=isolate({input$DEGmethod}),
                     MAstopmeth=isolate({input$MAstopmeth}),
                     MAcutp=isolate({input$MAcutp}),
                     MAfc=isolate({input$MAfc}),
                     Topgene=isolate({input$Topgene}),
                     MAgenesym=isolate({input$MAgenesym}),
                     MAXlab=isolate({input$MAXlab}),
                     MAYlab=isolate({input$MAYlab}),
                     MAlegendPosition=isolate({input$MAlegendPosition})
                   )
                 })
  })
  
  adjplot<-eventReactive(input$DEGvisbt,{
    input$DEGvisbt
    DEG<-isolate({DEGana()$DEG})
    col<-isolate({input$padjcol})
    withProgress(message = "Drawing adjusted P plot",
                 detail = "This may take a while...",
                 value = 3,
                 {
                   padjplt(DEG,col)
                 })
  }
  )
  
  observeEvent(DEGana(), {
    data <- DEGana()[[1]]
    updateSelectizeInput(
      session,
      "MAgenesym",
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
      selected = "TP53"
    )
  })
  
  
  
  
  observeEvent(input$DEGvisbt, {
    output$DEGvising  <- renderPlot({
      closeAlert(session, "DEGvismess")
      if(isolate({input$DEGvismeth})=="Heatmap"){
        print(heatmap())
        
      } else if(isolate({input$DEGvismeth})=="Volcano plot"){
        print(vocanaplt())
        
      } else if(isolate({input$DEGvismeth})=="MA plot"){
        MAplot()
        
      } else {
        adjplot()
      }
    })
  })
  
  observeEvent(input$DEGvisbt, {
    # updateCollapse(session, "collapsesurvivalplot", open = "Survival plot")
    output$DEGvis<- renderUI({
      plotOutput("DEGvising",
                 width = paste0(isolate({input$padjwidth}), "%"),
                 height = isolate({input$padjheight}))
    })})
  
  
  output$downloadDEGvis <- downloadHandler(
    
    filename = function() {
      paste0("DEG-visualization-plot-",Sys.Date(),'.pdf')
    },
    content = function(file) {
      pdf(file,width = isolate({input$DEGviswidthdl}),height = isolate({input$DEGvisheightdl})) # open the pdf device
      if(isolate({input$DEGvismeth})=="Heatmap"){
        print(heatmap())
        
      } else if(isolate({input$DEGvismeth})=="Volcano plot"){
        print(vocanaplt())
        
      } else if(isolate({input$DEGvismeth})=="MA plot"){
        print(MAplot())
        
      } else {
        print(adjplot())
      }
      dev.off()
    }
  )
  
  observeEvent(input$page_before_DEG, {
    newtab <- switch(input$tabs, "DEG" = "msurv","msurv" = "DEG")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  observeEvent(input$page_after_DEG, {
    newtab <- switch(input$tabs, "DEG" = "nestr","nestr" = "DEG")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  #WGCNA.R
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
  
  #nestr.R
  observe({
    data <- rawdata()$clinical
    updateSelectizeInput(session,
                         'nrtime',
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
                         'nrstatus',
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
  
  value <- reactiveValues(data = NULL)
  observeEvent(input$MTRbt, {
    value$data <- isolate({netsummary()})
    
  })
  
  observeEvent(input$DEGvisbt, {
    value$data <- isolate({DEGop()})
    
  })
  
  observeEvent(input$msurvopbt, {
    value$data <- isolate({msurvsig()})
  })
  
  nestrun<-eventReactive(input$nrbt,{
    input$nrbt
    req(value$data)
    
    data <- isolate(value$data)
    
    withProgress(
      message = "Starting resample experiment",
      detail = "According to your parameter settings, it may take a long time, please be patient",
      value = 5,{
        seed<-as.numeric(isolate({input$seed}))
        
        OS.time<-isolate({input$nrtime})
        OS<-isolate({input$nrstatus})
        clinical <- data$clinical
        expres <- data$expres
        
        clinical<-clinical[complete.cases(clinical[[OS.time]]) & complete.cases(clinical[[OS]]) & clinical[[OS.time]]>0,]
        
        os.time<-clinical[,OS.time]
        os<-clinical[,OS]
        if(mode(os.time)!="numeric" || mode(os)!="numeric" ){
          createAlert(
            session,
            "nestrmess",
            "exampleAlert",
            title = "Please note!",
            style =  "danger",
            content = "Either the survival time or survival status column is not numeric data, please check!",
            append = T,
            dismiss=T
          )
          return(NULL)
          
        }else if(any(levels(factor(os))== c("0", "1"))!=T){
          createAlert(
            session,
            "nestrmess",
            "exampleAlert",
            title = "Please note!",
            style =  "danger",
            content = "The survival status column is not numerically in 1, and 0, please check!",
            append = T,
            dismiss=T
          )
          
          return(NULL)
        }else{
          
          index <- intersect(names(expres), row.names(clinical))
          clinical <- clinical[index, ]
          expres <- expres[, index]
          expres<-as.data.frame(t(expres))
          genesel<-row.names(data$Output)
          expres<-expres[,genesel]
          names(expres)<-make.names(names(expres))
          
          optimization<-isolate({input$optimization})
          if(optimization=="Grid search"){
            ctrl = makeTuneControlGrid(resolution=isolate({input$resolution}))
          } else{
            ctrl = makeTuneControlRandom(maxit=isolate({input$maxit}))
          }
          
          inner = makeResampleDesc("CV", iters = isolate({input$ifold}))
          outer = makeResampleDesc("CV", iters = isolate({input$ofold}))
          
          Learner<- isolate({input$Learner})
          if("RandomForestSRC" %in% Learner){
            
            createAlert(
              session,
              "nestrmess",
              "exampleAlert",
              title = "Please note!",
              style =  "info",
              content = "Random forest is a machine learning method based on ensemble learning. Its calculation process will consume a lot of time and computing resources, especially when integrated into nested cross validation. Considering the limited computing power of the server, we recommend that users download and install the standalone application of CBioExplorer to perform related calculations locally.",
              append = T,
              dismiss=T
            )
          }
          
          split<-isolate({input$split})
          ratio<-isolate({input$sratio})
          biter<-isolate({input$biter})
          
          rdesc<-makeResampleDesc("Bootstrap", iters = biter)
          
          if(isolate({input$validat})=="Nested cross-validation"){
            res<-nestml(
              expres = expres,
              clinical = clinical,
              endpoint = c(OS.time,OS),
              ratio = ratio,
              ctrl = ctrl,
              inner = inner,
              outer = outer,
              lnid = Learner,
              split=split,
              seed=seed
            )
            
          }else {
            res<-cv(expres=expres,
                    clinical=clinical,
                    split=split,
                    endpoint =c(OS.time,OS),
                    ratio=ratio,
                    lnid=Learner,
                    inner=inner,
                    ctrl=ctrl,
                    rdesc=rdesc,
                    seed=seed
            )
            
          }
          if(is.null(res)){
            createAlert(
              session,
              "nestrmess",
              "exampleAlert",
              title = "Please note!",
              style =  "danger",
              content = paste(
                "At least one iteration of the model failed to converge, please reset the data or parameters."),
              append = T,
              dismiss=T
            )
            return(NULL)
          }else if(length(res[[2]])==0){
            createAlert(
              session,
              "nestrmess",
              "exampleAlert",
              title = "Please note!",
              style =  "danger",
              content = paste(
                "No biomarkers identified based on the benchmark experiment you specified."),
              append = T,
              dismiss=T
            )
            return(NULL)
            
          }else{
            res
          } }
      }
    )
  })
  
  observe({
    if(
      is.null(value$data)||
      is.null(input$Learner) ||
      input$Learner == "" ||
      is.null(input$validat) ||
      input$validat == "" ||
      is.null(input$nrtime) ||
      input$nrtime == "" ||
      is.null(input$nrstatus) ||
      input$nrstatus == ""
    ){
      disable("nrbt")
    }
    else{
      enable("nrbt")
    }
  })
  
  observeEvent(input$nrbt, {
    disable("nrbt")
  })
  
  observeEvent(nestrun(), {
    shinyjs::show("modelcomp_wrapper")
    shinyjs::show("biomarkout_wrapper.table")
    enable("nrbt")
  })
  
  observeEvent(is.null(nestrun()), {
    enable("nrbt")
  })
  
  
  observeEvent(input$nrbt,{
    output$modeldescrip <- renderText(
      {
        paste("Features selected base on", nestrun()$lrnid,":",paste(nestrun()[[2]],collapse = ", "),sep=" ")
      }
    )
  })
  
  observeEvent(input$nrbt, {
    output$biomarkout <-  DT::renderDT({
      nestrun()$listidx
    })
  })
  
  output$savebiomarkout <- downloadHandler(
    filename = function() {
      paste0("Benchmark-experiment-biomark-output-",Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(nestrun()$listidx, file)
    }
  )
  
  
  observeEvent(input$nrbt,{
    output$modelcomp <- renderPlot(
      width = isolate({input$nrwidth}),
      height = isolate({input$nrheight}),
      res = 96,
      {
        df<-as.data.frame(nestrun()$`Benchmark result`)
        df$learner.id<-gsub(".tuned","",df$learner.id)
        ggboxplot(df, "learner.id", "cindex",
                  fill = "learner.id",ylab="C-index",xlab=NULL
        )+theme(axis.title.x = element_blank(),legend.position = 'none')
      }
    )
  })
  
  
  observeEvent(input$nrbt, {
    output$modelcomping  <- renderPlot({
      closeAlert(session, "nestrmess")
      df<-as.data.frame(nestrun()$`Benchmark result`)
      df$learner.id<-gsub(".tuned","",df$learner.id)
      ggboxplot(df, "learner.id", "cindex",
                fill = "learner.id",ylab="C-index",xlab=NULL
      )+theme(axis.title.x = element_blank(),legend.position = 'none')
    })
  })
  
  observeEvent(input$nrbt, {
    output$modelcomp <- renderUI({
      plotOutput("modelcomping",
                 width = paste0(isolate({input$nrwidth}), "%"),
                 height = isolate({input$nrheight}))
    })})
  
  
  output$downloadmodelcomp <- downloadHandler(
    filename = function() {
      paste0("C-index-comparison",Sys.Date(),'.pdf')
    },
    content = function(file) {
      pdf(file,width = isolate({input$modelcompwidthdl}),height = isolate({input$modelcompheightdl})) # open the pdf device
      
      df<-as.data.frame(nestrun()$`Benchmark result`)
      df$learner.id<-gsub(".tuned","",df$learner.id)
      print(ggboxplot(df, "learner.id", "cindex",
                      fill = "learner.id",ylab="C-index",xlab=NULL
      )+theme(axis.title.x = element_blank(),legend.position = 'none'))
      dev.off()
    }
  )
  
  
  
  
  observeEvent(input$page_before_nestr, {
    newtab <- switch(input$tabs, "DEG" = "nestr","nestr" = "DEG")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  observeEvent(input$page_after_nestr, {
    newtab <- switch(input$tabs, "bmpm" = "nestr","nestr" = "bmpm")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  #predmo.R
  observeEvent(rawdata(), {
    data <- rawdata()$clinical
    updateSelectInput(session, "bmtime", choices = {
      if (!is.null(data)) {
        if (class(data) ==  class(data.frame()))
          choices <- as.character(colnames(data))
        if (class(data) !=  class(data.frame()))
          choices <-
            as.character(colnames(as.data.frame(data)))
        choices
      }
    },
    selected = "OS.time"
    # server = TRUE
    )
  })
  
  observeEvent(rawdata(), {
    data <- rawdata()$clinical
    updateSelectInput(session, "bmstatus", choices = {
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
    selected = "OS"
    )
  })
  
  yearlabel <- eventReactive(nestrun(),{
    req(input$bmtime)
    req(nestrun())
    res<-isolate({nestrun()})
    split<-isolate({input$split})
    if(split==F){
      data =res$clinical
    }else{
      data<-res$Testsurv
    }
    bmtime<-isolate({input$bmtime})
    # data <- rawdata()$clinical
    time<-data[,bmtime]
    if(is.numeric(time)==F){
      createAlert(
        session,
        "bmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival time you selected is not numeric data, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    } else {
      if(max(na.omit(data[,bmtime]))>600){
        by<- 365
      } else{
        by<-12
      }
      years<-seq(from=0,to=quantile(na.omit(data[,bmtime]),0.95),by=by)
      years<-years[-1]
      yearlabel<-c()
      for(i in 1:length(years)){
        yearlabel[i]<- paste(i,"year",sep="-")
      }
      yearlabel
    }
  })
  
  observeEvent(yearlabel(), {
    choices <- yearlabel()
    updateSelectInput(session, "bmpredictyear", choices = choices)
  })
  
  rocrun<-eventReactive(input$pmbt,{
    input$pmbt
    res<-isolate({nestrun()})
    
    predyear = isolate({input$bmpredictyear})
    clinical<-res$clinical
    OS.time<-clinical[,isolate({input$bmtime})]
    OS<-clinical[,isolate({input$bmstatus})]
    if(is.null(predyear)){
      createAlert(
        session,
        "bmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste(
          "Prediction years not provided, please provide the prediction years"),
        append = T,
        dismiss=T
      )
      return(NULL)
    }else if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
      createAlert(
        session,
        "bmmess",
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
        "bmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival status column is not numerically in 1, and 0, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    }else{
      
      split<-isolate({input$split})
      withProgress(
        message = "Constructing the prediction model",
        detail = "It may take a while, please be patient",
        value = 5,{
          if(split==F){
            survROC(
              data =res$clinical,
              bmtime = isolate({input$bmtime}),
              bmstatus = isolate({input$bmstatus}),
              marker = "Risk",
              method = isolate({input$bmSurvROCmethod}),
              predyear = predyear,
              cutpoint=isolate(input$bmcutoff)
            )
            
          }else{
            rocp <-
              lapply(
                list(res$Trainsurv, res$Testsurv),
                FUN = survROC,
                bmtime = isolate({input$bmtime}),
                bmstatus = isolate({input$bmstatus}),
                marker = "Risk",
                method = isolate({input$bmSurvROCmethod}),
                predyear = predyear,
                cutpoint = isolate(input$bmcutoff)
              )
            ggarrange(plotlist = rocp,ncol=2,nrow=1,align="hv",labels="AUTO")
          }
        })
    }
  })
  
  observe({
    if (is.null(input$bmtime) ||
        input$bmtime == "" ||
        is.null(input$bmstatus) ||
        input$bmstatus == "" ||
        is.null(input$bmpredictyear) ||
        input$bmpredictyear == "" || is.null(input$bmcoxclinvar) ||
        input$bmcoxclinvar == "")
    {
      disable("pmbt")
    }
    else{
      enable("pmbt")
    }
  })
  
  
  observeEvent(input$pmbt,{
    output$bmsurvROCing  <- renderPlot({
      closeAlert(session, "bmmess")
      rocrun()
    })
  })
  
  observeEvent(input$pmbt, {
    output$bmsurvROC <- renderUI({
      plotOutput("bmsurvROCing",
                 width = paste0(isolate({input$bmSurvROCwidth}), "%"),
                 height = isolate({input$bmSurvROCheight}))
    })})
  
  observeEvent(input$pmbt, {
    shinyjs::show("bmsurvROC_wrapper")
    shinyjs::show("bmKMplot_wrapper")
    shinyjs::show("bmCoxforest_wrapper")
    shinyjs::show("bmCoxtable_wrapper")
    disable("pmbt")
  })
  observeEvent(rocrun(), {
    enable("pmbt")
  })
  
  
  output$savebmsurvROC <- downloadHandler(
    filename = function(){
      paste0("Survival-ROC-validation-of-benchmark-experiment-",Sys.Date(),'.pdf')
    },
    content = function(file){
      pdf(file,width = isolate({input$bmsurvROCwidthdl}),height = isolate({input$bmsurvROCheightdl})) # open the pdf device
      print(rocrun())
      dev.off()
    }
  )
  
  
  kmrun<-eventReactive(input$pmbt,{
    input$pmbt
    res<-isolate({nestrun()})
    split<-isolate({input$split})
    clinical<-res$clinical
    OS.time<-clinical[,isolate({input$bmtime})]
    OS<-clinical[,isolate({input$bmstatus})]
    if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
      createAlert(
        session,
        "bmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "Eihter the survival time or survival status column is not numeric data, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
      
    }else if(any(levels(factor(OS))== c("0", "1"))!=T){
      createAlert(
        session,
        "bmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival status column is not numerically in 1, and 0, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    }else {
      if(split==F){
        km(
          data=res$clinical,
          time=isolate({input$bmtime}),
          status=isolate({input$bmstatus}),
          marker="Risk",
          groupby=isolate({input$bmkmgroupby}),
          ratio= isolate({input$riskcut}),
          value=isolate({input$bmkmgpvalue}),
          high="High risk group",
          low="Low risk group",
          survxlab=isolate({input$bmsurvxlab}),
          survP=isolate({input$bmsurvP}),
          survRT=isolate({input$bmsurvRT}),
          survCI=isolate({input$bmsurvCI}),
          color1=isolate({input$bmkmcolor1}),
          color2=isolate({input$bmkmcolor2})
        )
      } else{
        kmp<-lapply(
          list(res$Trainsurv,res$Testsurv),
          FUN=km,
          time=isolate({input$bmtime}),
          status=isolate({input$bmstatus}),
          marker="Risk",
          groupby=isolate({input$bmkmgroupby}),
          ratio= isolate({input$riskcut}),
          value=isolate({input$bmkmgpvalue}),
          high="High risk group",
          low="Low risk group",
          survxlab=isolate({input$bmsurvxlab}),
          survP=isolate({input$bmsurvP}),
          survRT=isolate({input$bmsurvRT}),
          survCI=isolate({input$bmsurvCI}),
          color1=isolate({input$bmkmcolor1}),
          color2=isolate({input$bmkmcolor2})
        )
        if(isolate({input$bmsurvRT})==T){
          ggarrange(
            ggarrange(
              kmp[[1]]$plot,
              kmp[[1]]$table,
              ncol = 1,
              align = "hv",
              heights = c(1.8, 0.7)
            ),
            ggarrange(
              kmp[[2]]$plot,
              kmp[[2]]$table,
              ncol = 1,
              align = "hv",
              heights = c(1.8, 0.7)
            ),
            ncol = 2,
            align = "hv",
            labels="AUTO"
          ) } else{
            ggarrange(kmp[[1]]$plot,kmp[[2]]$plot,ncol = 2,align = "v",labels="AUTO")
          }
      }
    }
  })
  
  
  observeEvent(input$pmbt,{
    output$bmKMploting  <- renderPlot({
      closeAlert(session, "bmmess")
      kmrun()
    })
  })
  
  observeEvent(input$pmbt, {
    output$bmKMplot <- renderUI({
      plotOutput("bmKMploting",
                 width = paste0(isolate({input$bmsurvwidth}), "%"),
                 height = isolate({input$bmsurvheight}))
    })})
  
  output$savebmKMplot <- downloadHandler(
    filename = function(){
      paste0("KM-curve-of-benchmark-experiment-",Sys.Date(),'.pdf')
    },
    content = function(file){
      pdf(file,width = isolate({input$bmKMplotwidthdl}),height = isolate({input$bmKMplotheightdl})) # open the pdf device
      print( kmrun())
      dev.off()
    }
  )
  
  
  
  observeEvent(rawdata(), {
    data <- rawdata()$clinical
    updateSelectInput(session, "bmcoxclinvar", choices = {
      if (!is.null(data)) {
        if (class(data) ==  class(data.frame()))
          choices <- as.character(colnames(data))
        if (class(data) !=  class(data.frame()))
          choices <-
            as.character(colnames(as.data.frame(data)))
        choices
      }
    },
    selected = "Age"
    )
  })
  
  
  coxrun<-eventReactive(input$pmbt,{
    input$pmbt
    res<-isolate({nestrun()})
    split<-isolate({input$split})
    varname<-isolate({input$bmcoxfeat})
    feature<-isolate({input$bmcoxclinvar})
    clinical<-res$clinical
    clinfeat<-subset(clinical,select=feature)
    name<-names(clinfeat)
    fc<-function(x){
      x<-factor(x)
      x<-length(levels(x))
      return(x)
    }
    clinfeature<-lapply(clinfeat,fc)
    marker<-"Risk"
    
    # clinical<-res$clinical
    OS.time<-clinical[,isolate({input$bmtime})]
    OS<-clinical[,isolate({input$bmstatus})]
    if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
      createAlert(
        session,
        "bmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "Eihter the survival time or survival status column is not numeric data, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
      
    }else if(any(levels(factor(OS))== c("0", "1"))!=T){
      createAlert(
        session,
        "bmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival status column is not numerically in 1, and 0, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    }else if(feature==""){
      createAlert(
        session,
        "bmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste(
          "No clinical covariates are provided, please provide the covarites you want to include in the CoxPH model"),
        append = T,
        dismiss=T
      )
      return(NULL)
      
    }else if(1 %in% clinfeature){
      createAlert(
        session,
        "bmmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste("Variable",paste(name[which(clinfeat==1)]),"you selected has only one level (only one value), which is not accepted by CoxPH model. Please check!"),
        append = T,
        dismiss = T
      )
      return(NULL)
    }else{
      var1 <- strsplit(varname, "|", fixed = T)[[1]]
      if(varname=="" || length(var1)==(length(feature)+1)){
        cox1(
          data = res,
          time = isolate({input$bmtime}),
          status = isolate({input$bmstatus}),
          feature = feature,
          marker = marker,
          maxtick = isolate({input$bmmaxtick}),
          split = split,
          varname = varname,
          legend.pos=isolate({input$bmlegend.pos})
        )
      }else{
        
        createAlert(
          session,
          "bmmess",
          "exampleAlert",
          title = "Please note!",
          style =  "danger",
          content = "The length of variable names does not equal to the clinical variables included in the CoxPH model",
          append = T,
          dismiss = T
        )
        return(NULL)
      }
    }
  })
  
  
  observeEvent(input$pmbt,{
    output$bmCoxforesting  <- renderPlot({
      closeAlert(session, "bmmess")
      if(isolate({input$split})==F){
        coxrun()$plotcox
      }else{
        ggarrange(
          coxrun()[[1]]$plotcox,
          coxrun()[[2]]$plotcox,
          nrow = 2,
          labels="AUTO",
          align ="hv"
        )
      }
      
    })
  })
  
  observeEvent(input$pmbt, {
    output$bmCoxforest <- renderUI({
      plotOutput("bmCoxforesting",
                 width = paste0(isolate({input$bmCoxwidth}), "%"),
                 height = isolate({input$bmCoxheight}))
    })})
  
  output$savebmCoxforest <- downloadHandler(
    filename = function(){
      paste0("Cox-forest-plot-of-benchmark-experiment-",Sys.Date(),'.pdf')
    },
    content = function(file){
      pdf(file,width = isolate({input$bmCoxforestwidthdl}),height = isolate({input$bmCoxforestheightdl})) # open the pdf device
      if(isolate({input$split})==F){
        print(coxrun()$plotcox)
      }else{
        print(ggarrange(
          coxrun()[[1]]$plotcox,
          coxrun()[[2]]$plotcox,
          nrow = 2,
          labels="AUTO",
          align ="hv"
        ))
      }
      dev.off()
    }
  )
  
  observeEvent(input$pmbt, {
    output$bmCoxtable <-  DT::renderDT(
      {
        if(isolate({input$split})==F){
          coxrun()$tablecox
        }else{
          tab<-list(TrainingSet=coxrun()[[1]]$tablecox, TestSet=coxrun()[[2]]$tablecox)
          do.call(rbind,tab)
        }
      }
    )
  })
  
  
  output$savebmCoxtable <- downloadHandler(
    filename = function() {
      paste0("CoxPH-table-of-benchmark-experiment-",Sys.Date(), ".csv")
    },
    content = function(file) {
      
      if(isolate({input$split})==F){
        write.csv(coxrun()$tablecox, file)
      }else{
        tab<-list(TrainingSet=coxrun()[[1]]$tablecox, TestSet=coxrun()[[2]]$tablecox)
        do.call(rbind,tab)
        write.csv(do.call(rbind,tab), file)
      }
      
    }
  )
  
  
  observeEvent(input$page_before_bmpm, {
    newtab <- switch(input$tabs, "bmpm" = "nestr","nestr" = "bmpm")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  observeEvent(input$page_after_bmpm, {
    newtab <- switch(input$tabs, "bmpm" = "valmo","valmo" = "bmpm")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  #valmo.R
  observe({
    data <- rawdata2()$clinical
    updateSelectizeInput(session,
                         'Valtime',
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
    data <- rawdata2()$clinical
    updateSelectizeInput(session,
                         'Valstatus',
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
  
  yearlabel1 <- reactive({
    req(input$Valtime)
    Valtime<-paste(input$Valtime)
    data <- rawdata2()$clinical
    if(Valtime %in% names(data)){
      time<-data[,Valtime]
      if(is.numeric(time)==F){
        createAlert(
          session,
          "valmess",
          "exampleAlert",
          title = "Please note!",
          style =  "danger",
          content = "The survival time you selected is not numeric data, please check!",
          append = T,
          dismiss=T
        )
        return(NULL)
      } else {
        if(max(na.omit(data[,Valtime]))>600){
          by<- 365
        } else{
          by<-12
        }
        years1<-seq(from=0,to=quantile(na.omit(data[,Valtime]),0.95),by=by)
        years1<-years1[-1]
        yearlabel1<-c()
        for(i in 1:length(years1)){
          yearlabel1[i]<- paste(i,"year",sep="-")
        }
        yearlabel1
      }
    }
  })
  
  observeEvent(yearlabel1(), {
    choices <- yearlabel1()
    updateSelectInput(session, "Valpredictyear", choices = choices)
  })
  
  
  
  
  predres<-eventReactive(input$Valbt,{
    input$Valbt
    withProgress(
      message = "Starting resample experiment",
      detail = "According to your parameter settings, it may take a long time, please be patient",
      value = 5,{
        
        data<-isolate({rawdata2()})
        res<-isolate({nestrun()})
        
        expval<-data$expres
        expval<-as.data.frame(expval)
        
        selvar=res$`Fitted model`$features
        clinical<-data$clinical
        OS.time<-clinical[,isolate({input$Valtime})]
        
        OS<-clinical[,isolate({input$Valstatus})]
        if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
          createAlert(
            session,
            "valmess",
            "exampleAlert",
            title = "Please note!",
            style =  "danger",
            content = "Either survival time or survival status you selected is not numeric data, please check!",
            append = T,
            dismiss=T
          )
          return(NULL)
          
        }else if(length(setdiff(selvar, row.names(expval)))!=0){
          createAlert(
            session,
            "valmess",
            "exampleAlert",
            title = "Please note!",
            style =  "danger",
            content = paste(
              paste(paste(setdiff(selvar, row.names(expval)),collapse =","), "were not found in the expression file you provided")
            ),
            append = T,
            dismiss=T
          )
          return(NULL)
        } else if(any(levels(factor(OS))== c("0", "1"))!=T){
          createAlert(
            session,
            "valmess",
            "exampleAlert",
            title = "Please note!",
            style =  "danger",
            content = "The survival status column is not numerically in 1, and 0, please check!",
            append = T,
            dismiss=T
          )
          return(NULL)
        }else{
          
          predm(expclin=data$clinical,
                expval=data$expres,
                time=isolate({input$Valtime}),
                status=isolate({input$Valstatus}),
                lnid=res$lrnid,
                modeltrain=res$`Fitted model`)
        }
      })
    
  })
  
  observe({
    if (
      # is.null(input$data1) ||
      # input$data1 == "" ||
      is.null(input$Valtime) ||
      input$Valtime == "" ||
      is.null(input$Valstatus) ||
      input$Valstatus == "" || is.null(input$Valpredictyear) ||
      input$Valpredictyear == ""||
      is.null(input$Valcoxclinvar)||
      input$Valcoxclinvar==""
    )
    {
      disable("Valbt")
    }
    else{
      enable("Valbt")
    }
  })
  
  
  observeEvent(predres(), {
    shinyjs::show("valsurvROC_wrapper")
    shinyjs::show("valKMplot_wrapper")
    shinyjs::show("valCoxforest_wrapper")
    shinyjs::show("valCoxtable_wrapper")
  })
  
  vrocrun<-eventReactive(input$Valbt,{
    input$Valbt
    withProgress(
      message = "Starting resample experiment",
      detail = "According to your parameter settings, it may take a long time, please be patient",
      value = 5,{
        predyear<-isolate({input$Valpredictyear})
        
        clinical = isolate({predres()})
        
        OS.time<-clinical[,isolate({input$Valtime})]
        OS<-clinical[,isolate({input$Valstatus})]
        if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
          createAlert(
            session,
            "valmess",
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
            "valmess",
            "exampleAlert",
            title = "Please note!",
            style =  "danger",
            content = "The survival status column is not numerically in 1, and 0, please check!",
            append = T,
            dismiss=T
          )
          return(NULL)
        }else if(is.null(predyear)){
          createAlert(
            session,
            "valmess",
            "exampleAlert",
            title = "Please note!",
            style =  "danger",
            content = paste(
              "Prediction years not provided, please provide the prediction years."),
            append = T,
            dismiss=T
          )
          return(NULL)
        }else{
          
          survROC(
            data = isolate({predres()}),
            bmtime = isolate({input$Valtime}),
            bmstatus = isolate({input$Valstatus}),
            marker = "Risk",
            method = isolate({input$ValSurvROCmethod}),
            predyear = predyear,
            cutpoint = isolate({input$Valcutoff})
          )
        }
      }
    )
  })
  
  
  observeEvent(input$Valbt,{
    output$valsurvROCing  <- renderPlot({
      closeAlert(session, "valmess")
      vrocrun()
    })
  })
  
  observeEvent(input$Valbt, {
    output$valsurvROC <- renderUI({
      plotOutput("valsurvROCing",
                 width = paste0(isolate({input$ValSurvROCwidth}), "%"),
                 height = isolate({input$ValSurvROCheight}))
    })})
  
  
  output$savevalsurvROC <- downloadHandler(
    filename = function(){
      paste0("Survival-ROC-in-validation-set-",Sys.Date(),'.pdf')
    },
    content = function(file){
      pdf(file,width = isolate({input$valsurvROCwidthdl}),height = isolate({input$valsurvROCheightdl})) # open the pdf device
      print(vrocrun())
      dev.off()
    }
  )
  
  
  vkmrun<-eventReactive(input$Valbt,{
    input$Valbt
    clinical = isolate({predres()})
    
    OS.time<-clinical[,isolate({input$Valtime})]
    OS<-clinical[,isolate({input$Valstatus})]
    if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
      createAlert(
        session,
        "valmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "Eihter the survival time or survival status column is not numeric data, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
      
    }else if(any(levels(factor(OS))== c("0", "1"))!=T){
      createAlert(
        session,
        "valmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival status column is not numerically in 1, and 0, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    }else{
      km(data=isolate({predres()}),
         time=isolate({input$Valtime}),
         status=isolate({input$Valstatus}),
         marker="Risk",
         groupby=isolate({input$Valkmgroupby}),
         ratio= isolate({input$Valriskcut}),
         value=isolate({input$Valkmgpvalue}),
         high="High risk group",
         low="Low risk group",
         survxlab=isolate({input$Valsurvxlab}),
         survP=isolate({input$ValsurvP}),
         survRT=isolate({input$ValsurvRT}),
         survCI=isolate({input$ValsurvCI}),
         color1=isolate({input$Valkmcolor1}),
         color2=isolate({input$Valkmcolor2})
      )
    }}
  )
  
  observeEvent(input$Valbt,{
    output$valKMploting  <- renderPlot({
      closeAlert(session, "valmess")
      vkmrun()
    })
  })
  
  observeEvent(input$Valbt, {
    
    output$valKMplot <- renderUI({
      plotOutput("valKMploting",
                 width = paste0(isolate({input$Valsurvwidth}), "%"),
                 height = isolate({input$Valsurvheight}))
    })})
  
  
  output$savevalvalKMplot <- downloadHandler(
    filename = function(){
      paste0("KM-curve-in-validation-set-",Sys.Date(),'.pdf')
    },
    content = function(file){
      pdf(file,width = isolate({input$valKMplotwidthdl}),height = isolate({input$valKMplotheightdl})) # open the pdf device
      print( vkmrun())
      dev.off()
    }
  )
  
  
  observeEvent(rawdata2(), {
    data <- rawdata2()$clinical
    updateSelectInput(session, "Valcoxclinvar", choices = {
      if (!is.null(data)) {
        if (class(data) ==  class(data.frame()))
          choices <- as.character(colnames(data))
        if (class(data) !=  class(data.frame()))
          choices <-
            as.character(colnames(as.data.frame(data)))
        choices
      }
    },
    selected = "Age"
    
    )
  })
  
  
  
  vcoxrun<-eventReactive(input$Valbt,{
    input$Valbt
    res<-isolate({predres()})
    
    varname<-isolate({input$Valcoxfeat})
    feature<-isolate({input$Valcoxclinvar})
    marker<-"Risk"
    clinical<-res
    clinfeat<-subset(clinical,select=feature)
    name<-names(clinfeat)
    fc<-function(x){
      x<-factor(x)
      x<-length(levels(x))
      return(x)
    }
    clinfeature<-lapply(clinfeat,fc)
    
    OS.time<-clinical[,isolate({input$Valtime})]
    OS<-clinical[,isolate({input$Valstatus})]
    if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
      createAlert(
        session,
        "valmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "Eihter the survival time or survival status column is not numeric data, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
      
    }else if(any(levels(factor(OS))== c("0", "1"))!=T){
      createAlert(
        session,
        "valmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival status column is not numerically in 1, and 0, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    } else if(feature==""){
      createAlert(
        session,
        "valmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste(
          "No clinical covariates are provided, please provide the covarites you want to include in the CoxPH model"),
        append = T,
        dismiss=T
      )
      return(NULL)
    }else if(1 %in% clinfeature){
      createAlert(
        session,
        "valmess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste("Variable",paste(name[which(clinfeat==1)]),"you selected has only one level (only one value), which is not accepted by CoxPH model. Please check!"),
        append = T,
        dismiss = T
      )
      return(NULL)
    }else{
      feature<-c(feature,marker)
      
      var1 <- strsplit(varname, "|", fixed = T)[[1]]
      
      if(varname=="" || length(var1)==length(feature)){
        
        if(varname==""){
          varname<-NULL
        }else{
          varname <- strsplit(varname, "|", fixed = T)[[1]]
        }
        
        cox(data=res,
            time=isolate({input$Valtime}),
            status=isolate({input$Valstatus}),
            feature=feature,
            maxtick=isolate({input$Valmaxtick}),
            varname=varname,
            legend.pos=isolate({input$valegend.pos})
        )
      }else{
        
        createAlert(
          session,
          "valmess",
          "exampleAlert",
          title = "Please note!",
          style =  "danger",
          content = "The length of variable names does not equal to the clinical variables included in the CoxPH model",
          append = T,
          dismiss = T
        )
        return(NULL)
      }
    }
  })
  
  
  
  observeEvent(input$Valbt,{
    output$valCoxforesting  <- renderPlot({
      closeAlert(session, "valmess")
      vcoxrun()$plotcox
    })
  })
  
  observeEvent(input$Valbt, {
    
    output$valCoxforest <- renderUI({
      plotOutput("valCoxforesting",
                 width = paste0(isolate({input$ValCoxwidth}), "%"),
                 height = isolate({input$ValCoxheight}))
    })})
  
  
  output$savevalvalCoxforest <- downloadHandler(
    filename = function(){
      paste0("Forestplot-in-validation-set-",Sys.Date(),'.pdf')
    },
    content = function(file){
      pdf(file,width = isolate({input$valCoxforestwidthdl}),height = isolate({input$valCoxforestheightdl})) # open the pdf device
      print(vcoxrun()$plotcox)
      dev.off()
    }
  )
  
  
  observeEvent(input$Valbt, {
    output$valCoxtable <-  DT::renderDT(
      {vcoxrun()$tablecox}
    )
  })
  
  
  output$savevalCoxtable <- downloadHandler(
    filename = function() {
      paste0("CoxPH-table-in-the-validation-set-",Sys.Date(), ".csv")
    },
    content = function(file) {
      
      write.csv(vcoxrun()$tablecox, file)
      
    }
  )
  
  
  observeEvent(input$page_before_valmo, {
    newtab <- switch(input$tabs, "bmpm" = "valmo","valmo" = "bmpm")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  observeEvent(input$page_after_valmo, {
    newtab <- switch(input$tabs, "nomo" = "valmo","valmo" = "nomo")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  #nomo.R
  observeEvent(rawdata(), {
    data <- rawdata()$clinical
    updateSelectInput(session, "nomotime", choices = {
      if (!is.null(data)) {
        if (class(data) ==  class(data.frame()))
          choices <- as.character(colnames(data))
        if (class(data) !=  class(data.frame()))
          choices <-
            as.character(colnames(as.data.frame(data)))
        choices
      }
    },
    selected = "OS.time"
    )
  })
  
  observeEvent(rawdata(), {
    data <- rawdata()$clinical
    updateSelectInput(session, "nomostatus", choices = {
      if (!is.null(data)) {
        if (class(data) ==  class(data.frame()))
          choices <- as.character(colnames(data))
        if (class(data) !=  class(data.frame()))
          choices <-
            as.character(colnames(as.data.frame(data)))
        choices
      }
    },
    selected = "OS"
    )
  })
  
  yearlabel2 <- eventReactive(nestrun(),
                              {
                                req(input$nomotime)
                                nomotime<-isolate({input$nomotime})
                                req(nestrun())
                                res<-isolate({nestrun()})
                                split<-isolate({input$split})
                                if(split==F){
                                  data =res$clinical
                                }else{
                                  data<-res$Trainsurv
                                }
                                # data <- rawdata()$clinical
                                
                                time<-data[,nomotime]
                                if(is.numeric(time)==F){
                                  createAlert(
                                    session,
                                    "nomomess",
                                    "exampleAlert",
                                    title = "Please note!",
                                    style =  "danger",
                                    content = "The survival time you selected is not numeric data, please check!",
                                    append = T,
                                    dismiss=T
                                  )
                                  return(NULL)
                                } else{
                                  
                                  if(max(na.omit(data[,nomotime]))>600){
                                    by<- 365
                                  } else{
                                    by<-12
                                  }
                                  years2<-seq(from=0,to=quantile(na.omit(data[,nomotime]),0.95),by=by)
                                  years2<-years2[-1]
                                  yearlabel2<-c()
                                  for(i in 1:length(years2)){
                                    yearlabel2[i]<- paste(i,"year",sep="-")
                                  }
                                  yearlabel2
                                }
                              })
  
  observeEvent(yearlabel2(), {
    choices <- yearlabel2()
    updateSelectInput(session, "nomoypoint", choices = choices)
  })
  
  observeEvent(nestrun(), {
    if(isolate({input$split})==T){
      data<-nestrun()$Trainsurv
    }else{
      data<-nestrun()$clinical
    }
    updateSelectInput(session, "nomovar", choices = {
      if (!is.null(data)) {
        if (class(data) ==  class(data.frame()))
          choices <- as.character(colnames(data))
        if (class(data) !=  class(data.frame()))
          choices <-
            as.character(colnames(as.data.frame(data)))
        choices
      }
    },
    selected = "Risk"
    )
  })
  
  observe({
    if(
      is.null(input$nomotime)||
      input$nomotime == "" ||
      is.null(input$nomostatus) ||
      input$nomostatus == "" ||
      is.null(input$nomovar) ||
      input$nomovar == "" ||
      is.null(input$nomoypoint) ||
      input$nomoypoint == ""
    )
    {
      disable("nomogrambt")
    }
    else{
      enable("nomogrambt")
    }
  })
  
  nomorun<-eventReactive(input$nomogrambt,{
    input$nomogrambt
    withProgress(
      message = "Constructing the nomogram",
      detail = "It may take a while, please be patient",
      value = 5,{
        split<-isolate({input$split})
        res<-isolate({nestrun()})
        if(split==T){
          train<-res$Trainsurv
        }else{
          train<-res$clinical
        }
        variable <-isolate({input$nomovar})
        varlabels<-isolate({input$nomovarlab})
        yearpoint<-isolate({input$nomoypoint})
        time<-isolate({input$nomotime})
        status<-isolate({input$nomostatus})
        
        clinfeat<-subset(train,select=variable)
        name<-names(clinfeat)
        fc<-function(x){
          x<-factor(x)
          x<-length(levels(x))
          return(x)
        }
        clinfeature<-lapply(clinfeat,fc)
        OS.time<-train[,time]
        OS<-train[,status]
        if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
          createAlert(
            session,
            "nomomess",
            "exampleAlert",
            title = "Please note!",
            style =  "danger",
            content = "Eihter the survival time or survival status column is not numeric data, please check!",
            append = T,
            dismiss=T
          )
          return(NULL)
          
        }else if(any(levels(factor(OS))== c("0", "1"))!=T){
          createAlert(
            session,
            "nomomess",
            "exampleAlert",
            title = "Please note!",
            style =  "danger",
            content = "The survival status column is not numerically in 1, and 0, please check!",
            append = T,
            dismiss=T
          )
          return(NULL)
        }else if(1 %in% clinfeature){
          createAlert(
            session,
            "valmess",
            "exampleAlert",
            title = "Please note!",
            style =  "danger",
            content = paste("Variable",paste(name[which(clinfeat==1)]),"you selected has only one level (only one value), which is not accepted by CoxPH model. Please check!"),
            append = T,
            dismiss = T
          )
          return(NULL)
        }else{
          var1 <- strsplit(varlabels, "|", fixed = T)[[1]]
          if(varlabels=="" || length(var1)==length(variable)){
            if(varlabels==""){
              varlabels<-NULL
            }else{
              varlabels <- strsplit(varlabels, "|", fixed = T)[[1]]
            }
            nomo(
              train = train,
              time = time,
              status = status,
              variable = variable,
              varlabels = varlabels,
              yearpoint = yearpoint
            )
            
          } else{
            createAlert(
              session,
              "nomomess",
              "exampleAlert",
              title = "Please note!",
              style =  "danger",
              content = "The length of variable names does not equal to the clinical variables included in the CoxPH model",
              append = T,
              dismiss = T
            )
            return(NULL)
          }
        }
      })
  })
  
  observeEvent(nomorun(), {
    shinyjs::show("nomogram_wrapper")
    shinyjs::show("invalcal_wrapper")
  })
  
  observeEvent(input$nomogrambt,{
    output$nomoploting  <- renderPlot({
      print(plot(nomorun()))
    })
  })
  
  observeEvent(input$nomogrambt, {
    output$nomoplot <- renderUI({
      plotOutput("nomoploting",
                 width = paste0(isolate({input$nomowidth}), "%"),
                 height = isolate({input$nomoheight}))
    })})
  
  output$savenomogram <- downloadHandler(
    filename = function(){
      paste0("Nomogram-",Sys.Date(),'.pdf')
    },
    content = function(file){
      pdf(file,width = isolate({input$nomogramwidthdl}),height = isolate({input$nomogramheightdl})) # open the pdf device
      print(plot(nomorun()))
      dev.off()
    }
  )
  
  yearlabel3 <- eventReactive(nestrun(),{
    req(input$nomotime)
    nomotime<-isolate({input$nomotime})
    split<-isolate({input$split})
    if(split==F){
      data <- nestrun()$clinical
    }else{
      data<-nestrun()$Testsurv
    }
    
    if(max(na.omit(data[,nomotime]))>600){
      by<- 365
    } else{
      by<-12
    }
    years3<-seq(from=0,to=quantile(na.omit(data[,nomotime]),0.95),by=by)
    years3<-years3[-1]
    yearlabel3<-c()
    for(i in 1:length(years3)){
      yearlabel3[i]<- paste(i,"year",sep="-")
    }
    yearlabel3
  })
  
  observeEvent(yearlabel3(), {
    choices <- yearlabel3()[-length(yearlabel3())]
    updateSelectInput(session, "invalidypoint", choices = choices)
  })
  
  
  invalidrun<-eventReactive(input$invalidbt,{
    input$invalidbt
    withProgress(
      message = "Internally validating and calibrating the CoxPH based nomogram with bootstrap",
      detail = "It may take a while, please be patient",
      value = 5,{
        split<-isolate({input$split})
        res<-isolate({nestrun()})
        internalvc(split=split,
                   time=isolate({input$nomotime}),
                   status=isolate({input$nomostatus}),
                   variable=isolate({input$nomovar}),
                   ratio=isolate({input$inratio}),
                   bootstrap=isolate({input$inreps}),
                   yearpoint=isolate({input$invalidypoint}),
                   res=res)
      })
  })
  
  
  observe({
    if(
      is.na(input$invalidypoint)||
      length(input$invalidypoint)==0||
      is.null(input$invalidypoint)||
      input$invalidypoint==""||
      
      is.na(input$inreps)||
      length(input$inreps)==0||
      is.null(input$inreps)||
      input$inreps == "" ||
      
      is.na(input$inratio)||
      length(input$inratio)==0||
      is.null(input$inratio) ||
      input$inratio == ""||
      input$inratio >1 ||
      input$inratio <0
    )
    {
      disable("invalidbt")
    }
    else{
      enable("invalidbt")
    }
  })
  
  
  observeEvent(invalidrun(), {
    shinyjs::show("invalid_wrapper")
    shinyjs::show("incalib_wrapper")
    shinyjs::show("exvalcal_wrapper")
    
  })
  
  observeEvent(input$invalidbt, {
    output$invaliddescrip <- renderText({
      if(isolate({input$split==T})){
        paste("Training set",invalidrun()$idex[[1]],"; Test set",invalidrun()$idex[[2]])
      }else{
        paste(invalidrun()$idex[[1]])
      }
    })
  })
  
  observeEvent(input$invalidbt,{
    output$invaliding  <- renderPlot({
      print(invalidrun()$p)
    })
  })
  
  observeEvent(input$invalidbt, {
    output$invalid <- renderUI({
      plotOutput("invaliding",
                 width = paste0(isolate({input$invalidwidth}), "%"),
                 height = isolate({input$invalidheight}))
    })})
  
  output$saveinvalid <- downloadHandler(
    filename = function(){
      paste0("Internal-validation-",Sys.Date(),'.pdf')
    },
    content = function(file){
      pdf(file,width = isolate({input$invalidwidthdl}),height = isolate({input$invalidheightdl})) # open the pdf device
      print(invalidrun()$p)
      dev.off()
    }
  )
  
  observeEvent(input$invalidbt,{
    output$incalibing  <- renderPlot({
      if(isolate({input$split})==T){
        calibtrain<-invalidrun()$validation$validtrain$calibration
        calibtest<-invalidrun()$validation$validtest$calibration
        yearpoint<-isolate({input$invalidypoint})
        calitration<-append(calibtrain,calibtest)
        yp<-c(yearpoint,yearpoint)
        par(mfrow = c(2, length(yearpoint)))
        for(i in 1:length(calitration)){
          plot(calitration[[i]],xlab=paste("Predicted", yp[i], "Survival Probability"),
               ylab=paste("Actual", yp[i], "survival Probability"),  subtitles=isolate({input$incalibsubtitle}))
          my.label <- paste(LETTERS[i])
          mtext(my.label, side = 3, adj = 0,font = 2,cex = 1,line = 1)
        }
      }else{
        calibration<-result1$validation$validtrain$calibration
        yearpoint<-isolate({input$invalidypoint})
        par(mfrow = c(1, length(yearpoint)))
        for(i in 1:length(calitration)){
          plot(calitration[[i]],xlab=paste("Predicted", yearpoint[i], "Survival Probability"),
               ylab=paste("Actual", yearpoint[i], "survival Probability"),  subtitles=isolate({input$incalibsubtitle}))
          my.label <- paste(LETTERS[i])
          mtext(my.label, side = 3, adj = 0,font = 2,cex = 1,line = 1)
        }
      }
    })
  })
  
  observeEvent(input$invalidbt, {
    output$incalib<- renderUI({
      plotOutput("incalibing",
                 width = paste0(isolate({input$incalibwidth}), "%"),
                 height = isolate({input$incalibheight}))
    })})
  
  output$saveincalib <- downloadHandler(
    filename = function(){
      paste0("Internal-calibration-",Sys.Date(),'.pdf')
    },
    content = function(file){
      pdf(file,width = isolate({input$incalibwidthdl}),height = isolate({input$incalibheightdl})) # open the pdf device
      if(isolate({input$split})==T){
        calibtrain<-invalidrun()$validation$validtrain$calibration
        calibtest<-invalidrun()$validation$validtest$calibration
        yearpoint<-isolate({input$invalidypoint})
        calitration<-append(calibtrain,calibtest)
        yp<-c(yearpoint,yearpoint)
        par(mfrow = c(2, length(yearpoint)))
        for(i in 1:length(calitration)){
          plot(calitration[[i]],xlab=paste("Predicted", yp[i], "Survival Probability"),
               ylab=paste("Actual", yp[i], "survival Probability"),  subtitles=isolate({input$incalibsubtitle}))
          my.label <- paste(LETTERS[i])
          mtext(my.label, side = 3, adj = 0,font = 2,cex = 1,line = 1)
        }
      }else{
        calibration<-result1$validation$validtrain$calibration
        yearpoint<-isolate({input$invalidypoint})
        par(mfrow = c(1, length(yearpoint)))
        for(i in 1:length(calitration)){
          plot(calitration[[i]],xlab=paste("Predicted", yearpoint[i], "Survival Probability"),
               ylab=paste("Actual", yearpoint[i], "survival Probability"),  subtitles=isolate({input$incalibsubtitle}))
          my.label <- paste(LETTERS[i])
          mtext(my.label, side = 3, adj = 0,font = 2,cex = 1,line = 1)
        }
      }
      dev.off()
    }
  )
  
  observeEvent(rawdata2(), {
    data <- rawdata2()$clinical
    updateSelectInput(session, "extvaltime", choices = {
      if (!is.null(data)) {
        if (class(data) ==  class(data.frame()))
          choices <- as.character(colnames(data))
        if (class(data) !=  class(data.frame()))
          choices <-
            as.character(colnames(as.data.frame(data)))
        choices
      }
    },
    selected = "OS.time"
    )
  })
  
  observeEvent(rawdata2(), {
    data <- rawdata2()$clinical
    updateSelectInput(session, "extvalstatus", choices = {
      if (!is.null(data)) {
        if (class(data) ==  class(data.frame()))
          choices <- as.character(colnames(data))
        if (class(data) !=  class(data.frame()))
          choices <-
            as.character(colnames(as.data.frame(data)))
        choices
      }
    },
    selected = "OS"
    )
  })
  
  observeEvent(rawdata2(), {
    data <- rawdata2()$clinical
    updateSelectInput(session, "extvalvar", choices = {
      if (!is.null(data)) {
        if (class(data) ==  class(data.frame()))
          choices <- as.character(colnames(data))
        if (class(data) !=  class(data.frame()))
          choices <-
            as.character(colnames(as.data.frame(data)))
        choices
      }
    },
    selected = "Age"
    )
  })
  ####################################################################################################################
  
  yearlabel4 <- reactive({
    req(input$extvaltime)
    nomotime<-isolate({input$extvaltime})
    data <- rawdata2()$clinical
    
    time<-data[,nomotime]
    if(is.numeric(time)==F){
      createAlert(
        session,
        "exnomomess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival time you selected is not numeric data, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    } else{
      if(max(na.omit(data[,nomotime]))>600){
        by<- 365
      } else{
        by<-12
      }
      years4<-seq(from=0,to=quantile(na.omit(data[,nomotime]),0.95),by=by)
      years4<-years4[-1]
      yearlabel4<-c()
      for(i in 1:length(years4)){
        yearlabel4[i]<- paste(i,"year",sep="-")
      }
      yearlabel4
    }
  })
  
  observeEvent(yearlabel4(), {
    choices <- yearlabel4()[-length(yearlabel4())]
    updateSelectInput(session, "exvalidypoint", choices = choices)
  })
  
  
  extvcrun<-eventReactive(input$exvalidbt,{
    input$exvalidbt
    data<- isolate({rawdata2()})
    time=isolate({input$extvaltime})
    status=isolate({input$extvalstatus})
    variable = isolate({input$extvalvar})
    
    clinical<-data$clinical
    OS.time<-clinical[,time]
    OS<-clinical[,status]
    
    clinfeat<-subset(clinical,select=variable)
    name<-names(clinfeat)
    fc<-function(x){
      x<-factor(x)
      x<-length(levels(x))
      return(x)
    }
    clinfeature<-lapply(clinfeat,fc)
    
    if(is.numeric(OS.time)==F || is.numeric(OS)==F ){
      createAlert(
        session,
        "exnomomess",
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
        "exnomomess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = "The survival status column is not numerically in 1, and 0, please check!",
        append = T,
        dismiss=T
      )
      return(NULL)
    } else if(1 %in% clinfeature){
      createAlert(
        session,
        "exnomomess",
        "exampleAlert",
        title = "Please note!",
        style =  "danger",
        content = paste("Variable",paste(name[which(clinfeat==1)]),"you selected has only one level (only one value), which is not accepted by CoxPH model. Please check!"),
        append = T,
        dismiss = T
      )
      return(NULL)
    } else{
      
      withProgress(
        message = "Externally validating and calibrating the CoxPH based nomogram with bootstrap",
        detail = "It may take a while, please be patient",
        value = 5,{
          res<-isolate({nestrun()})
          extvalid(
            extdata = isolate({rawdata2()}),
            variable = variable,
            res = res,
            ratio = isolate({input$exratio}),
            boostrap = isolate({input$exreps}),
            time=time,
            status=status,
            yearpoint=isolate({input$exvalidypoint})
          )
        })
      
    }
  })
  
  observe({
    if(
      is.na(input$extvaltime)||
      is.null(input$extvaltime)||
      input$extvaltime == "" ||
      
      is.null(input$extvalstatus)||
      is.na(input$extvalstatus)||
      input$extvalstatus== "" ||
      
      is.na(input$extvalvar)||
      is.null(input$extvalvar)||
      input$extvalvar== "" ||
      
      is.na(input$exvalidypoint) ||
      is.null(input$exvalidypoint) ||
      input$exvalidypoint == "" ||
      
      is.null(input$exreps) ||
      is.na(input$exreps) ||
      input$exreps == "" ||
      input$exratio > 1 ||
      input$exratio< 0||
      is.null(input$exratio) ||
      is.na(input$exratio) ||
      input$exratio == ""
    )
    {
      disable("exvalidbt")
    }
    else{
      enable("exvalidbt")
    }
  })
  
  
  observeEvent(input$exvalidbt, {
    disable("exvalidbt")
  })
  
  
  observeEvent(extvcrun(), {
    shinyjs::show("exvalid_wrapper")
    shinyjs::show("excalib_wrapper")
    enable("exvalidbt")
  })
  
  
  observeEvent(input$exvalidbt, {
    output$exvaliddescrip <- renderText({
      paste(extvcrun()$idex[[1]])
    })
  })
  
  observeEvent(input$exvalidbt,{
    output$exvaliding  <- renderPlot({
      print(extvcrun()$p)
    })
  })
  
  observeEvent(input$exvalidbt, {
    output$exvalid <- renderUI({
      plotOutput("exvaliding",
                 width = paste0(isolate({input$exvalidwidth}), "%"),
                 height = isolate({input$exvalidheight}))
    })})
  
  output$saveexvalid <- downloadHandler(
    filename = function(){
      paste0("External-validation-",Sys.Date(),'.pdf')
    },
    content = function(file){
      pdf(file,width = isolate({input$exvalidwidthdl}),height = isolate({input$exvalidheightdl})) # open the pdf device
      print(extvcrun()$p)
      dev.off()
    }
  )
  
  observeEvent(input$exvalidbt,{
    output$extcalibing  <- renderPlot({
      
      calibration<-extvcrun()$validation$validexternal$calibration
      yearpoint<-isolate({input$exvalidypoint})
      par(mfrow = c(1, length(yearpoint)))
      for(i in 1:length(calibration)){
        plot(calibration[[i]],xlab=paste("Predicted", yearpoint[i], "Survival Probability"),
             ylab=paste("Actual", yearpoint[i], "survival Probability"),  subtitles=isolate({input$excalibsubtitle}))
        my.label <- paste(LETTERS[i])
        mtext(my.label, side = 3, adj = 0,font = 2,cex = 1,line = 1)
      }
    })
  })
  
  observeEvent(input$exvalidbt, {
    output$extcalib<- renderUI({
      plotOutput("extcalibing",
                 width = paste0(isolate({input$excalibwidth}), "%"),
                 height = isolate({input$excalibheight}))
    })})
  
  output$saveexcalib <- downloadHandler(
    filename = function(){
      paste0("External-calibration-",Sys.Date(),'.pdf')
    },
    content = function(file){
      pdf(file,width = isolate({input$excalibwidthdl}),height = isolate({input$excalibheightdl})) # open the pdf device
      
      calibration<-extvcrun()$validation$validexternal$calibration
      yearpoint<-isolate({input$exvalidypoint})
      par(mfrow = c(1, length(yearpoint)))
      for(i in 1:length(calibration)){
        plot(calibration[[i]],xlab=paste("Predicted", yearpoint[i], "Survival Probability"),
             ylab=paste("Actual", yearpoint[i], "survival Probability"),  subtitles=isolate({input$excalibsubtitle}))
        my.label <- paste(LETTERS[i])
        mtext(my.label, side = 3, adj = 0,font = 2,cex = 1,line = 1)
      }
      dev.off()
    }
  )
  
  observeEvent(input$page_before_nomo, {
    newtab <- switch(input$tabs, "valmo" = "nomo","nomo" = "valmo")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  observeEvent(input$page_after_nomo, {
    newtab <- switch(input$tabs, "nomo" = "clinical","clinical" = "nomo")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  #Bioann.R
  
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
    newtab <- switch(input$tabs, "stemness" = "bioan","bioan" = "stemness")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  #immune.R
  observe({
    data <- rawdata()$expres
    updateSelectizeInput(session,
                         'immgene',
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
  
  immrun<-eventReactive(input$immunebt,{
    input$immunebt
    data <- isolate({
      rawdata()
    })
    expres<-data$expres
    withProgress(message = "Calculating immune cell composition",
                 detail = "This may take a while...",
                 value = 3,
                 {
                   geneset<-isolate({input$geneset})
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
                   # print(idx)
                   if(length(idx)==0){
                     createAlert(
                       session,
                       "immunemess",
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
                   }else
                   {
                     score<-ims(exp=as.matrix(expres),geneset=geneset)
                     res <-
                       plotcor(score = score,
                               exp = as.data.frame(expres),
                               gene=isolate({input$immgene}),
                               select=isolate({input$corselect}),
                               normmethod = isolate({input$immnorm}),
                               fdrcutoff=isolate({input$immcutoff}),
                               type=isolate({input$immmethod}),
                               plottype=isolate({input$immpltype}))
                     
                     if(nrow(res$cordata)==0){
                       createAlert(
                         session,
                         "immunemess",
                         "exampleAlert",
                         title = "Please note!",
                         style =  "danger",
                         content = paste(
                           "No immune cell infiltration score meet the significance cutoff (",
                           isolate({input$immcutoff}),
                           ") at the correction method",
                           paste("'", isolate({input$immnorm}), "'"),
                           "you specified",
                           sep = " "
                         ),
                         append = T,
                         dismiss=T
                       )
                       return(NULL)
                     }else{
                       res
                     }
                   }
                 })
  })
  
  observe({
    if(
      is.null(input$geneset) ||
      input$geneset == "" ||
      is.null(input$immgene) ||
      input$immgene == "" ||
      is.null(input$immmethod) ||
      input$immmethod == "" ||
      is.null(input$immnorm) ||
      input$immnorm == "" ||
      is.null(input$immpltype)||
      input$immpltype==""
    ){
      disable("immunebt")
    }
    else{
      enable("immunebt")
    }
  })
  
  observeEvent(input$immunebt, {
    disable("immunebt")
    
  })
  
  observeEvent(immrun(), {
    shinyjs::show("immuneplot_wrapper")
    shinyjs::show("immunecell_wrapper")
    shinyjs::show("immcor_wrapper")
    enable("immunebt")
  })
  
  observeEvent(is.null(immrun()), {
    enable("immunebt")
  })
  
  
  observeEvent(input$immunebt, {
    output$immcelltable <-  DT::renderDT({
      DT::datatable(immrun()$score,options=list(scrollX=TRUE))
    })
  })
  
  output$saveimmcelltable <- downloadHandler(
    filename = function() {
      paste0("Immune-cell-inflitration-composition-table-",Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(immrun()$score, file)
    }
  )
  
  
  observeEvent(input$immunebt, {
    output$immcortable <-  DT::renderDT({
      DT::datatable( immrun()$cordata,options=list(scrollX=TRUE))
      
    })
  })
  
  output$saveimmcortable <- downloadHandler(
    filename = function() {
      paste0("Immune-cell-inflitration-composition-table-",Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(immrun()$cordata, file)
    }
  )
  
  observeEvent(input$immunebt,{
    output$immuneploting  <- renderPlot({
      closeAlert(session, "immunemess")
      immrun()$corplot
    })
  })
  
  observeEvent(input$immunebt, {
    output$immuneplot<- renderUI({
      plotOutput("immuneploting",
                 width = paste0(isolate({input$immwidth}), "%"),
                 height = isolate({input$immheight}))
    })})
  
  
  
  output$saveimune <- downloadHandler(
    
    filename = function() {
      paste0("Immune-correlation-plot-",Sys.Date(),'.pdf')
    },
    content = function(file) {
      pdf(file,width = isolate({input$immpltwidth}),height = isolate({input$immpltheight}))
      print(immrun()$corplot)
      dev.off()
      
    }
  )
  
  observeEvent(input$page_before_immune, {
    newtab <- switch(input$tabs, "genediff" = "immune","immune" = "genediff")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  observeEvent(input$page_after_immune, {
    newtab <- switch(input$tabs, "immune" = "stemness","stemness" = "immune")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  #stemness.R
  observe({
    data <- rawdata()$expres
    updateSelectizeInput(session,
                         'stemgene',
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
  
  stemrun<-eventReactive(input$stembt,{
    input$stembt
    data <- isolate({
      rawdata()
    })
    expres<-data$expres
    withProgress(message = "Correlating with stemness score",
                 detail = "This may take a while...",
                 value = 3,
                 {
                   stem(exp=expres,gene=isolate({input$stemgene}),method=isolate({input$stemmethod}))
                 })
  })
  
  observe({
    if(
      is.null(input$stemgene) ||
      input$stemgene == "" ||
      is.null(input$stemmethod) ||
      input$stemmethod == ""
    ){
      disable("stembt")
    }
    else{
      enable("stembt")
    }
  })
  
  
  observeEvent(stemrun(), {
    shinyjs::show("stemnessplot_wrapper")
    shinyjs::show("stemtable_wrapper")
  })
  
  
  observeEvent(input$stembt, {
    output$stemtable <-  DT::renderDT({
      stemrun()$stemnesstable
    })
  })
  
  output$savestemtable <- downloadHandler(
    filename = function() {
      paste0("Stemness-score-table-",Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(stemrun()$stemnesstable, file)
    }
  )
  
  
  observeEvent(input$stembt,{
    output$stemnessploting  <- renderPlot({
      closeAlert(session, "stemnessmess")
      stemrun()$scatter
    })
  })
  
  observeEvent(input$stembt, {
    output$stemnessplot<- renderUI({
      plotOutput("stemnessploting",
                 width = paste0(isolate({input$stemwidth}), "%"),
                 height = isolate({input$stemheight}))
    })})
  
  output$savestemness <- downloadHandler(
    filename = function() {
      paste0("Stemness-score-correlation-plot-",Sys.Date(),'.pdf')
    },
    content = function(file) {
      pdf(file,width = isolate({input$stempltwidth}),height = isolate({input$stempltheight}))
      print(stemrun()$scatter)
      dev.off()
    }
  )
  
  observeEvent(input$page_before_stemness, {
    newtab <- switch(input$tabs, "stemness" = "immune","immune" = "stemness")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  observeEvent(input$page_after_stemness, {
    newtab <- switch(input$tabs, "bioan" = "stemness","stemness" = "bioan")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  
  observe({
    if(!is.null(input$test)) stopApp()  # stop shiny
  })
  
}
shinyApp(ui = ui, server = server)
