compdrug<-function(cohort,res,drug,tissueType,batchCorrect,paired,plottype,method){
  if(cohort=="Training set"){
    expres<-as.data.frame(res$traindata$expres)
    clinical<-res$traindata$clinical
  }else{
    expres<-as.data.frame(res$validata$expres)
    clinical<-res$validata$clinical
  }
  
  index<-intersect(row.names(clinical),names(expres))
  clinical<-clinical[index,]
  expres<-expres[,index]
  
  response <- pRRopheticPredict(testMatrix=as.matrix(expres), 
                                drug=drug,
                                tissueType = tissueType, 
                                batchCorrect = batchCorrect,
                                selection=1,
                                dataset = "cgp2016")
  
  if(identical(names(response),row.names(clinical))){
    dat<-matrix(ncol=2,nrow=length(clinical$Subtype))
    dat<-as.data.frame(dat)
    dat[,1]<-as.numeric(response)
    dat[,2]<-clinical$Subtype
  }
  names(dat)<-c("Response","Subtype")
  row.names(dat)<-row.names(clinical)
  
  if(identical(row.names(clinical),row.names(dat))){
    compare<-as.data.frame(combn(unique(as.vector(clinical$Subtype)),2))
    comp<-list()
    for(i in 1:length(compare)){
      comp[[i]]<-compare[[i]]
    }
  }
  name<-paste("Estimated IC50 of",drug)
  pp<-barp(dat,paired=paired,plottype=plottype,method=method,comp=comp,name=name)
  res<-list(plot=pp,table=dat)
  return(res)
}



coroncopath<-function(expres,gene,geneset,method,name){
  genesetscore<-gsva(as.matrix(expres),  geneset, annotation,method = "ssgsea", parallel.sz = 4, verbose = TRUE)
  genesetscore<-as.data.frame(t(genesetscore))
  names(genesetscore)<-name
  if(identical(names(expres),row.names(genesetscore))){
    genesetscore$gene<-as.numeric(expres[gene,])
  }else{
    index<-intersect(names(expres),row.names(genesetscore))
    genesetscore<-genesetscore[index,]
    expres<-expres[,index]
    genesetscore$gene<-as.numeric(expres[gene,])
  }
  
  # library(ggpubr)
  scatter<-list()
  for(i in 1:(ncol(genesetscore)-1)){
    scatter[[i]]<-ggscatter(data=genesetscore, y = paste(names(genesetscore)[i]),x = "gene",
                            size = 1,
                            add = "reg.line", conf.int = TRUE,
                            add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                            cor.coef = TRUE, cor.method = method,cor.coef.size=4.5,
                            cor.coeff.args = list(label.y=max(genesetscore[,i])+se(genesetscore[,i]),label.sep = ","),
                            xlab =paste("Normalized expression of",gene,sep=" ") , ylab = paste(names(genesetscore)[i]))
    
  }
  res<-list(plot=scatter,table=genesetscore)
  return(res)
}


compgeneset<-function(res,cohort,plottype,method,geneset,name){
  if(cohort=="Training set"){
    expres<-as.data.frame(res$traindata$expres)
    clinical<-res$traindata$clinical
  }else{
    expres<-as.data.frame(res$validata$expres)
    clinical<-res$validata$clinical
  }
  
  index<-intersect(row.names(clinical),names(expres))
  clinical<-clinical[index,]
  expres<-expres[,index]
  
  pathscore<-gsva(as.matrix(expres),geneset, annotation,method = "ssgsea", parallel.sz = 4, verbose = TRUE)
  
  pathscore<-as.data.frame(t(pathscore))
  # name<-names(pathscore)
  
  if(identical(row.names(clinical),row.names(pathscore))){
    pathscore$Subtype<-clinical$Subtype
  }else{
    index<-intersect(row.names(clinical),row.names(pathscore))
    clinical<-clinical[index,]
    pathscore<-pathscore[,index]
    pathscore$Subtype<-clinical$Subtype
  }
  compare<-as.data.frame(combn(unique(as.vector(clinical$Subtype)),2))
  comp<-list()
  for(i in 1:length(compare)){
    comp[[i]]<-compare[[i]]
  }
  res<-barp(pathscore,paired=F,plottype=plottype,method=method,comp=comp,name=name)
  res<-list(plot=res,table=pathscore)
  return(res)
}



geo_mean <- function(data){
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

corCYA<-function(expres,gene, method){
  expres<-as.data.frame(t(expres))
  GZMA<-expres[,"GZMA"]
  PRF1<-expres[,"PRF1"]
  cytotoxic.data<-as.data.frame(cbind(GZMA,PRF1))
  row.names(cytotoxic.data)<-row.names(expres)
  GM<-c()
  for(i in 1:nrow(cytotoxic.data)){
    data<-as.numeric(cytotoxic.data[i,])
    GM[i]<-geo_mean(data)
  }
  GM<-as.data.frame(GM)
  row.names(GM)<-row.names(cytotoxic.data)
  if(identical(row.names(expres),row.names(GM))){
    GM$gene<-as.numeric(expres[,gene])
  }
  names(GM)<-c("Cytotoxic activity","gene")
  scatter<-ggscatter(data=GM, y = paste(names(GM)[1]),x = "gene",
                     size = 1,
                     add = "reg.line", conf.int = TRUE,
                     add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                     cor.coef = TRUE, cor.method = method,cor.coef.size=4.5,
                     cor.coeff.args = list(label.y=max(GM[,1])+se(GM[,1]),label.sep = ","),
                     xlab =paste("Normalized expression of",gene,sep=" ") , ylab = paste(names(GM)[1]))
  res<-list(plot=scatter,table=GM)
  return(res)
  
}
compca<-function(expres,clinical,paired,plottype,method,name){
  GZMA<-expres[,"GZMA"]
  PRF1<-expres[,"PRF1"]
  cytotoxic.data<-as.data.frame(cbind(GZMA,PRF1))
  row.names(cytotoxic.data)<-row.names(expres)
  GM<-c()
  for(i in 1:nrow(cytotoxic.data)){
    data<-as.numeric(cytotoxic.data[i,])
    GM[i]<-geo_mean(data)
  }
  GM<-as.data.frame(GM)
  row.names(GM)<-row.names(cytotoxic.data)
  if(identical(row.names(clinical),row.names(GM))){
    compare<-as.data.frame(combn(unique(as.vector(clinical$Subtype)),2))
    comp<-list()
    for(i in 1:length(compare)){
      comp[[i]]<-compare[[i]]
    }
    GM$Subtype<-clinical$Subtype
  }
  pp<-barp(GM,paired=paired,plottype=plottype,method=method,comp=comp,name=name)
  res<-list(plot=pp,table=GM)
  return(res)
}

compIFN<-function(expres,clinical,IFN,plottype,method,name){
  res<-list(IFN)
  IFNscore<-gsva(as.matrix(expres),res, annotation,method = "ssgsea", parallel.sz = 4, verbose = TRUE)
  IFNscore<-as.data.frame(t(IFNscore))
  names(IFNscore)<-"Interferon-gamma score"
  if(identical(row.names(IFNscore),row.names(clinical))){
    IFNscore$Subtype<-clinical$Subtype
  }else{
    index<-intersect(row.names(clinical),row.names(IFNscore))
    IFNscore<-IFNscore[index,]
    clinical<-clinical[index,]
    IFNscore$Subtype<-clinical$Subtype
  }
  compare<-as.data.frame(combn(unique(as.vector(clinical$Subtype)),2))
  comp<-list()
  for(i in 1:length(compare)){
    comp[[i]]<-compare[[i]]
  }
  plot<-barp(IFNscore,paired=F,plottype=plottype,method=method,comp=comp,name=name)
  res<-list(plot=plot,IFN=IFNscore)
  return(res)
}

corIFN<-function(expres,gene, method,IFN ){
  res<-list(IFN)
  IFNscore<-gsva(as.matrix(expres),res, annotation,method = "ssgsea", parallel.sz = 4, verbose = TRUE)
  IFNscore<-as.data.frame(t(IFNscore))
  names(IFNscore)<-"Interferon-gamma score"
  if(identical(row.names(IFNscore),names(expres))){
    IFNscore$gene<-as.numeric(expres[gene,])
  }else{
    index<-intersect(row.names(IFNscore),names(expres))
    IFNscore<-IFNscore[index,]
    expres<-expres[,index]
    IFNscore$gene<-as.numeric(expres[gene,])
  }
  scatter<-ggscatter(data=IFNscore, y = paste(names(IFNscore)[1]),x = "gene",
                     size = 1,
                     add = "reg.line", conf.int = TRUE,
                     add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                     cor.coef = TRUE, cor.method = method,cor.coef.size=4.5,
                     cor.coeff.args = list(label.y=max(IFNscore[,1])+se(IFNscore[,1]),label.sep = ","),
                     xlab =paste("Normalized expression of",gene,sep=" ") , ylab = paste(names(IFNscore)[1]))
  res<-list(plot=scatter,IFNscore=IFNscore)
  return(res)
}

corICB<-function(expres,ICB,gene,method){
  expres<-expres[,c(ICB,gene)]
  names(expres)[ncol(expres)]<-"gene"
  scatter<-list()
  for(i in 1:(ncol(expres)-1)){
    scatter[[i]]<-ggscatter(data=expres, y = paste(names(expres)[i]),x = "gene",
                            size = 1,
                            add = "reg.line", conf.int = TRUE,
                            add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                            cor.coef = TRUE, cor.method = method,cor.coef.size=4.5,
                            cor.coeff.args = list(label.y=max(expres[,i])+se(expres[,i]),label.sep = ","),
                            xlab =paste("Normalized expression of",gene,sep=" ") , ylab = paste(names(expres)[i]))
    
  }
  res<-list(expres=expres,plot=scatter)
  return(res)
}

compicb<-function(ICB,expres,clinical,plottype,method, paired){
  expres<-expres[,c(ICB,"Subtype")]
  name<-names(expres)[-which(names(expres)=="Subtype")]
  compare<-as.data.frame(combn(unique(as.vector(clinical$Subtype)),2))
  comp<-list()
  for(i in 1:length(compare)){
    comp[[i]]<-compare[[i]]
  }
  pp<-barp(expres,paired=paired,plottype=plottype,method=method,comp=comp,name=name)
  res<-list(expres=expres,plot=pp)
  return(res)
}

se <- function(x) sqrt(var(x) / length(x))

outputgct<-function (input.f){
  if (is.data.frame(input.f) == TRUE) {
    exp.data <- input.f
  }
  else {
    exp.data <- read.table(input.f, header = TRUE, row.names = 1, 
                           sep = "\t", quote = "")
  }
  exp.data1 <- data.frame(NAME = rownames(exp.data), Description = rownames(exp.data), 
                          exp.data)
  column1 <- colnames(exp.data1)
  column1[1] <- "NAME"
  column1[2] <- "Description"
  exp.data1$NAME <- factor(exp.data1$NAME)
  exp.data1$Description <- factor(exp.data1$Description)
  levels(exp.data1[, 1]) <- c(levels(exp.data1[, 1]), "NAME")
  levels(exp.data1[, 2]) <- c(levels(exp.data1[, 2]), "Description")
  exp.data2 <- rbind(column1, exp.data1)
  row1 <- rep("", length(1:ncol(exp.data)))
  row1_2 <- data.frame(row1, row1)
  row1_2 <- t(row1_2)
  No_gene <- nrow(exp.data1)
  No_sample <- (ncol(exp.data1) - 2)
  GCT <- matrix(c("#1.2", No_gene, "", No_sample), 
                nrow = 2, ncol = 2)
  gct <- cbind(GCT, row1_2)
  colnames(gct) <- colnames(exp.data2)
  tmp <- rbind(gct, exp.data2)
  return(tmp)
}

EstimateScore<-function (ds, platform = c("affymetrix","agilent", "illumina")){
  require(estimate)   
  descs <- ds[, 1]
  ds <- ds[-1]
  row.names <- row.names(ds)
  names <- names(ds)
  dataset <- list(ds = ds, row.names = row.names, descs = descs, names = names)
  m <- data.matrix(dataset$ds)
  gene.names <- dataset$row.names
  sample.names <- dataset$names
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])
  for (j in 1:Ns) {
    m[, j] <- rank(m[, j], ties.method = "average")
  }
  m <- 10000 * m/Ng
  gs <- as.matrix(SI_geneset[, -1], dimnames = NULL)
  N.gs <- 2
  gs.names <- row.names(SI_geneset)
  score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
  for (gs.i in 1:N.gs) {
    gene.set <- gs[gs.i, ]
    gene.overlap <- intersect(gene.set, gene.names)
    print(paste(gs.i, "gene set:", gs.names[gs.i], 
                " overlap=", length(gene.overlap)))
    if (length(gene.overlap) == 0) {
      score.matrix[gs.i, ] <- rep(NA, Ns)
      next
    }
    else {
      ES.vector <- vector(length = Ns)
      for (S.index in 1:Ns) {
        gene.list <- order(m[, S.index], decreasing = TRUE)
        gene.set2 <- match(gene.overlap, gene.names)
        correl.vector <- m[gene.list, S.index]
        TAG <- sign(match(gene.list, gene.set2, nomatch = 0))
        no.TAG <- 1 - TAG
        N <- length(gene.list)
        Nh <- length(gene.set2)
        Nm <- N - Nh
        correl.vector <- abs(correl.vector)^0.25
        sum.correl <- sum(correl.vector[TAG == 1])
        P0 <- no.TAG/Nm
        F0 <- cumsum(P0)
        Pn <- TAG * correl.vector/sum.correl
        Fn <- cumsum(Pn)
        RES <- Fn - F0
        max.ES <- max(RES)
        min.ES <- min(RES)
        if (max.ES > -min.ES) {
          arg.ES <- which.max(RES)
        }
        else {
          arg.ES <- which.min(RES)
        }
        ES <- sum(RES)
        EnrichmentScore <- list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = TAG)
        ES.vector[S.index] <- EnrichmentScore$ES
      }
      score.matrix[gs.i, ] <- ES.vector
    }
  }
  score.data <- data.frame(score.matrix)
  names(score.data) <- sample.names
  row.names(score.data) <- gs.names
  estimate.score <- apply(score.data, 2, sum)
  if (platform != "affymetrix") {
    score.data <- rbind(score.data, estimate.score)
    rownames(score.data) <- c("StromalScore", "ImmuneScore","ESTIMATEScore")
  }
  else {
    convert_row_estimate_score_to_tumor_purity <- function(x) {
      stopifnot(is.numeric(x))
      cos(0.6049872018 + 0.0001467884 * x)
    }
    est.new <- NULL
    for (i in 1:length(estimate.score)) {
      est_i <- convert_row_estimate_score_to_tumor_purity(estimate.score[i])
      est.new <- rbind(est.new, est_i)
      if (est_i >= 0) {
        next
      }
      else {
        message(paste(sample.names[i], ": out of bounds",sep = ""))
      }
    }
    colnames(est.new) <- c("TumorPurity")
    estimate.t1 <- cbind(estimate.score, est.new)
    x.bad.tumor.purities <- estimate.t1[, "TumorPurity"] < 
      0
    estimate.t1[x.bad.tumor.purities, "TumorPurity"] <- NA
    score.data <- rbind(score.data, t(estimate.t1))
    rownames(score.data) <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")
  }
  score.data<-as.data.frame(t(score.data))
  return(score.data)
}
# PCBC_stemSig<-readRDS("PCBC_stemSig.rds")
stemgroup<-function(exp,clinical,method,plottype,paired){
  #require(TCGAbiolinks)
  exp<-as.data.frame(exp)
  stemness <- TCGAanalyze_Stemness(stemSig = PCBC_stemSig,dataGE = exp)
  index<-intersect(row.names(clinical),row.names(stemness))
  stemness<-stemness[index,]
  clinical<-clinical[index,]
  stemness$Subtype<-clinical$Subtype
  
  if(plottype=="Bar plot"){
    p <-
      ggbarplot(
        stemness,
        x = "Subtype",
        y = "StemnessScore",
        fill = "Subtype",
        xlab = FALSE,
        ylab = "Stemness score",
        font.y = 12,
        font.xtickslab = 12,
        font.ytickslab = 12,
        add = "mean_se"
      ) +theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +theme(legend.position = "top")
  }else if(plottype=="Box plot"){
    p <-
      ggboxplot(
        stemness,
        x = "Subtype",
        y = "StemnessScore",
        fill = "Subtype",
        xlab = FALSE,
        ylab = "Stemness score",
        font.y = 12,
        font.xtickslab = 12,
        font.ytickslab = 12,
        add = "mean_se"
      ) +theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +theme(legend.position = "top")
  }else {
    p <-ggviolin(
      stemness,
      x = "Subtype",
      y = "StemnessScore",
      fill = "Subtype",
      xlab = FALSE,
      ylab = "Stemness score",
      font.y = 12,
      font.xtickslab = 12,
      font.ytickslab = 12,
      add = "boxplot", add.params = list(fill = "white")
    ) +theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +theme(legend.position = "top")
  }
  if(method=="Nonparametric test"){
    meth<-"kruskal.test"
    meth1<-"wilcox.test"
  }else{
    meth<-"anova"
    meth1<-"t.test"
  }
  
  compare<-as.data.frame(combn(unique(as.vector(clinical$Subtype)),2))
  comp<-list()
  for(i in 1:length(compare)){
    comp[[i]]<-compare[[i]]
  }
  
  data<-clinical
  
  if(length(unique(data$Subtype))>2 && paired==T){
    p<-p+stat_compare_means(comparisons = comp,method=meth1, label = "p.signif")+
      stat_compare_means(vjust=T,method=meth)+theme(legend.title=element_blank())
  }else if(length(unique(data$Subtype))>2 && paired==F){
    # p<-p+ stat_compare_means(label.x.npc="center",label.y.npc="top",method=meth1 )
    p<-p+stat_compare_means(vjust=T,method=meth) + theme(legend.title=element_blank())
  } else {
    p<-p+ stat_compare_means(label.x.npc="center",label.y.npc="top",method=meth1)+
      theme(legend.title=element_blank())
  }
  return(p)
}

barp<-function(data,paired,plottype,method,comp,name){
  p<-list()
  p1<-list()
  for(i in 1:(length(names(data))-1)){
    if(plottype=="Bar plot"){
      
      p[[i]] <-
        ggbarplot(
          data,
          x = "Subtype",
          y = names(data)[i],
          fill = "Subtype",
          xlab = FALSE,
          ylab = name[i],
          font.y = 12,
          font.xtickslab = 12,
          font.ytickslab = 12,
          add = "mean_se"
        ) +theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ) +theme(legend.position = "top")
    }else if(plottype=="Box plot"){
      p[[i]] <-
        ggboxplot(
          data,
          x = "Subtype",
          y = names(data)[i],
          fill = "Subtype",
          xlab = FALSE,
          ylab = name[i],
          font.y = 12,
          font.xtickslab = 12,
          font.ytickslab = 12,
          add = "mean_se",
          
        ) +theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ) +theme(legend.position = "top")
    }else {
      p[[i]] <-ggviolin(
        data,
        x = "Subtype",
        y = names(data)[i],
        fill = "Subtype",
        xlab = FALSE,
        ylab = name[i],
        font.y = 12,
        font.xtickslab = 12,
        font.ytickslab = 12,
        add = "boxplot", add.params = list(fill = "white")
      ) +theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +theme(legend.position = "top")
    }
    if(method=="Nonparametric test"){
      meth<-"kruskal.test"
      meth1<-"wilcox.test"
    }else{
      meth<-"anova"
      meth1<-"t.test"
    }
    
    if(length(unique(data$Subtype))>2 && paired==T){
      p1[[i]]<-p[[i]]+stat_compare_means(comparisons = comp,method=meth1, label = "p.signif")+
        stat_compare_means(vjust=T,method=meth)+
        theme(legend.title=element_blank())
    }else if(length(unique(data$Subtype))>2 && paired==F){
      # p1[[i]]<-p[[i]]+ stat_compare_means(label.x.npc="center",label.y.npc="top",method=meth1 )+
      #   theme(legend.title=element_blank())
      p1[[i]]<-p[[i]]+stat_compare_means(vjust=T,method=meth) + theme(legend.title=element_blank())
      # p1[[i]]<-p[[i]]+ stat_compare_means(vjust=T,method=meth)
    } else {
      p1[[i]]<-p[[i]]+ stat_compare_means(label.x.npc="center",label.y.npc="top",method=meth1 )+
        theme(legend.title=element_blank())
    }
  }
  return(p1)
}

assignSubtype <- function(ImmuneMW, dataGE,method){
  if(method=="Spearman"){
    method<-"spearman"
  }else if(method=="Pearson"){
    method=="pearson"
  }
  dataImmuneGroups_merged <- matrix(0,ncol(dataGE),ncol(ImmuneMW))
  rownames(dataImmuneGroups_merged) <- colnames(dataGE)
  colnames(dataImmuneGroups_merged) <- colnames(ImmuneMW)
  dataImmuneGroups_merged <- as.data.frame(dataImmuneGroups_merged)
  for( i in 1: ncol(ImmuneMW)){
    
    cursubtype <- colnames(ImmuneMW)[i]
    print(cursubtype)
    ImmuneMW_cur <- ImmuneMW[,cursubtype,drop=F]
    
    reads <- dataGE
    X <- reads
    w <- ImmuneMW_cur
    commonStemsigGenes <- intersect(row.names(w),rownames(X))
    
    X <- X[commonStemsigGenes,]
    w <- w[ rownames(X), ]
    
    # Score the Matrix X using Spearman correlation.
    
    s <- apply( X, 2, function(z) {cor( z, w, method = method, use = "complete.obs" )} )
    
    ## Scale the scores to be between 0 and 1
    s <- s - min(s)
    s <- s / max(s)
    
    dataSce_immuneSubtypes <- cbind(s)
    dataSce_immuneSubtypes <- as.data.frame(dataSce_immuneSubtypes)
    colnames(dataSce_immuneSubtypes) <- cursubtype
    
    dataImmuneGroups_merged[rownames(dataSce_immuneSubtypes),cursubtype] <- dataSce_immuneSubtypes[,1]
  }
  
  dataImmuneSubtypes <- matrix(0,nrow(dataImmuneGroups_merged),2)
  rownames(dataImmuneSubtypes) <- rownames(dataImmuneGroups_merged)
  colnames(dataImmuneSubtypes) <- c("Sample","Subtype")
  dataImmuneSubtypes <- as.data.frame(dataImmuneSubtypes)
  
  dataImmuneSubtypes$Sample <- rownames(dataImmuneSubtypes)
  
  for( j in 1:nrow(dataImmuneSubtypes)){
    cursample <- dataImmuneSubtypes$Sample[j]
    
    idx <- which(dataImmuneGroups_merged[cursample,] == max(dataImmuneGroups_merged[cursample,]))
    
    if( length(idx)!=1) {
      idx <- idx[1]
    }
    dataImmuneSubtypes[cursample,"Subtype"] <- colnames(dataImmuneGroups_merged)[idx]
    
  }
  
  return(dataImmuneSubtypes)
}


FSVar<-function (Data, cut.type = "TopK", value) {
  vars = apply(Data, 1, var)
  feature_num = length(vars)
  # hist(vars, breaks = feature_num * 0.1, col = "red", 
  #      main = "Expression (Variance) distribution", xlab = "The Variance of feature")
  if (cut.type == "TopK") {
    index = sort(vars, decreasing = TRUE, index.return = TRUE)
    if (value > nrow(Data)) {
      value = nrow(Data)
      cat("Warning: the feature selection number is beyond the original feature numnber")
    }
    cutoff = index$x[value]
    # abline(v = cutoff, col = "blue", lty = 5, lwd = 1.5)
    index = index$ix[1:value]
    selectData = Data[index, ]
  }
  if (cut.type == "Cutoff") {
    # abline(v = value, col = "blue", lty = 5, lwd = 1.5)
    index = which(vars > value)
    selectData = Data[index, ]
  }
  return(selectData)
}

FSMAD<-function (Data, cut.type = "TopK", value) {
  mads = apply(Data, 1, mad)
  feature_num = length(mads)
  # hist(mads, breaks = feature_num * 0.1, col = "red", 
  #      main = "Expression (MAD) distribution", xlab = "The MAD of feature")
  if (cut.type == "TopK") {
    index = sort(mads, decreasing = TRUE, index.return = TRUE)
    if (value > nrow(Data)) {
      value = nrow(Data)
      cat("Warning: the feature selection number is beyond the original feature numnber")
    }
    cutoff = index$x[value]
    # abline(v = cutoff, col = "blue", lty = 5, lwd = 1.5)
    index = index$ix[1:value]
    selectData = Data[index, ]
  }
  if (cut.type == "Cutoff") {
    # abline(v = value, col = "blue", lty = 5, lwd = 1.5)
    index = which(mads > value)
    selectData = Data[index, ]
  }
  return(selectData)
}

FSPCA<-function(Data, PC_percent = 1, scale = TRUE){
  newData = t(Data)
  aa = prcomp(newData, scale = scale)
  vars <- apply(aa$x, 2, var)
  props <- vars/sum(vars)
  xx = as.vector(cumsum(props))
  num = which(xx > PC_percent)[1]
  coeff = aa$ro
  score = (newData %*% coeff)[, 1:num]
  result = t(score)
  return(result)
}

FSCox<-function (Data, time, status, cutoff = cutoff){
  feature_cox_test = NULL
  for (i in 1:nrow(Data)) {
    dataset = list(time, status, x <- Data[i, ])
    temp <- coxph(Surv(time, status) ~ x, dataset)
    feature_cox_test = rbind(feature_cox_test, summary(temp)$coefficients)
  }
  indexNum = which(feature_cox_test[, 5] < cutoff)
  selectData = Data[indexNum, ]
  return(selectData)
}


SBpreproc <- function(df, fsm, FScutmethod, FSvaluet,PC_percent, scale, clin,time,status,cutoff) {
  clin <- clin[complete.cases(clin[, time]) & clin[, time] > 0, ]
  df <- as.data.frame(df)
  idx <- intersect(names(df), row.names(clin))
  df <- df[, idx]
  clin <- clin[idx, ]
  time <- clin[, time]
  status <- clin[, status]
  if (fsm == "Variance") {
    data = FSVar(df, cut.type = FScutmethod, value = FSvaluet)
  } else if (fsm == "MAD") {
    data = FSMAD(df, cut.type = FScutmethod, value = FSvaluet)
  } else if (fsm == "PCA") {
    data = FSPCA(df, PC_percent = PC_percent, scale = scale)
  } else{
    data = FSCox(as.matrix(df), time, status, cutoff = cutoff)
  }
  return(data)
}


sampleCols <- function( d,
                        pSamp=NULL,
                        pRow=NULL,
                        weightsItem=NULL,
                        weightsFeature=NULL ){
  ## returns a list with the sample columns, as well as the sub-matrix & sample features (if necessary)
  ##  if no sampling over the features is performed, the submatrix & sample features are returned as NAs
  ##  to reduce memory overhead
  
  space <- ifelse( inherits( d, "dist" ), ncol( as.matrix(d) ), ncol(d) )
  sampleN <- floor(space*pSamp)
  sampCols <- sort( sample(space, sampleN, replace = FALSE, prob = weightsItem) )
  
  this_sample <- sampRows <- NA
  if ( inherits( d, "matrix" ) ) {
    if ( (! is.null( pRow ) ) &&
         ( (pRow < 1 ) || (! is.null( weightsFeature ) ) ) ) {
      ## only sample the rows and generate a sub-matrix if we're sampling over the row/gene/features
      space = nrow(d)
      sampleN = floor(space*pRow)
      sampRows = sort( sample(space, sampleN, replace = FALSE, prob = weightsFeature) )
      this_sample <- d[sampRows,sampCols]
      dimnames(this_sample) <- NULL
    } else {
      ## do nothing
    }
  }
  return( list( submat=this_sample,
                subrows=sampRows,
                subcols=sampCols ) )
}
connectivityMatrix <- function( clusterAssignments, m, sampleKey){
  ##input: named vector of cluster assignments, matrix to add connectivities
  ##output: connectivity matrix
  names( clusterAssignments ) <- sampleKey 
  cls <- lapply( unique( clusterAssignments ), function(i) as.numeric( names( clusterAssignments[ clusterAssignments %in% i ] ) ) )  #list samples by clusterId
  
  for ( i in 1:length( cls ) ) {
    nelts <- 1:ncol( m )
    cl <- as.numeric( nelts %in% cls[[i]] ) ## produces a binary vector
    updt <- outer( cl, cl ) #product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster;
    m <- m + updt
  }
  return(m)
}

triangle = function(m,mode=1){
  #mode=1 for CDF, vector of lower triangle.
  #mode==3 for full matrix.
  #mode==2 for calcICL; nonredundant half matrix coun
  #mode!=1 for summary 
  n=dim(m)[1]
  nm = matrix(0,ncol=n,nrow=n)
  fm = m
  nm[upper.tri(nm)] = m[upper.tri(m)] #only upper half
  fm = t(nm)+nm
  diag(fm) = diag(m)
  nm=fm
  nm[upper.tri(nm)] = NA
  diag(nm) = NA
  vm = m[lower.tri(nm)]
  if(mode==1){
    return(vm) #vector 		
  }else if(mode==3){
    return(fm) #return full matrix
  }else if(mode == 2){
    return(nm) #returns lower triangle and no diagonal. no double counts.
  }
}

myPal = function(n=10){
  #returns n colors
  seq = rev(seq(0,255,by=255/(n)))
  palRGB = cbind(seq,seq,255)
  rgb(palRGB,maxColorValue=255)
} 

setClusterColors = function(past_ct,ct,colorU,colorList){
  #description: sets common color of clusters between different K
  newColors = c()
  if(length(colorList)==0){
    #k==2
    newColors = colorU[ct]
    colori=2
  }else{
    newColors = rep(NULL,length(ct))
    colori = colorList[[2]]
    mo=table(past_ct,ct)
    m=mo/apply(mo,1,sum)
    for(tci in 1:ncol(m)){ # for each cluster
      maxC = max(m[,tci])
      pci = which(m[,tci] == maxC)				
      if( sum(m[,tci]==maxC)==1 & max(m[pci,])==maxC & sum(m[pci,]==maxC)==1  )  {
        #if new column maximum is unique, same cell is row maximum and is also unique
        ##Note: the greatest of the prior clusters' members are the greatest in a current cluster's members.
        newColors[which(ct==tci)] = unique(colorList[[1]][which(past_ct==pci)]) # one value
      }else{ #add new color.
        colori=colori+1
        newColors[which(ct==tci)] = colorU[colori]
      }
    }
  }
  return(list(newColors,colori,unique(newColors) ))
}

CDF=function(ml,breaks=100){
  #plot CDF distribution
  par(mfrow=c(1,2))
  plot(c(0),xlim=c(0,1),ylim=c(0,1),col="white",bg="white",xlab="Consensus index",ylab="CDF",main="Consensus CDF", las=2)
  k=length(ml)
  this_colors = rainbow(k-1)
  areaK = c()
  for (i in 2:length(ml)){
    v=triangle(ml[[i]],mode=1)
    
    #empirical CDF distribution. default number of breaks is 100    
    h = hist(v, plot=FALSE, breaks=seq(0,1,by=1/breaks))
    h$counts = cumsum(h$counts)/sum(h$counts)
    
    #calculate area under CDF curve, by histogram method.
    thisArea=0
    for (bi in 1:(length(h$breaks)-1)){
      thisArea = thisArea + h$counts[bi]*(h$breaks[bi+1]-h$breaks[bi]) #increment by height by width
      bi = bi + 1
    }
    areaK = c(areaK,thisArea)
    lines(h$mids,h$counts,col=this_colors[i-1],lwd=2,type='l')
  }
  legend(0.005,1,legend=paste(rep("",k-1),seq(2,k,by=1),sep=""),fill=this_colors)
  
  #plot area under CDF change.
  deltaK=areaK[1] #initial auc at k=2
  for(i in 2:(length(areaK))){
    #proportional increase relative to prior K.
    deltaK = c(deltaK,( areaK[i] - areaK[i-1])/areaK[i-1])
  }
  plot(1+(1:length(deltaK)),y=deltaK,xlab="k",ylab="relative change in area under CDF curve",main="Delta area",type="b")
}

ccRun <- function( d=d,
                   maxK=NULL,
                   repCount=NULL,
                   diss=inherits( d, "dist" ),
                   pItem=NULL,
                   pFeature=NULL,
                   innerLinkage=NULL,
                   distance=NULL, #ifelse( inherits(d,"dist"), attr( d, "method" ), "euclidean" ),@@@@@
                   clusterAlg=NULL,
                   weightsItem=NULL,
                   weightsFeature=NULL,
                   verbose=NULL,
                   corUse=NULL) {
  m = vector(mode='list', repCount)
  ml = vector(mode="list",maxK)
  n <- ifelse( diss, ncol( as.matrix(d) ), ncol(d) )
  mCount = mConsist = matrix(c(0),ncol=n,nrow=n)
  ml[[1]] = c(0);
  
  if (is.null( distance ) ) distance <- 'euclidean'  ## necessary if d is a dist object and attr( d, "method" ) == NULL
  
  acceptable.distance <- c( "euclidean", "maximum", "manhattan", "canberra", "binary","minkowski",
                            "pearson", "spearman" )
  
  main.dist.obj <- NULL
  if ( diss ){
    main.dist.obj <- d
    ## reset the pFeature & weightsFeature params if they've been set (irrelevant if d is a dist matrix)
    if ( ( !is.null(pFeature) ) &&
         ( pFeature < 1 ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified pFeature parameter\n" )
      pFeature <- 1 # set it to 1 to avoid problems with sampleCols
    }
    if ( ! is.null( weightsFeature ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified weightsFeature parameter\n" )
      weightsFeature <- NULL  # set it to NULL to avoid problems with sampleCols
    }
  } else { ## d is a data matrix
    ## we're not sampling over the features
    if ( ( clusterAlg != "km" ) &&
         ( is.null( pFeature ) ||
           ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) ) {
      ## only generate a main.dist.object IFF 1) d is a matrix, 2) we're not sampling the features, and 3) the algorithm isn't 'km'
      if ( inherits( distance, "character" ) ) {
        if ( ! distance %in%  acceptable.distance  &  ( class(try(get(distance),silent=T))!="function") ) stop("unsupported distance.")
        
        if(distance=="pearson" | distance=="spearman"){
          main.dist.obj <- as.dist( 1-cor(d,method=distance,use=corUse ))
        }else if( class(try(get(distance),silent=T))=="function"){
          main.dist.obj <- get(distance)( t( d )   )
        }else{
          main.dist.obj <- dist( t(d), method=distance )
        }
        attr( main.dist.obj, "method" ) <- distance  
      } else stop("unsupported distance specified.")
    } else {
      ## pFeature < 1 or a weightsFeature != NULL
      ## since d is a data matrix, the user wants to sample over the gene features, so main.dist.obj is left as NULL
    }
  }
  
  for (i in 1:repCount){
    if(verbose){
      message(paste("random subsample",i));
    }
    ## take expression matrix sample, samples and genes
    sample_x = sampleCols( d, pItem, pFeature, weightsItem, weightsFeature )
    
    this_dist = NA
    if ( ! is.null( main.dist.obj ) ) {
      boot.cols <- sample_x$subcols
      this_dist <- as.matrix( main.dist.obj )[ boot.cols, boot.cols ]
      if ( clusterAlg != "km" ) {
        ## if this isn't kmeans, then convert to a distance object
        this_dist <- as.dist( this_dist )
        attr( this_dist, "method" ) <- attr( main.dist.obj, "method" )
      }
    } else {

      if ( clusterAlg != "km" )  {
        if ( ! distance %in% acceptable.distance &  ( class(try(get(distance),silent=T))!="function")  ) stop("unsupported distance.")
        if( ( class(try(get(distance),silent=T))=="function") ){
          this_dist <- get(distance)( t( sample_x$submat ) )
        }else{
          if( distance == "pearson" | distance == "spearman"){
            this_dist <- as.dist( 1-cor(sample_x$submat,use=corUse,method=distance) )
          }else{
            this_dist <- dist( t( sample_x$submat ), method= distance  )
          }
        }
        attr( this_dist, "method" ) <- distance  
      } else {
        ## if we're not sampling the features, then grab the colslice
        if ( is.null( pFeature ) ||
             ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) {
          this_dist <- d[, sample_x$subcols ]
        } else {
          if ( is.na( sample_x$submat ) ) {
            stop( "error submat is NA" )
          }
          
          this_dist <- sample_x$submat
        } 
      }
    }
    
    ## cluster samples for HC.
    this_cluster=NA
    if(clusterAlg=="hc"){
      this_cluster = hclust( this_dist, method=innerLinkage)
    }
    ##mCount is possible number of times that two sample occur in same random sample, independent of k
    ##mCount stores number of times a sample pair was sampled together.
    mCount <- connectivityMatrix( rep( 1,length(sample_x[[3]])),
                                  mCount,
                                  sample_x[[3]] ) 
    
    ##use samples for each k		
    for (k in 2:maxK){
      if(verbose){
        message(paste("  k =",k))
      }
      if (i==1){
        ml[[k]] = mConsist #initialize
      }
      this_assignment=NA
      if(clusterAlg=="hc"){
        ##prune to k for hc
        this_assignment = cutree(this_cluster,k)
        
        
      }else if(clusterAlg=="km"){
        ##this_dist should now be a matrix corresponding to the result from sampleCols
        this_assignment <- kmeans( t( this_dist ),
                                   k,
                                   iter.max = 10,
                                   nstart = 1,
                                   algorithm = c("Hartigan-Wong") )$cluster
      }else if ( clusterAlg == "pam" ) {
        this_assignment <- pam( x=this_dist,
                                k,
                                diss=TRUE,
                                metric=distance, 
                                cluster.only=TRUE )
      } else{
        ##optional cluterArg Hook.
        this_assignment <- get(clusterAlg)(this_dist, k)
      }
      ##add to tally				
      ml[[k]] <- connectivityMatrix( this_assignment,
                                     ml[[k]],
                                     sample_x[[3]] )
    }
  }
  
  
  ##consensus fraction
  res = vector(mode="list",maxK)
  for (k in 2:maxK){
    ##fill in other half of matrix for tally and count.
    tmp = triangle(ml[[k]],mode=3)
    tmpCount = triangle(mCount,mode=3)
    res[[k]] = tmp / tmpCount
    res[[k]][which(tmpCount==0)] = 0
  }
  # message("end fraction")
  return(res)
}

clusterTrackingPlot = function(m){
  #description: plots cluster tracking plot
  #input: m - matrix where rows are k, columns are samples, and values are cluster assignments.
  plot(NULL,xlim=c(-0.1,1),ylim=c(0,1),axes=FALSE,xlab="samples",ylab="k",main="tracking plot")
  for(i in 1:nrow(m)){
    rect(  xleft=seq(0,1-1/ncol(m),by=1/ncol(m)),  ybottom=rep(1-i/nrow(m),ncol(m)) , xright=seq(1/ncol(m),1,by=1/ncol(m)), ytop=rep(1-(i-1)/nrow(m),ncol(m)), col=m[i,],border=NA)   
  }
  #hatch lines to indicate samples
  xl = seq(0,1-1/ncol(m),by=1/ncol(m))
  segments(  xl, rep(-0.1,ncol(m)) , xl, rep(0,ncol(m)), col="black")    #** alt white and black color?
  ypos = seq(1,0,by=-1/nrow(m))-1/(2*nrow(m))
  text(x=-0.1,y=ypos[-length(ypos)],labels=seq(2,nrow(m)+1,by=1))
}


ccp <- function (d = dc,
                 maxK = 4,
                 reps = 100,
                 pItem = 0.8,
                 pFeature = 1,
                 # title = "example",
                 distance = "pearson",
                 clusterAlg = "hc",
                 ml = NULL,
                 weightsFeature = NULL,
                 corUse = "everything",
                 verbose = F,
                 weightsItem = NULL,
                 innerLinkage = "average",
                 finalLinkage = "average",
                 tmyPal = NULL,
                 seed = NULL
) {
  if (is.null(seed) == TRUE) {
    seed = timeSeed = as.numeric(Sys.time())
  }
  set.seed(seed)
  if (is.null(ml) == TRUE) {
    if (!any(class(d) %in% c("dist", "matrix",
                             "ExpressionSet"))) {
      stop("d must be a matrix, distance object or ExpressionSet (eset object)")
    }
    if (inherits(d, "dist")) {
      if (is.null(attr(d, "method"))) {
        attr(d, "method") <- distance <- "unknown - user-specified"
      }
      if (is.null(distance) || (distance != attr(d, "method"))) {
        distance <- attr(d, "method")
      }
      if ((!is.null(pFeature)) && (pFeature < 1)) {
        message(
          "Cannot use the pFeatures parameter when specifying a distance matrix as the data object, setting pFeature to 1.\n"
        )
        pFeature <- 1
      }
      if (!is.null(weightsFeature)) {
        message(
          "Cannot use the weightsFeature parameter when specifying a distance matrix as the data object\n"
        )
        weightsFeature <- NULL
      }
    }
    else {
      if (is.null(distance)) {
        message("no specified distance, setting distance to Pearson")
        distance <- "pearson"
      }
    }
    if ((clusterAlg == "km") && inherits(distance,
                                         "character") &&
        (distance != "euclidean")) {
      message(
        "Note: The km (kmeans) option only supports a euclidean distance metric when supplying a data matrix.  If you want to cluster a distance matrix, use a different algorithm such as 'hc' or 'pam'.  Changing distance to euclidean"
      )
      distance <- "euclidean"
    }
    if (inherits(d, "ExpressionSet")) {
      d <- exprs(d)
    }
    ml <-
      ccRun(
        d = d,
        maxK = maxK,
        repCount = reps,
        diss = inherits(d, "dist"),
        pItem = pItem,
        pFeature = pFeature,
        innerLinkage = innerLinkage,
        clusterAlg = clusterAlg,
        weightsFeature = weightsFeature,
        weightsItem = weightsItem,
        distance = distance,
        verbose = verbose,
        corUse = corUse
      )
  }
  res = list()
  colorList = list()
  colorM = rbind()
  thisPal <- c(
    "#A6CEE3",
    "#1F78B4",
    "#B2DF8A",
    "#33A02C",
    "#FB9A99",
    "#E31A1C",
    "#FDBF6F",
    "#FF7F00",
    "#CAB2D6",
    "#6A3D9A",
    "#FFFF99",
    "#B15928",
    "#bd18ea",
    "#2ef4ca",
    "#f4cced",
    "#f4cc03",
    "#05188a",
    "#e5a25a",
    "#06f106",
    "#85848f",
    "#000000",
    "#076f25",
    "#93cd7f",
    "#4d0776",
    "#ffffff"
  )
  colBreaks = NA
  if(is.null(tmyPal) == TRUE) {
    colBreaks = 10
    tmyPal = myPal(colBreaks)
  } else {
    colBreaks = length(tmyPal)
  }
  sc = cbind(seq(0, 1, by = 1 / (colBreaks)))
  rownames(sc) = sc[, 1]
  sc = cbind(sc, sc)
  for (tk in 2:maxK) {
    fm = ml[[tk]]
    hc = hclust(as.dist(1 - fm), method = finalLinkage)
    ct = cutree(hc, tk)
    names(ct) = colnames(d)
    if (any(class(d) == "dist")) {
      names(ct) = colnames(as.matrix(d))
    }
    c = fm
    colorList = setClusterColors(res[[tk - 1]][[3]], ct, thisPal, colorList)
    pc = c
    pc = pc[hc$order,]
    pc = rbind(pc, 0)
    colorM = rbind(colorM, colorList[[1]])
    res[[tk]] = list(
      consensusMatrix = c,
      consensusTree = hc,
      consensusClass = ct,
      ml = ml[[tk]],
      clrs = colorList,
      colorM = colorM,
      pc = pc,
      tmyPal = tmyPal,
      sc=sc
    )
  }
  return(res)
}

clusterun<-function(clustermethod, data,maxK,pItem,clusterAlg,seed,reps,
                    pFeature, distance,corUse,
                    iters,
                    ref_method,repsref,repsreal,
                    pacx1,pacx2,objective,
                    # silent,
                    method,
                    lambdadefault,
                    tunelambda
                    # ,lseq
){
  if(clustermethod=="Monti consensus clustering"){
    
    if(mode(data)=='list'){
      data<-as.matrix(data)
    }
    res<-ccp(
      d = data,
      maxK = maxK,
      reps = reps,
      pItem = pItem,
      pFeature = pFeature,
      distance = distance,
      clusterAlg = clusterAlg,
      ml = NULL,
      weightsFeature = NULL,
      corUse = corUse,
      verbose = F,
      weightsItem = NULL,
      innerLinkage = "average",
      finalLinkage = "average",
      tmyPal = NULL,
      seed = NULL
    )
    ml<-list()
    maxK = 4
    for(i in 2:maxK){
      ml[[i]]<-res[[i]]$ml
    }
    Kvec = 2:maxK
    x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
    PAC = rep(NA,length(Kvec))
    names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
    for(i in Kvec){
      M = res[[i]]$consensusMatrix
      Fn = ecdf(M[lower.tri(M)])
      PAC[i-1] = Fn(x2) - Fn(x1)
    }#end for i
    # The optimal K
    optK = Kvec[which.min(PAC)]
    res<-res[[optK]]
    res<-list(res=res,ml=ml)
  }else {
    require(M3C)
    print(1)
    print(clustermethod)
    res <- M3C(
      mydata=as.data.frame(data),
      ref_method=ref_method,
      repsref=reps,
      iters =iters,
      objective = objective,
      repsreal=repsreal,
      clusteralg=clusterAlg,
      pacx1=pacx1,
      pacx2=pacx2,
      method=method,
      lambdadefault=lambdadefault,
      tunelambda=tunelambda,
      lseq=seq(0.02, 0.1, by =0.02),
      silent=FALSE,
      cores=4,
      des = NULL,
      removeplots = TRUE,
      fsize = 8,
      lthick = 1,
      dotsize = 1.25
    )
  }
  return(res)
}  


silhouette_SimilarityMatrix<-function (group, similarity_matrix){
  similarity_matrix = as.matrix(similarity_matrix)
  similarity_matrix <- (similarity_matrix + t(similarity_matrix))/2
  diag(similarity_matrix) = 0
  normalize <- function(X) X/rowSums(X)
  similarity_matrix <- normalize(similarity_matrix)
  n <- length(group)
  if (!all(group == round(group))) 
    stop("'group' must only have integer codes")
  cluster_id <- sort(unique(group <- as.integer(group)))
  k <- length(cluster_id)
  if (k <= 1 || k >= n) 
    return(NA)
  doRecode <- (any(cluster_id < 1) || any(cluster_id > k))
  if (doRecode) 
    group <- as.integer(fgroup <- factor(group))
  cluster_id <- sort(unique(group))
  wds <- matrix(NA, n, 3, dimnames = list(names(group), c("cluster", 
                                                          "neighbor", "sil_width")))
  for (j in 1:k) {
    index <- (group == cluster_id[j])
    Nj <- sum(index)
    wds[index, "cluster"] <- cluster_id[j]
    dindex <- rbind(apply(similarity_matrix[!index, index, 
                                            drop = FALSE], 2, function(r) tapply(r, group[!index], 
                                                                                 mean)))
    maxC <- apply(dindex, 2, which.max)
    wds[index, "neighbor"] <- cluster_id[-j][maxC]
    s.i <- if (Nj > 1) {
      a.i <- colSums(similarity_matrix[index, index])/(Nj - 
                                                         1)
      b.i <- dindex[cbind(maxC, seq(along = maxC))]
      ifelse(a.i != b.i, (a.i - b.i)/pmax(b.i, a.i), 0)
    }
    else 0
    wds[index, "sil_width"] <- s.i
  }
  attr(wds, "Ordered") <- FALSE
  class(wds) <- "silhouette"
  wds
}



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

TCGAanalyze_Stemness<-function (stemSig, dataGE, annotation = FALSE) 
{
  reads <- dataGE
  X <- reads
  w <- stemSig
  commonStemsigGenes <- intersect(names(w), rownames(X))
  X <- X[commonStemsigGenes, ]
  w <- w[rownames(X)]
  s <- apply(X, 2, function(z) {
    cor(z, w, method = "sp", use = "complete.obs")
  })
  s <- s - min(s)
  s <- s/max(s)
  dataSce_stemness <- cbind(s)
  dataAnnotationSC <- matrix(0, ncol(reads), 2)
  colnames(dataAnnotationSC) <- c("Sample", "Annotation")
  dataAnnotationSC <- as.data.frame(dataAnnotationSC)
  dataAnnotationSC$Sample <- colnames(reads)
  rownames(dataAnnotationSC) <- colnames(reads)
  dataAnnotationSC <- cbind(dataAnnotationSC, StemnessScore = rep(0, 
                                                                  nrow(dataAnnotationSC)))
  dataAnnotationSC[rownames(dataSce_stemness), "StemnessScore"] <- as.numeric(dataSce_stemness)
  colnames(dataAnnotationSC)[1] <- "Sample"
  if (annotation == "sampleType") {
    sampleTP <- TCGAquery_SampleTypes(barcode = dataAnnotationSC$Sample, 
                                      typesample = "TP")
    sampleNT <- TCGAquery_SampleTypes(barcode = dataAnnotationSC$Sample, 
                                      typesample = "NT")
    dataAnnotationSC[sampleTP, "Annotation"] <- "TP"
    dataAnnotationSC[sampleNT, "Annotation"] <- "NT"
  }
  if (annotation == "subtype") {
    dataSubt <- TCGAquery_subtype(tumor = "BRCA")
    for (i in 1:nrow(dataAnnotationSC)) {
      curSample <- dataAnnotationSC$Sample[i]
      dataAnnotationSC$Annotation <- dataSubt[dataSubt$patient %in% 
                                                substr(curSample, 1, 12), "BRCA_Subtype_PAM50"]
    }
  }
  return(dataAnnotationSC)
}

stem<-function(exp,gene,method){
  #require(TCGAbiolinks)
  exp<-as.data.frame(exp)
  stemness <- TCGAanalyze_Stemness(stemSig = PCBC_stemSig, dataGE = exp)
  
  stemness$gene<-as.numeric(exp[gene,])
  stemness<-subset(stemness,select=c(StemnessScore,gene))
  
  scatter<-ggscatter(data=stemness, y = "StemnessScore",x = "gene",
                     size = 1,
                     add = "reg.line", conf.int = TRUE,
                     add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                     cor.coef = TRUE, cor.method = method,cor.coef.size=4.5,
                     cor.coeff.args = list(label.y=max(stemness[,1])+se(stemness[,1]),label.sep = ","),
                     xlab =paste("Normalized expression of",gene,sep=" ") , ylab = "Stemness score")
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
  require(rms)
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
  require(rms)
  require(survival)
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
  require(rms)
  require(survival)
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
  
  ggplotify::as.ggplot(function() hist(
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
  require(clusterProfiler)
  require(ReactomePA)
  require(msigdbr)
  require(enrichplot)
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
  require(clusterProfiler)
  require(ReactomePA)
  require(msigdbr)
  require(enrichplot)
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
    require(clusterProfiler)
    require(ReactomePA)
    require(msigdbr)
    require(enrichplot)

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
    require(clusterProfiler)
    require(ReactomePA)
    require(msigdbr)
    require(enrichplot)
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
    require(clusterProfiler)
    require(ReactomePA)
    require(msigdbr)
    require(enrichplot)
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
  require(WGCNA)
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
    require(WGCNA)
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
  require(WGCNA)
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
  require(WGCNA)
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

sample.split<-function (Y, SplitRatio = 2/3, group = NULL) 
{
  nSamp = length(Y)
  nGroup = length(group)
  if (nGroup > 0 && nGroup != nSamp) 
    stop("Error in sample.split: Vectors 'Y' and 'group' have to have the same length")
  BinOne = logical(nSamp)
  SplitRatio = abs(SplitRatio)
  if (SplitRatio >= nSamp) 
    stop("Error in sample.split: 'SplitRatio' parameter has to be i [0, 1] range or [1, length(Y)] range")
  U = unique(Y)
  nU = length(U)
  if (2 * nU > nSamp | nU == 1) {
    n = if (SplitRatio >= 1) 
      SplitRatio
    else SplitRatio * nSamp
    rnd = runif(nSamp)
    if (nGroup) 
      split(rnd, group) <- lapply(split(rnd, group), mean)
    ord = order(rnd)
    BinOne[ord[1:n]] = TRUE
  }
  else {
    rat = if (SplitRatio >= 1) 
      SplitRatio/nSamp
    else SplitRatio
    for (iU in 1:nU) {
      idx = which(Y == U[iU])
      n = round(length(idx) * rat)
      rnd = runif(length(idx))
      if (nGroup) {
        grp = group[idx]
        split(rnd, grp) <- lapply(split(rnd, grp), mean)
      }
      ord = order(rnd)
      BinOne[idx[ord[1:n]]] = TRUE
    }
  }
  if (SplitRatio >= 1) {
    n = sum(BinOne) - SplitRatio
    if (n > 0) 
      BinOne[sample(which(BinOne), n)] = FALSE
    else if (n < 0) 
      BinOne[sample(which(!BinOne), -n)] = TRUE
  }
  return(BinOne)
}
#######################################################################################################################
#######################################################################################################################
nestml<-function(expres,clinical,split,endpoint,ratio,ctrl,inner,outer,lnid,seed){
  # require(caret)
  require(mlr)
  require(randomForestSRC)
  require(survival)
  expres$time<-clinical[,endpoint[1]]
  expres$status<-clinical[,endpoint[2]]
  if(split==T){
    set.seed(seed)
    # index <- createDataPartition(expres$status, p = ratio, list = FALSE)
    index = sample.split(expres$status, SplitRatio=ratio)
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
  # require(caret)
  require(mlr)
  require(randomForestSRC)
  require(survival)

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
    # index <- createDataPartition(expres$status, p = ratio, list = FALSE)
    index = sample.split(expres$status, SplitRatio=ratio)
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


# survROC<-function(data,bmtime,bmstatus,marker,method,predyear,cutpoint){
#   require(survivalROC)
# 
#   data[complete.cases(data[,bmtime]) & data[,bmtime]>0,]
#   if(max(data[,bmtime])>600){
#     by<- 365
#   } else{
#     by<-12
#   }
#   years<-seq(from=0,to=max(data[,bmtime]),by=by)
#   years<-years[-1]
#   yearlabel<-c()
#   for(i in 1:length(years)){
#     yearlabel[i]<- paste(i,"year",sep="-")
#   }
#   sumROC<-list()
#   for (i in 1:length(years)){
#     sumROC[[i]] <- survivalROC(Stime = data[,bmtime],status = data[,bmstatus],marker = data[[marker]],
#                                predict.time =years[i],method = method,span = 0.25*nrow(data)^(-0.20))
#   }
#   sumAUC<-list()
#   for (i in 1:length(sumROC)){
#     sumAUC[[i]]<-sumROC[[i]]$AUC
#   }
#   ROCdata<-c()
#   for(i in 1:length(sumROC)){
#     predict.time<-sumROC[[i]]$predict.time
#     TP<-sumROC[[i]]$TP
#     FP<-sumROC[[i]]$FP
#     auc<-sumROC[[i]]$AUC
#     tmp<-c(predict.time,TP,FP,auc)
#     ROCdata<-rbind(ROCdata,tmp)
#   }
#   survivalROC_helper <- function(t) {
#     survivalROC(Stime = data[,bmtime],status = data[,bmstatus],marker = data[[marker]],
#                 predict.time =t,method = method,span = 0.25*nrow(data)^(-0.20))
#   }
# 
#   yearid<-as.numeric(substr(predyear,1,1))
#   time<-years[yearid]
#   survivalROC_data <- data_frame(t = time) %>%
#     mutate(survivalROC = map(t, survivalROC_helper),
# 
#            auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
# 
#            df_survivalROC = map(survivalROC, function(obj) {
#              as_data_frame(obj[c("cut.values","TP","FP")])
#            })) %>%
#     dplyr::select(-survivalROC) %>%
#     unnest() %>%
#     arrange(t, FP, TP)
#   survivalROC_data1<-mutate(survivalROC_data,auc = sprintf("%.3f",auc))
#   survivalROC_data1$years<-survivalROC_data1$t/by
#   survivalROC_data1<-unite(survivalROC_data1,year, years,auc,sep = " year AUC: " )
#   AUC =factor(survivalROC_data1$year)
#   sumAUC1<-list()
#   for(i in 1:length(yearid)){
#     sumAUC1[[i]]<-sumAUC[[yearid[i]]]
#   }
#   sumROC1<-list()
#   for(i in 1:length(yearid)){
#     sumROC1[[i]]<-sumROC[[yearid[i]]]
#   }
#   ROC.1<-sumROC1[[which.max(sumAUC1)]]
# 
#   dot <- data.frame(TP = ROC.1$TP[with(ROC.1, which.min((1-TP)^2+ FP^2))],
#                     FP = ROC.1$FP[with(ROC.1, which.min((1-TP)^2+ FP^2))])
#   dot <- rbind(c(1,0),dot)
#   if(cutpoint==T){
#     cutoff.imm <- ROC.1$cut.values[with(ROC.1, which.min((1-TP)^2+ FP^2))]
#     ROC.plot<-ggplot(survivalROC_data1,mapping = aes(x = FP, y = TP)) +
#       geom_path(aes(color= AUC))+
#       geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
#       theme_bw() +scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+xlab("False positive rate")+
#       ylab("True positive rate") +
#       theme(legend.position = c(0.7,0.2))+
#       geom_path(mapping = aes(x = FP,y = TP),data = dot)+
#       annotate("text",x = dot$FP[2] + 0.05,y = dot$TP[2],label = paste0("Cutoff: ",round(cutoff.imm,3)))
#   } else{
#     ROC.plot<-ggplot(survivalROC_data1,mapping = aes(x = FP, y = TP)) +
#       geom_path(aes(color= AUC))+
#       geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
#       theme_bw() +scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+xlab("False positive rate")+
#       ylab("True positive rate") +
#       theme(legend.position = c(0.7,0.2))
#   }
#   return(ROC.plot)
# }

survROC<-function(data,bmtime,bmstatus,marker,method,predyear,cutpoint){
  require(survivalROC)
  require(survival)
  require(DT)
  require(tidyverse)
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















