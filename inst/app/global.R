  library(shiny)
  library(shinydashboard)
  library(shinyBS)
  # library(ggpubr)
  library(shinyjs)
  #library(CBPhelper)
  
source("funcs.R")
dataset<-read.csv("dataset.csv",header=T,row.names=1)

PCBC_stemSig<-readRDS("PCBC_stemSig.rds")
common_genes<-readRDS("common_genes.rds")
oncopath<-readRDS("oncopath.rds")
metabolism<-readRDS("metabolismpath.rds") 
