#!/bin/env R
library(plyr) #count function conflicts with the one in seqinR. This library must be loaded before seqinR
library(seqinr)
library(parallel)
library(doMC)
library(pROC)
library(caret)
library(protr)
load("nicemodels_ntree250_07AUC.RData") # LOAD TRAINED MODELS list # need to be changed
ge_cripr<-read.fasta("../data/drosophila/GE.fasta")
gne_cripr<-read.fasta("../data/drosophila/GNE.fasta")

#TcGNE<-read.fasta("Tribolium_NonEssential-less_than0.30lethal_from_ibeetle.fasta")
#TcGE<-read.fasta("Tribolium_top100_RNAiLethal_FROM_Ulrich2015_allCDSs_from_genes.fasta")

LAMBDA=length(grep( "Xc2.lambda" , colnames( modellist[[1]]$trainingData ))); #set like the training data in model
OMEGA=0.05;#set like the training data in model
source("functions/Serial_gene_features_extraction.R",echo=T)

cl<-makeCluster(6,type="FORK")
fcGE <-rbind.fill(parSapply(cl,ge_cripr ,function (x) Calc_feats(x,LAMBDA,OMEGA)))
fcGNE<-rbind.fill(parSapply(cl,gne_cripr,function (x) Calc_feats(x,LAMBDA,OMEGA)))

fcGE<- data.frame(fcGE,Class="E")   #Set when outcome is known
fcGNE<- data.frame(fcGNE,Class="NE")#Set when outcome is known

stopCluster(cl)
Crispr_set<-rbind.fill(list(fcGE,fcGNE))

ntr<-setdiff(names(modellist[[1]]$trainingData),names(Crispr_set))
ntr<-ntr[which(ntr!=".outcome")]
nte<-setdiff(names(Crispr_set),names(modellist[[1]]$trainingData))
nte<-nte[which(nte!="Class")]


Crispr_set[,nte]<-NULL
Crispr_set[,ntr]<-FALSE
Crispr_set[is.na(Crispr_set)] <- F

res<-predict(modellist[[1]],Crispr_set,type="prob");
result.roc.crispr<-roc(Crispr_set$Class,res$E)


#cl<-makeCluster(6,type="FORK")
#fcGE_Tc <-rbind.fill(parSapply(cl,TcGE ,function (x) Calc_feats(x,LAMBDA,OMEGA))) #tribolium Features Calc
#fcGNE_Tc<-rbind.fill(parSapply(cl,TcGNE,function (x) Calc_feats(x,LAMBDA,OMEGA)))
## 
#fcGE_Tc<- data.frame(fcGE_Tc,Class="E")  #Set when outcome is known
#fcGNE_Tc<- data.frame(fcGNE_Tc,Class="NE")#Set when outcome is known
# 
#stopCluster(cl)
#Tribolium_set<-rbind.fill(list(fcGE_Tc,fcGNE_Tc))
# 
#ntr<-setdiff(names(modellist[[1]]$trainingData),names(Tribolium_set))
#ntr<-ntr[which(ntr!=".outcome")]
#nte<-setdiff(names(Tribolium_set),names(modellist[[1]]$trainingData))
#nte<-nte[which(nte!="Class")]
#
# 
#Tribolium_set[,nte]<-NULL
#Tribolium_set[,ntr]<-FALSE
#Tribolium_set[is.na(Tribolium_set)] <- F
#
