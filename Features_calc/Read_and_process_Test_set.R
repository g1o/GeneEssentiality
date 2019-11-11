#!/bin/env R
library(data.table) #count function conflicts with the one in seqinR. This library must be loaded before seqinR
library(seqinr)
library(parallel)
library(doMC)
library(pROC)
library(caret)
library(protr)
script.dir<-dirname(sys.frame(1)$ofile) 
FEATURES_func<-paste0(script.dir,"/functions/Serial_gene_features_extraction.R") #when sourcing, get path to function relative to this script
source(FEATURES_func,echo=T)

CPU=30
load("nicemodels_ntree250_07AUC.RData") # LOAD TRAINED MODELS list # need to be changed
PFAM_path="/mnt/Databases/PFAM/Pfam-A.hmm"; #change the path to your hmmpressed pfam database
GE<-read.fasta("../data/drosophila/GE.fasta")
GNE<-read.fasta("../data/drosophila/GNE.fasta")

LAMBDA=length(grep( "Xc2.lambda" , colnames( modellist[[1]]$trainingData ))); #set like the training data in model
OMEGA=0.05;#set like the training data in model

cl<-makeCluster(CPU,type="FORK")
#-------------------------------------------------
Features_essential<-as.data.frame(rbindlist(
                                  parSapply(cl,GE,function(gene)
                                            Calc_feats(gene,PFAM_PATH=PFAM_path,LAMBDA=LAMBDA,OMEGA=OMEGA)),
                    fill=T,idcol=T))

rownames(Features_essential)<-t(Features_essential[,1])
Features_essential[is.na(Features_essential)] <- F #binnary NA to False
Features_essential<- data.frame(Features_essential,Class="E")

Features_notessential<-as.data.frame(rbindlist(
                                     parSapply(cl,GNE,function(gene)
                                               Calc_feats(gene,PFAM_PATH=PFAM_path,LAMBDA=LAMBDA,OMEGA=OMEGA)),
                       fill=T,idcol=T))

rownames(Features_notessential)<-t(Features_notessential[,1])
Features_notessential[is.na(Features_notessential)] <- F
Features_notessential<-data.frame(Features_notessential,Class="NE")

Complete_set<-rbindlist(list(Features_essential,Features_notessential),fill=T,idcol=T)
Complete_set<-as.data.frame(Complete_set)
rownames(Complete_set)<-t(Complete_set[,2])
Complete_set<-Complete_set[,-c(1,2)]
Complete_set[is.na(Complete_set)] <- F

#-------------------------------------------------
stopCluster(cl)
ntr<-setdiff(names(modellist[[1]]$trainingData),names(Complete_set))
ntr<-ntr[which(ntr!=".outcome")]
nte<-setdiff(names(Complete_set),names(modellist[[1]]$trainingData))
nte<-nte[which(nte!="Class")]

Complete_set[,nte]<-NULL
Complete_set[,ntr]<-FALSE
Complete_set[is.na(Complete_set)] <- F

res<-predict(modellist[[1]],Complete_set,type="prob");
result.roc.crispr<-roc(Complete_set$Class,res$E)

