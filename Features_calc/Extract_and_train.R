#!/bin/R
library(plyr) #count function conflicts with the one in seqinR. This library must be loaded before seqinR
library(seqinr)
library(parallel)
library(doMC)
library(pROC)
library(caret)
library(protr)
LAMBDA=50; #If sequence length is less than LAMBDA, then it will be skipped; 
OMEGA=0.05;
CPU=50

source("functions/Serial_gene_features_calcs.R",echo=T)
ESSENTIAL_GENES_PATH<-("GE_train-ZeroCrispr.fasta")
NONESSENTIAL_GENES_PATH<-("GNE_train-ZeroCrispr.fasta")

GE<-read.fasta(ESSENTIAL_GENES_PATH)
GNE<-read.fasta(NONESSENTIAL_GENES_PATH)

GNE["FBgn0013675"]<-NULL # removing as it is mitochondrial and it has another genetic code, it is the only mito gene here. 
start.time <- Sys.time() #timing code

cl<-makeCluster(CPU,type="FORK")

Features_essential<-rbind.fill(parSapply(cl,GE,function(gene) Calc_feats(gene,LAMBDA,OMEGA)))
rownames(Features_essential)<-getName(GE)
Features_essential[is.na(Features_essential)] <- F #binnary NA to False
Features_essential<- data.frame(Features_essential,Class="E")

Features_notessential<-rbind.fill(parSapply(cl,GNE,function(gene) Calc_feats(gene,LAMBDA,OMEGA)))
rownames(Features_notessential)<-getName(GNE)
Features_notessential[is.na(Features_notessential)] <- F
Features_notessential<-data.frame(Features_notessential,Class="NE")

stopCluster(cl)

#-------------------------------------------------

end.time <- Sys.time() #timing feature extraction
time.taken <- end.time - start.time
time.taken

Complete_set<-rbind.fill(Features_essential,Features_notessential)
Complete_set[is.na(Complete_set)] <- F

train<-Complete_set[-sampleindex,];
train<-Complete_set;
tunegrid  <- expand.grid(.mtry = round(sqrt(length(train))) )

modellist <- list()
predictionlist <- list()
registerDoMC(CPU) #--- needed for random forest parallel
for (ntree in c(500)){
for (CV in c(10)) {
	control  <- trainControl(method="repeatedcv", number=CV, repeats=5,classProbs = TRUE,summaryFunction=twoClassSummary)

	#TRAIN By maximizing the ROC METRIC
	fit<-train(Class ~ .,data=train,metric="ROC",method="rf",trControl=control,tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=ntree,importance=T) 
	key <- paste(toString(CV),"cv-metric_ROC-ntree",toString(ntree)) 
	modellist[[key]] <- fit              

	#TRAIN By maximizing Specificity 
	fit<-train(Class ~ .,data=train,metric="Spec",method="rf",trControl=control,tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=ntree,importance=T)
	key <- paste(toString(CV),"cv-metric_Spec-ntree",toString(ntree))
	modellist[[key]] <- fit                    

}
}

namet<-format(Sys.time(),"Drosophila_MODELS.%d_%H-%M-%S.RData")
save(modellist,file=namet)
source("Read_and_process_Test_set.R")

res<-predict(modellist[[1]],Crispr_set,type="prob");result.roc.crispr<-roc(Crispr_set$Class,res$E)
res<-predict(modellist[[1]],Tribolium_set,type="prob");result.roc.tribolium<-roc(Tribolium_set$Class,res$E)
