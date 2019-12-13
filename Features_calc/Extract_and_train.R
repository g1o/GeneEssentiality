#!/bin/R
library(data.table)
library(seqinr)
library(parallel)
library(doMC)
library(pROC)
library(caret)
library(protr)
#HMMSCAN and PFAM are needed in Calc_feats
LAMBDA=50; #If sequence length is less than LAMBDA, then it will be skipped; 
OMEGA=0.05;
CPU=20;
trees=750;#number of trees in rf
cross_validation=10;

script.dir<-dirname(sys.frame(1)$ofile)
FEATURES_func<-paste0(script.dir,"/functions/Serial_gene_features_extraction.R")
source(FEATURES_func,echo=T)

PFAM_path="/home/programs/DATABASES/PFAM/Pfam-A.hmm"; #change the path to your hmmpressed pfam database
ESSENTIAL_GENES_PATH<-("../data/drosophila/GE.fasta") #CDS
NONESSENTIAL_GENES_PATH<-("../data/drosophila/GNE.fasta") #CDS

GE<-read.fasta(ESSENTIAL_GENES_PATH)
GNE<-read.fasta(NONESSENTIAL_GENES_PATH)
GNE["FBgn0013675"]<-NULL # removing as it is mitochondrial and it has another genetic code, it is the only mito gene here. 

#=============CALCULATE FEATURES=============
start.time <- Sys.time() #timing code
cl<-makeCluster(CPU,type="FORK")

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

stopCluster(cl)
end.time <- Sys.time() #timing feature extraction
time.taken <- end.time - start.time
time.taken
#-------------------------------------------------

Complete_set<-rbindlist(list(Features_essential,Features_notessential),fill=T,idcol=T)
Complete_set<-as.data.frame(Complete_set)
rownames(Complete_set)<-t(Complete_set[,2])
Complete_set<-Complete_set[,-c(1,2)]
Complete_set[is.na(Complete_set)] <- F

train<-Complete_set;
mtry<-round(sqrt(length(train)))
tunegrid  <- expand.grid(.mtry = c(mtry,mtry*2) )

#==========TRAIN MODELS==========================
modellist <- list()
predictionlist <- list()
registerDoMC(CPU) #--- needed for random forest parallelization
for (ntree in trees){
for (CV in c(cross_validation)) {
	control  <- trainControl(method="repeatedcv", number=CV, repeats=5,classProbs = TRUE,summaryFunction=twoClassSummary,savePredictions = TRUE)

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

namet<-format(Sys.time(),"MODELS.%d_%b_%Y_%X.RData")
save(modellist,file=namet)
				     
anot<-paste0('\n',paste0(label[1,1:2],collapse=" AUC="),'\n',paste0(label[2,1:2],collapse=" AUC="),'\n',paste0(label[3,1:2],collapse=" AUC="),'\n',paste0(label[4,1:2],collapse=" AUC="),'\n')				     

graphic<-ggplot(modellist[[1]]$pred,aes(m = NE, d = obs, color=factor(mtry))) + geom_roc(n.cut=0) + style_roc(theme=theme_classic,ylab="Sensibilidade (Taxa de verdadeiros positivos)",xlab="Taxa de falsos positivos (Especificidade - 1)") + annotate(geom="text",x = 0.7, y = 0.2, label=c(anot))+ geom_abline(intercept=0, slope=1, colour="orange") 
pdf("ROC.pdf");
graphic;dev.off()
