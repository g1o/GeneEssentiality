#!/bin/R
library(seqinr)
library(parallel)
library(doMC)
library(pROC)
library(caret)
library(protr)
CPU=60

start_time <- Sys.time()

EXTRACT_PAAC<- function (x, props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"),
    lambda = 30, w = 0.05, customprops = NULL)
{
    if (protcheck(x) == FALSE)
        stop("x has unrecognized amino acid type")
    if (length(x) <= lambda)  #edited to use with SeqFastaAA class from seqinR
        stop("Length of the protein sequence must be greater than \"lambda\"")
    AAidx = read.csv(system.file("sysdata/AAidx.csv", package = "protr"),
        header = TRUE)
    tmp = data.frame(AccNo = c("Hydrophobicity", "Hydrophilicity",
        "SideChainMass"), A = c(0.62, -0.5, 15), R = c(-2.53,
        3, 101), N = c(-0.78, 0.2, 58), D = c(-0.9, 3, 59), C = c(0.29,
        -1, 47), E = c(-0.74, 3, 73), Q = c(-0.85, 0.2, 72),
        G = c(0.48, 0, 1), H = c(-0.4, -0.5, 82), I = c(1.38,
            -1.8, 57), L = c(1.06, -1.8, 57), K = c(-1.5, 3,
            73), M = c(0.64, -1.3, 75), F = c(1.19, -2.5, 91),
        P = c(0.12, 0, 42), S = c(-0.18, 0.3, 31), T = c(-0.05,
            -0.4, 45), W = c(0.81, -3.4, 130), Y = c(0.26, -2.3,
            107), V = c(1.08, -1.5, 43))
    AAidx = rbind(AAidx, tmp)
    if (!is.null(customprops))
        AAidx = rbind(AAidx, customprops)
    aaidx = AAidx[, -1]
    row.names(aaidx) = AAidx[, 1]
    n = length(props)
    H0 = as.matrix(aaidx[props, ])
    H = matrix(ncol = 20, nrow = n)
    for (i in 1:n) H[i, ] = (H0[i, ] - mean(H0[i, ]))/(sqrt(sum((H0[i,
        ] - mean(H0[i, ]))^2)/20))
    AADict = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
    dimnames(H) = list(props, AADict)
    Theta = vector("list", lambda)
    xSplitted = x #edited to use with SeqFastaAA class from seqinR
    N = length(xSplitted)
    for (i in 1:lambda) {
        for (j in 1:(N - i)) {
            Theta[[i]][j] = mean((H[, xSplitted[j]] - H[, xSplitted[j +
                i]])^2)
        }
    }
    theta = sapply(Theta, mean)
    fc = summary(factor(xSplitted, levels = AADict), maxsum = 21)
    Xc1 = fc/(1 + (w * sum(theta)))
    names(Xc1) = paste("Xc1.", names(Xc1), sep = "")
    Xc2 = (w * theta)/(1 + (w * sum(theta)))
    names(Xc2) = paste("Xc2.lambda.", 1:lambda, sep = "")
    Xc = c(Xc1, Xc2)
    Xc
}

#LAMBDA=22;
#for (OMEGA in seq(0.01,1,by=0.04)){
Calc_feats<-function(seqs,longest_orf){
        cl <- makeCluster(5,type="FORK") # limit to 5, over it it decreases speed
        #----------------parallel
        H2<-parSapply(cl,seqs,function(x) sum(-count(x,2,freq=T)*log2(count(x,2,freq=T)),na.rm=T))
        H3<-parSapply(cl,seqs,function(x) sum(-count(x,3,freq=T)*log2(count(x,3,freq=T)),na.rm=T))
        f1<-parSapply(cl,seqs, function(x) count(x,1,freq=T) )
        f2<-parSapply(cl,seqs, function(x) count(x,2,freq=T) )
        f3<-parSapply(cl,seqs, function(x) count(x,3,freq=T) )

        MI<-parSapply(cl,names(f1[1,]), function(g) sapply (names(f1[,1]) , function(x) sapply (names(f1[,1]) , function(y) f2[paste(x,y,sep=""),g]*log2(f2[paste(x,y,sep=""),g]/(f1[x,g]*f1[y,g])) ) ) )
        row.names(MI)<-names(f2[,1])
        SomaMI<-apply(MI,2,sum,na.rm=T)
        CMI<-parSapply (cl,names(f1[1,]), function(g) sapply (names(f1[,1]) ,function(x) sapply (names(f1[,1]) , function(y) sapply (names(f1[,1]) , function(z)  f3[paste(x,z,y,sep=""),g]*log2(f1[z,g]*f3[paste(x,z,y,sep=""),g]/(f2[paste(x,z,sep=""),g]*f2[paste(z,y,sep=""),g]))) ) ) )
        row.names(CMI)<-names(f3[,1])
        SomaCMI<-apply(CMI,2,sum,na.rm=T)

        f1<-parSapply(cl,longest_orf, function(x) count(x,1,freq=T,alphabet=s2c("ACDEFGHIKLMNPQRSTVWY")))
        f2<-parSapply(cl,longest_orf, function(x) count(x,2,freq=T,alphabet=s2c("ACDEFGHIKLMNPQRSTVWY")))
        #f2<-f2[-(1:21),] #remove dimers que começam com stop, pois nao tem numa orf
        protMI <-parSapply(cl,names(f1[1,]), function(g) sapply (names(f2[,1]) , function(x) f2[x,g]*log2(f2[x,g]/(f1[strsplit(x,"")[[1]][1],g]*f1[strsplit(x,"")[[1]][2],g])) ) )
        propPEP<-parLapply(cl,longest_orf,function(x) AAstat(x,plot=F))
        Soma_protMI<-apply(protMI,2,sum,na.rm=T)
	pepFeatures<-sapply(propPEP,function (x) {
                                    pep<-c(as.numeric(x[[2]]),as.numeric(x[3]))
                                    names(pep)<-c(names(x[[2]]),names(x[3]))
                                    pep} )
	aalength<-sapply(longest_orf,function(x) length(x))
	PseAA<-parSapply(cl,longest_orf, function(x) {EXTRACT_PAAC(x,lambda=LAMBDA,w=OMEGA)}) #Chou's pseudoAA #22 min for AB
        Features<-rbind(SomaMI,MI,SomaCMI,CMI,H2,H3,pepFeatures,protMI,Soma_protMI,PseAA,aalength)
        Features<-data.frame(t(Features))
stopCluster(cl)
        return(Features)
}
seqs<-read.fasta("Acinetobacter_baylyi_ADP1-EG.fasta")
longest_orf<-read.fasta("Acinetobacter_baylyi_ADP1-EG.faa",seqtype="AA")
Features_essential<-Calc_feats(seqs,longest_orf)
Features_essential[is.na(Features_essential)] <- 0
#Features_essential<-cbind(Features_essential,read.table("signalP/EG.tsv",skip=1,row.names=1)[1:nrow(Features_essential),8:9])
Features_essential<-cbind(Features_essential,Class="E")

seqs<-read.fasta("back/Acinetobacter_baylyi_ADP1-NEG.fasta")
longest_orf<-read.fasta("back/Acinetobacter_baylyi_ADP1-NEG.faa",seqtype="AA")
Features_notessential<-Calc_feats(seqs,longest_orf)
Features_notessential[is.na(Features_notessential)] <- 0
#Features_notessential<-cbind(Features_notessential,read.table("signalP/NEG.tsv",skip=1,row.names=1)[1:nrow(Features_notessential),8:9])
Features_notessential<-cbind(Features_notessential,Class="NE")


seqs<-read.fasta("AB_notFINE.fasta")
longest_orf<-read.fasta("AB_notFINE.faa",seqtype="AA")
Feat_negCDS<-Calc_feats(seqs,longest_orf)
Feat_negCDS[is.na(Feat_negCDS)] <- 0
Feat_negCDS<-cbind(Feat_negCDS,Class="NE")

seqs<-read.fasta("back/ABunaltered.fasta")
longest_orf<-read.fasta("back/ABunaltered.faa",seqtype="AA")
Feat_Unaltered_DNEG<-Calc_feats(seqs,longest_orf)
Feat_Unaltered_DNEG[is.na(Feat_Unaltered_DNEG)] <- 0
Feat_Unaltered_DNEG<-cbind(Feat_Unaltered_DNEG,Class="NE")

teste<- rbind(Features_essential[61:120,],Features_notessential[61:600,]);
negteste<-rbind(Features_essential[61:120,],Feat_negCDS[61:600,]);
unaltered_test<-rbind(Features_essential[61:120,],Feat_Unaltered_DNEG[61:600,]);

train<-rbind(Features_essential[c(1:60,121:421),],Features_notessential[c(1:60,601:901),])
NEG_train<-rbind(Features_essential[c(1:60,121:421),],Feat_negCDS[c(1:60,601:901),])
unaltered_train<-rbind(Features_essential[c(1:60,121:421),],Feat_Unaltered_DNEG[c(1:60,601:901),])

end_time <- Sys.time();
end_time - start_time

#set.seed(7)
tunegrid <- expand.grid(.mtry=c(1:11,21,31,41,51,61,71,81,91,101))
control  <- trainControl(method="repeatedcv", number=10, repeats=5,classProbs = TRUE,summaryFunction=twoClassSummary)
#control<-trainControl(method="cv", number=5)

#--- needed for random forest parallel
registerDoMC(CPU)
start_time <- Sys.time()
ntree=1000
modellist <- list()
#for (ntree in c(2000))
	fit<-train(Class ~ .,data=train,metric="ROC",method="rf",trControl=control,tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=ntree,importance=T) 
	fit_NEG<-train(Class ~ .,data=NEG_train,metric="ROC",method="rf",trControl=control,tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=ntree,importance=T) 
	fit_unaltered<-train(Class ~ .,data=unaltered_train,metric="ROC",method="rf",trControl=control,tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=ntree,importance=T) 
#	key <- toString(OMEGA)
#	modellist[[key]] <- fit
	modellist[[1]] <- fit
	modellist[[2]] <- fit_NEG
	modellist[[3]] <- fit_unaltered
#}
#results <- resamples(modellist)
#summary(results)
#dotplot(results)
#dev.off()
end_time <- Sys.time()
end_time - start_time


#layout(matrix(c(1,2), ncol=2, byrow=TRUE), heights=c(4, 1))

pdf(file="DEG-CDS_TRAIN_TEST.pdf")
for (model in modellist){
	if(model$call$data == "train"){
		mainn<-"A: Treinado com CDSs verdadeiras"
	}
	if(model$call$data == "NEG_train"){
		mainn<-"B: treinado com CDSs falsas"
	}else{
		mainn<-"C: treinado com CDSs não filtradas"
	}
	res<-predict(model,teste,type="prob");	result.roc<-roc(teste$Class,res$E)
	#plot(result.roc,print.auc=T,main = "RandomForest(rf) ROC: True CDS ")
	plot(result.roc,print.auc=T,axes=F,ylab='',xlab='',main = mainn,col="blue")
	axis(1, pos=0);	axis(2, pos=1);	title(ylab="Sensitividade",xlab="Especificidade", mgp=c(2,3,2),cex.lab=1.2)
        
	res<-predict(model,negteste,type="prob")
        result.roc<-roc(negteste$Class,res$E)
        plot(result.roc,print.auc=T,axes=F,ylab='',xlab='',col="red",print.auc.y=0.95,print.auc.x=0.65)

        res<-predict(model,unaltered_test,type="prob")
        result.roc<-roc(unaltered_test$Class,res$E)        
	plot(result.roc,print.auc=T,axes=F,ylab='',xlab='',col="green",print.auc.y=0.85,print.auc.x=0.55)

	par(new=TRUE)
	legend(0.5,0.3,title="Parameters",legend=c(paste("ntree: ",model$dots$ntree),paste("mtry: ",model$finalModel$mtry)));	dim(model$resample)[1]
}#	legend(0.93, 0.9, legend=c("Teste com CDS verdadeiras","Teste com CDS falsas"), col=c("blue","red"), lty=1, cex=1,lwd=2)
dev.off()
##############################


#plot(result.roc,print.thres="best", print.thres.best.method="closest.topleft")


