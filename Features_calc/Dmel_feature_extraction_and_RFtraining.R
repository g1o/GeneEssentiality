#!/bin/R
library(seqinr)
library(parallel)
library(doMC)
library(pROC)
library(caret)
library(protr)
CPU=50

ESSENTIAL_GENES_PATH<-("example/GE-cds.fasta") #CDS region of a positive gene group
NONESSENTIAL_GENES_PATH<-("example/GNE-cds.fasta") #CDS region of a negative gene group

GE<-read.fasta(ESSENTIAL_GENES_PATH)
GNE<-read.fasta(NONESSENTIAL_GENES_PATH)

start.time <- Sys.time() #timing code

LAMBDA=22; #parameter for pseudo-aa composition, lambda is limited to the minimum peptide length minus 1 
OMEGA=0.05;#parameter for pseudo-aa composition
#for (OMEGA in seq(0.01,1,by=0.04)){  #future loop to optmize parameter
Calc_feats<-function(seqs){
        cl <- makeCluster(3,type="FORK") # limit to 5 threads, more might slow it 

##Translate in 3 frames and get the largest ORF. It is assumed that the gene was extracted in the correct strand. 
        frame0<-parLapply(cl,seqs,function (x) translate(x,frame=0))
        frame1<-parLapply(cl,seqs,function (x) translate(x,frame=1))
        frame2<-parLapply(cl,seqs,function (x) translate(x,frame=2))

        longest_orf<-parSapply(cl,names(frame0),function (x){
                                       start0<-which( frame0[[x]] %in% "M")    #vetor posicao de codon de metionina (possivel inicio)
                                       stop0<-which( frame0[[x]] %in% "*")     #vetor posicao de codon de parada
                                       max0<-as.numeric(sapply(seq_along(start0),function (x){ stop0[stop0>start0[x]][1]-start0[x]})) #vetor de comprimento de orfs (Pfim - Pini) 
                                       start1<-which( frame1[[x]] %in% "M")
                                       stop1<-which( frame1[[x]] %in% "*")
                                       max1<-as.numeric(sapply(seq_along(start1),function (x){ stop1[stop1>start1[x]][1]-start1[x]}))
                                       start2<-which( frame2[[x]] %in% "M")
                                       stop2<-which( frame2[[x]] %in% "*")
                                       max2<-as.numeric(sapply(seq_along(start2),function (x){ stop2[stop2>start2[x]][1]-start2[x]}))
                                       maxi<-max(max0,max1,max2,na.rm=T)
                                       if( max(c(max0,0),na.rm=T) == maxi){
                                               maxi<-which.max(max0)
                                               ret<-frame0[[x]][start0[maxi]:stop0[stop0>start0[maxi]][1]]
                                               ret[-max(NROW(ret))] # REMOVE * (STOP CODON) and return
                                       }else if (max(c(max1,0),na.rm=T) == maxi){
                                               maxi<-which.max(max1)
                                               ret<-frame1[[x]][start1[maxi]:stop1[stop1>start1[maxi]][1]]
                                               ret[-max(NROW(ret))]
                                       }else if (max(c(max2,0),na.rm=T) == maxi){
                                               maxi<-which.max(max2)
                                               ret<-frame2[[x]][start2[maxi]:stop2[stop2>start2[maxi]][1]]
                                               ret[-max(NROW(ret))]
                                       }
        })

        ##Counting words frequencies
        
        f1<-parSapply(cl,seqs, function(x) count(x,1,freq=T) )
        f2<-parSapply(cl,seqs, function(x) count(x,2,freq=T) )
        f3<-parSapply(cl,seqs, function(x) count(x,3,freq=T) )

        ##Shannon Entropy
        
        H2<-apply(f2,2,function(x) sum(-x*log2(x),na.rm=T)) #shannon entropy for word size 2
        H3<-apply(f3,2,function(x) sum(-x*log2(x),na.rm=T)) #shannon entropy for word size 3

        ##Nucleotide Mutual Information
        
        MI<-parSapply(cl,names(f1[1,]), function(g) sapply (names(f1[,1]) , function(x) sapply (names(f1[,1]) , function(y) f2[paste(x,y,sep=""),g]*log2(f2[paste(x,y,sep=""),g]/(f1[x,g]*f1[y,g])) ) ) )
        row.names(MI)<-names(f2[,1])
        SomaMI<-apply(MI,2,sum,na.rm=T)

        ##Nucleotide Conditional Mutual Information
        
        CMI<-parSapply (cl,names(f1[1,]), function(g) sapply (names(f1[,1]) ,function(x) sapply (names(f1[,1]) , function(y) sapply (names(f1[,1]) , function(z)  f3[paste(x,z,y,sep=""),g]*log2(f1[z,g]*f3[paste(x,z,y,sep=""),g]/(f2[paste(x,z,sep=""),g]*f2[paste(z,y,sep=""),g]))) ) ) )
        row.names(CMI)<-names(f3[,1])
        SomaCMI<-apply(CMI,2,sum,na.rm=T)

        ##Counting aminoacid word frequencies
        
        f1<-parSapply(cl,longest_orf, function(x) count(x,1,freq=T,alphabet=s2c("ACDEFGHIKLMNPQRSTVWY")))
        f2<-parSapply(cl,longest_orf, function(x) count(x,2,freq=T,alphabet=s2c("ACDEFGHIKLMNPQRSTVWY"))) #Do not count the STOP "*"

        ##Shannon Entropy for aminoacids
        
        pH2<-apply(f2,2,function(x) sum(-x*log2(x),na.rm=T))

        ##Mutual Information for aminoacids
        
        protMI <-parSapply(cl,names(f1[1,]), function(g) sapply (names(f2[,1]) , function(x) f2[x,g]*log2(f2[x,g]/(f1[strsplit(x,"")[[1]][1],g]*f1[strsplit(x,"")[[1]][2],g])) ) )
        Soma_protMI<-apply(protMI,2,sum,na.rm=T)

        ##Peptide properties
        
        propPEP<-parLapply(cl,longest_orf,function(x) AAstat(x,plot=F))
        pepFeatures<-sapply(propPEP,function (x) {
                                    pep<-c(as.numeric(x[[2]]),as.numeric(x[3]))
                                    names(pep)<-c(names(x[[2]]),names(x[3]))
                                    pep} )
        ##peptide length
        
        aalength<-sapply(longest_orf,function(x) length(x))

        ##Pseudo aminoacid
        
        PseAA<-parSapply(cl,longest_orf, function(x) {EXTRACT_PAAC(x,lambda=LAMBDA,w=OMEGA)}) #Chou's pseudoAA #22 min for AB

        Features<-rbind(SomaMI,MI,SomaCMI,CMI,H2,H3,pH2,pepFeatures,protMI,Soma_protMI,aalength,PseAA)
        Features<-data.frame(t(Features))
stopCluster(cl)
        return(Features)
}

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


Features_essential<-Calc_feats(GE)
Features_essential[is.na(Features_essential)] <- 0
#Features_essential<-cbind(Features_essential,read.table("signalP/EG.tsv",skip=1,row.names=1)[1:nrow(Features_essential),8:9])
#bed<-read.table("/home/giovannimc/projects/essential/drosophila/data/essential.bed",header=F,row.names=4)
#exon<-bed[,9]
# coluna 10 do bed tem numero de exons, -1 pq virou rowname
#exon<-t(exon)
#colnames(exon)<-row.names(bed)
#row.names(exon)<-"Nºexon"
#Features_essential<-cbind(Features_essential,exon)
Features_essential<-cbind(Features_essential,Class="E")

#Features_notessential<-Calc_feats(GNE[sample(length(GNE),600)])
Features_notessential<-Calc_feats(GNE)
Features_notessential[is.na(Features_notessential)] <- 0
#Features_notessential<-cbind(Features_notessential,read.table("signalP/NEG.tsv",skip=1,row.names=1)[1:nrow(Features_notessential),8:9])
#bed<-read.table("/home/giovannimc/projects/essential/drosophila/data/notessential.bed",header=F,row.names=5)
#exon<-bed[,9]
#exon<-t(exon)
#colnames(exon)<-row.names(bed)
#row.names(exon)<-"Nºexon"
#Features_notessential<-cbind(Features_essential,exon)
Features_notessential<-cbind(Features_notessential,Class="NE")
#MUST DEAL WITH EXON NUMBER  FEATURE
#-------------------------------------------------

end.time <- Sys.time() #timing code
time.taken <- end.time - start.time
time.taken

teste<- rbind(Features_essential[1:40,],Features_notessential[1:60,]); # sem overlap com treino
train<- rbind(Features_essential[c(41:nrow(Features_essential)),],Features_notessential[c(61:nrow(Features_notessential)),]) 

tunegrid <- expand.grid(.mtry=c(9,21,61,111,151,181,201))
#control  <- trainControl(method="repeatedcv", number=5, repeats=10,classProbs = TRUE,summaryFunction=twoClassSummary)

NE_count<-dim(train[train$Class=="NE",])[1]
balanced_Trainingset<-rbind(train[1:NE_count,],train[train$Class=="NE",]) 

ntree<-20
modellist <- list()
predictionlist <- list()
registerDoMC(CPU) #--- needed for random forest parallel
for (CV in c(10)) {
	control  <- trainControl(method="repeatedcv", number=CV, repeats=5,classProbs = TRUE,summaryFunction=twoClassSummary)

	#TRAIN By maximizing the ROC METRIC
	fit<-train(Class ~ .,data=train,metric="ROC",method="rf",trControl=control,tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=ntree,importance=T) 
	key <- paste(toString(CV),"cv-metric_ROC-ntree",toString(ntree)) 
	modellist[[key]] <- fit                   
	#Get prediction
	res<-predict(fit,teste,type="prob");      
	result.roc<-roc(teste$Class,res$E)        
	predictionlist[[key]] <- result.roc

	#TRAIN By maximizing Specificity 
	fit<-train(Class ~ .,data=train,metric="Spec",method="rf",trControl=control,tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=ntree,importance=T)
	key <- paste(toString(CV),"cv-metric_Spec-ntree",toString(ntree))
	modellist[[key]] <- fit                    
	#Get prediction
	res<-predict(fit,teste,type="prob");       
	result.roc<-roc(teste$Class,res$E)        
	predictionlist[[key]] <- result.roc

	#Train using a balanced set (392GE-392E)
	        #TRAIN By maximizing the ROC METRIC
        fit<-train(Class ~ .,data=balanced_Trainingset,metric="ROC",method="rf",trControl=control,tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=ntree,importance=T)
        key <- paste(toString(CV),"cv-metric_ROC-balanced-ntree",toString(ntree))
        modellist[[key]] <- fit
        #Get prediction
        res<-predict(fit,teste,type="prob");
        result.roc<-roc(teste$Class,res$E)
        predictionlist[[key]] <- result.roc

        #TRAIN By maximizing Specificity
        fit<-train(Class ~ .,data=balanced_Trainingset,metric="Spec",method="rf",trControl=control,tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=ntree,importance=T)
        key <- paste(toString(CV),"cv-metric_Spec-balanced-ntree",toString(ntree))
        modellist[[key]] <- fit
        #Get prediction
        res<-predict(fit,teste,type="prob");
        result.roc<-roc(teste$Class,res$E)
        predictionlist[[key]] <- result.roc

}
