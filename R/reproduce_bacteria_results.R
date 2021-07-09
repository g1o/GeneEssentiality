#' Extract features using the procariote data from the package, then do a complete train and test to reproduces the results from the chapter 1.
#' 
#' @param CPU number of threads to use, more threads uses more RAM and is faster. 
#' @export
reproduce_bacteria_results<-function(CPU=10){

extracttimes<-c(0,0,0,0,0,0);
starttime<- Sys.time()
AcinetoNE_pass<-Extract(FASTA_PATH=system.file("extdata", "Acinetobacter_baylyi_ADP1_NE_PASS.fasta.gz", package = "GeneEssentiality"),CPU=CPU)
gc()
endtime  <- Sys.time()
extracttimes[1] <- endtime - starttime
starttime<- Sys.time()
AcinetoNE_fail<-Extract(FASTA_PATH=system.file("extdata", "Acinetobacter_baylyi_ADP1_NE_Fail.fasta.gz", package = "GeneEssentiality"), AAfile =system.file("extdata","Acinetobacter_baylyi_ADP1_NE_Fail.faa.gz", package = "GeneEssentiality"),CPU=CPU);
gc();
endtime  <- Sys.time()
extracttimes[2] <- endtime - starttime
starttime<- Sys.time()
AcinetoE<-Extract(FASTA_PATH=system.file("extdata", "Acinetobacter_baylyi_ADP1_E.fasta.gz", package = "GeneEssentiality"),CPU=CPU)
gc()
endtime  <- Sys.time()
extracttimes[3] <- endtime - starttime
starttime<- Sys.time()
StaphyloNE_pass<-Extract(FASTA_PATH=system.file("extdata", "Staphylococcus_aureus_NCTC_8325_NE_PASS.fasta.gz", package = "GeneEssentiality"),CPU=CPU)
gc()
endtime  <- Sys.time()
extracttimes[4] <- endtime - starttime
starttime<- Sys.time()
StaphyloNE_fail<-Extract(FASTA_PATH=system.file("extdata", "Staphylococcus_aureus_NCTC_8325_NE_Fail.fasta.gz", package = "GeneEssentiality"),AAfile =system.file("extdata","Staphylococcus_aureus_NCTC_8325_NE_Fail.faa.gz", package = "GeneEssentiality"),CPU=CPU);
gc();
endtime  <- Sys.time()
extracttimes[5] <- endtime - starttime
starttime<- Sys.time()
StaphyloE<-Extract(FASTA_PATH=system.file("extdata", "Staphylococcus_aureus_NCTC_8325_E.fasta.gz", package = "GeneEssentiality"),CPU=CPU)
endtime  <- Sys.time()
extracttimes[6] <- endtime - starttime

StaphyloNE_fail$Class<-"NE"
StaphyloNE_pass$Class<-"NE"
Staphylo_fakeNEcds<-rbind(StaphyloE,StaphyloNE_fail)
Staphylo_trueNEcds<-rbind(StaphyloE,StaphyloNE_pass)
Staphylo<-rbind(StaphyloE,StaphyloNE_pass,StaphyloNE_fail)

AcinetoNE_fail$Class<-"NE"
AcinetoNE_pass$Class<-"NE"
Acineto_fakeNEcds<-rbind(AcinetoE,AcinetoNE_fail)
Acineto_trueNEcds<-rbind(AcinetoE,AcinetoNE_pass)
Acineto<-rbind(AcinetoE,AcinetoNE_pass,AcinetoNE_fail)
gc()
save(AcinetoNE_pass,AcinetoNE_fail,AcinetoE,StaphyloNE_pass,StaphyloNE_fail,StaphyloE,Acineto,Staphylo,file="2Bacterias_features.Rdata")

starttime<- Sys.time()
RF_Acineto_fakeNE<-train_rf(Acineto_fakeNEcds,CPU=CPU);gc()
RF_Acineto_trueNE<-train_rf(Acineto_trueNEcds,CPU=CPU);gc()
RF_Staphylo_fakeNE<-train_rf(Staphylo_fakeNEcds,CPU=CPU);gc()
RF_Staphylo_trueNE<-train_rf(Staphylo_trueNEcds,CPU=CPU);gc()
RF_Staphylo_mix<-train_rf(Staphylo,CPU=CPU);gc()
RF_Acineto_mix<-train_rf(Acineto,CPU=CPU);gc()
endtime  <- Sys.time()
Training_times <- endtime - starttime

#A models
Acinetotrue_Model_Staphylo_fake_Pred<-predict(RF_Acineto_trueNE$original_fit,Staphylo_fakeNEcds,type="prob")
Acinetotrue_Model_Staphylo_true_Pred<-predict(RF_Acineto_trueNE$original_fit,Staphylo_trueNEcds,type="prob")
Acinetotrue_Model_Staphylo_mix_Pred <-predict(RF_Acineto_trueNE$original_fit,Staphylo,type="prob")
ATSF_roc<-pROC::roc(Staphylo_fakeNEcds$Class, Acinetotrue_Model_Staphylo_fake_Pred$E,direction=">")
ATST_roc<-pROC::roc(Staphylo_trueNEcds$Class, Acinetotrue_Model_Staphylo_true_Pred$E,direction=">")
ATSM_roc<-pROC::roc(Staphylo$Class, Acinetotrue_Model_Staphylo_mix_Pred$E,direction=">")
                                                   
Staphylotrue_Model_Acineto_fake_Pred<-predict(RF_Staphylo_trueNE$original_fit,Acineto_fakeNEcds,type="prob")
Staphylotrue_Model_Acineto_true_Pred<-predict(RF_Staphylo_trueNE$original_fit,Acineto_trueNEcds,type="prob")
Staphylotrue_Model_Acineto_mix_Pred <-predict(RF_Staphylo_trueNE$original_fit,Acineto,type="prob")
STAF_roc<-pROC::roc(Acineto_fakeNEcds$Class,Staphylotrue_Model_Acineto_fake_Pred$E,direction=">")
STAT_roc<-pROC::roc(Acineto_trueNEcds$Class,Staphylotrue_Model_Acineto_true_Pred$E,direction=">")
STAM_roc<-pROC::roc(Acineto$Class,Staphylotrue_Model_Acineto_mix_Pred$E,direction=">")

#B models
Acinetofake_Model_Staphylo_fake_Pred<-predict(RF_Acineto_fakeNE$original_fit,Staphylo_fakeNEcds,type="prob")
Acinetofake_Model_Staphylo_true_Pred<-predict(RF_Acineto_fakeNE$original_fit,Staphylo_trueNEcds,type="prob")
Acinetofake_Model_Staphylo_mix_Pred <-predict(RF_Acineto_fakeNE$original_fit,Staphylo,type="prob")
AFSF_roc<-pROC::roc(Staphylo_fakeNEcds$Class, Acinetofake_Model_Staphylo_fake_Pred$E,direction=">")
AFST_roc<-pROC::roc(Staphylo_trueNEcds$Class, Acinetofake_Model_Staphylo_true_Pred$E,direction=">")
AFSM_roc<-pROC::roc(Staphylo$Class, Acinetofake_Model_Staphylo_mix_Pred$E,direction=">")

Staphylofake_Model_Acineto_fake_Pred<-predict(RF_Staphylo_fakeNE$original_fit,Acineto_fakeNEcds,type="prob")
Staphylofake_Model_Acineto_true_Pred<-predict(RF_Staphylo_fakeNE$original_fit,Acineto_trueNEcds,type="prob")
Staphylofake_Model_Acineto_mix_Pred <-predict(RF_Staphylo_fakeNE$original_fit,Acineto,type="prob")
SFAF_roc<-pROC::roc(Acineto_fakeNEcds$Class,Staphylofake_Model_Acineto_fake_Pred$E,direction=">")
SFAT_roc<-pROC::roc(Acineto_trueNEcds$Class,Staphylofake_Model_Acineto_true_Pred$E,direction=">")
SFAM_roc<-pROC::roc(Acineto$Class,Staphylofake_Model_Acineto_mix_Pred$E,direction=">")
#C models
Acinetomix_Model_Staphylo_fake_Pred<-predict(RF_Acineto_mix$original_fit,Staphylo_fakeNEcds,type="prob")
Acinetomix_Model_Staphylo_true_Pred<-predict(RF_Acineto_mix$original_fit,Staphylo_trueNEcds,type="prob")
Acinetomix_Model_Staphylo_mix_Pred <-predict(RF_Acineto_mix$original_fit,Staphylo,type="prob")
AMSF_roc<-pROC::roc(Staphylo_fakeNEcds$Class, Acinetomix_Model_Staphylo_fake_Pred$E,direction=">")
AMST_roc<-pROC::roc(Staphylo_trueNEcds$Class, Acinetomix_Model_Staphylo_true_Pred$E,direction=">")
AMSM_roc<-pROC::roc(Staphylo$Class, Acinetomix_Model_Staphylo_mix_Pred$E,direction=">")

Staphylomix_Model_Acineto_fake_Pred<-predict(RF_Staphylo_mix$original_fit,Acineto_fakeNEcds,type="prob")
Staphylomix_Model_Acineto_true_Pred<-predict(RF_Staphylo_mix$original_fit,Acineto_trueNEcds,type="prob")
Staphylomix_Model_Acineto_mix_Pred <-predict(RF_Staphylo_mix$original_fit,Acineto,type="prob")
SMAF_roc<-pROC::roc(Acineto_fakeNEcds$Class,Staphylomix_Model_Acineto_fake_Pred$E,direction=">")
SMAT_roc<-pROC::roc(Acineto_trueNEcds$Class,Staphylomix_Model_Acineto_true_Pred$E,direction=">")
SMAM_roc<-pROC::roc(Acineto$Class,Staphylomix_Model_Acineto_mix_Pred$E,direction=">")

#Zero Rule Model
ZRMSF<-pROC::roc(Staphylo_fakeNEcds$Class, rep(1,length(Acinetotrue_Model_Staphylo_fake_Pred$E)),direction=">")
ZRMST<-pROC::roc(Staphylo_trueNEcds$Class, rep(1,length(Acinetotrue_Model_Staphylo_true_Pred$E)),direction=">")
ZRMSM<-pROC::roc(Staphylo$Class, rep(1,length(Acinetotrue_Model_Staphylo_mix_Pred$E)),direction=">")
ZRMAF<-pROC::roc(Acineto_fakeNEcds$Class, rep(1,length(Staphylotrue_Model_Acineto_fake_Pred$E)),direction=">")
ZRMAT<-pROC::roc(Acineto_trueNEcds$Class, rep(1,length(Staphylotrue_Model_Acineto_true_Pred$E)),direction=">")
ZRMAM<-pROC::roc(Acineto$Class, rep(1,length(Staphylotrue_Model_Acineto_mix_Pred$E)),direction=">")


save(SMAF_roc,SMAT_roc,SMAM_roc,AMSF_roc,AMST_roc,AMSM_roc,
     ATSF_roc,ATST_roc,ATSM_roc,STAF_roc,STAT_roc,STAM_roc,
     AFSF_roc,AFSM_roc,AFST_roc,SFAF_roc,SFAT_roc,SFAM_roc,
     ZRMSF,ZRMST,ZRMSM,ZRMAF,ZRMAT,ZRMAM,file="Bacteria_ROCKs.Rdata")


############### -PLOTS- ###############
library(pROC)


##############A-true model
svg("A_trueCDS_rfmodels.svg",height=6,width=12)
size=0.9
par(mfrow=c(1,2),oma=c(0,0,1,0)) 
plot(ATST_roc,legend=F,color="black",auc.main=F,print.auc=T,print.auc.cex=size,print.auc.x=0.2,print.auc.y=0.7,
        main="Modelo A. baylyi: teste em S. aureus")
plot(ATSF_roc,legend=F,col="red",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.65,lty=2)
plot(ATSM_roc,legend=F,col="blue",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.60,lty=3)
pvalue1<-roc.test(ATST_roc,ATSF_roc)
pvalue2<-roc.test(ATSF_roc,ATSM_roc)
pvalue3<-roc.test(ATSM_roc,ATST_roc)
legend(x=0.35,y=0.55,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue2$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue3$p.value,digits=3,scientific = T))),col=c("red","blue","black"), lty=c(2,3,1),xpd = TRUE,lwd=2, bty = "n")
legend(x=0.5 ,y=0.55,legend=c("vs","vs","vs"),x.intersp=0.5,col=c("black","red","blue"), lty=1:3,xpd = TRUE,lwd=2, bty = "n")
legend("bottomright",legend=c("CDSs verdadeiras","CDSs falsas","CDSs misturadas"),col=c("black","red","blue"), lty=1:3,lwd=2,xpd = TRUE)
pvalue1<-roc.test(ATST_roc,ZRMST)
pvalue2<-roc.test(ATSF_roc,ZRMSF)
pvalue3<-roc.test(ATSM_roc,ZRMSM)
legend(x=0.35,y=0.38,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue2$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue3$p.value,digits=3,scientific = T))),col=c("GRAY","GRAY","GRAY"), lty=1,xpd = TRUE,lwd=2, bty = "n")
legend(x=0.5 ,y=0.38,legend=c("vs","vs","vs"),x.intersp=0.5,col=c("black","red","blue"), lty=1:3,xpd = TRUE,lwd=2, bty = "n")

plot(STAT_roc,legend=F,col="black",auc.main=F,print.auc=T,print.auc.cex=size,print.auc.x=0.2,print.auc.y=0.7,
        main="Modelo S. aureus: teste em A. baylyi")
plot(STAF_roc,legend=F,col="red",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.65,lty=2)
plot(STAM_roc,legend=F,col="blue",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.60,lty=3)
mtext("A) Modelos treinados com CDSs verdadeiras", side = 3,line=-0.5, outer = TRUE,cex=1.5)
pvalue1<-roc.test(STAT_roc,STAF_roc)
pvalue2<-roc.test(STAF_roc,STAM_roc)
pvalue3<-roc.test(STAM_roc,STAT_roc)
legend(x=0.35,y=0.55,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue2$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue3$p.value,digits=3,scientific = T))),col=c("red","blue","black"), lty=c(2,3,1),xpd = TRUE,lwd=2, bty = "n")
legend(x=0.5 ,y=0.55,legend=c("vs","vs","vs"),x.intersp=0.5,inset =c(-1.5,-1.2) ,col=c("black","red","blue"), lty=1:3,xpd = TRUE,lwd=2, bty = "n")
legend("bottomright",legend=c("CDSs verdadeiras","CDSs falsas","CDSs misturadas"),col=c("black","red","blue"), lty=1:3,lwd=2,xpd = TRUE)
pvalue1<-roc.test(STAT_roc,ZRMAT)
pvalue2<-roc.test(STAF_roc,ZRMAF)
pvalue3<-roc.test(STAM_roc,ZRMAM)
legend(x=0.35,y=0.38,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue2$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue3$p.value,digits=3,scientific = T))),col=c("GRAY","GRAY","GRAY"), lty=1,xpd = TRUE,lwd=2, bty = "n")
legend(x=0.5 ,y=0.38,legend=c("vs","vs","vs"),x.intersp=0.5,col=c("black","red","blue"), lty=1:3,xpd = TRUE,lwd=2, bty = "n")
dev.off()

###############
###############B-fake model
svg("B_fakeCDS_rfmodels.svg",height=6,width=12)
par(mfrow=c(1,2),oma=c(0,0,1,0)) 
plot(AFST_roc,legend=F,col="black",auc.main=F,print.auc=T,print.auc.cex=size,print.auc.x=0.2,print.auc.y=0.7,
        main="Modelo A. baylyi: teste em S. aureus")
plot(AFSF_roc,legend=F,col="red",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.65,lty=2)
plot(AFSM_roc,legend=F,col="blue",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.6,lty=3)
pvalue1<-roc.test(AFST_roc,AFSF_roc)
pvalue2<-roc.test(AFSF_roc,AFSM_roc)
pvalue3<-roc.test(AFSM_roc,AFST_roc)
legend(x=0.35,y=0.55,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue2$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue3$p.value,digits=3,scientific = T))),col=c("red","blue","black"), lty=c(2,3,1),xpd = TRUE,lwd=2, bty = "n")
legend(x=0.5 ,y=0.55,legend=c("vs","vs","vs"),x.intersp=0.5,inset =c(-1.5,-1.2) ,col=c("black","red","blue"), lty=1:3,xpd = TRUE,lwd=2, bty = "n")
legend("bottomright",legend=c("CDSs verdadeiras","CDSs falsas","CDSs misturadas"),col=c("black","red","blue"), lty=1:3,lwd=2,xpd = TRUE)
pvalue1<-roc.test(AFST_roc,ZRMST)
pvalue2<-roc.test(AFSF_roc,ZRMSF)
pvalue3<-roc.test(AFSM_roc,ZRMSM)
legend(x=0.35,y=0.38,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue2$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue3$p.value,digits=3,scientific = T))),col=c("GRAY","GRAY","GRAY"), lty=1,xpd = TRUE,lwd=2, bty = "n")
legend(x=0.5 ,y=0.38,legend=c("vs","vs","vs"),x.intersp=0.5,col=c("black","red","blue"), lty=1:3,xpd = TRUE,lwd=2, bty = "n")

plot(SFAT_roc,legend=F,col="black",auc.main=F,print.auc=T,print.auc.cex=size,print.auc.x=0.2,print.auc.y=0.7,
        main="Modelo S. aureus: teste em A. baylyi")
plot(SFAF_roc,legend=F,col="red",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.65,lty=2)
plot(SFAM_roc,legend=F,col="blue",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.6,lty=3)
mtext("B) Modelos treinados com CDSs falsas", side = 3,line=-0.5, outer = TRUE,cex=1.5)
pvalue1<-roc.test(SFAT_roc,SFAF_roc)
pvalue2<-roc.test(SFAF_roc,SFAM_roc)
pvalue3<-roc.test(SFAM_roc,SFAT_roc)
legend(x=0.35,y=0.55,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue2$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue3$p.value,digits=3,scientific = T))),col=c("red","blue","black"), lty=c(2,3,1),xpd = TRUE,lwd=2, bty = "n")
legend(x=0.5 ,y=0.55,legend=c("vs","vs","vs"),x.intersp=0.5,inset =c(-1.5,-1.2) ,col=c("black","red","blue"), lty=1:3,xpd = TRUE,lwd=2, bty = "n")
legend("bottomright",legend=c("CDSs verdadeiras","CDSs falsas","CDSs misturadas"),col=c("black","red","blue"), lty=1:3,lwd=2,xpd = TRUE)
pvalue1<-roc.test(SFAT_roc,ZRMAT)
pvalue2<-roc.test(SFAF_roc,ZRMAF)
pvalue3<-roc.test(SFAM_roc,ZRMAM)
legend(x=0.35,y=0.38,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue2$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue3$p.value,digits=3,scientific = T))),col=c("GRAY","GRAY","GRAY"), lty=1,xpd = TRUE,lwd=2, bty = "n")
legend(x=0.5 ,y=0.38,legend=c("vs","vs","vs"),x.intersp=0.5,col=c("black","red","blue"), lty=1:3,xpd = TRUE,lwd=2, bty = "n")

dev.off()
###############C-mix model
svg("C_mixCDS_rfmodels.svg",height=6,width=12)
par(mfrow=c(1,2),oma=c(0,0,1,0) )
plot(AMST_roc,legend=F,col="black",auc.main=F, print.auc=T,print.auc.cex=size,print.auc.x=0.2,print.auc.y=0.7,
        main="Modelo A. baylyi: teste em S. aureus")
plot(AMSF_roc,legend=F,col="red",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.65,lty=2)
plot(AMSM_roc,legend=F,col="blue",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.6,lty=3)
pvalue1<-roc.test(AMST_roc,AMSF_roc)
pvalue2<-roc.test(AMSF_roc,AMSM_roc)
pvalue3<-roc.test(AMSM_roc,AMST_roc)
legend(x=0.35,y=0.55,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue2$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue3$p.value,digits=3,scientific = T))),col=c("red","blue","black"), lty=c(2,3,1),xpd = TRUE,lwd=2, bty = "n")
legend(x=0.5 ,y=0.55,legend=c("vs","vs","vs"),x.intersp=0.5,inset =c(-1.5,-1.2) ,col=c("black","red","blue"), lty=1:3,xpd = TRUE,lwd=2, bty = "n")
legend("bottomright",legend=c("CDSs verdadeiras","CDSs falsas","CDSs misturadas"),col=c("black","red","blue"), lty=1:3,lwd=2,xpd = TRUE)
pvalue1<-roc.test(AMST_roc,ZRMST)
pvalue2<-roc.test(AMSF_roc,ZRMSF)
pvalue3<-roc.test(AMSM_roc,ZRMSM)
legend(x=0.35,y=0.38,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue2$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue3$p.value,digits=3,scientific = T))),col=c("GRAY","GRAY","GRAY"), lty=1,xpd = TRUE,lwd=2, bty = "n")
legend(x=0.5 ,y=0.38,legend=c("vs","vs","vs"),x.intersp=0.5,col=c("black","red","blue"), lty=1:3,xpd = TRUE,lwd=2, bty = "n")

plot(SMAT_roc,legend=F,col="black",auc.main=F,print.auc=T,print.auc.cex=size,print.auc.x=0.2,print.auc.y=0.7,
        main="Modelo S. aureus: teste em A. baylyi")
plot(SMAF_roc,legend=F,col="red",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.65,lty=2)
plot(SMAM_roc,legend=F,col="blue",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.6,lty=3)
mtext("C) Modelos treinados com CDSs misturadas", side = 3,line=-0.5, outer = T,cex=1.5)
pvalue1<-roc.test(SMAT_roc,SMAF_roc)
pvalue2<-roc.test(SMAF_roc,SMAM_roc)
pvalue3<-roc.test(SMAM_roc,SMAT_roc)
legend(x=0.35,y=0.55,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue2$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue3$p.value,digits=3,scientific = T))),col=c("red","blue","black"), lty=c(2,3,1),xpd = TRUE,lwd=2, bty = "n")
legend(x=0.5 ,y=0.55,legend=c("vs","vs","vs"),x.intersp=0.5,inset =c(-1.5,-1.2) ,col=c("black","red","blue"), lty=1:3,xpd = TRUE,lwd=2, bty = "n")
legend("bottomright",legend=c("CDSs verdadeiras","CDSs falsas","CDSs misturadas"),col=c("black","red","blue"), lty=1:3,lwd=2,xpd = TRUE)
pvalue1<-roc.test(SMAT_roc,ZRMAT)
pvalue2<-roc.test(SMAF_roc,ZRMAF)
pvalue3<-roc.test(SMAM_roc,ZRMAM)
legend(x=0.35,y=0.38,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue2$p.value,digits=3,scientific = T)),paste0("P = ",format(pvalue3$p.value,digits=3,scientific = T))),col=c("GRAY","GRAY","GRAY"), lty=1,xpd = TRUE,lwd=2, bty = "n")
legend(x=0.5 ,y=0.38,legend=c("vs","vs","vs"),x.intersp=0.5,col=c("black","red","blue"), lty=1:3,xpd = TRUE,lwd=2, bty = "n")

dev.off()
###############


svg("Figura8.svg",height=10,width=10)
par(mfrow=c(2,2),oma=c(0,0,1,0))
#preta 5 vs 7 (A vs C)
size=0.9
plot(ATST_roc,legend=F,color="black",auc.main=F,print.auc=T,print.auc.cex=size,print.auc.x=0.2,print.auc.y=0.7,
     main="a) Modelo A. baylyi: teste em S. aureus")
plot(AMST_roc,legend=F,col="black",auc.main=F, add=T,print.auc=T,print.auc.cex=size,print.auc.x=0.2,print.auc.y=0.65,lty=3)
pvalue1<-roc.test(ATST_roc,AMST_roc)
legend(x=0.32,y=0.735,legend=c("",""),col=c("black","BLACK"), lty=c(1,3),xpd = TRUE,lwd=2, bty = "n")
legend(x=0.35,y=0.38,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T))),col=c("BLACK"), bty = "n")
legend("bottomright",legend=c("A: teste CDSs verdadeiras","C: teste CDSs verdadeiras"),col=c("black","black"), lty=c(1,3),lwd=2,xpd = TRUE)

plot(STAT_roc,legend=F,col="black",auc.main=F,print.auc=T,print.auc.cex=size,print.auc.x=0.2,print.auc.y=0.7,
     main="b) Modelo S. aureus: teste em A. baylyi")
plot(SMAT_roc,legend=F,col="black",auc.main=F, add=T,print.auc=T,print.auc.cex=size,print.auc.x=0.2,print.auc.y=0.65,lty=3)
pvalue1<-roc.test(STAT_roc,SMAT_roc)
legend(x=0.35,y=0.38,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T))),col=c("BLACK"), bty = "n")
legend(x=0.32,y=0.735,legend=c("",""),col=c("black","BLACK"), lty=c(1,3),xpd = TRUE,lwd=2, bty = "n")
legend("bottomright",legend=c("A: teste CDSs verdadeiras","C: teste CDSs verdadeiras"),col=c("black","black"), lty=c(1,3),lwd=2,xpd = TRUE)
#preta 5 vs azul 7
plot(ATST_roc,legend=F,color="black",auc.main=F,print.auc=T,print.auc.cex=size,print.auc.x=0.2,print.auc.y=0.7,
     main="c) Modelo A. baylyi: teste em S. aureus")
plot(AMSM_roc,legend=F,col="blue",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.65,lty=3)
pvalue1<-roc.test(ATST_roc,AMSM_roc)
legend(x=0.35,y=0.38,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T))),col=c("BLACK"), bty = "n")
legend(x=0.32,y=0.735,legend=c("",""),col=c("black","blue"), lty=c(1,3),xpd = TRUE,lwd=2, bty = "n")
legend("bottomright",legend=c("A: teste CDSs verdadeiras","C: teste CDSs misturadas"),col=c("black","blue"), lty=c(1,3),lwd=2,xpd = TRUE)

plot(STAT_roc,legend=F,col="black",auc.main=F,print.auc=T,print.auc.cex=size,print.auc.x=0.2,print.auc.y=0.7,
     main="d) Modelo S. aureus: teste em A. baylyi")
plot(SMAM_roc,legend=F,col="blue",auc.main=F,print.auc=T,print.auc.cex=size,add=T,print.auc.x=0.2,print.auc.y=0.65,lty=3)
pvalue1<-roc.test(STAT_roc,SMAM_roc)
legend(x=0.35,y=0.38,legend=c(paste0("P = ",format(pvalue1$p.value,digits=3,scientific = T))),col=c("BLACK"), bty = "n")
legend(x=0.32,y=0.735,legend=c("",""),col=c("black","blue"), lty=c(1,3),xpd = TRUE,lwd=2, bty = "n")
legend("bottomright",legend=c("A: teste CDSs verdadeiras","C: teste CDSs misturadas"),col=c("black","blue"), lty=c(1,3),lwd=2,xpd = TRUE)
dev.off()

resultados<-(list(extracttimes,Training_times,RF_Acineto_fakeNE,RF_Acineto_trueNE,RF_Staphylo_fakeNE,RF_Staphylo_trueNE,RF_Staphylo_mix,RF_Acineto_mix))
names(resultados)<-c("extracttimes","Training_times","RF_Acineto_fakeNE","RF_Acineto_trueNE","RF_Staphylo_fakeNE","RF_Staphylo_trueNE","RF_Staphylo_mix","RF_Acineto_mix")
return (resultados)

}
