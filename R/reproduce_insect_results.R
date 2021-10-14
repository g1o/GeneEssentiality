#' Reproduce Drosophila melanogaster and Tribolium castaneum train and tests 
#'
#' Run a complete train and test using the data from the package. Reproduces the paper results 
#' Can also be used to compare two datasets with the same comparsion made for the insects. The required column is the Clas, and its values must be "E" or "NE".
#' Reditect to an object to save the models, ROCs results and Pvalues from the De longs test.  Plots are outputed to a file in the working dir.
#'
#' @param CPU Number of threads to use
#' @param set1 Data.frame with the features of the first set. This will be used to train a model and as a test for the model trained from the second set. Default is the Drosophila melanogaster features
#' @param set2 Data.frame with the features of the second set. This will be used to train a model and as a test for the model trained from the first set. Default is the Tribolium castaneum features
#' @param seeds List of vectors, 30 lists of vectors with 6 elements and the last list with a single number, for the final model
#' @param plot_prefix Input a prefix for the filename of the plots. Default is "Insects_"
#' @param set1.name Name to be used for the title of the plots for the first set. Default is "Dmel"
#' @param set2.name Name to be used for the title of the plots for the first set. Default is "Trib"
#' @export
reproduce_insect_results<-function(set1=GeneEssentiality::drosophila_features,set2=GeneEssentiality::tribolium_features,CPU=20,seeds=GeneEssentiality::seed,plot_prefix="Insects_",set1.name="Dmel",set2.name="Trib"){
	trees=1000;
	CV=10;
	repeats=3;
	seed<-seeds;
	models<-lapply(list(set1,set2), function(Complete_set){
			Select_Complete_set<- Select(data=Complete_set);
			rfmodels<-train_rf(features=Select_Complete_set,CPU=CPU,trees=trees,CV=CV,nrepeats=repeats,seeds=seed);
			xgbtmodels<-train_xgbt(features=Select_Complete_set,CPU=CPU,CV=CV,nrepeats=repeats,seeds=seed);
			models<-list(rfmodels[[1]],xgbtmodels)
			return(	models)	} )

#dmel_importances<-lapply(models[[1]],function(x){varImp(x)});

DMELM<-models[[1]];
TRIBM<-models[[2]];

models<-DMELM;
	roc_noh_trib<-lapply(models,function(x){ res<-predict(x, GeneEssentiality::noh_trib, type="prob");
		pROC::roc(GeneEssentiality::noh_trib$Class,res$E,direction=">") } );
	alle<-as.data.frame(roc_noh_trib[[1]]$predictor/roc_noh_trib[[1]]$predictor)
	alle$NE<-0
	colnames(alle)<-c("E","NE")
	ZRroc<-pROC::roc(GeneEssentiality::noh_trib$Class,alle$E,direction=">")
	pvalues_noh_trib<-lapply(roc_noh_trib,function(rroc){ pROC::roc.test(rroc,ZRroc) } )

	roc_trib<-lapply(models,function(x){ res<-predict(x,set2,type="prob");
		pROC::roc(set2$Class,res$E,direction=">") } );
	alle<-as.data.frame(roc_trib[[1]]$predictor/roc_trib[[1]]$predictor)
	alle$NE<-0
	colnames(alle)<-c("E","NE")
	ZRroc<-pROC::roc(set2$Class,alle$E,direction=">")
	pvalues_trib<-lapply(roc_trib,function(rroc){ pROC::roc.test(rroc,ZRroc) } )

	finaldmel<-list(DMELM,roc_noh_trib,pvalues_noh_trib,roc_trib,pvalues_trib)
	names(finaldmel)<-c("Dmel_models","roc_noh_trib","pvalues_noh_trib","roc_trib","pvalues_trib")
models<-TRIBM;

        alle<-as.data.frame(roc_dmel[[1]]$predictor/roc_dmel[[1]]$predictor)
        alle$NE<-0
        colnames(alle)<-c("E","NE")
        ZRroc<-pROC::roc(set1$Class,alle$E,direction=">")
        roc_dmel<-lapply(models,function(x){ res<-predict(x,set1,type="prob");
                pROC::roc(set1$Class,res$E,direction=">") } );
        pvalues_dmel<-lapply(roc_dmel,function(rroc){ pROC::roc.test(rroc,ZRroc) } )
        
	finaltrib<-list(TRIBM,roc_dmel,pvalues_dmel)
	names(finaltrib)<-c("Trib_models","roc_dmel","pvalues_dmel")
	names(finaltrib$pvalues_dmel)<-c("RF","XGBT")
	names(finaltrib$roc_dmel)<-c("RF","XGBT")

	names(finaldmel$pvalues_noh_trib)<-c("RF","XGBT")
	names(finaldmel$pvalues_trib)<-c("RF","XGBT")
	names(finaldmel$roc_noh_trib)<-c("RF","XGBT")
	names(finaldmel$roc_trib)<-c("RF","XGBT")

#Ploting results
library("pROC")
library("PRROC")
plotname<-paste0(plot_prefix,"ROCs.svg")
svg(plotname,height=4,width=12)
par(mfrow=c(1,3))
size=0.9;
######################  Dmel models vs Trib
pvalue<-  finaldmel$pvalues_trib$RF$p.value
ROC_title<-paste0("ROC: ",set1.name," models vs ",set2.name)
plot(finaldmel$roc_trib$RF,legend=F,color="black",print.auc=T,print.auc.cex=size,main=ROC_title)
text(x=0.25, y = 0.35, labels = paste0("P-value= ",signif(pvalue,digits=3)),cex=size)
pvalue<- finaldmel$pvalues_trib$XGBT$p.value

plot(finaldmel$roc_trib$XGBT,legend=F,print.auc=T,print.auc.y=0.455 , add=T,print.auc.cex=size,col="red")
text(x=0.25, y = 0.30, labels = paste0("P-value= ",signif(pvalue,digits=3)),cex=size,col="red")
legend(x="topleft",c("RF","XGBT"),col=c(1,2),lty=1)
title(adj=0,line=2.5,main="A)",cex.main=2)
#####################  Trib models vs dmel
pvalue<- finaltrib$pvalues_dmel$RF$p.value
ROC_title<-paste0("ROC: ",set2.name," models vs ",set1.name)
plot(finaltrib$roc_dmel$RF ,legend=F,color="black",print.auc=T,print.auc.cex=size,main=ROC_title)
text(x=0.25, y = 0.35, labels = paste0("P-value= ",signif(pvalue,digits=3)),cex=size)
pvalue<- finaltrib$pvalues_dmel$XGBT$p.value 
plot(finaltrib$roc_dmel$XGBT,legend=F,print.auc=T,print.auc.y=0.455 , add=T,print.auc.cex=size,col="red")
text(x=0.25, y = 0.30, labels = paste0("P-value= ",signif(pvalue,digits=3)),cex=size,col="red")
legend(x="topleft",c("RF","XGBT"),col=c(1,2),lty=1)
title(adj=0,line=2.5,main="B)",cex.main=2)

#####################  Dmel models vs noh-Trib
pvalue<- finaldmel$pvalues_noh_trib$RF$p.value 
ROC_title<-paste0("ROC: ",set1.name," models vs noh-Trib")
plot(finaldmel$roc_noh_trib$RF,legend=F,color="black",print.auc=T,print.auc.cex=size,main=ROC_title)
text(x=0.25, y = 0.35, labels = paste0("P-value= ",signif(pvalue,digits=3)),cex=size)
pvalue<- finaldmel$pvalues_noh_trib$XGBT$p.value
plot(finaldmel$roc_noh_trib$XGBT,legend=F,print.auc=T,print.auc.y=0.455 , add=T,print.auc.cex=size,col="red")
text(x=0.25, y = 0.30, labels = paste0("P-value= ",signif(pvalue,digits=3)),cex=size,col="red")
legend(x="topleft",c("RF","XGBT"),col=c(1,2),lty=1)
title(adj=0,line=2.5,main="C)",cex.main=2)
dev.off()


####################  Precision Recall curves 
#svg("Insects_PRCs.svg",height=4,width=12)
plotname<-paste0(plot_prefix,"PRCs.svg")
svg(plotname,height=4,width=12)
par(mfrow=c(1,3))
prcurve_noh1 = pr.curve(c(finaldmel$roc_trib$RF$original.predictor),weights.class0=(finaldmel$roc_trib$RF$original.response=='E')*1,curve=T) #rf
prcurve_noh2 = pr.curve(c(finaldmel$roc_trib$XGBT$original.predictor),weights.class0=(finaldmel$roc_trib$XGBT$original.response=='E')*1,curve=T) #xgbt
prcurve_noh_ZR = pr.curve(c(finaldmel$pvalues_trib$RF$roc2$original.predictor[-1],0),weights.class0=(finaldmel$pvalues_trib$RF$roc2$original.response=='E')*1,curve=T)
PRC_title<-paste0("PRC: ",set1.name," models vs ",set2.name)
plot(prcurve_noh1,auc.main=F,main=PRC_title,color="black")
plot(prcurve_noh2,auc.main=F,add=T,color="red")
plot(prcurve_noh_ZR,auc.main=F,add=T,color="gray",lty=2)
text(x=0.35, y = 0.95, labels =paste0(" RF  AUC= ", signif(prcurve_noh1$auc.integral,digits=3)))
text(x=0.35, y = 0.90, labels =paste0("XGBT AUC= ",  signif(prcurve_noh2$auc.integral,digits=3)),col="red")
text(x=0.35, y = 0.85, labels =paste0("ZR AUC= ", signif(prcurve_noh_ZR$auc.integral,digits=3)),col="dark gray")
legend(x="topright",c("RF","XGBT"),col=c(1,2),lty=1)
title(adj=0,line=1.5,main="A)",cex.main=2)

##########################################################################################################
prcurve_noh1 = pr.curve(c(finaltrib$roc_dmel$RF$original.predictor),weights.class0=(finaltrib$roc_dmel$RF$original.response=='E')*1,curve=T) #rf 
prcurve_noh2 = pr.curve(c(finaltrib$roc_dmel$XGBT$original.predictor),weights.class0=(finaltrib$roc_dmel$XGBT$original.response=='E')*1,curve=T) #xgbt
prcurve_noh_ZR = pr.curve(c(finaltrib$pvalues_dmel$RF$roc2$original.predictor[-1],0),weights.class0=(finaltrib$pvalues_dmel$RF$roc2$original.response=='E')*1,curve=T)

PRC_title<-paste0("PRC: ",set2.name," models vs ",set1.name)
plot(prcurve_noh1,auc.main=F,main=PRC_title,color="black")
plot(prcurve_noh2,auc.main=F,add=T,color="red")
plot(prcurve_noh_ZR,auc.main=F,add=T,color="gray",lty=2)
text(x=0.35, y = 0.95, labels =paste0(" RF  AUC= ", signif(prcurve_noh1$auc.integral,digits=3)))
text(x=0.35, y = 0.90, labels =paste0("XGBT AUC= ", signif(prcurve_noh2$auc.integral,digits=3)),col="red")
text(x=0.35, y = 0.85, labels =paste0("ZR AUC= ", signif(prcurve_noh_ZR$auc.integral,digits=3)),col="dark gray")
legend(x="topright",c("RF","XGBT"),col=c(1,2),lty=1)
title(adj=0,line=1.5,main="B)",cex.main=2)

#############################################################################################################
prcurve_noh1 = pr.curve(c(finaldmel$roc_noh_trib$RF$original.predictor),weights.class0=(finaldmel$roc_noh_trib$RF$original.response=='E')*1,curve=T) #rf 
prcurve_noh2 = pr.curve(c(finaldmel$roc_noh_trib$XGBT$original.predictor),weights.class0=(finaldmel$roc_noh_trib$XGBT$original.response=='E')*1,curve=T) #xgbt
prcurve_noh_ZR = pr.curve(c(finaldmel$pvalues_noh_trib$RF$roc2$original.predictor[-1],0),weights.class0=(finaldmel$pvalues_noh_trib$RF$roc2$original.response=='E')*1,curve=T)

PRC_title<-paste0("PRC: ",set1.name," models vs noh-Trib")
plot(prcurve_noh1,auc.main=F,main=PRC_title,color="black")
plot(prcurve_noh2,auc.main=F,add=T,color="red")
plot(prcurve_noh_ZR,auc.main=F,add=T,color="gray",lty=2)
text(x=0.35, y = 0.95, labels =paste0(" RF  AUC= ", prcurve_noh1$auc.integral))
text(x=0.35, y = 0.90, labels =paste0("XGBT AUC= ", prcurve_noh2$auc.integral),col="red")
text(x=0.35, y = 0.85, labels =paste0("ZR AUC= ", prcurve_noh_ZR$auc.integral),col="dark gray")
legend(x="topright",c("RF","XGBT"),col=c(1,2),lty=1)
title(adj=0,line=1.5,main="C)",cex.main=2)
dev.off()


#Cross validation Plots

#svg("Insects_CV_ROCs.svg",height=6,width=12)
plotname<-paste0(plot_prefix,"CV_ROCs.svg")
svg(plotname,height=6,width=12)

size=0.9;
par(mfrow=c(1,2))
#,2=> Obs labels; ,3=> probabilities for "E"
Dmtry<-DMELM[[1]]$bestTune[1,1]
Tmtry<-TRIBM[[1]]$bestTune[1,1]
roc_dmelm244 <-roc (DMELM[[1]]$pred[DMELM[[1]]$pred$mtry == Dmtry,2], DMELM[[1]]$pred[DMELM[[1]]$pred$mtry == Dmtry,3],direction=">")
roc_tribm244 <-roc (TRIBM[[1]]$pred[TRIBM[[1]]$pred$mtry == Tmtry,2], TRIBM[[1]]$pred[TRIBM[[1]]$pred$mtry == Tmtry,3],direction=">")
#,2=> Obs labels; ,4=> probabilities for "E"
roc_xg_dmelmd10<-roc (DMELM[[2]]$pred[DMELM[[2]]$pred$max_depth == 10,2], DMELM[[2]]$pred[DMELM[[2]]$pred$max_depth == 10,4],direction=">") 
roc_xg_tribmd10<-roc (TRIBM[[2]]$pred[TRIBM[[2]]$pred$max_depth == 10,2], TRIBM[[2]]$pred[TRIBM[[2]]$pred$max_depth == 10,4],direction=">")

alle<-as.data.frame(roc_dmel[[1]]$predictor/roc_dmel[[1]]$predictor)
alle$NE<-0
colnames(alle)<-c("E","NE")
ZRroc<-pROC::roc(set1$Class,alle$E,direction=">")
pvalue_cv_rf<- pROC::roc.test(roc_dmelm244,ZRroc)
pvalue_cv_xg<- pROC::roc.test(roc_xg_dmelmd10,ZRroc)
ROC_title<-paste0("ROC: ",set1.name,"  cross-validation, XGBT vs RF")
plot(roc_dmelm244,legend=F,color="black",auc.main=F,main=ROC_title, print.auc=T,print.auc.cex=size)
plot(roc_xg_dmelmd10,legend=F,color="black",auc.main=F, print.auc=T,print.auc.cex=size,add=T,col="red",print.auc.y=0.45)
pvalue<-roc.test(roc_xg_dmelmd10,roc_dmelm244)
text(x=0.315, y = 0.4, labels = paste0("XGBT vs RF: P= ",signif(pvalue$p.value,digits=2)),cex=size)
text(x=0.315, y = 0.35, labels = paste0("XGBT vs ZR: P= ",signif(pvalue_cv_xg$p.value,digits=2)),cex=size)
text(x=0.315, y = 0.3, labels = paste0("RF vs ZR: P= ",signif(pvalue_cv_rf$p.value,digits=2)),cex=size)
legend(x="bottomright",c("RF","XGBT"),col=c(1,2),lty=1)

alle<-as.data.frame(roc_trib[[1]]$predictor/roc_trib[[1]]$predictor)
alle$NE<-0
colnames(alle)<-c("E","NE")
ZRroc<-pROC::roc(set1$Class,alle$E,direction=">")
pvalue_cv_rf<- pROC::roc.test(roc_dmelm244,ZRroc)
pvalue_cv_xg<- pROC::roc.test(roc_xg_dmelmd10,ZRroc)
ROC_title<-paste0("ROC: ",set2.name,"  cross-validation, XGBT vs RF")
plot(roc_tribm244,legend=F,color="black",auc.main=F,main=ROC_title, print.auc=T,print.auc.cex=size)
plot(roc_xg_tribmd10,legend=F,color="black",auc.main=F, print.auc=T,print.auc.cex=size,add=T,col="red",print.auc.y=0.45)
pvalue<-roc.test(roc_xg_tribmd10,roc_tribm244)
legend(x="bottomright",c("RF","XGBT"),col=c(1,2),lty=1)
text(x=0.315, y = 0.4, labels = paste0("XGBT vs RF: P= ",signif(pvalue$p.value,digits=2)),cex=size)
text(x=0.315, y = 0.35, labels = paste0("XGBT vs ZR: P= ",signif(pvalue_cv_xg$p.value,digits=2)),cex=size)
text(x=0.315, y = 0.3, labels = paste0("RF vs ZR: P= ",signif(pvalue_cv_rf$p.value,digits=2)),cex=size)
dev.off()

	final_results<-(list(finaldmel,finaltrib))
	names(final_results)<-c("DMEL","TRIB")
	return(final_results)
}
