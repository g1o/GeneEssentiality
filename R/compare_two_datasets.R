#' Compare two datasets
#'
#' Run a complete train and test using the provided features, both dataframes must have the same features.
#' Used to compare two datasets with the same comparsion made for the insects, leaving one dataset out. The required column is the Class, and its values must be "E" or "NE".
#' Reditect to an object to save the models, ROCs results and Pvalues.  Plots are outputed to a file in the working dir.
#'
#' @param CPU Number of threads to use
#' @param set1 Data.frame with the features of the first set. This will be used to train a model and as a test for the model trained from the second set. Default is the Drosophila melanogaster features
#' @param set2 Data.frame with the features of the second set. This will be used to train a model and as a test for the model trained from the first set. Default is the Tribolium castaneum features
#' @param seeds List of vectors, 30 lists of vectors with 6 elements, and the last list with a single number for the final model
#' @param plot_prefix Input a prefix for the filename of the plots. Default is "Insects_"
#' @param set1.name Name to be used for the title of the plots for the first set. Default is "Set1"
#' @param set2.name Name to be used for the title of the plots for the first set. Default is "Set2"
#' @param XGBT Use eXtreme Gradient Boosting trees model #slow
#' @param RF Use Randon Forest model
#' @export
compare_two_datasets<-function(set1=GeneEssentiality::drosophila_features,set2=GeneEssentiality::tribolium_features,CPU=20,seeds=GeneEssentiality::seed,plot_prefix="Compared_",set1.name="Set1",set2.name="Set2",XGBT=T,RF=T){
	trees=1000;
	CV=10;
	repeats=3;
	seed<-seeds;
library("PRROC")
	models<-lapply(list(set1,set2), function(Complete_set){
			models<-list();
			i<-0;
			if(RF==T){
			  rfmodels<-train_rf(features=Complete_set,CPU=CPU,trees=trees,CV=CV,nrepeats=repeats,seeds=seed);
			  i<-i+1;
			  models[[i]]<-rfmodels[[1]];
			}else{
			 #place holder to avoid downstream problems
			  i<-i+1;
			  models[[i]] <- train(Class ~ ., data = Complete_set , metric = "ROC", method = "null" ,trControl=trainControl(classProbs = TRUE,summaryFunction = twoClassSummary,savePredictions = T ,method = "repeatedcv", number = CV,repeats = 1))
			}
			if(XGBT==T){
			  xgbtmodels<-train_xgbt(features=Complete_set,CPU=CPU,CV=CV,nrepeats=repeats,seeds=seed);
			  i<-i+1;
			  models[[i]]<-xgbtmodels;
			}else{
			 #place holder to avoid downstream problems
			  i<-i+1;
			  models[[i]] <- train(Class ~ ., data = Complete_set , metric = "ROC", method = "null" ,trControl=trainControl(classProbs = TRUE,summaryFunction = twoClassSummary,savePredictions = T, method = "repeatedcv", number = CV,repeats = 1))
			}
			return(	models)	}
			)


#modelss<-models;
DMELM<-models[[1]];
TRIBM<-models[[2]];

models<-DMELM;
	roc_trib<-lapply(models,function(x){ res<-predict(x,set2,type="prob");
		pROC::roc(set2$Class,res$E,direction=">") } );
	alle<-as.data.frame(roc_trib[[1]]$predictor/roc_trib[[1]]$predictor)
	alle$NE<-0
	colnames(alle)<-c("E","NE")
	ZRroc<-pROC::roc(set2$Class,alle$E,direction=">")
	pvalues_trib<-lapply(roc_trib,function(rroc){ pROC::roc.test(rroc,ZRroc) } )

	finaldmel<-list(DMELM,roc_trib,pvalues_trib)
	names(finaldmel)<-c("Dmel_models","roc_trib","pvalues_trib")
models<-TRIBM;

        roc_dmel<-lapply(models,function(x){ res<-predict(x,set1,type="prob");
                pROC::roc(set1$Class,res$E,direction=">") } );
        alle<-as.data.frame(roc_dmel[[1]]$predictor/roc_dmel[[1]]$predictor)
        alle$NE<-0
        colnames(alle)<-c("E","NE")
        ZRroc<-pROC::roc(set1$Class,alle$E,direction=">")
        pvalues_dmel<-lapply(roc_dmel,function(rroc){ pROC::roc.test(rroc,ZRroc) } )
        
	finaltrib<-list(TRIBM,roc_dmel,pvalues_dmel)
	
	rf.name="RF";
	xgbt.name="XGBT";
	names(finaltrib)<-c("Trib_models","roc_dmel","pvalues_dmel")
	names(finaltrib$pvalues_dmel)<-c(rf.name,xgbt.name)
	names(finaltrib$roc_dmel)<-c(rf.name,xgbt.name)
	names(finaldmel$pvalues_trib)<-c(rf.name,xgbt.name)
	names(finaldmel$roc_trib)<-c(rf.name,xgbt.name)
	if(XGBT==F){xgbt.name="null_model"};
	if(RF==F){rf.name="null_model"};

#Ploting results
#library("pROC")
plotname<-paste0(plot_prefix,"ROCs.svg")
svg(plotname,height=6,width=12)
par(mfrow=c(1,2))
size=0.9;
######################  Dmel models vs Trib
pvalue<-  finaldmel$pvalues_trib$RF$p.value
ROC_title<-paste0("ROC: ",set1.name," models vs ",set2.name)
pROC::plot.roc(finaldmel$roc_trib$RF,legend=F,color="black",print.auc=T,print.auc.cex=size,main=ROC_title)
text(x=0.25, y = 0.35, labels = paste0("P-value= ",signif(pvalue,digits=3)),cex=size)
pvalue<- finaldmel$pvalues_trib$XGBT$p.value

pROC::plot.roc(finaldmel$roc_trib$XGBT,legend=F,print.auc=T,print.auc.y=0.455 , add=T,print.auc.cex=size,col="red")
text(x=0.25, y = 0.30, labels = paste0("P-value= ",signif(pvalue,digits=3)),cex=size,col="red")
legend(x="topleft",c(rf.name,xgbt.name),col=c(1,2),lty=1)
title(adj=0,line=2.5,main="A)",cex.main=2)
#####################  Trib models vs dmel
pvalue<- finaltrib$pvalues_dmel$RF$p.value
ROC_title<-paste0("ROC: ",set2.name," models vs ",set1.name)
pROC::plot.roc(finaltrib$roc_dmel$RF ,legend=F,color="black",print.auc=T,print.auc.cex=size,main=ROC_title)
text(x=0.25, y = 0.35, labels = paste0("P-value= ",signif(pvalue,digits=3)),cex=size)
pvalue<- finaltrib$pvalues_dmel$XGBT$p.value 
pROC::plot.roc(finaltrib$roc_dmel$XGBT,legend=F,print.auc=T,print.auc.y=0.455 , add=T,print.auc.cex=size,col="red")
text(x=0.25, y = 0.30, labels = paste0("P-value= ",signif(pvalue,digits=3)),cex=size,col="red")
legend(x="topleft",c(rf.name,xgbt.name),col=c(1,2),lty=1)
title(adj=0,line=2.5,main="B)",cex.main=2)

dev.off()


####################  Precision Recall curves #### Set 1 vs Set 2
#svg("Insects_PRCs.svg",height=4,width=12)
plotname<-paste0(plot_prefix,"PRCs.svg")
svg(plotname,height=4,width=8)
par(mfrow=c(1,2))
if(var(finaldmel$roc_trib$RF$original.predictor)==0 || var(finaldmel$roc_trib$XGBT$original.predictor) ==0 ) {
  print("Skipping set1 vs set 2 PR curve, a set1 model has no variance for the prediction score") 
}else{
prcurve__1 = pr.curve(c(finaldmel$roc_trib$RF$original.predictor),weights.class0=(finaldmel$roc_trib$RF$original.response=='E')*1,curve=T) #rf
prcurve__2 = pr.curve(c(finaldmel$roc_trib$XGBT$original.predictor),weights.class0=(finaldmel$roc_trib$XGBT$original.response=='E')*1,curve=T) #xgbt
prcurve___ZR = pr.curve(c(finaldmel$pvalues_trib$RF$roc2$original.predictor[-1],0),weights.class0=(finaldmel$pvalues_trib$RF$roc2$original.response=='E')*1,curve=T)
PRC_title<-paste0("PRC: ",set1.name," models vs ",set2.name)
plot(prcurve__1,auc.main=F,main=PRC_title,color="black")
plot(prcurve__2,auc.main=F,add=T,color="red")
plot(prcurve___ZR,auc.main=F,add=T,color="gray",lty=2)
text(x=0.35, y = 0.95, labels =paste0(" RF  AUC= ", signif(prcurve__1$auc.integral,digits=3)))
text(x=0.35, y = 0.90, labels =paste0("XGBT AUC= ",  signif(prcurve__2$auc.integral,digits=3)),col="red")
text(x=0.35, y = 0.85, labels =paste0("ZR AUC= ", signif(prcurve___ZR$auc.integral,digits=3)),col="dark gray")
legend(x="topright",c(rf.name,xgbt.name),col=c(1,2),lty=1)
title(adj=0,line=1.5,main="A)",cex.main=2)
}
###################### Set 2 vs Set 1 #######################################
if(var(finaltrib$roc_dmel$RF$original.predictor)==0 || var(finaltrib$roc_dmel$XGBT$original.predictor) ==0 ) {
  print("Skipping set 2 vs set 1 PR curve, a set 2 model has no variance for the prediction score") 
}else{
prcurve__1 = pr.curve(c(finaltrib$roc_dmel$RF$original.predictor),weights.class0=(finaltrib$roc_dmel$RF$original.response=='E')*1,curve=T) #rf 
prcurve__2 = pr.curve(c(finaltrib$roc_dmel$XGBT$original.predictor),weights.class0=(finaltrib$roc_dmel$XGBT$original.response=='E')*1,curve=T) #xgbt
prcurve___ZR = pr.curve(c(finaltrib$pvalues_dmel$RF$roc2$original.predictor[-1],0),weights.class0=(finaltrib$pvalues_dmel$RF$roc2$original.response=='E')*1,curve=T)
PRC_title<-paste0("PRC: ",set2.name," models vs ",set1.name)
plot(prcurve__1,auc.main=F,main=PRC_title,color="black")
plot(prcurve__2,auc.main=F,add=T,color="red")
plot(prcurve___ZR,auc.main=F,add=T,color="gray",lty=2)
text(x=0.35, y = 0.95, labels =paste0(" RF  AUC= ", signif(prcurve__1$auc.integral,digits=3)))
text(x=0.35, y = 0.90, labels =paste0("XGBT AUC= ", signif(prcurve__2$auc.integral,digits=3)),col="red")
text(x=0.35, y = 0.85, labels =paste0("ZR AUC= ", signif(prcurve___ZR$auc.integral,digits=3)),col="dark gray")
legend(x="topright",c(rf.name,xgbt.name),col=c(1,2),lty=1)
title(adj=0,line=1.5,main="B)",cex.main=2)
}
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
if(xgbt.name=="null_model"){
  roc_xg_dmelmd10<-roc(DMELM[[2]]$pred[,"obs"], DMELM[[2]]$pred[,"E"],direction=">")
  roc_xg_tribmd10<-roc(TRIBM[[2]]$pred[,"obs"], TRIBM[[2]]$pred[,"E"],direction=">")
}else{
  roc_xg_dmelmd10<-roc (DMELM[[2]]$pred[DMELM[[2]]$pred$max_depth == 10,2], DMELM[[2]]$pred[DMELM[[2]]$pred$max_depth == 10,4],direction=">") 
  roc_xg_tribmd10<-roc (TRIBM[[2]]$pred[TRIBM[[2]]$pred$max_depth == 10,2], TRIBM[[2]]$pred[TRIBM[[2]]$pred$max_depth == 10,4],direction=">")
}
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
legend(x="bottomright",c(rf.name,xgbt.name),col=c(1,2),lty=1)

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
legend(x="bottomright",c(rf.name,xgbt.name),col=c(1,2),lty=1)
text(x=0.315, y = 0.4, labels = paste0("XGBT vs RF: P= ",signif(pvalue$p.value,digits=2)),cex=size)
text(x=0.315, y = 0.35, labels = paste0("XGBT vs ZR: P= ",signif(pvalue_cv_xg$p.value,digits=2)),cex=size)
text(x=0.315, y = 0.3, labels = paste0("RF vs ZR: P= ",signif(pvalue_cv_rf$p.value,digits=2)),cex=size)
dev.off()

	final_results<-(list(finaldmel,finaltrib))
	names(final_results)<-c("DMEL","TRIB")
}
