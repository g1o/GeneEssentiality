#w' Compare two datasets
#'
#' Test using the provided models and features then plot the ROC and PRROC.
#' Used to compare two datasets with the same comparsion made for the insects, leaving one dataset out. The required column is the Class, and its values must be "E" or "NE".
#' Reditect to an object to save the models, ROCs results and Pvalues.  Plots are outputed to a file in the working dir.
#'
#' @param set1 list of models. Each model is an object given by the train function of the caret package.
#' @param set2 list of models. Each model is an object given by the train function of the caret package.
#' @param file_prefix Input a prefix for the filename of the plots. Default is "Insects_"
#' @param set1.name Name to be used for the title of the plots for the first set. Default is "Set1"
#' @param set2.name Name to be used for the title of the plots for the first set. Default is "Set2"
#' @param test_vs_ZR Use pROC::roc.test of the tested ROC-AUC against a model that classifies every sample as the same class with the same probability (Zero Rule: ROC-AUC of 1.5). If set to FALSE, then test the results in a all vs all. If there is a single model and test_vs_ZR = FALSE, it will revert to TRUE. If set to "both", show both. 
#' @return The AUCs of the tests
#' @export
VS_models_plot<-function(set1_model=list() , set2_model=list() , set1_data="", set2_data="", file_prefix="Compared_",set1.dataname="Set1", set2.dataname="Set2", set1.modelname="set1_model", set2.modelname="set2_model", test_vs_ZR=T){



	test_mem=0;
	if( length(set1_model)==1 & test_vs_ZR==F ){
		test_vs_ZR=T ;
			test_mem=1
	}

############Set 1 model Predictions 
	rocs<-lapply(set1_model,function(x){ 
			res <- predict(x,set2_data,type="prob");
			pROC::roc(set2_data$Class,res$E,direction=">") ;
			} );
	if(test_vs_ZR==T || test_vs_ZR=="both"){
#positive curve
		positive_null<-as.data.frame(rocs[[1]]$predictor/rocs[[1]]$predictor) ;
		positive_null$NE<-0;
		colnames(positive_null)<-c("E","NE");
		ZRroc<-pROC::roc(set2_data$Class,positive_null$E,direction=">");
		pvalues<-lapply(rocs,function(rroc){ pROC::roc.test(rroc,ZRroc) } );
	}else{
		pvalues<-list() #preallocate list
			for (i in 1:length(rocs) ) { 
				max=i;
				if(i==1){ next	}
				while (i>1){
					i=i-1;
					names_set1m<-paste0(max,"_vs_",i);
					pvalues[[names_set1m]]<-pROC::roc.test(rocs[[max]],rocs[[i]]) ;
				}
				i=max;
			}
	}
		finaldmel<-list(rocs,pvalues);
		names(finaldmel)<-c("roc_set2_prediction","pvalues_set2_prediction");
##set2 model

	if( length(set2_model)==1 & (test_vs_ZR==F || test_mem==1)){
		test_vs_ZR=T
	}
############Set 2 model Predictions 
	rocs<-lapply(set2_model,function(x){ 
			res<-predict(x,set1_data,type="prob");
			pROC::roc(set1_data$Class,res$E,direction=">") 
			} );
	if(test_vs_ZR==T){
#positive curve
		positive_null<-as.data.frame(rocs[[1]]$predictor/rocs[[1]]$predictor);
		positive_null$NE<-0;
		colnames(positive_null)<-c("E","NE");
		ZRroc<-pROC::roc(set1_data$Class,positive_null$E,direction=">");
		pvalues<-lapply(rocs,function(rroc){ pROC::roc.test(rroc,ZRroc) } );
	}else{
		pvalues<-list() #preallocate list
			for (i in 1:length(rocs) ) {
				max=i;
				if(i==1){ next  }
				while (i>1){
					i=i-1;
					names_set2m<-paste0(max,"_vs_",i);
					pvalues[[names_set2m]]<-pROC::roc.test(rocs[[max]],rocs[[i]]) ;
				}
				i=max;
			}
	}
		finaltrib<-list(rocs,pvalues);
		names(finaltrib)<-c("roc_set1_prediction","pvalues_set1_prediction");

#Ploting results
	plotname<-paste0(file_prefix,"ROCs.svg");
	svg(plotname,height=6,width=12);
	par(mfrow=c(1,2));
	size=0.9;
######################  Dmel models vs Trib
	i=1;
	ROC_title <-paste0("ROC: ",set1.modelname," models vs ",set2.dataname);
	pROC::plot.roc(finaldmel[[1]][[i]],legend=F,color="black",print.auc=F,main=ROC_title );
	if( length(set1_model)==1 & test_vs_ZR==F ){
		test_vs_ZR=T
	}
	AUCS_1 <- data.frame( ROC = 1:(length(set1_model) ) , PRC= 1:(length(set1_model)) );
	AUCS_1$ROC[i] <- signif(finaldmel[[1]][[i]]$auc,digits=3) ;
	AUCS_1$model_names <-c("");
	AUCS_1$model_names[i] <- names(finaldmel[[1]][i]) ;
	if (test_vs_ZR==T){
		pvalue <-  finaldmel[[2]][[i]]$p.value;
		legend(x=0.6, y =0.5 , legend = paste0("AUC: ", signif(finaldmel[[1]][[1]]$auc,digits=3)," | P= ",signif(pvalue,digits=2)),cex=size,col=1,lty=i,bty="n");
		AUCS_1$pvalue<-0;
		AUCS_1$pvalue[i]<-signif(pvalue,digits=3);
	}

####### Plot ROCs 
	if(length(finaldmel[[1]]) > 1){ 
		for (i in 2:length(finaldmel[[1]])){
			textY=0.5-(0.05*(i-1));
			pROC::plot.roc(finaldmel[[1]][[i]],legend=F,print.auc=F,add=T,lty=i,col=1);
			AUCS_1$ROC[i] <- signif(finaldmel[[1]][[i]]$auc,digits=3) ;
			AUCS_1$model_names[i] <- names(finaldmel[[1]][i]) ;

			if (test_vs_ZR==T){
				pvalue<-  finaldmel[[2]][[i]]$p.value;
				legend(x=0.6, y = textY , legend = paste0("AUC: ",AUCS_1$ROC[i]," | P= ",signif(pvalue,digits=2)),cex=size,col=1,lty=i,bty="n");
		                AUCS_1$pvalue[i]<-signif(pvalue,digits=3);
			}
		}
	}
####### Text: vs models p-value ##
        if (test_vs_ZR==F){
                textY=0.5;
                pvalues_vector<-signif ( sapply(finaldmel[[2]],'[[',"p.value") , digits=2) ;
                i=0;
                for (char in (strsplit(names(finaldmel[[2]]),'_'))){
                        textY=0.5-(0.05*i);
                        i=i+1;
                        legend(x=0.6, y = textY , legend = c("VS",paste0("P= ", pvalues_vector[i])),cex=size,col=1,lty=as.numeric(char[c(3,1)]),bty="n",ncol=2);
                }
	}
	
	if(length(names(finaldmel[[1]])) == 0){
		names(finaldmel[[1]])<-1:length(finaldmel[[1]])
	}
	legend(x="bottomright",legend=paste0(names(finaldmel[[1]]),": AUC = ", AUCS_1$ROC ),lty=1:length(finaldmel[[1]])) ;
	title(adj=0,line=2.5,main="A)",cex.main=2);

#####################  Trib models vs dmel
	keepi <- i ;
	i=1;
	if( length(set2_model)==1 & (test_vs_ZR==F || test_mem==1)){
		test_vs_ZR=T
	}
	ROC_title<-paste0("ROC: ",set2.modelname," models vs ",set1.dataname);
	pROC::plot.roc(finaltrib[[1]][[i]],legend=F,color="black",print.auc=F,main=ROC_title ,lwd=2) ;

	AUCS_2 <- data.frame( ROC = 1:(length(set2_model) ) , PRC= 1:(length(set2_model) ) ) ;
	AUCS_2$ROC[i ] <- signif(finaltrib[[1]][[i]]$auc,digits=3) ;
	AUCS_2$model_names[i ] <- names(finaltrib[[1]][i]) ;

	if (test_vs_ZR==T){
		pvalue<-  finaltrib[[2]][[i]]$p.value;
		legend(x=0.6, y =0.5 , legend = paste0("AUC: ", signif(finaltrib[[1]][[1]]$auc,digits=3)," | P= ",signif(pvalue,digits=2)),cex=size,col=1,lty=i,bty="n")
                AUCS_2$pvalue<-0;
                AUCS_2$pvalue[i]<-signif(pvalue,digits=3);
	}

####### Plot ROCs 
	if(length(finaltrib[[1]]) > 1){
		for (i in 2:length(finaltrib[[1]])){
			aucY=0.5-(0.05*(i-1))-.15;
			textY=0.5-(0.05*(i-1));
			pROC::plot.roc(finaltrib[[1]][[i]],legend=F,print.auc=F,add=T,lty=i,col=1,lwd=2);
			AUCS_2$ROC[i ] <-  signif(finaltrib[[1]][[i]]$auc,digits=3) ;
			AUCS_2$model_names[i ] <- names(finaltrib[[1]][i]) ;

			if (test_vs_ZR==T){
				pvalue<-  finaltrib[[2]][[i]]$p.value;
				legend(x=0.6, y = textY , legend = paste0("AUC: ", AUCS_2$ROC[i ] ," | P= ",signif(pvalue,digits=2)),cex=size,col=1,lty=i,bty="n");
	               		AUCS_2$pvalue[i]<-signif(pvalue,digits=3);
			}
		}
	}
####### Text: vs models p-value ##
	if (test_vs_ZR==F){
		textY=0.5;
		pvalues_vector<-signif ( sapply(finaltrib[[2]],'[[',"p.value") , digits=2) ;
		i=0;
		for (char in (strsplit(names(finaltrib[[2]]),'_'))){ 
			textY=0.5-(0.05*i);
			i=i+1;
			legend(x=0.6, y = textY , legend = c("VS",paste0("P= ", pvalues_vector[i])),cex=size,col=1,lty=as.numeric(char[c(3,1)]),bty="n",ncol=2);
		}
	}

	if(length(names(finaltrib[[1]])) == 0){
		names(finaltrib[[1]])<-1:length(finaltrib[[1]])
	}
	legend(x="bottomright",legend=paste0(names(finaltrib[[1]]),": AUC = ", AUCS_2$ROC ),lty=1:length(finaltrib[[1]]));
	title(adj=0,line=2.5,main="B)",cex.main=2);
	dev.off();

####################  Precision Recall curves #### Set 1 vs Set 2

  plotname<-paste0(file_prefix,"PRCs.svg") ;
    svg(plotname,height=6,width=12) ;
    par(mfrow=c(1,2)) ;
    if(var(finaldmel[[1]][[1]]$original.predictor)==0 )  {
      print("Skipping set1 vs set 2 PR curve, a set1 model has no variance for the prediction score") 
    }else{
      PRC_title<-paste0("PRC: ",set1.modelname," models vs ",set2.dataname);
      i=1;
      prcurve = PRROC::pr.curve(c(finaldmel[[1]][[i]]$original.predictor),weights.class0=(finaldmel[[1]][[i]]$original.response=='E')*1,curve=T) ;
      plot(prcurve,auc.main=F,main=PRC_title,color=1,lwd=2);
  
        AUCS_1$PRC[i] <- signif(prcurve$auc.integral,digits=3)

      if(length(finaldmel[[1]]) > 1){
        for (i in 2:length(finaldmel[[1]])){
          textY=0.95-(0.05*(i-1));
          prcurve = PRROC::pr.curve(c(finaldmel[[1]][[i]]$original.predictor),weights.class0=(finaldmel[[1]][[i]]$original.response=='E')*1,curve=T) ;
          plot(prcurve,auc.main=F,add=T,color=1,lty=i,lwd=2);
     	  AUCS_1$PRC[i] <- signif(prcurve$auc.integral,digits=3);
        };
      }else{
        i=i+1;
        textY=0.95-(0.05*(i-1))
      };
##########  ZR PRC  ############
      prcurve = PRROC::pr.curve(c(finaldmel[[1]][[1]]$original.predictor[-1]/finaldmel[[1]][[1]]$original.predictor[-1],0),weights.class0=(finaldmel[[1]][[1]]$original.response=='E')*1,curve=T);
      plot(prcurve,auc.main=F,add=T,color="gray",lty=2,lwd=2);
      legend(x="bottomright", legend= paste0(c(names(finaldmel[[1]]),"ZR"), ": AUC = ", c(AUCS_1$PRC, signif(prcurve$auc.integral,digits=3) ) ),col=c(rep(1,length(finaldmel[[1]])),"gray"),lty=c(1:length(finaldmel[[1]]),2) );
      title(adj=0,line=1.5,main="A)",cex.main=2);
    }
###################### Set 2 vs Set 1 #######################################
  if(var(finaltrib[[1]][[1]]$original.predictor)==0 ) {
    print("Skipping set1 vs set 2 PR curve, a set1 model has no variance for the prediction score")
  }else{
    PRC_title<-paste0("PRC: ",set2.modelname," models vs ",set1.dataname);
    i=1;
    prcurve = PRROC::pr.curve(c(finaltrib[[1]][[i]]$original.predictor),weights.class0=(finaltrib[[1]][[i]]$original.response=='E')*1,curve=T) ;
    plot(prcurve,auc.main=F,main=PRC_title,color=1,lwd=2);

       AUCS_2$PRC[i] <- signif(prcurve$auc.integral,digits=3)

    if(length(finaltrib[[1]]) > 1){
      for (i in 2:length(finaltrib[[1]])){
        textY=0.95-(0.05*(i-1));
        prcurve = PRROC::pr.curve(c(finaltrib[[1]][[i]]$original.predictor),weights.class0=(finaltrib[[1]][[i]]$original.response=='E')*1,curve=T) ;
        plot(prcurve,auc.main=F,add=T,color=1,lty=i,lwd=2);

        AUCS_2$PRC[i] <- signif(prcurve$auc.integral,digits=3);
      };
    }else{
      i=i+1;
      textY=0.95-(0.05*(i-1))
    };
    prcurve = PRROC::pr.curve(c(finaltrib[[1]][[1]]$original.predictor[-1]/finaltrib[[1]][[1]]$original.predictor[-1],0),weights.class0=(finaltrib[[1]][[1]]$original.response=='E')*1,curve=T);
    plot(prcurve,auc.main=F,add=T,color="gray",lty=2,lwd=2);
    legend(x="bottomright",y= textY ,legend= paste0(c(names(finaltrib[[1]]),"ZR"), ": AUC = ", c(AUCS_2$PRC,signif(prcurve$auc.integral,digits=3)) ),col=c(rep(1,length(finaltrib[[1]])),"gray"),lty=c(1:length(finaltrib[[1]]),2) );
    title(adj=0,line=1.5,main="B)",cex.main=2);
  }
  dev.off()

  AUCS_1$testSet<-set2.dataname;
  AUCS_2$testSet<-set1.dataname;
  AUCS<-rbind(AUCS_1,AUCS_2);
  return(AUCS)
}
