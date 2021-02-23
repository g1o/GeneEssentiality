#' Run a complete train and test using the data from the package. Reproduces the paper results 
#' @export
reproduce_results<-function(CPU=20,trees=1000,CV=10,repeats=3){
	seeds<-GeneEssentiality::seeds;
	models<-lapply(list(GeneEssentiality::drosophila_features, GeneEssentiality::tribolium_features), function(Complete_set){
			Select_Complete_set<- Select(data=Complete_set);
			dmel_rfmodels<-train_rf(features=Select_Complete_set,CPU=CPU,trees=trees,CV=CV,repeats=repeats,seeds=seeds);
			dmel_xmodels<-train_xgbt(features=Select_Complete_set,CPU=CPU,CV=CV,repeats=repeats,seeds=seeds);
			models<-list(dmel_rfmodels[[1]],dmel_rfmodels[[2]],dmel_xmodels)
			return(	models)	} )

dmel_importances<-lapply(models[[1]],function(x){varImp(x)});
trib_importances<-lapply(models[[2]],function(x){varImp(x)});

modelss<-models;
DMELM<-models[[1]];
TRIBM<-models[[2]];

models<-c(DMELM,TRIBM);
	result_noh_trib<-lapply(models,function(x){ res<-predict(x, GeneEssentiality::noh_trib, type="prob");
		pROC::roc(GeneEssentiality::noh_trib$Class,res$E,direction=">") } );
	alle<-as.data.frame(result_noh_trib[[1]]$predictor/result_noh_trib[[1]]$predictor)
	alle$NE<-0
	colnames(alle)<-c("E","NE")
	nullroc<-pROC::roc(GeneEssentiality::noh_trib$Class,alle$E,direction=">")
	pvalues_noh_trib<-lapply(result_noh_trib,function(rroc){ pROC::roc.test(rroc,nullroc) } )

	result_dmel<-lapply(models,function(x){ res<-predict(x,GeneEssentiality::drosophila_features,type="prob");
		pROC::roc(GeneEssentiality::drosophila_features$Class,res$E,direction=">") } );
	alle<-as.data.frame(result_dmel[[1]]$predictor/result_dmel[[1]]$predictor)
	alle$NE<-0
	colnames(alle)<-c("E","NE")
	nullroc<-pROC::roc(GeneEssentiality::drosophila_features$Class,alle$E,direction=">")
	pvalues_dmel<-lapply(result_dmel,function(rroc){ pROC::roc.test(rroc,nullroc) } )

	result_trib<-lapply(models,function(x){ res<-predict(x,GeneEssentiality::tribolium_features,type="prob");
		pROC::roc(GeneEssentiality::tribolium_features$Class,res$E,direction=">") } );
	alle<-as.data.frame(result_trib[[1]]$predictor/result_trib[[1]]$predictor)
	alle$NE<-0
	colnames(alle)<-c("E","NE")
	nullroc<-pROC::roc(GeneEssentiality::tribolium_features$Class,alle$E,direction=">")
	pvalues_trib<-lapply(result_trib,function(rroc){ pROC::roc.test(rroc,nullroc) } )

models<-TRIBM;
        result_noh_trib<-lapply(models,function(x){ res<-predict(x, GeneEssentiality::noh_trib, type="prob");
                pROC::roc(GeneEssentiality::noh_trib$Class,res$E,direction=">") } );
        alle<-as.data.frame(result_noh_trib[[1]]$predictor/result_noh_trib[[1]]$predictor)
        alle$NE<-0
        colnames(alle)<-c("E","NE")
        nullroc<-pROC::roc(GeneEssentiality::noh_trib$Class,alle$E,direction=">")
        pvalues_noh_trib<-lapply(result_noh_trib,function(rroc){ pROC::roc.test(rroc,nullroc) } )

        alle<-as.data.frame(result_dmel[[1]]$predictor/result_dmel[[1]]$predictor)
        alle$NE<-0
        colnames(alle)<-c("E","NE")
        nullroc<-pROC::roc(dmel$Class,alle$E,direction=">")
        result_dmel<-lapply(models,function(x){ res<-predict(x,GeneEssentiality::drosophila_features,type="prob");
                pROC::roc(GeneEssentiality::drosophila_features$Class,res$E,direction=">") } );
        pvalues_dmel<-lapply(result_dmel,function(rroc){ pROC::roc.test(rroc,nullroc) } )
        
        alle<-as.data.frame(result_trib[[1]]$predictor/result_trib[[1]]$predictor)
        alle$NE<-0
        colnames(alle)<-c("E","NE")
        nullroc<-pROC::roc(trib$Class,alle$E,direction=">")
        result_trib<-lapply(models,function(x){ res<-predict(x,GeneEssentiality::tribolium_features,type="prob");
                pROC::roc(GeneEssentiality::tribolium_features$Class,res$E,direction=">") } );
        pvalues_trib<-lapply(result_trib,function(rroc){ pROC::roc.test(rroc,nullroc) } )


	result_noh_trib<-lapply(models,function(x){ res<-predict(x, GeneEssentiality::noh_trib, type="prob");
		pROC::roc(GeneEssentiality::noh_trib$Class,res$E,direction=">") } );
	pvalues_noh_trib<-list(	pROC::roc.test(result_noh_trib[[1]],result_noh_trib[[2]]),
		pROC::roc.test(result_noh_trib[[3]],result_noh_trib[[4]]))
	result_dmel<-lapply(models,function(x){ res<-predict(x,GeneEssentiality::drosophila_features,type="prob");
		pROC::roc(GeneEssentiality::drosophila_features$Class,res$E,direction=">") } );
	pvalues_dmel<-list(	pROC::roc.test(result_dmel[[1]],result_dmel[[2]]),
		pROC::roc.test(result_dmel[[3]],result_dmel[[4]]) )
	result_trib<-lapply(models,function(x){ res<-predict(x,GeneEssentiality::tribolium_features,type="prob");
		pROC::roc(GeneEssentiality::tribolium_features$Class,res$E,direction=">") } );
	pvalues_trib<-list(	pROC::roc.test(result_trib[[1]],result_trib[[2]]),
			pROC::roc.test(result_trib[[3]],result_trib[[4]]) )

		final<-list(DMELM,TRIBM,result_noh_trib,pvalues_noh_trib,result_dmel,pvalues_dmel,result_trib,pvalues_trib)
		names(final)<-c("Dmel_models","Trib_models","result_noh_trib","pvalues_noh_trib","result_dmel","pvalues_dmel","result_trib","pvalues_trib")
		return(final)
}
