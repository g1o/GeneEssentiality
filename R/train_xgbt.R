#' Generate 1 SVM and 1 XGB models 
#' 
#' @param data A feature data frame including the label of two Classes 
#' @param CPU Number of threads to use when training the model
#' @param CV Number of sets for the cross-validation
#' @param repeats 
#' @return A list of models
#' @importFrom caret train trainControl twoClassSummary downSample upSample
#' @export

train_xgbt_svm<-function(features=data, seeds=seeds,CPU=2, CV=10, repeats=3){
set.seed(111);
control <- trainControl(method="repeatedcv", number=CV, repeats=repeats,
                         classProbs = TRUE,summaryFunction=twoClassSummary,
                         savePredictions = FALSE,seeds=seeds)

	xgbt_m<-train(Class ~ .,data= features,metric="ROC",
		method="xgbTree", trControl=control,
		tuneGrid=expand.grid(nrounds = c(100,200,500),
		max_depth = c(4,10),colsample_bytree = 1,eta = 0.1,  gamma=1, 
		min_child_weight = 1,   subsample = 1),nthread=CPU)

#       cls = parallel::makeCluster(CPU)
#       doParallel::registerDoParallel(cls) #--- needed for random forest parallelization
#	svm_m<-train(Class ~ .,data= features,metric="ROC",method="svmRadial",
#		preProc = c("center", "scale"), trControl=control,
#		tuneGrid=expand.grid(.sigma=c(2^-8,2^-9,2^-10),.C=c(0.25,0.5)) )
#       parallel::stopCluster(cls)

#	return(list(xgbt_m,svm_m))
	return(xgbt_m)
}
