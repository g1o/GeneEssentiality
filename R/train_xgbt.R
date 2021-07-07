#' Generate 1 XGBT model 
#' 
#' @param data A feature data frame including the label of two Classes 
#' @param CPU Number of threads to use when training the model
#' @param CV Number of sets for the cross-validation
#' @param repeats 
#' @return A list of models
#' @importFrom caret train trainControl twoClassSummary downSample upSample
#' @export

train_xgbt<-function(features=data, seeds=seed,CPU=4, CV=10, nrepeats=3,saveprediction=T){
set.seed(111);
control <- trainControl(method="repeatedcv", number=CV, repeats=nrepeats,
                         classProbs = TRUE,summaryFunction=twoClassSummary,
                         savePredictions = saveprediction,seeds=seed)

	xgbt_m<-train(Class ~ .,data= features,metric="ROC",
		method="xgbTree", trControl=control,
		tuneGrid=expand.grid(nrounds = c(100,200,500),
		max_depth = c(4,10),colsample_bytree = 1,eta = 0.1,  gamma=1, 
		min_child_weight = 1,   subsample = 1),nthread=CPU)

	return(xgbt_m)
}
