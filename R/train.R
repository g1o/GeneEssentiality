#' Generate six RF models from two classes
#' 
#' It will return a list with 3 normal rf models and 3 rf models made with randomized labels.
#' The models 1 and 2 are trained without any subsampling.
#' The models 3 and 4 are trained using upsampling during the cross-validation (Subsampling During Resampling).
#' The models 5 and 6 are trained using downsampling during the cross-validation (Subsampling During Resampling).
#' 
#' @param data A feature data frame including the label of two Classes 
#' @param CPU Number of threads to use when training the model
#' @param trees Number of trees for the random forest (rf) model
#' @param CV Number of sets for the cross-validation
#' @param repeats 
#' @return A list of models
#' @importFrom caret train trainControl twoClassSummary downSample upSample
#' @export

train_rf<-function(features=Complete_set,CPU=2,trees=1000,CV=10,nrepeats=3 ,seeds=seed ,saveprediction=T){
set.seed(111)
	#Original training data
	mtries<-round(sqrt(length(features)))
	mtries<-c( mtries, mtries*2) 
control <- trainControl(method="repeatedcv", number=CV, repeats=nrepeats,
                         classProbs = TRUE,summaryFunction=twoClassSummary,
                         savePredictions = saveprediction ,seeds=seed)

	#TRAIN By maximizing the ROC METRIC
	original_fit<-train(Class ~ .,data=features, metric="ROC",method="ranger",
		tuneGrid=expand.grid(.mtry=mtries , .splitrule="gini",.min.node.size=1),
		trControl=control,num.trees=trees,num.threads=CPU,importance = 'impurity')

	##Cross validation with  Subsampling During Resampling # no diff
#	control$sampling<-"down"
#	down_fit<-train(Class ~ .,data=features, metric="ROC",method="ranger",
#		tuneGrid=expand.grid(.mtry=mtries , .splitrule="gini",.min.node.size=1),
#		trControl=control,num.trees=trees,num.threads=CPU,importance = 'impurity')

	models<-list(original_fit)
	names(models)<-c("original_fit")
	return(models)
}
