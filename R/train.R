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

train_rfmodels<-function(features=Complete_set,CPU=2,trees=1000,CV=10,repeats=3){
	#Original training data
	mtry<-round(sqrt(length(features)))
	tunegrid  <- expand.grid(.mtry = c(mtry*2,mtry*4) )
	control  <- trainControl(method="repeatedcv", number=CV, repeats=repeats,
		       classProbs = TRUE,summaryFunction=twoClassSummary,
		       savePredictions = TRUE)

	cls = parallel::makeCluster(CPU)
	doParallel::registerDoParallel(cls) #--- needed for random forest parallelization

	original_fit<-train(Class ~ .,data=features,metric="ROC",method="rf",trControl=control,tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=trees,importance=T)

	##Cross validation with  Subsampling During Resampling
	control <- trainControl(method="repeatedcv", number=CV, repeats=repeats,
			classProbs = TRUE,summaryFunction=twoClassSummary,
			savePredictions = TRUE,sampling="down")
	#TRAIN By maximizing the ROC METRIC
	down_fit<- train(Class ~ .,data=features,metric="ROC",method="rf",trControl=control,
			tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=trees,importance=T)
#
	control$sampling<-"up"
	up_fit  <- train(Class ~ .,data=features,metric="ROC",method="rf",trControl=control,
			tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=trees,importance=T)
#	#Randomize train labels to generate random models
	features$Class<-sample(features$Class); 
	up_randomfit  <- train(Class ~ .,data=features,metric="ROC",method="rf",trControl=control,
			tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=trees,importance=T)
#
	control$sampling<-"down"
	down_randomfit<- train(Class ~ .,data=features,metric="ROC",method="rf",trControl=control,
			tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=trees,importance=T)
#	##random 
	control  <- trainControl(method="repeatedcv", number=CV, repeats=repeats,
		       classProbs = TRUE,summaryFunction=twoClassSummary,
		       savePredictions = TRUE)
	random_fit<-train(Class ~ .,data=features,metric="ROC",method="rf",trControl=control,
			tuneGrid=tunegrid,prox=T,allowParallel=TRUE,ntree=trees,importance=T)
#
#	#
	models<-list(original_fit,random_fit,up_fit,up_randomfit,down_randomfit,down_fit)
	names(models)<-c("original_fit","random_fit","up_fit","up_randomfit","down_randomfit","down_fit")
	parallel::stopCluster(cls)
	return(models)
}
