
model_test<-function(models=models,test=features){
	res<-predict(models,test, type="prob")
	rocs<-lapply(res,function(x){pROC::roc(test$Class,x$E,direction=">")})
	pvalues<-lapply(c(1,3),function(x){pROC::roc.test(rocs[[x]],rocs[[x+1]])})
	results<-list(rocs,pvalues,res)
	names(results)<-c("rocs","pvalues","res")
	return(results)
}
