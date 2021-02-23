#' Remove features of little importance
#'
#' Remove features with zero variance and logical values that occur at a single point.
#'
#' @param data A feature data frame including the Class label
#' @return A reduced set of features 
#' @export
Select<-function(data=Complete_set,remove_logical=T){
`%notin%`<-Negate(`%in%`)
train<-data
train$Class<-NULL
logica<-sapply(train,is.logical)
if (remove_logical){
	train<-train[,!logica]
}
else{
#Remove boolean values that exists in only one gene
	s<-sapply(train[,logica]*1,sum)
	train<-train[,names(train) %notin% names(s[s==1])]
}
##Removes variables that have Zero variance
	logica<-sapply(train,is.logical)
	zerovar<-sapply(train[,!logica],var)
	zerovar<-sapply(train,var)
	train<-train[,names(train) %notin% names(zerovar[zerovar==0])]
	train$Class<-data$Class
return(train)
}
