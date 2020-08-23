#' Join features from two datasets, changing the Class to E and NE
#' @export
Join_datasets<-function(essential=data.frame(),nonessential=data.frame()){
	nonessential$Class<-factor("NE")
	essential$Class<-factor("E")
	Complete_set<-data.table::rbindlist(list(essential,nonessential),fill=T,idcol=T)
	Complete_set<-as.data.frame(Complete_set)
	rownames(Complete_set)<-t(Complete_set[,2])
	Complete_set<-Complete_set[,-c(1,2)]
	Complete_set[is.na(Complete_set)] <- F
	return(Complete_set)
}
