#' A feature extractor function.
#'
#' Input a DNA fasta file and it will output a vector of features for each sequence.
#'
#' @param PFAM_path Path to the PfamA database.
#' @param FASTA_PATH Path to the DNA fasta file.
#' @param CPU Number of threads to use. Each sequence will use one thread (2). 
#' @param LAMBDA Pseudo amino acid composition parameter. Sequences with less amino acids than LAMBDA will be removed (50). 
#' @param OMEGA Pseudo amino acid composition parameter (0.05). 
#' @return A data frame of features
#' @export

Extract<-function(FASTA_PATH="",AAfile="",PFAM_path="",LAMBDA=50,OMEGA=0.05,CPU=2,nuc_only=F,varGibbs_Model_path="/mnt/DATABASES/bin/VarGibbs-2.2/data/AOP-CMB.par"){
	#HMMSCAN and PFAM are needed in Calc_feats
	#If sequence length is less than LAMBDA, then it will be skipped; 
	SEQS<-seqinr::read.fasta(FASTA_PATH)
	if(AAfile==""){
			cl<-parallel::makeCluster(CPU,type="FORK",timeout=12000)
			on.exit(parallel::stopCluster(cl))
	#the parallel is importing global variables that are not being used by the loop. tried a lot of things, and still no success. 
			features_list<-parallel::parLapply(cl,SEQS,function(SEQ){ 
							Calc_feats(SEQ,PFAM_PATH=PFAM_path,LAMBDA=LAMBDA,OMEGA=OMEGA,nuc_only=nuc_only,varGibbs_Model_PATH=varGibbs_Model_path)})
			Features<-as.data.frame(data.table::rbindlist(features_list),fill=T,idcol=T)
			rm(features_list)
	}else{
		AAs<-seqinr::read.fasta(AAfile,seqtype="AA") # must have the same number of sequencies as SEQS
			n<-length(SEQS)
			cl<-parallel::makeCluster(CPU,type="FORK",timeout=12000) #200 minutes unlikely to happen, only in case of error. 
			on.exit(parallel::stopCluster(cl)) #garantee that the cluster will be stoped
			Features<-as.data.frame(data.table::rbindlist(
						parallel::parLapply(cl,1:n,function(N)
							Calc_feats(SEQS[[N]],AAs[[N]],PFAM_PATH=PFAM_path,LAMBDA=LAMBDA,OMEGA=OMEGA,nuc_only=nuc_only,varGibbs_Model_PATH=varGibbs_Model_path)),
						fill=T,idcol=T))
	}
	#=============CALCULATE FEATURES=============
	rownames(Features)<-Features$rownames
	Features$rownames<- NULL
#	Features<-Features[,-1]
	Features[is.na(Features)] <- F #binnary NA to False
	Features$Class<-"E" #change Class later
	gc();
	return(Features)
}
