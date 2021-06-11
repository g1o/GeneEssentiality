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
##example Extract<-function(FASTA_PATH="/mnt/DATABASES/essenciais/drosophila/fastas/essential_complete-curated.fasta",PFAM_path="/home/programs/DATABASES/PFAM/Pfam-A.hmm",LAMBDA=50,OMEGA=0.05,CPU=2){
Extract<-function(FASTA_PATH="",AAfile="",PFAM_path="",LAMBDA=50,OMEGA=0.05,CPU=2,nuc_only=F){
	#HMMSCAN and PFAM are needed in Calc_feats
	#If sequence length is less than LAMBDA, then it will be skipped; 
	SEQS<-seqinr::read.fasta(FASTA_PATH)
	n<-length(SEQS)
	if(AAfile==""){
			cl<-parallel::makeCluster(CPU,type="FORK")
			Features<-as.data.frame(data.table::rbindlist(
						parallel::parLapply(cl,1:n,function(N)
							Calc_feats(SEQS[[N]],PFAM_PATH=PFAM_path,LAMBDA=LAMBDA,OMEGA=OMEGA,nuc_only=nuc_only)),

						fill=T,idcol=T))
	}else{
		AAs<-seqinr::read.fasta(AAfile,seqtype="AA") # must have the same number of sequencies as SEQS
			n<-length(SEQS)
			cl<-parallel::makeCluster(CPU,type="FORK")
			Features<-as.data.frame(data.table::rbindlist(
						parallel::parLapply(cl,1:n,function(N)
							Calc_feats(SEQS[[N]],AAs[[N]],PFAM_PATH=PFAM_path,LAMBDA=LAMBDA,OMEGA=OMEGA,nuc_only=nuc_only)),
						fill=T,idcol=T))
	}
	#=============CALCULATE FEATURES=============
#	rownames(Features)<-t(Features[,1]) #given when joining both classes
	Features[is.na(Features)] <- F #binnary NA to False
	Features$Class<-"E" #change Class later
	parallel::stopCluster(cl)
	gc();
	return(Features)
}
