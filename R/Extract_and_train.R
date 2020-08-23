#' Run a complete train using data
RUN_all<-function(PFAM_path="/mnt/DATABASES/PFAM/Pfam-A.hmm",GE_PATH=system.file("extdata", "GE.fasta", package = "GeneEssentiality"),GNE_PATH=system.file("extdata", "GNE.fasta", package = "GeneEssentiality"),LAMBDA=50,OMEGA=0.05,CPU=1,trees=10,CV=10,repeats=3){
	ge <-Extract( PFAM_path=PFAM_path, FASTA_PATH=GE_PATH , LAMBDA=LAMBDA, OMEGA=OMEGA, CPU=CPU );
	gne<-Extract( PFAM_path=PFAM_path, FASTA_PATH=GNE_PATH, LAMBDA=LAMBDA, OMEGA=OMEGA, CPU=CPU );
	Complete_set       <- Join_datasets(essential=ge,nonessential=gne);
	Select_Complete_set<- Select(data=Complete_set);
	models<-train_rfmodels(features=Select_Complete_set,CPU=CPU,trees=trees,CV=CV,repeats=repeats);
	models
}
