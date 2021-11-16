library(GeneEssentiality)
library(caret)

load("/mnt/DATABASES/essenciais/MAJOR_REVISION/WORK/importance_united.rda")
if(!(exists("importance_united"))){
	load("/mnt/DATABASES/essenciais/MAJOR_REVISION/WORK/United_model_nzv.rda")
		if(!(exists("United_model_nzv0"))){
#### colocar lá 
noh_dmel_ids <-read.table("/mnt/DATABASES/essenciais/MAJOR_REVISION/Orthofinder_no_homology/dmel_ids")
noh_trib_ids <-read.table("/mnt/DATABASES/essenciais/MAJOR_REVISION/Orthofinder_no_homology/trib_ids")

noh_trib_features <- tribolium_features[  ( rownames(tribolium_features ) %in% noh_trib_ids[,1] ) , ]
noh_dmel_features <- drosophila_features[ ( rownames(drosophila_features) %in% noh_dmel_ids[,1] ) , ]

#### prepare data ### dmel load all genes ##
function(Experiment_1_labels){
##this will plot results and return a list with the models and results from comparisons
	ES_features   <- Extract(ES.fasta.gz , ES.faa.gz, LAMBDA=50) #likely wont work as too many sequences
	Cell_features <- Extract(Cell_level.fasta.gz, Cell_level.faa.gz , LAMBDA=50) 
##nearZeroVar metrics (all_nzv) calculated after using the features from all longest CDS from both Drosophila melanogaster (6.32) and Tribolium castaneum genes (OGS3) 
	ES_features   <-   ES_features[, c(!(all_nzv$percentUnique < 10 & all_nzv$freqRatio > 95/5)) ] 
	Cell_features <- Cell_features[, c(!(all_nzv$percentUnique < 10 & all_nzv$freqRatio > 95/5)) ]   
	Experiment_1$Cell_level <- (train_rf(features= Cell_features ,seeds=seed,nrepeats=3,CV=10,trees=1000,CPU=20))[[1]] ;
	Experiment_1$ES         <- (train_rf(features= ES_features ,seeds=seed,nrepeats=3,CV=10,trees=1000,CPU=20))[[1]]
	Experiment_1$DMEL       <- (train_rf(features= drosophila_features_nzv ,seeds=seed,nrepeats=3,CV=10,trees=1000,CPU=20))[[1]]

	names(experimento_1_modelos)<-c("Cell-level","ES","DMEL")

	AUCS<-VS_models_plot(set1_model=  Experiment_1 , set2_model= Experiment_1 , set2_data= TRIB_nzv ,set1_data= NOH_TRIB_nzv , set1.dataname= "NOH TRIB", set2.dataname="TRIB", set1.modelname="D. melanogaster", set2.modelname="D. melanogaster", file_prefix= "ES_cell_full_comparision_" ,test_vs_ZR=F)
	VS_models_plot(set1_model=  Experiment_1 , set2_model= Experiment_1 , set2_data= TRIB_nzv ,set1_data= NOH_TRIB_nzv , set1.dataname= "NOH TRIB", set2.dataname="TRIB", set1.modelname="D. melanogaster", set2.modelname="D. melanogaster", file_prefix= "ES_cell_full_comparision_" ,test_vs_ZR=T)
	return(list(Experiment_1 , AUCS)) 
}

function(Experiment_2_feature_selection){
	United_model_nzv0<-train_rf(features=rbind(drosophila_features_nzv[,1:8199],tribolium_features_nzv[,1:8199]),seeds=seed,nrepeats=3,CV=10,trees=1000,CPU=15);
	importance_united <- varImp(United_model_nzv0[[1]]) ;
	models_norna<-list() ;
	for(IMPORTANCE in c(0,5,10:20,25,30) ){
		imp<-c(importance_united$importance$Overall>IMPORTANCE,TRUE) ;
		name1<-paste0(IMPORTANCE,"_Dmel") ;
		name2<-paste0(IMPORTANCE,"_Trib") ;

		models_norna[[paste0("rf_norna_",name1)]] <-(train_rf(features= drosophila_features_nzv[,c(imp,rep(F,13))],seeds=seed,nrepeats=3,CV=10,trees=1000,CPU=15))[[1]] ;
		models_norna[[paste0("rf_norna_",name2)]] <- (train_rf(features= tribolium_features_nzv[,c(imp,rep(F,13))],seeds=seed,nrepeats=3,CV=10,trees=1000,CPU=15))[[1]] ;
		models_norna[[paste0("xgbt_norna_",name1)]] <-(train_xgbt(features= drosophila_features_nzv[,c(imp,rep(F,13))],seeds=seed,nrepeats=3,CV=10,CPU=15))[[1]] ; 
		models_norna[[paste0("xgbt_norna_",name2)]] <- (train_xgbt(features= tribolium_features_nzv[,c(imp,rep(F,13))],seeds=seed,nrepeats=3,CV=10,CPU=15))[[1]] ; 
	}

	AUCS_noRNA<-VS_models_plot (set1_model=  models_norna[grep("Dmel", names(models_norna))] ,
				set2_model=  models_norna[grep("Trib", names(models_norna))] , 
				set2_data= tribolium_features ,set1_data= drosophila_features , 
				set1.dataname= "dmel", set2.dataname="trib", 
				set1.modelname="Dmel_selection", set2.modelname="Trib_selection", 
				file_prefix="Selection" ) ;

	Selecion_AUCS <- cbind(AUCS_noRNA,do.call(rbind, strsplit(AUCS_noRNA$model_names, '_') )) ;
	names(Selecion_AUCS)<-c("ROC" ,"PRC" , "model_names" , "pvalue", "testSet" , "Algorithm" , "norna", "Importance" , "Trained in") ;
	Selecion_AUCS$Importance<-as.numeric(as.character(Selecion_AUCS$Importance)) ;
	Selecion_AUCS <- Selecion_AUCS[grep("[.]",Selecion_AUCS$model_names ,invert=T),];
	Selecion_AUCS$group<-paste(Selecion_AUCS$`Trained in`,Selecion_AUCS$Algorithm) ;
	DT<-data.table::data.table(Selecion_AUCS);
	DT[, .SD[which.max(ROC)], by = group] ; # max in rocs
	DT[, .SD[which.max(PRC)], by = group] ; # max in prs
	svg("AUCS_ROC_PRC_selection.svg" );
	p1<-ggplot(Selecion_AUCS ,aes(x= Importance, y=ROC, col=`Trained in`,linetype=Algorithm ) ) + geom_line() +
		scale_x_continuous(breaks=c(0,5,10:20,25,30,(DT[, .SD[which.max(ROC)], by = group]$Importance)) ,name="Importance cut-off", expand = c(0, 0)) +
		ylab(label="AUC-ROC") +
		ggtitle("AUCs-ROCs from features selected by importance over each cut-off") +
		geom_vline(xintercept= (DT[, .SD[which.max(ROC)], by = group]$Importance),color="gray",linetype=1,alpha=1/4)  +
		scale_colour_grey() +  theme_classic();
	p2<-ggplot(Selecion_AUCS,aes(x= Importance, y=PRC, col=`Trained in`,linetype=Algorithm ) ) + geom_line() +
		scale_x_continuous(breaks=c(0,5,10:20,25,30,DT[, .SD[which.max(PRC)], by = group]$Importance) ,name="Importance cut-off", expand = c(0, 0)) + 
		geom_vline(xintercept= (DT[, .SD[which.max(PRC)], by = group]$Importance),color="gray",linetype=1,alpha=1/4)  +
		ggtitle("AUCs-PRCs from features selected by importance over each cut-off") +
		ylab(label="AUC-PRC") + scale_colour_grey() +  theme_classic();
	legend <- cowplot::get_legend(	p1 + guides(color = guide_legend(nrow = 1)) +
			theme(legend.position = "bottom") );
	cowplot::plot_grid(p1+ theme(legend.position="none"),legend,
			p2+ theme(legend.position="none"),
			labels = c('A','' ,'B'),
			label_size = 12,
			hjust = -1 , nrow=3,rel_heights = c(1, .1, 1) );
	dev.off();

}

function (experiment3_AA_vs_NT){
	tribolium_features_nzv_nt <- (tribolium_features_nzv[,c(rep(T,5093),rep(F,dim(tribolium_features_nzv)[2]-5093-14),T,rep(F,13))])
	tribolium_features_nzv_aa <- (tribolium_features_nzv[,c(rep(F,5093),rep(T,dim(tribolium_features_nzv)[2]-5093-14),T,rep(F,13))])
	drosophila_features_nzv_nt<-(drosophila_features_nzv[,c(rep(T,5093),rep(F,dim(tribolium_features_nzv)[2]-5093-14),T,rep(F,13))])
	drosophila_features_nzv_aa<-(drosophila_features_nzv[,c(rep(F,5093),rep(T,dim(tribolium_features_nzv)[2]-5093-14),T,rep(F,13))])
#models_campos_our_RF_noselection<-list(Campos2019_OGEEv2 = Campos2019_nzv$RF_model, Campos2020 = Campos2020_nzv$RF_model, Our =  models_norna$rf_norna_0_Dmel)
	aa_nt_models<-list()
	aa_nt_models[["rf_nt_dmel"]] <-(train_rf(features= drosophila_features_nzv_nt,seeds=seed,nrepeats=3,CV=10,trees=1000,CPU=15))[[1]]
	aa_nt_models[["rf_aa_dmel"]] <-(train_rf(features= drosophila_features_nzv_aa,seeds=seed,nrepeats=3,CV=10,trees=1000,CPU=15))[[1]]
	aa_nt_models[["rf_nt_trib"]] <- (train_rf(features= tribolium_features_nzv_nt,seeds=seed,nrepeats=3,CV=10,trees=1000,CPU=15))[[1]]
	aa_nt_models[["rf_aa_trib"]] <- (train_rf(features= tribolium_features_nzv_aa,seeds=seed,nrepeats=3,CV=10,trees=1000,CPU=15))[[1]]
	aa_nt_models[["xgbt_nt_dmel"]] <-(train_xgbt(features= drosophila_features_nzv_nt,seeds=seed,nrepeats=3,CV=10,CPU=15))[[1]]
	aa_nt_models[["xgbt_aa_dmel"]] <-(train_xgbt(features= drosophila_features_nzv_aa,seeds=seed,nrepeats=3,CV=10,CPU=15))[[1]]
	aa_nt_models[["xgbt_nt_trib"]] <- (train_xgbt(features= tribolium_features_nzv_nt,seeds=seed,nrepeats=3,CV=10,CPU=15))[[1]]
	aa_nt_models[["xgbt_aa_trib"]] <- (train_xgbt(features= tribolium_features_nzv_aa,seeds=seed,nrepeats=3,CV=10,CPU=15))[[1]]

## AA and NT ## 
	experimento_3_modelos<-list() ;
	experimento_3_modelos$dmel<-( c(aa_nt_models[grep("dmel", names(aa_nt_models))],models_norna[grep("_norna_0_Dmel",names(models_norna))]) );
	experimento_3_modelos$trib<-( c(aa_nt_models[grep("trib", names(aa_nt_models))],models_norna[grep("_norna_1_Trib",names(models_norna))]) );
	names(experimento_3_modelos$trib);
	names(experimento_3_modelos$trib)<-c("NT rf","AA rf","NT xgbt","AA xgbt","NT+AA rf","NT+AA xgbt");
	names(experimento_3_modelos$dmel)<-c("NT rf","AA rf","NT xgbt","AA xgbt","NT+AA rf","NT+AA xgbt");

	VS_models_plot(set1_model=  experimento_3_modelos$dmel  ,set2_model=  experimento_3_modelos$trib , set2_data= tribolium_features ,set1_data= drosophila_features , set1.dataname= "dmel", set2.dataname="trib", set1.modelname="Dmel", set2.modelname="Trib", file_prefix="AA_nt" ) ;
	AUCS_aa_nt_full <- VS_models_plot(set1_model=  experimento_3_modelos$dmel , set2_model=  experimento_3_modelos$trib , set2_data= tribolium_features ,set1_data= drosophila_features , set1.dataname= "dmel", set2.dataname="trib", set1.modelname="Dmel", set2.modelname="Trib", file_prefix="AA_nt2" , test_vs_ZR=F ) ;
# NOH ##
	VS_models_plot(set1_model=  experimento_3_modelos$dmel ,set2_model= experimento_3_modelos$trib ,  set2_data= noh_trib_features ,set1_data= noh_dmel_features , set1.dataname= "noh dmel", set2.dataname="noh trib", set1.modelname="Dmel", set2.modelname="Trib", file_prefix="AA_nt_noh"  );
	AUCS_aa_nt_full_noh<-	     VS_models_plot(set1_model=  experimento_3_modelos$dmel ,set2_model= experimento_3_modelos$trib ,  set2_data= noh_trib_features ,set1_data= noh_dmel_features , set1.dataname= "noh dmel", set2.dataname="noh trib", set1.modelname="Dmel", set2.modelname="Trib", file_prefix="AA_nt2_noh" , test_vs_ZR=F );

###plot aAUCS AA NT

	AUCS <- cbind(AUCS_aa_nt_full[[1]] ,do.call(rbind, strsplit( AUCS_aa_nt_full[[1]]$model_names, ' ') ));
	names(AUCS)<-c("ROC" ,"PRC" , "model_names" , "testSet","Trained in","algo_long","ZR","Type","Algorithm");

	AUC_noh <- cbind(AUCS_aa_nt_full_noh[[1]] ,do.call(rbind, strsplit(AUCS_aa_nt_full_noh[[1]]$model_names, ' ') ));
	names(AUC_noh)<-c("ROC" ,"PRC" , "model_names" , "testSet","Trained in","algo_long","ZR","Type","Algorithm");

	p1 <-ggplot(AUCS ,aes(x= `Trained in` , y=ROC, col= Algorithm, shape= Type ) ) + geom_point(size=4, position=position_dodge(width=0.5)) +
		ylab(label="AUC-ROC") +
		ggtitle("ROC- AA vs NT: Complete") +
		scale_y_continuous(limits=c(0.5,0.8)) +
		scale_colour_grey() +  theme_classic();
	p2 <-ggplot(AUCS ,aes(x= `Trained in` , y=PRC, col= Algorithm, shape= Type) ) + geom_point(size=4, position=position_dodge(width=0.5)) +
		ggtitle("PR- AA vs NT: Complete") +
		geom_segment(aes(x=1.5, xend=2.5,y = ZR[10] ,yend=ZR[10]) , linetype=2 )+
		geom_segment(aes(x=0.5, xend=1.5,y = ZR[1] ,yend=ZR[1]  ) , linetype=2 )+
		ylab(label="AUC-PRC ") + scale_colour_grey() +  theme_classic();
	legend <- cowplot::get_legend(
			p1 + guides(color = guide_legend(nrow = 1)) +
			theme(legend.position = "bottom"));

	prow<-cowplot::plot_grid(p1+ theme(legend.position="none"),
			p2+ theme(legend.position="none"),
			labels = c('A', 'B'),
			label_size = 12,
			align = 'vh',  hjust = -1 , nrow=1 );
	pab<-cowplot::plot_grid(prow, legend, rel_heights = c(1, .1)  , ncol =1 );

	p1 <-ggplot(AUC_noh ,aes(x= `Trained in` , y=ROC, col= Algorithm, shape= Type) ) + geom_point(size=4,position=position_dodge(width=0.5) ) +
		ggtitle("ROC- AA vs NT: NOH") +
		ylab(label="AUC-ROC") +
		scale_y_continuous(limits=c(0.5,0.8)) +
		scale_colour_grey() +  theme_classic() ;

	p2 <-ggplot(AUC_noh ,aes(x= `Trained in` , y=PRC, col= Algorithm, shape= Type ) ) + geom_point(size=4 ,position=position_dodge(width=0.5)) +
		ggtitle("PR- AA vs NT: NOH") +
		geom_segment(aes(x=1.5, xend=2.5,y = ZR[10], yend=ZR[10]), linetype=2  )+
		geom_segment(aes(x=0.5, xend=1.5,y = ZR[1] , yend=ZR[1] ), linetype=2  )+
		ylab(label="AUC-PRC") +
		scale_colour_grey() +  theme_classic();

	prow<-cowplot::plot_grid(p1+ theme(legend.position="none"),
			p2+ theme(legend.position="none"),
			labels = c('C', 'D'),
			label_size = 12,
			align = 'vh',  hjust = -1 , nrow=1);
	pcd<-cowplot::plot_grid(prow, rel_heights = c(1, .1)  , ncol =1 );
	svg("AUCS_ROC_PRC_AA_NT.svg");
	cowplot::plot_grid(pab,pcd,nrow=2);
	dev.off();

}

function(experiment4_classifiers_and_extrinsic) {
#../../GeneEssentiality/data/SVM_seeds.rda
	set1<-drosophila_features_nzv[,c(imp,rep(F,13))];
	set2<-tribolium_features_nzv[,c(imp,rep(F,13))];
	IMPORTANCE=12;
	imp<-c(importance_united$importance$Overall>IMPORTANCE,TRUE);
	name1<-paste0(IMPORTANCE,"_Dmel");
	name2<-paste0(IMPORTANCE,"_Trib");
	set1_rna<-drosophila_features_nzv[,c(imp,rep(T,13))];
	set2_rna<-tribolium_features_nzv[,c(imp,rep(T,13))];

	control <- trainControl(method = "repeatedcv", number = 10, repeats = 3,
			classProbs = TRUE, summaryFunction = twoClassSummary,
			savePredictions = "final", seeds = SVM_seeds, preProcOptions = NULL);

	library(doParallel);
	cl <- makeCluster(10, type = "FORK");
	registerDoParallel(cl);

	#selected_models_norna$polynomial_norna_12_Trib



	Dmel_svm_12 <- list();
	Trib_svm_12 <- list();

	selected_models_norna$radial_norna_12_Trib

	Dmel_svm_12$radial_13 <- train(form = Class ~ ., data = set1, metric = "ROC", method = "ranger",
			trControl = control,
			tuneGrid = expand.grid(
				sigma = 2^c(-25,-20, -15 ),
				C = 2^c(0, 0.5, 0.75, 1, 1.5)),
			prox = T, allowParallel = TRUE, preProcess = c("center",  "scale"));

	Dmel_svm_12$radial_13 <- train(form = Class ~ ., data = set1, metric = "ROC", method = "svmRadial",
			trControl = control,
			tuneGrid = expand.grid(
				sigma = 2^c(-25,-20, -15 ),
				C = 2^c(0, 0.5, 0.75, 1, 1.5)),
			prox = T, allowParallel = TRUE, preProcess = c("center",  "scale"));

	Trib_svm_12$radial_13 <- train(form = Class ~ ., data = set2, metric = "ROC", method = "svmRadial",
			trControl = control,
			tuneGrid = expand.grid(
				sigma = 2^c(-25,-20, -15 ),
				C = 2^c(0, 0.5, 0.75, 1, 1.5)),
			prox = T, allowParallel = TRUE, preProcess = c("center",  "scale"));


	Dmel_svm_12$polynomial_13 <- train(form = Class ~ ., data = set1, metric = "ROC", method = "svmPoly",
			trControl = control, tuneGrid = expand.grid(degree = (2:4),
				C = 2^c(0, 1,  2), scale = c(0.1,0.2)),
			prox = T, allowParallel = TRUE, preProcess = c("center",      "scale"));
	Trib_svm_12$polynomial_13 <- train(form = Class ~ ., data = set2, metric = "ROC", method = "svmPoly",
			trControl = control, tuneGrid = expand.grid(degree = (2:4),
				C = 2^c(0, 1,  2), scale = c(0.1, 0.2)),
			prox = T, allowParallel = TRUE, preProcess = c("center",      "scale"));

		stopCluster(cl)



#after looking at the plot, over 12 was decided
	models_rna<-list() ;
	models_rna[[paste0("rf_rna_",name1)]] <-(train_rf(features= set1_rna ,seeds=seed,nrepeats=3,CV=10,trees=1000,CPU=5))[[1]];
	models_rna[[paste0("rf_rna_",name2)]] <- (train_rf(features=set2_rna ,seeds=seed,nrepeats=3,CV=10,trees=1000,CPU=5))[[1]];
	models_rna[[paste0("xgbt_rna_",name1)]] <-(train_xgbt(features= set1_rna ,seeds=seed,nrepeats=3,CV=10,CPU=5))[[1]];
	models_rna[[paste0("xgbt_rna_",name2)]] <- (train_xgbt(features=set2_rna ,seeds=seed,nrepeats=3,CV=10,CPU=5))[[1]];


#SVM seeds... doparralel and make cluster done
	models_rna[[paste0("radial_rna_",name1)]] <- train(form = Class ~ ., data = set1_rna, metric = "ROC", method = "svmRadial",
			trControl = control,
			tuneGrid = expand.grid(
				sigma = 2^c(-25,-20, -15 ),
				C = 2^c(0, 0.5, 0.75, 1, 1.5)),
			prox = T, allowParallel = TRUE, preProcess = c("center",  "scale"));

	models_rna[[paste0("radial_rna_",name2)]] <- train(form = Class ~ ., data = set2_rna, metric = "ROC", method = "svmRadial",
			trControl = control,
			tuneGrid = expand.grid(
				sigma = 2^c(-25,-20, -15 ),
				C = 2^c(0, 0.5, 0.75, 1, 1.5)),
			prox = T, allowParallel = TRUE, preProcess = c("center",  "scale"));

	models_rna[[paste0("polynomial_rna_",name1)]]  <- train(form = Class ~ ., data = set1_rna, metric = "ROC", method = "svmPoly",
			trControl = control, tuneGrid = expand.grid(degree = (2:4),
				C = 2^c(0, 1,  2), scale = c(0.1,0.2)),
			prox = T, allowParallel = TRUE, preProcess = c("center",      "scale"));
	models_rna[[paste0("polynomial_rna_",name2)]]  <- train(form = Class ~ ., data = set2_rna, metric = "ROC", method = "svmPoly",
			trControl = control, tuneGrid = expand.grid(degree = (2:4),
				C = 2^c(0, 1,  2), scale = c(0.1, 0.2)),
			prox = T, allowParallel = TRUE, preProcess = c("center",      "scale"));

	stopCluster(cl);

	RNA_AUCs<-list();
	VS_models_plot(set1_model= c(models_rna[grep ("12_Dmel" ,names(models_rna))] ,
				models_norna[grep ("12_Dmel" ,	names(models_norna))])  ,
			set2_model=  c(models_rna[grep ("12_Trib" ,	names(models_rna))] ,
				models_norna[grep ("12_Trib" ,	names(models_norna))])  ,
			set2_data= tribolium_features_nzv  ,
			set1_data= drosophila_features_nzv  ,
			set1.dataname= "Dmel" ,
			set2.dataname="Trib" ,
			set1.modelname="Dmel_12" ,
			set2.modelname="Trib_12" ,
			file_prefix="RNA_12_" ,
			test_vs_ZR=T )
		RNA_AUCs <- VS_models_plot(set1_model= c(models_rna[grep ("12_Dmel" ,	names(models_rna))] ,
					models_norna[grep ("12_Dmel" ,	names(models_norna))])  ,
				set2_model=  c(models_rna[grep ("12_Trib" ,	names(models_rna))] ,
					models_norna[grep ("12_Trib" ,	names(models_norna))])  ,
				set2_data= tribolium_features_nzv  ,
				set1_data= drosophila_features_nzv  ,
				set1.dataname= "Dmel" ,
				set2.dataname="Trib" ,
				set1.modelname="Dmel_12" ,
				set2.modelname="Trib_12" ,
				file_prefix="RNA_12_" ,
				test_vs_ZR=F );
	AUCS <- cbind(RNA_AUCs[[1]],do.call(rbind, strsplit(RNA_AUCs[[1]]$model_names, '_') ));
	names(AUCS)<-c("ROC" ,"PRC" , "model_names" ,"testSet" , "model", "Algorith","ZR" ,"Algorithm", "Feature", "Importance","Trained in");
	AUCS$Feature<-gsub("noextrinsic","intrinsic",gsub("rna","extrinsic",AUCS$Feature));

	VS_models_plot(set1_model= c(models_rna[grep ("12_Dmel" , names(models_rna))] ,
				models_norna[grep ("12_Dmel" ,	names(models_norna))])  ,
			set2_model=  c(models_rna[grep ("12_Trib" , names(models_rna))] ,
				models_norna[grep ("12_Trib" ,	names(models_norna))])  ,
			set2_data= noh_trib_features  ,
			set1_data= noh_dmel_features  ,
			set1.dataname= "noh Dmel" ,
			set2.dataname="noh Trib" ,
			set1.modelname="Dmel_12" ,
			set2.modelname="Trib_12" ,
			file_prefix="RNA_12_" ,
			test_vs_ZR=T );
	RNA_AUC <- VS_models_plot(set1_model= c(models_rna[grep ("12_Dmel" , names(models_rna))] ,
				models_norna[grep ("12_Dmel" , 	names(models_norna))])  ,
			set2_model=  c(models_rna[grep ("12_Trib" , names(models_rna))] ,
				models_norna[grep ("12_Trib" ,	names(models_norna))])  ,
			set2_data= noh_trib_features  ,
			set1_data= noh_dmel_features  ,
			set1.dataname= "noh Dmel" ,
			set2.dataname="noh Trib" ,
			set1.modelname="Dmel_12" ,
			set2.modelname="Trib_12" ,
			file_prefix="RNA_12_" ,
			test_vs_ZR=F );
	AUC_noh <- cbind(RNA_AUC[[1]],do.call(rbind, strsplit(RNA_AUC[[1]]$model_names, '_') ));
	names(AUC_noh)<-c("ROC" ,"PRC" , "model_names" , "testSet" , "model", "Algorith","ZR" ,"Algorithm", "Feature", "Importance","Trained in");
	AUC_noh$Feature<-gsub("noextrinsic","intrinsic",gsub("rna","extrinsic",AUC_noh$Feature));

	p1 <-ggplot(AUCS ,aes(x= Algorithm, y=ROC, col=`Trained in`,shape=Feature) ) + geom_point(size=4, position=position_dodge(width=0.5)) +
		ylab(label="AUC-ROC") +
		ggtitle("Algorithms AUC-ROC: Complete ") +
		scale_y_continuous(limits=c(0.5,0.8)) +
		scale_colour_grey() +  theme_classic();
	p2 <-ggplot(AUCS ,aes(x= Algorithm, y=PRC, col=`Trained in`,shape=Feature) ) + geom_point(size=4 , position=position_dodge(width=0.5)) +
		ggtitle("Algorithms AUC-PR: Complete ") +
		geom_hline ( aes( yintercept =ZR, col=`Trained in`) , linetype=2 ) +
		ylab(label="AUC-PRC ") + scale_colour_grey() +  theme_classic();
	legend <- cowplot::get_legend(
			p1 + guides(color = guide_legend(nrow = 1)) +
			theme(legend.position = "bottom") );
	prow<-cowplot::plot_grid(p1+ theme(legend.position="none"), 
			p2+ theme(legend.position="none"), 
			labels = c('A', 'B'), 
			label_size = 12, 
			align = 'vh',  hjust = -1 , nrow=1);
	pab<-cowplot::plot_grid(prow, legend, rel_heights = c(1, .1)  , ncol =1 );

	p1 <-ggplot(AUC_noh ,aes(x= Algorithm, y=ROC, col=`Trained in`,shape=Feature) ) + geom_point(size=4 , position=position_dodge(width=0.5)) +
		ggtitle("Algorithms AUC-ROC: NOH") +
		ylab(label="AUC-ROC") +
		scale_y_continuous(limits=c(0.5,0.8)) +
		scale_colour_grey() +  theme_classic() ;
	p2 <-ggplot(AUC_noh ,aes(x= Algorithm, y=PRC, col=`Trained in`,shape=Feature) ) + geom_point(size=4, position=position_dodge(width=0.5) ) +
		ggtitle("Algorithms AUC-PR: NOH") +
		ylab(label="AUC-PRC") + 
		geom_hline ( aes( yintercept=ZR, col=`Trained in`) , linetype=2  ) +
		scale_colour_grey() +  theme_classic();
	prow<-cowplot::plot_grid(p1+ theme(legend.position="none"),
			p2+ theme(legend.position="none"),
			labels = c('C', 'D'),
			label_size = 12,
			align = 'vh',  hjust = -1 , nrow=1);
	pcd<-cowplot::plot_grid(prow, rel_heights = c(1, .1)  , ncol =1 );
	svg("AUCS_ROC_PRC_RNA.svg",height=8,width=8);
	cowplot::plot_grid(pab,pcd,nrow=2);
	dev.off();
}



experiment5_extrinsic_add<-function(CPU=5){
  library(caret)
  THREADS=CPU;
  IMPORTANCE=12;
    imp<-c(importance_united$importance$Overall>IMPORTANCE,TRUE);
  name1<-paste0(IMPORTANCE,"_Dmel");
  name2<-paste0(IMPORTANCE,"_Trib");
  set1_rna<-drosophila_features_nzv[,c(imp,rep(T,13))];
  set2_rna<-tribolium_features_nzv[,c(imp,rep(T,13))];
  experimento_5_models<-list();
  experimento_5_models$dmel<-list();
  experimento_5_models$trib<-list();
  
  
  noh_dmel_ids <-read.csv(system.file("extdata", "noh_dmel_ids.csv", package="GeneEssentiality"),row.names = 1)
  noh_trib_ids <-read.csv(system.file("extdata", "noh_trib_ids.csv", package="GeneEssentiality"),row.names = 1)
  noh_trib_features <- tribolium_features_nzv[  ( rownames(tribolium_features_nzv ) %in% noh_trib_ids[,1] ) , ]
  noh_dmel_features <- drosophila_features_nzv[ ( rownames(drosophila_features_nzv) %in% noh_dmel_ids[,1] ) , ]
  
  
	mtries<-round(sqrt(length(set1_rna)))
		mtries<-c( mtries, mtries*2)
		control <- trainControl(method = "repeatedcv", number = 10, repeats = 3,
				classProbs = TRUE, summaryFunction = twoClassSummary,
				savePredictions = "final", seeds = seed, preProcOptions = NULL, allowParallel=F);
		rf_tests<-list();
## training this way generates much more compact models when saving (24MB to 2,7MB), for some reason using doparallel reduce compression power
		features <- set1_rna[,1:632]
		experimento_5_models$dmel$dmel_adult<-train(Class ~ .,data=features, metric="ROC",method="ranger",
				tuneGrid=expand.grid(.mtry=mtries , .splitrule="gini",.min.node.size=1),
				trControl=control,num.trees= 1000,num.threads=THREADS,importance = 'impurity')
		features <- set2_rna[,1:632]
		experimento_5_models$trib$trib_adult<-train(Class ~ .,data=features, metric="ROC",method="ranger",
				tuneGrid=expand.grid(.mtry=mtries , .splitrule="gini",.min.node.size=1),
				trControl=control,num.trees= 1000,num.threads=THREADS,importance = 'impurity')

#only larva
		features <- set1_rna[,c(1:631,633)]
		experimento_5_models$dmel$dmel_larva <-train(Class ~ .,data=features, metric="ROC",method="ranger",
				tuneGrid=expand.grid(.mtry=mtries , .splitrule="gini",.min.node.size=1),
				trControl=control,num.trees= 1000,num.threads=THREADS,importance = 'impurity')
		features <- set2_rna[,c(1:631,633)]
		experimento_5_models$trib$trib_larva <-train(Class ~ .,data=features, metric="ROC",method="ranger",
				tuneGrid=expand.grid(.mtry=mtries , .splitrule="gini",.min.node.size=1),
				trControl=control,num.trees= 1000,num.threads=THREADS,importance = 'impurity')
#only subcell
		features <- set1_rna[,c(1:631)]
		experimento_5_models$dmel$dmel_subcel<-train(Class ~ .,data=features, metric="ROC",method="ranger",
				tuneGrid=expand.grid(.mtry=mtries , .splitrule="gini",.min.node.size=1),
				trControl=control,num.trees= 1000,num.threads=THREADS,importance = 'impurity')
		features <- set2_rna[,c(1:631)]
		experimento_5_models$trib$trib_subcel <-train(Class ~ .,data=features, metric="ROC",method="ranger",
				tuneGrid=expand.grid(.mtry=mtries , .splitrule="gini",.min.node.size=1),
				trControl=control,num.trees= 1000,num.threads=THREADS,importance = 'impurity')
#both
		features <- set1_rna
		experimento_5_models$dmel$dmel_2rnaseq<-train(Class ~ .,data=features, metric="ROC",method="ranger",
				tuneGrid=expand.grid(.mtry=mtries , .splitrule="gini",.min.node.size=1),
				trControl=control,num.trees= 1000,num.threads=THREADS,importance = 'impurity')
		features <- set2_rna
		experimento_5_models$trib$trib_2rnaseq <-train(Class ~ .,data=features, metric="ROC",method="ranger",
				tuneGrid=expand.grid(.mtry=mtries , .splitrule="gini",.min.node.size=1),
				trControl=control,num.trees= 1000,num.threads=THREADS,importance = 'impurity')
		
		features <-drosophila_features_nzv[,c(imp,rep(F,13))];
		experimento_5_models$dmel$intrinsic_Dmel <-train(Class ~ .,data=features, metric="ROC",method="ranger",
		                                                   tuneGrid=expand.grid(.mtry=mtries , .splitrule="gini",.min.node.size=1),
		                                                   trControl=control,num.trees= 1000,num.threads=THREADS,importance = 'impurity')
		features <-tribolium_features_nzv[,c(imp,rep(F,13))];
		experimento_5_models$trib$intrinsic_Trib <-train(Class ~ .,data=features, metric="ROC",method="ranger",
		                                                   tuneGrid=expand.grid(.mtry=mtries , .splitrule="gini",.min.node.size=1),
		                                                   trControl=control,num.trees= 1000,num.threads=THREADS,importance = 'impurity')

		experimento_5_models$dmel<- c(models_norna[grep ("rf_norna_12_Dmel",names(models_norna) )] , models_rna[grep ("rf_rna_12_Dmel",names(models_rna) )] , rf_tests[grep ("dmel",names(rf_tests))])
		experimento_5_models$trib<- c(models_norna[grep ("rf_norna_12_Trib",names(models_norna) )] , models_rna[grep ("rf_rna_12_Trib",names(models_rna) )] , rf_tests[grep ("trib",names(rf_tests))])
		
		
		experimento5 <- VS_models_plot(set1_model= experimento_5_models$dmel , set2_model= experimento_5_models$trib ,
		                               set2_data= tribolium_features_nzv ,set1_data= drosophila_features_nzv ,
		                               set1.dataname= "Dmel", set2.dataname="Trib", set1.modelname="Dmel", set2.modelname="Trib", file_prefix="experimento5_",test_vs_ZR=F )
		experimento5[[1]]$Features<-c( "I + S + Adult", "I + S + Larva","I + Subcellular (S)","I + S + Adult+Larva","Intrinsic (I)")
		
		experimento5_noh <- VS_models_plot(set1_model= experimento_5_models$dmel , set2_model= experimento_5_models$trib , set2_data= noh_trib_features ,set1_data= noh_dmel_features , set1.dataname= "NOH Dmel", set2.dataname="NOH Trib", set1.modelname="Dmel", set2.modelname="Trib", file_prefix="experimento5_noh_",test_vs_ZR=F )
		experimento5_noh[[1]]$Features<-c( "I + S + Adult", "I + S + Larva","I + Subcellular (S)","I + S + Adult+Larva","Intrinsic (I)")
		
		
		p1 <-ggplot(experimento5[[1]] ,aes(x= Features , y=ROC, col= `Trained in`  )) + geom_point(size=4, position=position_dodge(width=0.5)) +
		  ylab(label="AUC-ROC") +
		  ggtitle("ROC- Extrinsic addition: Complete") +
		  scale_y_continuous(limits=c(0.5,0.8)) +
		  scale_x_discrete(limits= c("Intrinsic (I)", "I + Subcellular (S)", "I + S + Larva", "I + S + Adult", "I + S + Adult+Larva")) +
		  scale_colour_grey() +  theme_classic() +
		  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1))
		p2 <-ggplot(experimento5[[1]] ,aes(x= Features , y=PRC, col= `Trained in`  )) + geom_point(size=4, position=position_dodge(width=0.5)) +
		  ggtitle("PR- Extrinsic addition: Complete") +
		  geom_hline ( aes( yintercept =PR_ZR, col=`Trained in`) , linetype=2 ) +
		  scale_x_discrete(limits= c("Intrinsic (I)", "I + Subcellular (S)", "I + S + Larva", "I + S + Adult", "I + S + Adult+Larva")) +
		  ylab(label="AUC-PRC ") + scale_colour_grey() +  theme_classic() +
		  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1))
		legend <- cowplot::get_legend(
		  p1 + guides(color = guide_legend(nrow = 1)) +
		    theme(legend.position = "bottom")
		)
		prow<-cowplot::plot_grid(p1+ theme(legend.position="none"),
		                         p2+ theme(legend.position="none"),
		                         labels = c('A', 'B'),
		                         label_size = 12,
		                         align = 'vh',  hjust = -1 , nrow=1 )
		p1 <-ggplot(experimento5_noh[[1]] ,aes(x= Features , y=ROC, col= `Trained in`  )) + geom_point(size=4, position=position_dodge(width=0.5)) +
		  ylab(label="AUC-ROC") +
		  ggtitle("ROC- Extrinsic addition: NOH") +
		  scale_y_continuous(limits=c(0.5,0.8)) +
		  scale_x_discrete(limits= c("Intrinsic (I)", "I + Subcellular (S)", "I + S + Larva", "I + S + Adult", "I + S + Adult+Larva") )+
		  scale_colour_grey() +  theme_classic() +
		  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1))
		p2 <-ggplot(experimento5_noh[[1]] ,aes(x= Features , y=PRC, col= `Trained in`  )) + geom_point(size=4, position=position_dodge(width=0.5)) +
		  ggtitle("PR- Extrinsic addition: NOH") +
		  geom_hline ( aes( yintercept =PR_ZR, col=`Trained in`) , linetype=2 ) +
		  scale_x_discrete(limits= c("Intrinsic (I)", "I + Subcellular (S)", "I + S + Larva", "I + S + Adult", "I + S + Adult+Larva") )+
		  ylab(label="AUC-PRC ") + scale_colour_grey() +  theme_classic() +
		  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1))
		#legend <- cowplot::get_legend(
		  #p1 + guides(color = guide_legend(nrow = 1)) +
		    #theme(legend.position = "bottom")
		#)
		prow2<-cowplot::plot_grid(p1+ theme(legend.position="none"),
		                         p2+ theme(legend.position="none"),
		                         labels = c('C', 'D'),
		                         label_size = 12,
		                         align = 'vh',  hjust = -1 , nrow=1)
		
		svg("Comparacao_rnaseq.svg",width=8,height = 9)
		cowplot::plot_grid(prow,legend,prow2, rel_heights = c(1, .1, 1)  , ncol =1 )
		
		dev.off()
		
}

