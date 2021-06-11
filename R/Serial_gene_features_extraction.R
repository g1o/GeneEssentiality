#' Calculate features of a sequence
#'
#' Extract the features from a single nucleotide sequence. Also translate in the 3+ frames seaching for the longest ORF, using it to extract the features for the amino acid sequence. 
#' 
#' @param seq nucleotide sequence, works using a SeqFastadna object from the seqinr package. 
#' @param PFAM_PATH Path to the Pfam-A database
#' @param LAMBDA Pseudo AA composition parameter
#' @param OMEGA Pseudo AA composition parameter
#' 
#protein features not working well with another genetic code
#just tested in protein coding regions and the standard genetic code (ncbi 1).
Calc_feats<-function(seq,aa="",PFAM_PATH="",LAMBDA=30,OMEGA=0.05,nuc_only=F){ #calculates the features of one sequence
#                        fs.time <- Sys.time() #timing code
	if(rDNAse::dnacheck(paste(seq,collapse=''))){return();}
	if(length(seq)<=LAMBDA*3){
	#	stop("AA length too short: length < LAMBDA"); 	#Can't calculate pseudo aa with it... 
		return();
	}
if(nuc_only==F){
	if(aa==""){
		frame0<-seqinr::translate(seq,frame=0,ambiguous=T)
			frame1<-seqinr::translate(seq,frame=1,ambiguous=T)
			frame2<-seqinr::translate(seq,frame=2,ambiguous=T)

			longest_orf<-{
				start0<-which( frame0 %in% "M")    #vetor posicao de codon de metionina (possivel inicio)
					stop0<-which( frame0 %in% "*")     #vetor posicao de codon de parada
					max0<-as.numeric(sapply(seq_along(start0),function (x){ stop0[stop0>start0[x]][1]-start0[x]})) #vetor de comprimento de orfs (Pfim - Pini)
					start1<-which( frame1 %in% "M")
					stop1<-which( frame1 %in% "*")
					max1<-as.numeric(sapply(seq_along(start1),function (x){ stop1[stop1>start1[x]][1]-start1[x]}))
					start2<-which( frame2 %in% "M")
					stop2<-which( frame2 %in% "*")
					max2<-as.numeric(sapply(seq_along(start2),function (x){ stop2[stop2>start2[x]][1]-start2[x]}))
					maxi<-max(max0,max1,max2,na.rm=T)
					if( max(c(max0,0),na.rm=T) == maxi){
						maxi<-which.max(max0)
							ret<-frame0[start0[maxi]:stop0[stop0>start0[maxi]][1]]
							ret[-max(NROW(ret))] # REMOVE * (STOP CODON) and return
					}else if (max(c(max1,0),na.rm=T) == maxi){
						maxi<-which.max(max1)
							ret<-frame1[start1[maxi]:stop1[stop1>start1[maxi]][1]]
							ret[-max(NROW(ret))]
					}else if (max(c(max2,0),na.rm=T) == maxi){
						maxi<-which.max(max2)
							ret<-frame2[start2[maxi]:stop2[stop2>start2[maxi]][1]]
							ret[-max(NROW(ret))]
					}
			}
	}else {
		longest_orf<-aa ;
		rm(aa)
	}
	#peptide length
        aalength<-length(longest_orf)
        names(aalength)<-c("aalength")
	if(aalength<=LAMBDA){
	#Can't calculate pseudo aa with it... 
		return();
	}
	#Resolve problems due to a ambiguos base
	Xaa<-seqinr::count(longest_orf,1,freq=F,alphabet=seqinr::s2c("X")) # Check number of surviving ambiguous bases 
	if(Xaa>=1){ #if there are ambiguos bases 
		return();
	}
	rm(Xaa);
}
#                        fs.time <- Sys.time() #timing code
#Counting words frequencies
                        f1<-seqinr::count(seq,1,freq=T)
                        f2<-seqinr::count(seq,2,freq=T)
                        f3<-seqinr::count(seq,3,freq=T)
	
#Gibbs Entropy and Entalpy from varGibbs
                        tmp_file<-Sys.getpid();
                        sequence<-toupper(paste(seq,collapse=''))
		
                        varGibbs_Model_PATH<-c("/mnt/DATABASES/bin/VarGibbs-2.2/data/AOP-CMB.par") #change to paramter
                        VarGibbsCMD<-sprintf("vargibbs -seq=%s -ct=1 -o='%s' -calc=prediction -par=%s >/dev/null;
                                        awk -F' ' -v OFS=';' '{print $10,$13}'  %s.dat | head -n2 ",
                                        sequence,tmp_file,varGibbs_Model_PATH,tmp_file) #build command
                        VarGibbs<-system(VarGibbsCMD,intern=T)
                        Gibbs<-as.numeric((strsplit(VarGibbs,';'))[[2]])
                        names(Gibbs)<-strsplit(VarGibbs,';')[[1]]
                        system(paste0('rm ',tmp_file,'*')) # remove tmp files
				
	
#Shannon Entropy
                        H2<-sum(-f2*log2(f2),na.rm=T) #shannon entropy for word size 2
                        H3<-sum(-f3*log2(f3),na.rm=T) #shannon entropy for word size 3
			names(H2)<-c("H2")
			names(H3)<-c("H3")

#Nucleotide Mutual Information
MI<-sapply (names(f2) , function(dinucleotide) f2[dinucleotide]*log2(f2[dinucleotide]/(f1[substring(dinucleotide,1,1)]*f1[substring(dinucleotide,2,2)])) ) 
 
                        names(MI)<-names(f2)
                        SomaMI<-sum(MI,na.rm=T)
			names(SomaMI)<-c("SomaMI")

#Nucleotide Conditional Mutual Information
			CMI<- sapply (names(f3) , function(trinucl)  f3[trinucl]*log2(f1[substring(trinucl,2,2)]*f3[trinucl]/(f2[substring(trinucl,1,2)]*f2[ substring(trinucl,2,3)])))
                        names(CMI)<-names(f3)
                        SomaCMI<-sum(CMI,na.rm=T)
			names(SomaCMI)<-c("SomaCMI")
## rDNAse features
	DACC<-	rDNAse::extrDACC(sequence,allprop=T,nlag=2)
	TACC<-	rDNAse::extrTACC(sequence,allprop=T,nlag=LAMBDA*3)
	PseDNC<-rDNAse::extrPseDNC(sequence,lambda=3) # 
	rm(sequence)
if(nuc_only==F){
 #Peptide properties
        PEP<-seqinr::AAstat(longest_orf,plot=F)
        pepFeatures<-c(as.numeric(PEP[[2]]),as.numeric(PEP[3]))
        names(pepFeatures)<-c(names(PEP[[2]]),names(PEP[3]))
        
        #Counting aminoacid word frequencies
        f1<-seqinr::count(longest_orf,1,freq=T,alphabet=seqinr::s2c("ACDEFGHIKLMNPQRSTVWY"))
        f2<-seqinr::count(longest_orf,2,freq=T,alphabet=seqinr::s2c("ACDEFGHIKLMNPQRSTVWY")) #Do not count the STOP "*"
        f3<-seqinr::count(longest_orf,3,freq=T,alphabet=seqinr::s2c("ACDEFGHIKLMNPQRSTVWY")) #Do not count the STOP "*"

        #Shannon Entropy for aminoacids
	pH2<-sum(-f2*log2(f2),na.rm=T)
        names(pH2)<-c("pH2")

                        pH3<-sum(-f3*log2(f3),na.rm=T) #shannon entropy for amino acid word size 3
			names(pH3)<-c("pH3")

	#Protein Mutual Information
	ProtMI<-sapply (names(f2) , function(dinucleotide) f2[dinucleotide]*log2(f2[dinucleotide]/(f1[substring(dinucleotide,1,1)]*f1[substring(dinucleotide,2,2)]) ) )
                        names(ProtMI)<-names(f2)
                        SomaProtMI<-sum(ProtMI,na.rm=T)
                        names(SomaProtMI)<-c("SomaProtMI")

	#Protein Conditional Mutual Information (close to 8000 features)
                        pCMI<- sapply (names(f3) , function(trinucl)  f3[trinucl]*log2(f1[substring(trinucl,2,2)]*f3[trinucl]/(f2[substring(trinucl,1,2)]*f2[ substring(trinucl,2,3)])))
                        names(pCMI)<-names(f3)
                        pSomaCMI<-sum(pCMI,na.rm=T)
			names(pSomaCMI)<-c("pSomaCMI")
	rm(f1,f2,f3);
##protR
        PseAA  <- protr::extractPAAC(paste(longest_orf,collapse=''),lambda=LAMBDA,w=OMEGA) # pseudoAA
	ApseAA <- protr::extractAPAAC(paste(longest_orf,collapse=''),lambda=LAMBDA,w=OMEGA)
	CTriad <- protr::extractCTriad(paste(longest_orf,collapse='')) #Conjoint Triad 
	MoreauBroto<-protr::extractMoreauBroto(paste(longest_orf,collapse=''))
	Moran<-	protr::extractMoran(paste(longest_orf,collapse=''))
	Geary<-	protr::extractGeary(paste(longest_orf,collapse=''))
	CTDC<-	protr::extractCTDC(paste(longest_orf,collapse='')) 
	CTDD<-	protr::extractCTDD(paste(longest_orf,collapse=''))
	CTDT<-	protr::extractCTDT(paste(longest_orf,collapse=''))

	Features<-c(SomaMI,MI,SomaCMI,CMI,H2,H3,Gibbs,DACC,TACC,PseDNC,pH2,pH3,pSomaCMI,pCMI,pepFeatures,aalength,ProtMI,SomaProtMI,PseAA,CTriad,ApseAA,MoreauBroto,Moran,Geary,CTDC,CTDD,CTDT)
	Features[is.na(Features)] <- 0 
	rm(SomaMI,MI,SomaCMI,CMI,H2,H3,Gibbs,DACC,TACC,PseDNC,pH2,pH3,pSomaCMI,pCMI,pepFeatures,aalength,ProtMI,SomaProtMI,PseAA,CTriad,ApseAA,MoreauBroto,Moran,Geary,CTDC,CTDD,CTDT);
}
	if(nuc_only==T){
	        Features<-c(SomaMI,MI,SomaCMI,CMI,H2,H3,Gibbs,DACC,TACC,PseDNC);
        	Features[is.na(Features)] <- 0;
		rm(SomaMI,MI,SomaCMI,CMI,H2,H3,Gibbs,DACC,TACC,PseDNC);
	}
#Pfam--- Uses external HMMSCAN and perl. Needs pfam database located in the working dir     ↓
	if(PFAM_PATH!="" && nuc_only==F){
		hmmCMD<-sprintf("perl -e \"while(<>){print '>1\n';print } \" -| hmmscan --acc --noali %s - |\
				grep -Po '>> PF\\S+'|sed 's/>> //'", PFAM_PATH) #build command
			PFAM<-system(hmmCMD,input=paste(longest_orf,collapse=''),intern=T)  #execute command inputing amino acid sequence from longest_orf
			names(PFAM)<-PFAM;
			PFAMa <- ifelse(grepl(".", PFAM), T, F)
			names(PFAMa)<-names(PFAM) 
# numeric vectors
			Features<-data.frame(t(Features),t(PFAMa)) #join with binary vector in a data.frame as they have different types
			row.names(Features)<-seqinr::getName(seq)                        
	}else{
		Features<-data.frame(t(Features))
			row.names(Features)<-seqinr::getName(seq)
	}
	gc();
	return(Features)
}


