#protein features not working well with another genetic code
#just tested with protein coding regions and the standard genetic code (ncbi 1).
Calc_feats<-function(seq){ #calculates the features of one sequence
#                        fs.time <- Sys.time() #timing code
                        frame0<-translate(seq,frame=0)
                        frame1<-translate(seq,frame=1)
                        frame2<-translate(seq,frame=2)

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

                        fs.time <- Sys.time() #timing code
#Counting words frequencies
                        f1<-count(seq,1,freq=T)
                        f2<-count(seq,2,freq=T)
                        f3<-count(seq,3,freq=T)

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


        #Peptide properties
        PEP<-AAstat(longest_orf,plot=F)
        pepFeatures<-c(as.numeric(PEP[[2]]),as.numeric(PEP[3]))
        names(pepFeatures)<-c(names(PEP[[2]]),names(PEP[3]))
        

	#peptide length
        aalength<-length(longest_orf)
        names(aalength)<-c("aalength")

        #Counting aminoacid word frequencies
        f1<-count(longest_orf,1,freq=T,alphabet=s2c("ACDEFGHIKLMNPQRSTVWY"))
        f2<-count(longest_orf,2,freq=T,alphabet=s2c("ACDEFGHIKLMNPQRSTVWY")) #Do not count the STOP "*"
        f3<-count(longest_orf,3,freq=T,alphabet=s2c("ACDEFGHIKLMNPQRSTVWY")) #Do not count the STOP "*"

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


        #Pseudo aminoacid
	LAMBDA=22;
	OMEGA=0.05;
        PseAA<-EXTRACT_PAAC(longest_orf,lambda=LAMBDA,w=OMEGA) #Chou's pseudoAA #22 min for AB


#       Features<-rbind(SomaMI,MI,SomaCMI,CMI,H2,H3,pH2,pepFeatures,protMI,Soma_protMI,aalength,PseAA)
#       Features<-data.frame(t(Features))
	Features<-c(SomaMI,MI,SomaCMI,CMI,H2,H3,pH2,pH3,pepFeatures,aalength,ProtMI,SomaProtMI,PseAA)
                        
#	fe.time <- Sys.time() #timing code
#        time.taken <- fe.time - fs.time
#        print(c("Time for features calc=>",time.taken))
        
	return(Features)

}


#EXTRACT_PAAC is an adaptation from the function inside protR package in order to accept the SeqFastaAA class from seqinR package
EXTRACT_PAAC<- function (x, props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"), 
    lambda = 30, w = 0.05, customprops = NULL)
{
    if (protcheck(x) == FALSE)
        stop("x has unrecognized amino acid type")
    if (length(x) <= lambda)  #edited to use with SeqFastaAA class from seqinR
        stop("Length of the protein sequence must be greater than \"lambda\"")
    AAidx = read.csv(system.file("sysdata/AAidx.csv", package = "protr"),
        header = TRUE)
    tmp = data.frame(AccNo = c("Hydrophobicity", "Hydrophilicity",
        "SideChainMass"), A = c(0.62, -0.5, 15), R = c(-2.53,
        3, 101), N = c(-0.78, 0.2, 58), D = c(-0.9, 3, 59), C = c(0.29,
        -1, 47), E = c(-0.74, 3, 73), Q = c(-0.85, 0.2, 72),
        G = c(0.48, 0, 1), H = c(-0.4, -0.5, 82), I = c(1.38,
            -1.8, 57), L = c(1.06, -1.8, 57), K = c(-1.5, 3,
            73), M = c(0.64, -1.3, 75), F = c(1.19, -2.5, 91),
        P = c(0.12, 0, 42), S = c(-0.18, 0.3, 31), T = c(-0.05,
            -0.4, 45), W = c(0.81, -3.4, 130), Y = c(0.26, -2.3,
            107), V = c(1.08, -1.5, 43))
    AAidx = rbind(AAidx, tmp)
    if (!is.null(customprops))
        AAidx = rbind(AAidx, customprops)
    aaidx = AAidx[, -1]
    row.names(aaidx) = AAidx[, 1]
    n = length(props)
    H0 = as.matrix(aaidx[props, ])
    H = matrix(ncol = 20, nrow = n)
    for (i in 1:n) H[i, ] = (H0[i, ] - mean(H0[i, ]))/(sqrt(sum((H0[i,
        ] - mean(H0[i, ]))^2)/20))
    AADict = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
    dimnames(H) = list(props, AADict)
    Theta = vector("list", lambda)
    xSplitted = x #edited to use with SeqFastaAA class from seqinR
    N = length(xSplitted)
    for (i in 1:lambda) {
        for (j in 1:(N - i)) {
            Theta[[i]][j] = mean((H[, xSplitted[j]] - H[, xSplitted[j +
                i]])^2)
        }
    }
    theta = sapply(Theta, mean)
    fc = summary(factor(xSplitted, levels = AADict), maxsum = 21)
    Xc1 = fc/(1 + (w * sum(theta)))
    names(Xc1) = paste("Xc1.", names(Xc1), sep = "")
    Xc2 = (w * theta)/(1 + (w * sum(theta)))
    names(Xc2) = paste("Xc2.lambda.", 1:lambda, sep = "")
    Xc = c(Xc1, Xc2)
    Xc
}


