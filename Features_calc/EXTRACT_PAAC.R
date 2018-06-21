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
