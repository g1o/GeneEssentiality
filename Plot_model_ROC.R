#!/bin/env R
#Plots the ROC from the model cross-validation 
library(data.table)
library(seqinr)
library(parallel)
library(doMC)
library(pROC)
library(caret)
library(protr)
library(ggplot2)
library(plotROC)
anot<-paste0('\n',paste0(label[1,1:2],collapse=" AUC="),'\n',paste0(label[2,1:2],collapse=" AUC="),'\n',paste0(label[3,1:2],collapse=" AUC="),'\n',paste0(label[4,1:2],collapse=" AUC="),'\n')

g<-ggplot(modellist[[1]]$pred,aes(m = NE, d = obs, color=factor(mtry))) + geom_roc(n.cut=0) + style_roc(theme=theme_classic,ylab="Sensibilidade (Taxa de verdadeiros positivos)",xlab="Taxa de falsos positivos (Especificidade - 1)") + annotate(geom="text",x = 0.7, y = 0.2, label=c(anot))+ geom_abline(intercept=0, slope=1, colour="orange")  ;g;dev.off()

g; dev.off()
