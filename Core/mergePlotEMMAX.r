#!/bin/R

#################################################################
# mergePlotEMMAX.r - Matt Ravenhall								#
# - Convert EMMAX output to Manhattan and QQ plots 				#
# - Includes the addition of appropriate info, such as lambda 	#
# - File prefix should be based at the command-line				#
# - This script is a component of EMMAX_Pipeline.sh				#
#################################################################

args <- commandArgs(trailingOnly=TRUE)			# args[1] = outprefix

library(qqman)

filePrefix <- args[1]
#filePrefix <- "FINAL_Cand_Merge_preImpute_wMasked"
indxFile <- paste(filePrefix,'_MergeIndex.tmp',sep='')

print('Parsing P value files')
addFile <- read.table(paste(filePrefix,"_Add.ps",sep=''),header=F)
hetFile <- read.table(paste(filePrefix,"_Het.ps",sep=''),header=F)
domFile <- read.table(paste(filePrefix,"_Dom.ps",sep=''),header=F)
recFile <- read.table(paste(filePrefix,"_Rec.ps",sep=''),header=F)

names(addFile) <- c('SNP', 'addB', 'ADD')
names(hetFile) <- c('SNP', 'hetB', 'HET')
names(domFile) <- c('SNP', 'domB', 'DOM')
names(recFile) <- c('SNP', 'recB', 'REC')

numModels <- 4 # Additive, Heterozygous, Recessive & Dominant

# Get Index.tmp
print('Parsing SNP Index')
indx <- read.table(indxFile,header=F)
names(indx) <- c('CHR', 'SNP', 'BP')

# Merge
print('Merging all models')
allP.AH <- merge(addFile, hetFile)
allP.DR <- merge(domFile, recFile)
allP <- merge(allP.AH, allP.DR)

onlyP <- data.frame(allP$SNP, allP$ADD, allP$HET, allP$DOM, allP$REC)
names(onlyP) <- c('SNP', 'ADD', 'HET', 'DOM', 'REC')
onlyP$P <- do.call(pmin,onlyP[,2:(numModels+1)]) 					# Pmin
onlyP$MODEL <- names(onlyP)[apply(onlyP[,2:(numModels+1)], 1, which.min)+1]		# Model of Pmin

toPlot <- merge(indx, onlyP)
### deal with missing values here?
# toPlot names should be: SNP, CHR, BP, ADD, HET, DOM, REC, P, MODEL
write.table(toPlot, file=paste(filePrefix,'.ps',sep=''), quote=F, row.names=F)

## Manhattan plots

# Chromosomes
for (i in unique(toPlot$CHR)) {
	print(paste('Plotting chr',i))
	png(paste(filePrefix,'_chr',i,'_Pmin_EMMAX_MPlot.png',sep=''))
	manhattan(toPlot[toPlot$CHR == i,])
	invisible(dev.off())
}

# Additive
print('Plotting Additive Ps')
png(paste(filePrefix,'_Add_EMMAX_MPlot.png',sep=''))
manhattan(toPlot, p='ADD')
invisible(dev.off())

# Heterozygous
print('Plotting Heterozygous Ps')
png(paste(filePrefix,'_Het_EMMAX_MPlot.png',sep=''))
manhattan(toPlot, p='HET')
invisible(dev.off())

# Dominant
print('Plotting Dominant Ps')
png(paste(filePrefix,'_Dom_EMMAX_MPlot.png',sep=''))
manhattan(toPlot, p='DOM')
invisible(dev.off())

# Recessive
print('Plotting Recessive Ps')
png(paste(filePrefix,'_Rec_EMMAX_MPlot.png',sep=''))
manhattan(toPlot, p='REC')
invisible(dev.off())

# Pmin
print('Plotting minimum Ps')
png(paste(filePrefix,'_Pmin_EMMAX_MPlot.png',sep=''))
manhattan(toPlot)
invisible(dev.off())

#calculate lambda
print('Processing Genomic Inflation Factors')
        
# QQ plots and Genomic Inflation Factors

# Additive
chisq <- qchisq(1-toPlot$ADD,1)
lambda <- median(chisq)/qchisq(0.5,1)

print(paste('Additive GIF: ',lambda,sep=''))

png(paste(filePrefix,'_Add_EMMAX_QQplot.png',sep=''))
qq(toPlot$ADD,main=paste('lambda: ',lambda,sep=''))
invisible(dev.off())

# Heterozygous
chisq <- qchisq(1-toPlot$HET,1)
lambda <- median(chisq)/qchisq(0.5,1)

print(paste('Heterozygous GIF: ',lambda,sep=''))

png(paste(filePrefix,'_Het_EMMAX_QQplot.png',sep=''))
qq(toPlot$HET,main=paste('lambda: ',lambda,sep=''))
invisible(dev.off())

# Dominant
chisq <- qchisq(1-toPlot$DOM,1)
lambda <- median(chisq)/qchisq(0.5,1)

print(paste('Dominant GIF: ',lambda,sep=''))

png(paste(filePrefix,'_Dom_EMMAX_QQplot.png',sep=''))
qq(toPlot$DOM,main=paste('lambda: ',lambda,sep=''))
invisible(dev.off())

# Recessive
chisq <- qchisq(1-toPlot$REC,1)
lambda <- median(chisq)/qchisq(0.5,1)

print(paste('Recessive GIF: ',lambda,sep=''))

png(paste(filePrefix,'_Rec_EMMAX_QQplot.png',sep=''))
qq(toPlot$REC,main=paste('lambda: ',lambda,sep=''))
invisible(dev.off())
