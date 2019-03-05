#!/bin/R

#####################################################################################
# plotSNPIntensity.r - Matt Ravenhall												#
# - Given an array of SNPs as shown below, produce SNP intensity plots for each. 	#
# - e.g. 'Rscript plotSNPIntensity.r CHR1:BP1 CHR2:BP2 CHR3:BP3'	 				#
#####################################################################################

# Note: Currently XYfile is subset for each SNP, this could be improved with pre-loop chromosome subsets

# Optional
xlim.fixed <- 2
ylim.fixed <- 2
sepChar <- ':'
XYfile <- './PLACEHOLDER.XY.FORMAT' # Requires POS (position), CHR (chromosome) and XY (intensity as X,Y) columns
GTfile <- './PLACEHOLDER.GT.FORMAT' # Requires POS (position), CHR (chromosome) and colID columns
CCfile <- './PLACEHOLDER.sample'	# Requires 'scanID' (IDs) and 'caseorcontrol' ('CASE'/'CONTROL') columns

# Grab Position (in BP) of SNP
args <- commandArgs(trailingOnly=TRUE)

# Bring in all SNPs
print('Parsing XY Data')
dat <- read.table(XYfile,header=T)
print('Parsing GT data (for point colours)')
cols <- read.table(GTfile,header=T)

# Tickers
skippedSNPs=0
totalSNPs=0

print('--------------------------------------------------------------------')

# Loop SNPs, currently probably not the most efficient approach but it'll work
for (eachsnp in args) {

	totalSNPs = totalSNPs + 1

	# Skip if separation character is not in a given argument
	if (!grepl(sepChar, eachsnp)) {
		print(paste('Skipping ',eachsnp,' due to missing separation character (',sepChar,').',sep=''))
		skippedSNPs = skippedSNPs + 1
		print('--------------------------------------------------------------------')
		next
	}

	# Use i for BP, j for CHR
	i <- strsplit(eachsnp, ':')[[1]][2]
	j <- strsplit(eachsnp, ':')[[1]][1]

	# Isolate SNP of interest
	print(paste('Isolating XY data for SNP on CHR ',j,' at BP ',i,sep=''))
	snp <- data.frame(t(dat[dat$POS == i & dat$CHR == j,]))
	# Check the SNP, and only one version of that SNP, exists in dataset
	if (dim(snp)[2] != 1) {
		print(paste('Skipping ',eachsnp,' due to missing or duplicate value in XY dataset.',sep=''))
		skippedSNPs = skippedSNPs + 1
		print('--------------------------------------------------------------------')
		next
	}
	print('Removing CHROM and POS data from isolated subset')
	snp <- subset(snp, row.names(snp) != 'CHROM' & row.names(snp) != 'POS')


	# Split XY columns
	names(snp) <- c('XY')
	print('Splitting XY data into appropriate columns')
	snp$X <- lapply(strsplit(as.character(snp$XY), ','), '[', 1)
	snp$Y <- lapply(strsplit(as.character(snp$XY), ','), '[', 2)

	# Add colours
	print(paste('Isolating GT data for SNP on CHR ',j,' at BP ',i,sep=''))
	snpcol <- data.frame(t(cols[cols$POS == i & cols$CHR == j,]))
	print('Removing CHROM and POS data from isolated subset')
	snpcol <- subset(snpcol, row.names(snpcol) != 'CHROM' & row.names(snpcol) != 'POS')

	print('Splitting GT data into appropriate columns')
	names(snpcol) <- c('XY')
	snpcol$X <- lapply(strsplit(as.character(snpcol$XY), '/'), '[', 1)
	snpcol$Y <- lapply(strsplit(as.character(snpcol$XY), '/'), '[', 2)
	print('Parsing GT to single numeric format')
	snpcol$colID <- as.numeric(snpcol$X) + as.numeric(snpcol$Y)
	print('Accounting for missing GT calls')
	snpcol$colID[is.na(snpcol$colID)] <- 4

	# Identify IDs as Case/Control
	print('Parsing phenotypic data')
	pheno <- read.table(CCfile,header=T)
	print('Identifying Case & Control individuals')
	CaseIDs <- pheno$scanID[toupper(pheno$caseorcontrol)=='CASE']
	ControlIDs <- pheno$scanID[toupper(pheno$caseorcontrol)=='CONTROL']

	# Plot to file
	print(paste('Plotting figure as SNPIntensity_',j,':',i,'_Combined.png',sep=''))
	png(paste('SNPIntensity_',j,'-',i,'_Combined.png',sep=''))
	plot(snp$X, snp$Y, xlim=c(0,xlim.fixed), ylim=c(0,ylim.fixed), col=rainbow(4)[1 + as.numeric(snpcol$colID)], pch=16, main=paste(j,':',i,' Intensity (Case+Control)',sep=''), xlab='X', ylab='Y')
	dev.off()

	print(paste('Plotting figure as SNPIntensity_',j,':',i,'_Case.png',sep=''))
	png(paste('SNPIntensity_',j,'-',i,'_Case.png',sep=''))
	plot(snp$X[CaseIDs], snp$Y[CaseIDs], xlim=c(0,xlim.fixed), ylim=c(0,ylim.fixed), col=rainbow(4)[1 + as.numeric(snpcol$colID[CaseIDs])], pch=16, main=paste(j,':',i,' Intensity (Case Only)',sep=''), xlab='X', ylab='Y')
	dev.off()

	print(paste('Plotting figure as SNPIntensity_',j,':',i,'_Control.png',sep=''))
	png(paste('SNPIntensity_',j,'-',i,'_Control.png',sep=''))
	plot(snp$X[ControlIDs], snp$Y[ControlIDs], xlim=c(0,xlim.fixed), ylim=c(0,ylim.fixed), col=rainbow(4)[1 + as.numeric(snpcol$colID[ControlIDs])], pch=16, main=paste(j,':',i,' Intensity (Control Only)',sep=''), xlab='X', ylab='Y')
	dev.off()

	print('--------------------------------------------------------------------')
}

print('All SNPs plotted successfully.')
print(paste(skippedSNPs,'of',totalSNPs, 'SNPs skipped.',sep=' '))
