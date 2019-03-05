library(GenABEL)

datList <- c('PLACEHOLDER_Add.ps', 'PLACEHOLDER_Het.ps', 'PLACEHOLDER_Rec.ps', 'PLACEHOLDER_Dom.ps')

for (i in datList) {
	print(paste('Processing',i))
	dat <- read.table(i, header=F)
	
	# Estimated GIF (fast)
	estL <- estlambda(dat$V3)
	print(paste('Estimated Lambda:', estL$estimate))

	# Calculated GIF (slow)
	chisq <- qchisq(1-dat$V3,1)
	lambda <- median(chisq)/qchisq(0.5,1)
	print(paste('Calculated Lambda:', lambda))
	print('')
}
