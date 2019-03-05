library(data.table)

modelID <- c('Additive','Heterozygous','Dominant','Recessive')
names(modelID) <- c('ADD','HET','DOM','REC')

args <- commandArgs(trailingOnly=T)
prefix <- args[1]
threshold <- 1e-6

dat <- fread(paste(prefix,'_ALL.ps',sep=''))

dat <- dat[order(dat$P),]
dat <- dat[dat$P < threshold,]

dat$Gene <- ''
dat$SNP_ID <- unlist(lapply(strsplit(as.character(dat$SNP),':'),'[',1))
dat$Location <- paste(dat$CHR,dat$BP,sep=':')
dat$MinP <- dat$P
dat$Model <- as.character(dat$MODEL)
dat$Model[dat$Model=='ADD'] <- 'Additive'
dat$Model[dat$Model=='HET'] <- 'Heterozygous'
dat$Model[dat$Model=='DOM'] <- 'Dominant'
dat$Model[dat$Model=='REC'] <- 'Recessive'

dat <- dat[,c('Gene','SNP_ID','Location','MinP','Model')]

write.csv(dat, file=paste(prefix,'top.csv',sep='_'), row.names=F)

system(paste("sed 's/e-/x10-/g' ",prefix,"_top.csv > remove.tmp",sep=''))
system(paste("mv remove.tmp ",prefix,"_top.csv",sep=''))
