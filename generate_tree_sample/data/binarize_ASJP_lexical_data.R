taxa.merged <- read.csv('ASJP_merged.csv',sep='\t')

ASJP.wide <- read.csv('asjp_multistate_coding.csv')

ASJP.merged <- merge(taxa.merged,ASJP.wide,by.x='ASJP_ID',by.y='ID')

rownames(ASJP.merged) <- ASJP.merged$Oxford_ID

ASJP.merged <- ASJP.merged[,c(8:47)]

require(phytools)

ASJP.binarized <- NULL

for (i in 1:ncol(ASJP.merged)) {
	binarized <- to.matrix(ASJP.merged[,i],seq=levels(ASJP.merged[,i]))
	ASJP.binarized <- cbind(ASJP.binarized,binarized)
}

rownames(ASJP.binarized) <- rownames(ASJP.merged)

ASJP.binarized <- ASJP.binarized[,colSums(ASJP.binarized)!=0]

write.csv(file='ASJP_romance_binarized.csv',ASJP.binarized,quote=F,sep='\t',row.names=T)