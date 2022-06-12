cog <- read.csv('ASJP_romance_binarized.csv',sep=',',row.names=1)

sc <- read.csv('soundcomp_romance_binarized.csv',sep='\t',row.names=1)

merged <- merge(cog,sc,by='row.names')

rownames(merged) <- merged$Row.names

merged <- merged[,c(2:ncol(merged))]

merged <- merged[,colSums(merged)!=nrow(merged)]

merged <- merged[,colSums(merged)!=0]

merged <- merged[,colSums(merged)!=1]

write.csv(file='all_data_binarized.csv',merged,quote=F,row.names=T)