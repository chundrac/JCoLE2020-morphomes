require(phytools)
require(rstan)
#require(stringdist)

get.ngrams <- function(w,k) {
  ngrams <- c()
  ww <- paste('#',lemma.data[lemma.data$Lemma==w,]$Stem,'#',sep=' ')
  ww <- unlist(strsplit(ww,' '))
  for (i in 1:(length(ww)-k+1)) {
    ngrams <- c(ngrams,paste(ww[i:(i+k-1)],collapse=' '))
  }
  labels <- rep(w,length(ngrams))
  return(cbind(labels,ngrams))
}

gen.ngrams <- function(x,k) {
  all.ngrams <- NULL
  for (w in x) {
    all.ngrams <- rbind(all.ngrams,get.ngrams(w,k))
  }
  ngram.table <- xtabs( ~ labels + ngrams,all.ngrams)
  #ngram.table <- ngram.table[,which(colSums(ngram.table)!=1)]
  inds.to.remove <- c()
  for (i in 1:(ncol(ngram.table)-1)) {
    for (j in (i+1):ncol(ngram.table)) {
      if (all(ngram.table[,i]==ngram.table[,j])) {
        inds.to.remove <- c(inds.to.remove,j)
      }
    }
  }
  ngram.table <- ngram.table[,-inds.to.remove]
  return(ngram.table)
}

lemma.freq <- read.csv('lemma_frequencies.csv')
lemma.stem <- read.csv('stemmed_verbs.csv')
lemma.data <- merge(lemma.stem,lemma.freq,by='Lemma')
rownames(lemma.data) <- lemma.data$Lemma

trees <- read.nexus('thinned_tree_sample.nex')

trees <- trees[floor(seq(1,length(trees),length.out=50))]

morphome.data <- read.csv('morphomesbylexemeICRomance.csv',row.names=1)

inds.to.exclude <- sort(unique(c(
  which(apply(morphome.data,2,function(x){length(xtabs( ~ x))})==1),
  which(apply(morphome.data,2,function(x){length(xtabs( ~ unlist(lapply(strsplit(rownames(morphome.data)[which(!is.na(x))],'_'),function(x){x[1]}))))})==1),
  which(apply(morphome.data,2,function(x){length(na.omit(x))}) < 5)
)))

morphome.data <- morphome.data[,-inds.to.exclude]

lemma.type = as.factor(gsub("[[:upper:]]+", "", colnames(morphome.data)))
morphome.data <- morphome.data[,which(lemma.type %in% lemma.data$Lemma)]

lemma.type = as.factor(gsub("[[:upper:]]+", "", colnames(morphome.data)))
morphome.type = as.factor(gsub("[[:lower:]]+", "", colnames(morphome.data)))
conj.type = lemma.data[lemma.type,]$Conjugation

log.lemma.freq <- log(lemma.data[sort(as.character(unique(lemma.type))),]$freqTFTL)

lemma.unique <- as.character(sort(unique(lemma.type)))
lemma.ngrams <- gen.ngrams(lemma.unique,3)

embed = read.delim('latin_verb_embeddings.txt',header=F,row.names=1,sep=' ')
embed <- as.matrix(embed)
sim <- embed / sqrt(rowSums(embed * embed))
sim <- sim %*% t(sim)
sem.dist <- as.dist(1 - sim)
sem.dist <- as.matrix(sem.dist)[lemma.unique,lemma.unique]

save.image('input_data.Rdata')
