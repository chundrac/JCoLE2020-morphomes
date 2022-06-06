require(phytools)
require(rstan)
require(cmdstanr)

t = as.integer(commandArgs(trailingOnly = TRUE)[1])

load('input_data.Rdata')

lemma.type <- factor(lemma.type,levels=rev(lemma.freq$Lemma))
morphome.type <- factor(morphome.type,levels=c('N','LU','P'))

mod <- cmdstan_model('model.stan')

fit.list <- list()

tree <- trees[[t]]
tree <- reorder.phylo(tree,'pruningwise')
morphome.data.curr <- morphome.data[tree$tip.label,]
bin.states <- NULL
for (d in 1:ncol(morphome.data)) {
  bin.states.d <- to.matrix(as.character(morphome.data.curr[,d]),seq=c('no','yes'))
  bin.states.d[rowSums(bin.states.d)==0,] <- c(1,1)
  bin.states <- cbind(bin.states,bin.states.d)
}
bin.states <- rbind(bin.states,matrix(1,nrow=tree$Nnode,ncol=ncol(bin.states)))
parent <- tree$edge[,1]
child <- tree$edge[,2]
b.lens <- tree$edge.length/1000
N <- length(unique(c(parent,child)))
T <- length(child[which(!child %in% parent)])
tip.lik <- bin.states
data.list <- list(N=N,
                  T=T,
                  B=length(parent),
                  brlen=b.lens,
                  child=child,
                  parent=parent,
                  tiplik=tip.lik,
                  D=ncol(tip.lik)/2,
                  L=max(as.numeric(lemma.type)),
                  J=max(as.numeric(morphome.type)),
                  lemma_id=as.numeric(lemma.type),
                  pattern_id=as.numeric(morphome.type))
fit <- mod$sample(data=data.list,chains=3,parallel_chains=3,refresh = 200)
stanfit <- rstan::read_stan_csv(fit$output_files())
fit.list <- append(fit.list,stanfit)
fit.full <- sflist2stanfit(fit.list)
save.image(paste('output/model_fit_',t,'.Rdata',sep=''))
