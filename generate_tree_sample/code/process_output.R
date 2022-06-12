require(phytools)
require(phangorn)
require(ape)
require(tikzDevice)

j = 1
tree.list <- list()
for (chain in 1:4) {
  
  trees.text <- readLines(paste('output/romance_run_',chain,'.t',sep=''))
  trees.text <- trees.text[-1]
  
  for (i in 1:length(trees.text)) {
    if (i > 5000) {
      tree <- read.newick(textConnection(strsplit(trees.text[[i]],'\t')[[1]][5]))
      tree$tip.label <- gsub("\\[[^\\]]+\\]",'',tree$tip.label,perl=T)
      #for (k in 1:length(tree$edge.length)) {
      #  if (tree$edge[k,2] <= length(tree$tip.label)) {
      #    tree$edge.length[k] <- tree$edge.length[k] + 200
      #  }
      #}
      if (length(tree$tip.label)==73) {
        tree$node.label <- NULL
        tree.list[[j]] <- collapse.singles(tree)
        j <- j + 1
      }
    }
  }
}

ntip <- max(unlist(lapply(tree.list,function(x) {length(x$tip.label)})))

tree.list.thinned <- list()
for (i in 1:length(tree.list)) {
  tree <- tree.list[[i]]
  if (length(tree$tip.label)==73) {
    tree$node.label <- NULL
    tree.list.thinned[[i]] <- tree
  }
}

tree.list <- tree.list.thinned

trees <- tree.list

tree.list.thinned <- list()
for (i in 1:length(trees)) {
  tree <- trees[[i]]
  if (length(tree$tip.label)!=0) {
    tree$node.label <- NULL
    tree.list.thinned[[i]] <- tree
  }
}

tree.list <- tree.list.thinned

trees <- tree.list

class(trees) <- 'multiPhylo'

plot.trees <- list()
for (i in 1:length(trees)) {
  t <- trees[[i]]
  t$tip.label <- gsub('_','',t$tip.label)
  plot.trees[[i]] <- t
}

class(plot.trees) <- 'multiPhylo'

thinned.trees <- trees[floor(seq(0,length(trees),length.out=500))]

write.nexus(thinned.trees,file='thinned_tree_sample.nex')

