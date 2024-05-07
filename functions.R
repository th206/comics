
##################################################################
# function: neighbors - extract all neighbors of a given node list
# input: a node list
# output: a node list
##################################################################
neighbors = function(nodes, network) {
  # test input:
  # network = fullnet
  # nodes = c("cg19693031", "3485-28_2", "BRAIN:Phospholipids in chylomicrons and extremely large VLDL  [mmol/l] [BRAIN]")
  # BRAIN:Phospholipids in chylomicrons and extremely large VLDL  [mmol/l] [BRAIN]
  ix = lapply(nodes, function(x){union(which(network$edges$from == x), which(network$edges$to == x))}) %>% unlist() %>% unique()
  d = union(network$edges$from[ix], network$edges$to[ix])
  d
}
##################################################################
# function: maxneighbors - extract all nodes connected to a given node list
# input: a node list
# output: a node list
# limit: stop if more than limit nodes were found
##################################################################
maxneighbors = function(nodes, network, limit = 0) {
  # test input:
  # network = fullnet  
  # nodes = c("cg19693031", "3485-28_2")
  if (limit == 0) {limit = 1E99}
  nnodes = length(nodes)
  nnodeslast = 0
  nodeslast = nodes
  while ((nnodes > nnodeslast) & (nnodes<=limit)) {
    nodeslast = nodes
    nodes = neighbors(nodes, network)
    nnodeslast = nnodes
    nnodes = length(nodes)
  }
  nodeslast
}
##################################################################
# function: maxneighbors_noSTAT - extract all nodes connected to a given node list
#           but stop growing STAT: nodes
# input: a node list
# output: a node list
# limit: stop if more than limit nodes were found
##################################################################
maxneighbors_noSTAT = function(nodes, network, limit = 0) {
  # test input:
  # network = fullnet  
  # nodes = c("cg19693031", "3485-28_2")
  if (limit == 0) {limit = 1E99}
  nodesin = nodes
  nnodes = length(nodes)
  nnodeslast = 0
  nodeslast = nodes
  while ((nnodes > nnodeslast) & (nnodes<=limit)) {
    nodeslast = nodes
    nnodeslast = nnodes
    nogrow_nodes = nodes[grep("STAT: |GWAS: ", nodes)]
    grow_nodes   = nodes[grep("STAT: |GWAS: ", nodes, invert = TRUE)]
    nodes = c(nogrow_nodes, neighbors(grow_nodes, network))
    nnodes = length(nodes)
  }
  unique(c(nodesin, nodeslast))
}

##################################################################
# function: nodes2network - connect a node list (extract all nodes between them)
# input: a node list
# output: a network
##################################################################
nodes2network = function(nodes, network) {
  
  # test input:
  # network = fullnet  
  # nodes = maxneighbors("cg19693031", fullnet, limit = 100)
  
  ixfrom = lapply(nodes, function(x){which(network$edges$from == x)}) %>% unlist() %>% unique()
  ixto = lapply(nodes, function(x){which(network$edges$to == x)}) %>% unlist() %>% unique()
  ix = intersect(ixfrom, ixto)
  
  iy = lapply(nodes, function(x){which(network$nodes$id == x)}) %>% unlist() %>% unique()
  list ( edges = network$edges[ix,], nodes = network$nodes[iy,])
}
##################################################################
# function: network2node - extract all nodes from a network
# input: a network
# output: a node list
##################################################################
network2nodes = function(network) {
  
  # test input:
  # network = maxneighbors(c("1","2"), fullnet) %>%  nodes2network(fullnet)
  
  nodes = union(network$edges$from, network$edges$to) %>% unlist() %>% unique()
  iy = lapply(nodes, function(x){which(network$nodes$id == x)}) %>% unlist() %>% unique()
  
  network$nodes$id[iy]
}
