require(tidyverse)
require(tidygraph)
require(rlang)
`%NT%` <- function(lhs, rhs) {
      rhs <- enexpr(rhs)
      lhs <- activate(lhs, 'nodes') %>% as.tibble

    # Magrittr does not support inlining so caller
    # _must_ have `%>%` in scope
    expr <- call('%>%', lhs, rhs)
    eval_bare(expr, caller_env())
}
`%ET%` <- function(lhs, rhs) {
      rhs <- enexpr(rhs)
      lhs <- activate(lhs, 'edges') %>% as.tibble

    # Magrittr does not support inlining so caller
    # _must_ have `%>%` in scope
    expr <- call('%>%', lhs, rhs)
    eval_bare(expr, caller_env())
}


get_primary = function(node, path, ...){
	nodes = .N()
	edges = .E()
	dg = nodes$degree[node]
	if(dg <= 1) return(node)
	for(i in 2:dg){
		my_edge = edges$to == node
		node = edges$from[my_edge]
	}
	return(node)
}

classify_axon <- function(tree){
        tree = tree %N>% mutate(id = 1:(tree %NT% nrow),
                        parent = map_bfs_int(node_is_root(), .f = get_primary),
			root_dist=node_distance_from(node_is_root(), weights=length)
			) %E>%
	mutate(primary = .N()$parent[to])

	called_axons = tree %ET% group_by(primary) %>% 
		summarize(axon_length=max(edist), 
			  axon_sum=sum(length),
			  n_dendrites=n()) %>% 
		mutate(axon= axon_length==max(axon_length))

	tree %E>% left_join(called_axons, by=c('primary'))
	#tree %E>% left_join(called_axons, by=c('length', 'neurite'))
}

