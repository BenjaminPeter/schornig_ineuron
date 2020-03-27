suppressPackageStartupMessages({
require(tidyverse)
require(tidygraph)
require(ggraph)
require(reshape2)
source("scripts/utils.R")
})


get_tree <- function(data){
    vertices <- data %>% select(PathID)
    tree <- data %>% select(from=StartsOnPath, to=PathID, length=PathLength) %>%
        mutate(from=coalesce(from, -1L)) %>% as_tbl_graph()

    tree <- activate(tree, nodes) %>% 
        mutate(degree=bfs_dist(root=node_is_root()),
               root_dist=node_distance_from(node_is_root(), weights=length),
               is_leaf=node_is_leaf(), root=node_is_root()) %E>%
        mutate(edist= .N()$root_dist[to], degree=.N()$degree[to]) %>%
        group_by(from) %>% mutate(local_rank = rank(-edist)) %>% ungroup %>%
        mutate(local_rank = ifelse(.N()$root[from], 1L, local_rank))

    tree <- tree %E>% mutate(is_secondary=as.integer(local_rank>1 ) ) %N>%
        mutate(n=as.integer(node_distance_from(node_is_root(), weights=is_secondary)))
}

make_tree_df <- function(data){
    data2 <- data %>% select(cell, ID, filament, staining, batch, day, ID, cellline, species,
                            from=parent, to=neurite, length=length) %>%
        mutate(from=coalesce(from, -1L)) %>% 
        mutate(from=ifelse(from>=0,filament*1000 + from, -1),
                            to = ifelse(to >= 0,filament * 1000 + to, -1)) %>%
        select(-filament)
            
}

get_tree2 <- function(data){
                    
    vertices <- data %>% select(from)

    tree <- as_tbl_graph(data)
    tree <- activate(tree, nodes) %>% 
        mutate(degree=bfs_dist(root=node_is_root()),
               root_dist=node_distance_from(node_is_root(), weights=length),
               is_leaf=node_is_leaf(), root=node_is_root()) %E>%
        mutate(edist= .N()$root_dist[to], degree=.N()$degree[to]) %>%
        group_by(from) %>% mutate(local_rank = rank(-edist)) %>% ungroup %>%
        mutate(local_rank = ifelse(.N()$root[from], 1L, local_rank))

    tree <- tree %E>% mutate(is_secondary=as.integer(local_rank>1 ) ) %N>%
        mutate(n=as.integer(node_distance_from(node_is_root(), weights=is_secondary)))
}

get_tree_neurites <- function(data){
    vertices <- data %>% select(from)

    tree <- as_tbl_graph(data)
    tree <- activate(tree, nodes) %>% 
        morph(to_components) %>%
            mutate(degree=bfs_dist(root=node_is_root()),
                   root_dist=node_distance_from(node_is_root(), weights=length),
                   is_leaf=node_is_leaf(), root=node_is_root()) %E>%
            mutate(edist= .N()$root_dist[to], degree=.N()$degree[to]) %>%
            group_by(from) %>% mutate(local_rank = rank(-edist)) %>% ungroup %>%
            mutate(local_rank = ifelse(.N()$root[from], 1L, local_rank)) %E>%
            mutate(is_secondary=as.integer(local_rank>1 ) ) %N>%
            mutate(n=as.integer(node_distance_from(node_is_root(), weights=is_secondary))) %>%
        unmorph()

}

get_branches <- function(tree){
    n_max <- tree %NT% select(n) %>% max
    branches <- tree %NT% filter(n==0, is_leaf) 

    if(n_max == 0) return(branches)
    for(i in 1:n_max){
        tree <-  tree %E>% 
            mutate(length = ifelse(.N()$n[to] < i, 0, length)) %N>%
            mutate(root_dist=node_distance_from(node_is_root(), weights=length))

        branches <- tree %NT% filter(n==i, is_leaf) %>% bind_rows(branches)
    }

    branches %>% select(name, root_dist, n)
}

write_tree <- function(tree, file, ...){
    tree %ET% select(from, to, length, degree) %>% write_csv(file, ...)
}

write_branches <- function(branches, file, ...){
    branches %>% rename(to=name, length=root_dist, degree=n)  %>%
        write_csv(file, ...)
}


get_neurites <- function(tree){
    tree  %>% 
        morph(to_components) %>% crystallize() %>% select(graph)
}

#read a csv file with all hoc files for all neurons
readtrees  <- function(fname="scripts/merged_ublind.csv"){
  data <- read.csv(fname)
  data2 <- make_tree_df(data) 
   trees <- data2 %>% group_by(cell, ID, staining, batch, day, cellline, species) %>%
	        do(tree=get_tree2(.))                                      
}
lookatdegrees <- function(trees){
  df_degree <- trees %>% 
    do(x=get_degrees(.$tree), y=get_tree_summaries(.$tree),
       z=get_tree_branch_lengths(.$tree)) %>% 
    bind_cols(trees) %>% 
    unnest(y) 
}

classify_polarity <- function(df){
    df=df %>%
          mutate(unipolar = polarity ==1) %>%
          mutate(bipolar = polarity==2) %>%
          mutate(multipolar = polarity>2) %>%
          mutate(polarity = factor(NA, levels = c("unipolar", "bipolar", "multipolar")))

      df$polarity[df$unipolar]="unipolar"
      df$polarity[df$bipolar]="bipolar"
      df$polarity[df$multipolar]="multipolar"
      df %>% select(-unipolar, -bipolar, -multipolar) 
}

getpolarity <- function(df_degree){
  df_polarity <- unnest(df_degree, x) %>% 
    filter(degree == 1) %>%
  mutate(unipolar = n ==1) %>%
  mutate(bipolar = n==2) %>%
  mutate(multipolar = n>2) %>%
  
  mutate(polarity = factor(NA, levels = c("unipolar", "bipolar", "multipolar")))
  df_polarity$polarity[df_polarity$unipolar]="unipolar"
  df_polarity$polarity[df_polarity$bipolar]="bipolar"
  df_polarity$polarity[df_polarity$multipolar]="multipolar"
  df_polarity %>% select(cell, ID, polarity) %>% 
    right_join(df_degree)
}

getpolarity_branches <- function(df_branches){
  df_polarity <-df_branches %>%
    filter(degree == 0) %>%
    group_by(cell, ID, staining, batch, day, cellline, species) %>%
    summarize(polarity= sum(degree==0)) %>%
  mutate(unipolar = polarity == 1) %>%
  mutate(bipolar = polarity == 2) %>%
  mutate(multipolar = polarity >=2)  
  
  df_polarity %>% 
  mutate(polarity = factor(NA, levels = c("unipolar", "bipolar", "multipolar")))
  df_polarity$polarity[df_polarity$unipolar]="unipolar"
  df_polarity$polarity[df_polarity$bipolar]="bipolar"
  df_polarity$polarity[df_polarity$multipolar]="multipolar"
  df_polarity %>% select(cell, ID, polarity)  %>%
      inner_join(df_branches)
}

get_axon_table <- function(tree){
    tree %ET% select(length, primary, axon_length, n_dendrites, axon, axon_sum, degree)
}
  
basic_analysis <- function(fname){
    defaults = list(n_edges=0, length=0, max_dist=0, max_degree=0, n=0, 
                    polarity='unipolar')
    trees <- readtrees(fname)

    trees2 = trees %>% rowwise %>% do(tree = classify_axon(.$tree))
    trees = bind_cols(trees %>% select(-tree), trees2)


    df_axon <- trees %>% do(axon=get_axon_table(.$tree)) %>% 
        bind_cols(trees, .) %>% 
        unnest(axon) %>% 
        group_by(cell, ID, staining, batch, day, cellline, species) %>% 
        mutate(polarity = n_distinct(primary)) %>%
        classify_polarity()
    saveRDS(df_axon, "df_axon.rds")
    
    df_primary <- trees %>% do(axon=get_axon_table(.$tree)) %>% 
      bind_cols(trees, .) %>% 
      unnest(axon) %>% 
      group_by(cell, ID, staining, batch, day, cellline, species) %>% 
      mutate(polarity = n_distinct(primary)) %>%
      classify_polarity()%>%
      mutate(number_primary = n_distinct(primary)) %>%
      mutate()
    saveRDS(df_primary, "df_primary.rds")


    branches <- trees %>% do(branches=get_branches(.$tree)) 
    df_branches <- branches %>% bind_cols(trees %>% select(-tree)) %>%
        unnest(branches) %>%
        select(cell, ID,  staining, batch, day, cellline, species, length=root_dist, degree=n) %>%
        getpolarity_branches()
    saveRDS(df_branches, "df_branches.rds")
    df_degree <- lookatdegrees(trees)
    saveRDS(df_degree, "df_degree.rds")

        
    df_neurites <- neurites %>% rowwise %>% 
        do(x=get_degrees(.$tree), y=get_tree_summaries(.$tree),
           z=get_tree_branch_lengths(.$tree)) %>% 
        bind_cols(neurites) %>% 
        unnest(y)  %>%
        group_by(cell, ID, staining, batch, day, cellline, species) %>% mutate(neurite=rank(-length))
   
    df_branches <- branches %>% bind_cols(trees %>% select(-tree)) %>%
        unnest(branches) %>%
        select(cell, ID, staining, batch, day, cellline, species, length=root_dist, degree=n) 

}

get_tree_summaries <- function(tree){
    tree %E>% as.tibble %>% 
        summarize(length=sum(length), n_edges=n(), 
                  max_dist=max(edist), max_degree=max(degree))
}

get_degrees <- function(tree){
    tree %ET% group_by(degree) %>% summarize(n=n())
}

get_tree_branch_lengths <- function(tree){
    tree %E>% as.tibble %>% select(length)
}



basic_analysis(fname = "scripts/merged_ublind.csv" )



