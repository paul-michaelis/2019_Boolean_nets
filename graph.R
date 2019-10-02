library(igraph)
library(dplyr)
library(RCy3)
library(gtools)


# load list of ensembl_id / gene names
load("/home/michaelp/Downloads/data_proc/data_proc_naproxen.Rd")


# function to get a list of names, ensemblIDs and probe names for all genes
get_namen = function(){
  load("/home/michaelp/Downloads/data_proc/data_proc_naproxen.Rd")
  namen = data_proc_naproxen$genes[c(6, 7, 9)]
  namen = na.omit(namen)
  namen = distinct(namen)
  rm(data_proc_naproxen)
  return(namen)
}

# takes an ARACNe graph and adds all kinds of vertex attributes
for (i in V(net_embryo)){
  gene_ensembl = V(net_embryo)[i]$name
  # set external_gene_name, toxnode and clustername as vertex attribute
  net_embryo = set_vertex_attr(net_embryo, 'external_gene_name', index = i, namen[namen$ensembl_gene_id == gene_ensembl,]$external_gene_name[1])
  net_embryo = set_vertex_attr(net_embryo, 'toxnode', index = i, cluster_table[cluster_table$ensembl == gene_ensembl,]$toxnode[1])
  
  clustername = toString(cluster_table[cluster_table$ensembl == gene_ensembl,]$clustername[1])
  net_embryo = set_vertex_attr(net_embryo, 'clustername', index = i, clustername)
  
  # set control, diclofenac and naproxen timeseries as attributes
  gene_df = get_gene_timeseries(gene_ensembl, data_list, 'control')
  timeseries = apply(gene_df, 2, mean)
  net_embryo = set_vertex_attr(net_embryo, 'c_t_1', index = i, timeseries[1])
  net_embryo = set_vertex_attr(net_embryo, 'c_t_2', index = i, timeseries[2])
  net_embryo = set_vertex_attr(net_embryo, 'c_t_3', index = i, timeseries[3])
  net_embryo = set_vertex_attr(net_embryo, 'c_t_4', index = i, timeseries[4])
  net_embryo = set_vertex_attr(net_embryo, 'c_t_5', index = i, timeseries[5])
  net_embryo = set_vertex_attr(net_embryo, 'c_t_6', index = i, timeseries[6])

  gene_df = get_gene_timeseries(gene_ensembl, list(data_list$difclofenac), 'treatment')
  timeseries = apply(gene_df, 2, mean)
  net_embryo = set_vertex_attr(net_embryo, 'd_t_1', index = i, timeseries[1])
  net_embryo = set_vertex_attr(net_embryo, 'd_t_2', index = i, timeseries[2])
  net_embryo = set_vertex_attr(net_embryo, 'd_t_3', index = i, timeseries[3])
  net_embryo = set_vertex_attr(net_embryo, 'd_t_4', index = i, timeseries[4])
  net_embryo = set_vertex_attr(net_embryo, 'd_t_5', index = i, timeseries[5])
  net_embryo = set_vertex_attr(net_embryo, 'd_t_6', index = i, timeseries[6])
  # 
  gene_df = get_gene_timeseries(gene_ensembl, list(data_list$naproxen), 'treatment')
  timeseries = apply(gene_df, 2, mean)
  net_embryo = set_vertex_attr(net_embryo, 'n_t_1', index = i, timeseries[1])
  net_embryo = set_vertex_attr(net_embryo, 'n_t_2', index = i, timeseries[2])
  net_embryo = set_vertex_attr(net_embryo, 'n_t_3', index = i, timeseries[3])
  net_embryo = set_vertex_attr(net_embryo, 'n_t_4', index = i, timeseries[4])
  net_embryo = set_vertex_attr(net_embryo, 'n_t_5', index = i, timeseries[5])
  net_embryo = set_vertex_attr(net_embryo, 'n_t_6', index = i, timeseries[6])
}
createNetworkFromIgraph(net_embryo, '125_graph')

# import genes from Wibkes hypotheses
wibke = read.csv('/home/michaelp/Downloads/aracne/wibke/wibke_gene_2.csv', header = TRUE, sep = '\t', colClasses = c(rep('character', 3)))
wibke_merged = merge(wibke, namen, by.x = 'gene_name', by.y = 'external_gene_name')
wibke_split = split(wibke_merged, wibke_merged$group)


# get subgraphs with shortest paths inside of boxes
for (x in wibke_split){
  path_lengths = c()
  node_list = c()
  for (i in x$ensembl_gene_id){
    for (j in x$ensembl_gene_id){
      if (!i == j){
        paths = all_shortest_paths(net_embryo, i, to = j)$res
        best_path = paths[[1]]
        for (path in paths){
          path_score = mean(E(induced_subgraph(net_embryo, path))$MI)
          if (path_score > mean(E(induced_subgraph(net_embryo, best_path))$MI)){
            best_path = path
          }
        }
        path_lengths = c(path_lengths, length(best_path -1 ))
        node_list = c(node_list, best_path)
      }
    }
  }
  subgraph = induced_subgraph(net_embryo, node_list)
  list_parameters(subgraph, x$group[1])
  
  # set custom attributes for wibkes genes
  for (node in V(subgraph)){
  for (i in 1:nrow(wibke_merged)){
    if (vertex_attr(subgraph, 'name', index = node) == wibke_merged[i,]$ensembl_gene_id){
      subgraph = set_vertex_attr(subgraph, 'up_down', index = node, wibke_merged[i,]$up_down)
      subgraph = set_vertex_attr(subgraph, 'group', index = node, wibke_merged[i,]$group)
      subgraph = set_vertex_attr(subgraph, 'wibke', index = node, TRUE)
      }
    }
  }

  sub_d = subset(dicnap_df, ensembl_gene_id %in% names(V(subgraph)))
  for (node in V(subgraph)){
    for (i in 1:nrow(sub_d)){
      if (vertex_attr(subgraph, 'name', index = node) == sub_d[i,]$ensembl_gene_id){
        for (col in names(sub_d)){
          subgraph = set_vertex_attr(subgraph, col, index = node, sub_d[[col]][i])
        }
      }
    }
  }

  # send to cytoscape 
  createNetworkFromIgraph(subgraph, paste('best_shrtst_pths_in_box_', toString(x$group[1]), '_n'))
}


# create subgraphs with shortest paths between boxes
# type: 'cytoscape' to send graph to cytoscape, 'graph' to return igraph object
combine_boxes = function(box_1, box_2, type){
  node_list = c()
  for (i in box_1$ensembl_gene_id){
    for (j in box_2$ensembl_gene_id){
      if (!i == j){
        paths = all_shortest_paths(net_embryo, i, to = j)$res
        best_path = paths[[1]]
        for (path in paths){
          path_score = mean(E(induced_subgraph(net_embryo, path))$MI)
          print(path_score)
          if (path_score > mean(E(induced_subgraph(net_embryo, best_path))$MI)){
            best_path = path
          }
        }
        node_list = c(node_list, best_path)
      }
    }
  }
  subgraph = induced_subgraph(net_embryo, node_list)
  
  for (node in V(subgraph)){
    for (i in 1:nrow(wibke_merged)){
      if (vertex_attr(subgraph, 'name', index = node) == wibke_merged[i,]$ensembl_gene_id){
        subgraph = set_vertex_attr(subgraph, 'up_down', index = node, wibke_merged[i,]$up_down)
        subgraph = set_vertex_attr(subgraph, 'group', index = node, wibke_merged[i,]$group)
        subgraph = set_vertex_attr(subgraph, 'wibke', index = node, TRUE)
      }
    }
  }
  
  if (type == 'cytoscape'){
    createNetworkFromIgraph(subgraph, paste('best_combined_boxes_', toString(box_1$group[1]), '_', toString(box_2$group[1])))
  }
  else if (type == 'graph'){
    return(subgraph)
  }
}

combine_boxes(wibke_split$'1', wibke_split$'1.1', 'cytoscape')
combine_boxes(wibke_split$'2', wibke_split$'2.1')
combine_boxes(wibke_split$'2', wibke_split$'2.2')
combine_boxes(wibke_split$'2', wibke_split$'2.3')
subgraph = combine_boxes(wibke_split$'2.3', wibke_split$'2.4', 'graph')


# combine one box with multiple other boxes
combine_box_with_boxes = function(box_1, boxes){
  node_list = c()
  for (i in box_1$ensembl_gene_id){
    for (box in boxes){
      #print(box)
      for (j in box$ensembl_gene_id){
        if (!i == j){
          paths = all_shortest_paths(net_embryo, i, to = j)$res
          best_path = paths[[1]]
          for (path in paths){
            path_score = mean(E(induced_subgraph(net_embryo, path))$MI)
            #print(path_score)
            if (path_score > mean(E(induced_subgraph(net_embryo, best_path))$MI)){
              best_path = path
            }
          }
          node_list = c(node_list, best_path)
        }
      }
    }
  }
  subgraph = induced_subgraph(net_embryo, node_list)
  
  return(subgraph)
}


# this function takes the nodes of a graph and plots their expression profiles
plot_subgraph = function(graph, data_file){
  namen = get_namen()
  ids_wibke = namen[namen$ensembl_gene_id %in% as_ids(V(subgraph))[as_ids(V(subgraph)) %in% wibke_merged$ensembl_gene_id],]$ProbeID
  ids_else = namen[namen$ensembl_gene_id %in% as_ids(V(subgraph))[!(as_ids(V(subgraph)) %in% wibke_merged$ensembl_gene_id)],]$ProbeID
  
  df = data_file[data_file$genes$ProbeID %in% ids_wibke, data_file$targets$concentration_level == 'C5' & data_file$targets$type == 'treatment']
  df_df = as.data.frame(df$E)
  new_df = data.frame(row.names(df_df))
  for (time in data_file$targets$time_hpe){
    time_str = toString(time)
    new_df[, time_str] = rowMeans(select(df_df, matches(time_str)))
  }
  row.names(new_df) = new_df[,1]
  new_df[,1] = NULL
  new_df = new_df[, mixedsort(names(new_df))]
  new_df = as.data.frame(t(new_df))
  new_df$id <- 1:nrow(new_df)
  plot_data_wibke <- melt(new_df,id.var="id")
  plot_data_wibke$type = 'wibke'
  
  df = data_file[data_file$genes$ProbeID %in% ids_else, data_file$targets$concentration_level == 'C5' & data_file$targets$type == 'treatment']
  df_df = as.data.frame(df$E)
  new_df = data.frame(row.names(df_df))
  for (time in data_file$targets$time_hpe){
    time_str = toString(time)
    new_df[, time_str] = rowMeans(select(df_df, matches(time_str)))
  }
  row.names(new_df) = new_df[,1]
  new_df[,1] = NULL
  new_df = new_df[, mixedsort(names(new_df))]
  new_df = as.data.frame(t(new_df))
  new_df$id <- 1:nrow(new_df)
  plot_data_else <- melt(new_df,id.var="id")
  plot_data_else$type = 'else'
  
  plot_data = merge(plot_data_else, plot_data_wibke, all = TRUE)
  
  return(plot_data)

  p = ggplot(plot_data, aes(x=(id), y=value, group=variable, color=type)) +
    geom_line(aes(lty=variable), linetype = 'solid', size = 0.1) +
    theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'),
          axis.text.y = element_blank(), axis.title=element_text(size=9), axis.text.x = element_text(size = 8),
          legend.position = 'none') +
    labs(x = 'time (hpe)', y = 'logFC expression') +
    scale_x_discrete(breaks = c('1', '2', '3', '4', '5', '6'), labels = c('3', '6', '12', '24', '48', '72'))
  return(p)
}


p = ggplot(plot_data, aes(x=factor(id), y=value, group=variable, color=type)) +
  geom_line(aes(lty=variable), linetype = 'solid', size = 0.1) +
  theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'),
        axis.text.y = element_blank(), axis.title=element_text(size=9), axis.text.x = element_text(size = 9),
        legend.position = 'none') +
  labs(y = 'logFC expression') +
  scale_x_discrete(name = 'time (hpe)', breaks = c('1', '2', '3', '4', '5', '6'), labels = c('3', '6', '12', '24', '48', '72'),
                   expand = c(0.03,0.03))
p

ggsave(plot = p, filename = '/home/michaelp/Downloads/latex_template (copy)/fig/aracne_vectors_2.png', width = 5, height = 4, units = 'cm')
