library(dplyr)
library(BoolNet)
library(igraph)

# load cluster information and remove lines whithout probe ids
get_cluster_info = function(){
  cluster_table = read.table('./TableS2_ClusterTable.txt', header=TRUE, sep='\t')
  cluster_table = na.omit(cluster_table)
  return(cluster_table)
}


# load data and get data frame of expression of logE data
load('./data_proc/data_proc_diclofenac.Rd')
load('./data_proc/data_proc_naproxen.Rd')
load('./data_proc/data_proc_diuron.Rd')
load('./data_proc/data_proc_mix3.Rd')
load('./data_proc/data_proc_mix13.Rd')
data_list = list(difclofenac = data_proc_diclofenac,
                 diuron = data_proc_diuron,
                 naproxen = data_proc_naproxen,
                 mix3 = data_proc_mix3,
                 mix13 = data_proc_mix13)
rm(data_proc_diclofenac, data_proc_diuron, data_proc_mix13, data_proc_mix3, data_proc_naproxen)

# load cluster information and remove lines whithout probe ids
cluster_table = read.table('./TableS2_ClusterTable.txt', header=TRUE, sep='\t')
cluster_table = na.omit(cluster_table)
# merge by probe IDs with expression data
cluster_table = merge(cluster_table, dic_df, by.x = 'ProbeID', by.y = 'row.names')


# function to generate pointwise and total variances for every toxnode
get_node_var_df = function(cluster_table, list_of_datasets){
  nodes = unique(cluster_table$toxnode)
  node_var_df = data.frame(toxnode = rep('', length(nodes)),
                           var_3 = rep(NA, length(nodes)),
                           var_6 = rep(NA, length(nodes)),
                           var_12 = rep(NA, length(nodes)),
                           var_24 = rep(NA, length(nodes)),
                           var_48 = rep(NA, length(nodes)),
                           var_72 = rep(NA, length(nodes)),
                           var_total = rep(NA, length(nodes)),
                           stringsAsFactors = FALSE)
  
  for (node_index in length(nodes)){
    print(node_index)
    node = nodes[node_index]
    node_df = get_node_timeseries(cluster_table, list_of_datasets, node)
    to_add = c()
    to_add = c(to_add, node)
    for (i in c(1:6)){
      to_add = c(to_add, var(node_df[, i]))
    }
    to_add = c(to_add, var(unlist(node_df)))
    
    node_var_df[node_index, ]  = to_add
    print(to_add)
  }
  return(node_var_df)
}

# calculates pointwise and total variances for all genes
get_gene_var_df = function(cluster_table, list_of_datasets){
  genes = unique(cluster_table$ensembl)
  gene_var_df = data.frame(ensembl = rep('', length(genes)),
                           var_3 = rep(NA, length(genes)),
                           var_6 = rep(NA, length(genes)),
                           var_12 = rep(NA, length(genes)),
                           var_24 = rep(NA, length(genes)),
                           var_48 = rep(NA, length(genes)),
                           var_72 = rep(NA, length(genes)),
                           var_total = rep(NA, length(genes)),
                           stringsAsFactors = FALSE)
  
  for (gene_index in c(1:length(genes))){
    print(gene_index)
    gene = genes[gene_index]
    gene_df = get_gene_timeseries(gene, list_of_datasets)
    to_add = c()
    to_add = c(to_add, toString(gene))
    for (i in c(1:6)){
      to_add = c(to_add, var(gene_df[, i]))
    }
    to_add = c(to_add, var(unlist(gene_df)))
    gene_var_df[gene_index, ]  = to_add
  }
  return(gene_var_df)
}


# given a list of genes, outputs a list of toxnodes which contain at least one of the genes
get_toxnodes = function(genes, cluster_table){
  genes_df = data.frame(ensembl = genes)
  df = merge(genes_df, cluster_table, by.x = 'ensembl', by.y = 'ensembl', all = FALSE)
  return(unique(df$toxnode))
}


# selection step for the genes
# 372 genes, 214 nodes
reduced = gene_var_df[gene_var_df$var_total > (3*sd(gene_var_df$var_total) + mean(gene_var_df$var_total)),]
# calculate max pointwise variance
reduced$pointwise_max_var = apply(reduced[,c("var_3", "var_6", "var_12", "var_24", "var_48", "var_72")], 1, max)
# 350 genes, 205 nodes
# reduced = reduced[reduced$var_total > reduced$pointwise_max_var,]
# calculate mean var over timepoints
reduced$pointwise_mean_var = apply(reduced[,c("var_3", "var_6", "var_12", "var_24", "var_48", "var_72")], 1, mean)
# 124 genes, 91 nodes
reduced = reduced[reduced$var_total / reduced$pointwise_mean_var > 30,]
# merge with cluster information
reduced = merge(cluster_table, reduced, by.x = 'ensembl', by.y = 'ensembl_id', all = FALSE)
# merge with external gene names
reduced = merge(namen, reduced, by.x = 'ensembl_gene_id', by.y = 'ensembl', all = FALSE)
write.table(test_2, file='./reduced_nodes.txt', sep='\t', row.names = FALSE)
reduced_all_genes = cluster_table[cluster_table$toxnode %in% reduced$toxnode, ]


# selection step for the toxnodes
reduce_nodes = function(node_var_df, cluster_info){
  # 2sd: 141 nodes
  reduced_nodes = node_var_df[node_var_df$var_total > (1*sd(node_var_df$var_total) + mean(node_var_df$var_total)),]
  # merge with information about how many genes each cluster contains
  reduced_nodes = merge(reduced_nodes, data.frame(table(cluster_info$toxnode)), by.x = 'toxnode', by.y = 'Var1', all = FALSE)
  
  reduced_nodes$pointwise_max_var = apply(reduced_nodes[, c("var_3", "var_6", "var_12", "var_24", "var_48", "var_72")], 1, max)
  
  reduced_nodes$pointwise_mean_var = apply(reduced_nodes[, c("var_3", "var_6", "var_12", "var_24", "var_48", "var_72")], 1, mean)
  
  reduced_nodes = reduced_nodes[reduced_nodes$var_total/ reduced_nodes$pointwise_mean_var > 2,]
  # merge with cluster information etc
  reduced_nodes = unique(merge(cluster_info[, c('clustername', 'toxnode', 'ensembl')], reduced_nodes, by.x = 'toxnode', by.y = 'toxnode', all = FALSE))
  print(length(unique(reduced_nodes$toxnode)))
  return(reduced_nodes)
}


# given a gene and a list of logE datasets returns a dataframe with time series data for all probes in the gene
get_gene_timeseries = function(ensembl_id, list_of_datasets, type){
  gene_df = data.frame(h_3 = rep(NA, 30),
                       h_6 = rep(NA, 30),
                       h_12 = rep(NA, 30),
                       h_24 = rep(NA, 30),
                       h_48 = rep(NA, 30),
                       h_72 = rep(NA, 30))
  index = 1
  
  if (type == 'control'){
    # datasets
    for (dataset in list_of_datasets){
      dataset_df = na.omit(data.frame(dataset[dataset$genes$ensembl_gene_id == ensembl_id, dataset$targets$type == 'control']$E))
      # replicates
      for (replicate in c('control$', 'control\\.1', 'control\\.2', 'control\\.3')){
        dataset_df_rep = select(dataset_df, matches(replicate))
        dataset_df_rep = dataset_df_rep[, mixedsort(names(dataset_df_rep))]
        # probes
        for (row in nrow(dataset_df_rep)){
          if (ncol(dataset_df_rep) == 6) {
            gene_df[index,] = dataset_df_rep[row, ]
            index = index + 1
          }
        }
      }
    }
    return(na.omit(gene_df))
  }
  else if (type == 'treatment'){
    # datasets
    for (dataset in list_of_datasets){
      dataset_df = na.omit(data.frame(dataset[dataset$genes$ensembl_gene_id == ensembl_id, dataset$targets$type == 'treatment' & dataset$targets$concentration_level == 'C5']$E))
      # replicates
      for (replicate in c('treatment$', 'treatment\\.1', 'treatment\\.2', 'treatment\\.3')){
        dataset_df_rep = select(dataset_df, matches(replicate))
        dataset_df_rep = dataset_df_rep[, mixedsort(names(dataset_df_rep))]
        # probes
        for (row in nrow(dataset_df_rep)){
          if (ncol(dataset_df_rep) == 6) {
            gene_df[index,] = dataset_df_rep[row, ]
            index = index + 1
          }
        }
      }
    }
    return(na.omit(gene_df))
  }
}

# function to get a dataframe containing all time series of a toxnode
get_node_timeseries = function(df, list_of_datasets, node, type) {
  genes_in_node = df[df$toxnode == node, ]$ensembl
  if (length(genes_in_node) == 1) {
    return(get_gene_timeseries(genes_in_node[1], list_of_datasets, type))
  } else {
    out = get_gene_timeseries(genes_in_node[1], list_of_datasets, type)
    for (gene in genes_in_node[2:length(genes_in_node)]) {
      x = get_gene_timeseries(gene, list_of_datasets, type)
      out = merge(out, x, all = TRUE)
    }
    return(out)
  }
}





