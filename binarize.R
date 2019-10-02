
# binarise a node dataframe, returns consens binary vector for the toxnode
binarize_node = function(df, node, bin_method){
  majority = function(col){
    if(sum(col) >= length(col) / 2){
      return(1)
    }
    else {
      return(0)
    }
  }
  
  df_bin = binarizeTimeSeries(df, method = bin_method, edge = 'maxEdge')
  bin = apply(data.frame(df_bin$binarizedMeasurements), 2, majority)
  return(c(node, bin))
}


# given a dataframe with toxnodes, creates a data frame containing a binary vector for each toxnode
binarize_df_nodes = function(df, list_of_datasets, type, method){
  out = data.frame(times = c(3,6,12,24,48,72))
  row.names(out) = out$times
  out$times = NULL
  for (node in unique(df$toxnode)){
    bin = binarize_node(get_node_timeseries(df, list_of_datasets, node, type), node, method)
    out[, toString(node)] = bin[2:7]
  }
  
  return(t(out))
}

# given a dataframe with genes, creates a data frame containing a binary vector for each gene
binarize_df_genes = function(df, list_of_datasets, type, method){
  out = data.frame(times = c(3,6,12,24,48,72))
  row.names(out) = out$times
  out$times = NULL
  for (gene in unique(df$ensembl)){
    bin = binarize_node(get_gene_timeseries(gene, list_of_datasets, type), gene, method)
    out[, toString(gene)] = bin[2:7]
  }
  
  return(t(out))
}


# creates a table that compares control, treatment and consistency binarised data frames to each other. cons_list should contain the dataframes of the consistency vectors
compare_bin_df = function(df_control, df_dic, df_nap, cons_list){
  out = data.frame(name = rep('', length(row.names(df_control))),
                   cons_dic = rep('', length(row.names(df_control))),
                   cons_nap = rep('', length(row.names(df_control))),
                   cons_diu = rep('', length(row.names(df_control))),
                   cons_mix3 = rep('', length(row.names(df_control))),
                   cons_mix13 = rep('', length(row.names(df_control))),
                   control_seq = rep('', length(row.names(df_control))),
                   dic_seq = rep('', length(row.names(df_control))),
                   change_dic = rep('', length(row.names(df_control))),
                   nap_seq = rep('', length(row.names(df_control))),
                   change_nap = rep('', length(row.names(df_control))),
                   same_change = rep('', length(row.names(df_control))),
                   stringsAsFactors = FALSE)
  
  i = 1
  for (name in row.names(df_control)){
    cons_dic = paste(cons_list$diclofenac[name,], collapse = '.')
    cons_nap = paste(cons_list$naproxen[name,], collapse = '.')
    cons_diu = paste(cons_list$diuron[name,], collapse = '.')
    cons_mix3 = paste(cons_list$mix3[name,], collapse = '.')
    cons_mix13 = paste(cons_list$mix13[name,], collapse = '.')
    seq_control = paste(df_control[name,], collapse = '.')
    seq_dic = paste(df_dic[name,], collapse = '.')
    seq_nap = paste(df_nap[name,], collapse = '.')
    if (seq_control == seq_dic){ type_dic = 'same' }
    else if (seq_control != seq_dic){ type_dic = 'change' }
    if (seq_control == seq_nap){ type_nap = 'same' }
    else if (seq_control != seq_nap){ type_nap = 'change' }
    if (seq_dic == seq_nap) { same_change = TRUE }
    else { same_change = FALSE}
    out[i, ] = c(name, cons_dic, cons_nap, cons_diu, cons_mix3, cons_mix13,
                 seq_control, seq_dic, type_dic, seq_nap, type_nap, same_change)
    i = i + 1
  }
  return(out)
}













