# consistency check

reduced_bin_list = list()
for (data in data_list){
  print(data$targets$substance[1])
  df = binarize_df_nodes(reduced_35, list(data), type = 'control')
  df = unique(df)
  rownames(df) = apply(df, 1, paste, collapse = '')
  reduced_bin_list[[toString(data$targets$substance[1])]] = df
}

for (i in reduced_bin_list){
  for (j in reduced_bin_list){
    print(table(i == j))
  }
}

# binarize toxnode time series
reduced_bin_control_35 = binarize_df_nodes(reduced_35, data_list, type = 'control', method = 'edgeDetector')
reduced_bin_treatment_dic_35 = binarize_df_nodes(reduced_35, list(data_list$difclofenac), type = 'treatment', method = 'edgeDetector')
reduced_bin_treatment_nap_35 = binarize_df_nodes(reduced_35, list(data_list$naproxen), type = 'treatment', method = 'edgeDetector')

# filter unique vectors
reduced_bin_control_unique_35 = unique(reduced_bin_control_35)
reduced_bin_treatment_unique_dic_35 = unique(reduced_bin_treatment_dic_35)
reduced_bin_treatment_unique_nap_35 = unique(reduced_bin_treatment_nap_35)

# name rownames as the binary vector
rownames(reduced_bin_control_unique_35) = apply(reduced_bin_control_unique_35, 1, paste, collapse = '')
rownames(reduced_bin_treatment_unique_dic_35) = apply(reduced_bin_treatment_unique_dic_35, 1, paste, collapse = '')
rownames(reduced_bin_treatment_unique_nap_35) = apply(reduced_bin_treatment_unique_nap_35, 1, paste, collapse = '')

# reconstruct PBN
net_control = reconstructNetwork(reduced_bin_control_unique_35, returnPBN = TRUE, readableFunctions = TRUE)
net_dic = reconstructNetwork(reduced_bin_treatment_unique_dic_35, returnPBN = TRUE, readableFunctions = TRUE)
net_nap = reconstructNetwork(reduced_bin_treatment_unique_nap_35, returnPBN = TRUE, readableFunctions = TRUE)

# does the same as above for the 124 genes
reduced_bin_control_124 = binarize_df_genes(reduced_124, data_list, type = 'control')
reduced_bin_treatment_dic_124 = binarize_df_genes(reduced_124, list(data_list$difclofenac), type = 'treatment')
reduced_bin_treatment_nap_124 = binarize_df_genes(reduced_124, list(data_list$naproxen), type = 'treatment')
reduced_bin_control_unique_124 = unique(reduced_bin_control_124)
reduced_bin_treatment_unique_dic_124 = unique(reduced_bin_treatment_dic_124)
reduced_bin_treatment_unique_nap_124 = unique(reduced_bin_treatment_nap_124)
rownames(reduced_bin_control_unique_124) = apply(reduced_bin_control_unique_124, 1, paste, collapse = '')
rownames(reduced_bin_treatment_unique_dic_124) = apply(reduced_bin_treatment_unique_dic_124, 1, paste, collapse = '')
rownames(reduced_bin_treatment_unique_nap_124) = apply(reduced_bin_treatment_unique_nap_124, 1, paste, collapse = '')
net_control_124 = reconstructNetwork(reduced_bin_control_unique_124, returnPBN = TRUE, readableFunctions = TRUE)
net_dic_124 = reconstructNetwork(reduced_bin_treatment_unique_dic_124, returnPBN = TRUE, readableFunctions = TRUE)
net_nap_124 = reconstructNetwork(reduced_bin_treatment_unique_nap_124, returnPBN = TRUE, readableFunctions = TRUE)

# make consistency dataframes
control_35_dic = binarize_df_nodes(reduced_35, data_list[c(2,3,4,5)], type = 'control', method = 'edgeDetector')
control_35_diu = binarize_df_nodes(reduced_35, data_list[c(1,3,4,5)], type = 'control', method = 'edgeDetector')
control_35_nap = binarize_df_nodes(reduced_35, data_list[c(1,2,4,5)], type = 'control', method = 'edgeDetector')
control_35_mix3 = binarize_df_nodes(reduced_35, data_list[c(1,2,3,5)], type = 'control', method = 'edgeDetector')
control_35_mix13 = binarize_df_nodes(reduced_35, data_list[c(1,2,3,4)], type = 'control', method = 'edgeDetector')
consistency_checklist_35 = list(diclofenac = control_35_dic, diuron = control_35_diu, naproxen = control_35_nap, mix3 = control_35_mix3, mix13 = control_35_mix13)

for (name in names(consistency_checklist)){
  df = unique(consistency_checklist[[name]])
  rownames(df) = apply(df, 1, paste, collapse = '')
  consistency_checklist[[name]] = df
}

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
