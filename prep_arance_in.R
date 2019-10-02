# creates an input file for ARACNe

load("./data_proc/data_proc_diclofenac.Rd")
file = data_proc_diclofenac
namen = data_list$difclofenac$genes[c(4, 7, 9)]

load("./data/single_substance_data/data_logfc_diclofenac.Rd")
load("./data/single_substance_data/data_logfc_naproxen.Rd")

both_merged = na.omit(namen)
for (file2 in list(data_list$naproxen, data_list$difclofenac)){
  namen2 = file2$genes[4]
  expr2 = file2[file2$genes$ensembl_gene_id %in% unique(reduced_124$ensembl_gene_id), file2$targets$type == 'control' |
                  (file2$targets$type == 'treatment' &
                     (file2$targets$concentration_level == 'C4' | file2$targets$concentration_level == 'C5'))]$E
  merged2 = merge(namen2, expr2, by.x = 'ProbeName', by.y = 'row.names')
  merged_filtered2 = na.omit(merged2)
  both_merged = merge(both_merged, merged_filtered2, by.x = 'ProbeName', by.y = 'ProbeName', all = FALSE)
}

both_merged = aggregate(both_merged, by = list(both_merged$ensembl_gene_id), FUN = mean)
both_merged$ProbeName = NULL
both_merged$ensembl_gene_id = NULL
both_merged$external_gene_name = NULL
both_merged = merge(both_merged, namen[,c('ensembl_gene_id', 'external_gene_name')],
                    by.x = 'Group.1', by.y = 'ensembl_gene_id', all = FALSE)
both_merged = unique(both_merged)
both_merged = both_merged %>% select(external_gene_name, everything())
both_merged = both_merged %>% select(Group.1, everything())
write.table(both_merged, '/home/michaelp/Downloads/aracne/aracne_in/124_aracne_all.txt', sep = '\t', row.names=FALSE)
