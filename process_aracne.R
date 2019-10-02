library(dplyr)

# This script was provided by Andreas Schuettler

# create network from aracne output -----------------------------------------------------------------------------------------
maxreg<-max(count.fields(file = "./125_aracne_all_out.txt",sep="\t",quote = "",comment.char = ">"))
adj_table<-read.table(file = "./125_aracne_all_out.txt",sep="\t",quote = "",dec = ",",fill = T,as.is=T,comment.char = ">",col.names = c("Regulator",rep(c("Target","MI"),((maxreg-1)/2))))



x<-adj_table[4,]

adj_list<-apply(adj_table,MARGIN = 1,FUN = function(x){
  first0<-which(x=="")[1]
  n<-(first0-2)/2
  if(!is.na(n)){
  y<-data.frame(Regulator=rep(as.character(x[1]),n),Target=as.character(x[seq(from = 2,to = (first0-2),by = 2)]),MI=as.numeric(x[seq(from = 3,to = (first0-1),by = 2)]))
  y}
  })

adj_dataframe<-do.call("rbind",adj_list)
adj_dataframe = as.data.frame(sapply(adj_dataframe, function(x) gsub("\"", "", x)))
adj_dataframe = transform(adj_dataframe, MI = as.numeric(levels(MI)[MI]))

rm(adj_table)
rm(adj_list)
gc()

save(adj_dataframe,file = "./fullnetwork_125_all.Rd")

net_embryo<-igraph::graph_from_data_frame(d = adj_dataframe, directed = F)
net_embryo<-igraph::simplify(graph = net_embryo, remove.multiple = T, edge.attr.comb = "median")

save(net_embryo,file = "./fullnetwork_125_all_igraph.Rd")

rm(data_proc_naproxen)
rm(adj_dataframe)

