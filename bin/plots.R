args = commandArgs(trailingOnly=TRUE)
label.file = args[1]

ngscheckmateoutput.file = "NCM_output/output_all.txt"

source("ngscheckmate2xgmml.R")
create.xgmml.from.ngscheckmateout(label.file,ngscheckmateoutput.file,"NCM_graph_wrongmatch.xgmml")
create.xgmml.from.ngscheckmateout(label.file,ngscheckmateoutput.file,"NCM_graph.xgmml",filter=FALSE)