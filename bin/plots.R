args = commandArgs(trailingOnly=TRUE)
NCM_folder = args[1]
label.file = args[2]
ngscheckmateoutput.file = "NCM_output/output_all.txt"

source(paste0(NCM_folder,"/graph/ngscheckmate2xgmml.R") )
create.xgmml.from.ngscheckmateout(label.file,ngscheckmateoutput.file,"NCM_graph.xgmml")