args = commandArgs(trailingOnly=TRUE)
NCM_folder = Sys.getenv("NCM_HOME") #args[1]
label.file = args[1]

print(c(NCM_folder,label.file) )

ngscheckmateoutput.file = "NCM_output/output_all.txt"

source(paste0(NCM_folder,"/graph/ngscheckmate2xgmml.R") )
create.xgmml.from.ngscheckmateout(label.file,ngscheckmateoutput.file,"NCM_graph_wrongmatch.xgmml")
create.xgmml.from.ngscheckmateout(label.file,ngscheckmateoutput.file,"NCM_graph.xgmml",filter=FALSE)