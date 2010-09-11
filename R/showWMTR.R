showWMTR <-
function(fdir = getwd())



{



# library(rggobi)



ReadyTR <- read.table(as.character(paste(fdir,"/","TR_vertices.dat",sep = "")))



ggobi(ReadyTR)



}

