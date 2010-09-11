loadWMTR <-
function(dim, fdir = getwd())

{

wmtreg <- matrix(scan(as.character(paste(fdir,"/","TRegion.dat",sep = ""))), ncol = dim+1)

wmtreg

}

