loadWMTR <-
function (fname = "TRegion.dat", fdir = getwd()) 
{
    wmtreg <- read.table(as.character(paste(fdir, "/", fname, 
        sep = "")))
    as.matrix(wmtreg)
}

