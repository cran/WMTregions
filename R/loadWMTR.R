loadWMTR <-
function (fdir = getwd()) 
{
    wmtreg <- read.table(as.character(paste(fdir, "/", "TRegion.dat", 
        sep = "")))
    as.matrix(wmtreg)
}

