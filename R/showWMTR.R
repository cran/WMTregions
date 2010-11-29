showWMTR <-
function (fdir = getwd()) 
{
    ReadyTR <- read.table(as.character(paste(fdir, "/", "tmp_vrtheap.dat", 
        sep = "")))
    ggobi(ReadyTR)
}

