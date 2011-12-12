showWMTR <-
function (fdir = getwd()) 
{
	if(!require('rggobi', quietly = TRUE)) return("rggobi is not available!")
    ReadyTR <- read.table(as.character(paste(fdir, "/", "tmp_vrtheap.dat", 
        sep = "")))
    rggobi::ggobi(ReadyTR)
}

