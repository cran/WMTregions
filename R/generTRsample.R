generTRsample <-
function (fname = "Cloud.dat", fdir = getwd(), dim = 3, num = 20, 
    alpha = 0.05, trtype = "zonoid") 
{
    ffullname = as.character(paste(fdir, "/", fname, sep = ""))
    write(trtype, ffullname)
    write(alpha, ffullname, append = TRUE)
    write(dim, ffullname, append = TRUE)
    write(num, ffullname, append = TRUE)
    write(array(runif(dim * num, 0, 1000)/100, dim = c(num, dim)), 
        ffullname, ncolumns = dim, append = TRUE)
}

