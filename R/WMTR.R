WMTR <-
function (fname = "Cloud.dat", fdir = getwd()) 
{
    .C("ComputeWMTR", as.character(fname), as.character(paste(fdir, 
        "/", sep = "")), PACKAGE = "WMTregions")
    "The trimmed region has been successfully calculated!"
}

