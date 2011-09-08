WMTR <-
function (fname = "Cloud.dat", fdir = getwd(), bound = 0) 
{
    .C("ComputeWMTR", as.character(fname), as.character(paste(fdir, 
        "/", sep = "")), as.integer(bound), PACKAGE = "WMTregions")
    "The trimmed region has been successfully calculated!"
}

