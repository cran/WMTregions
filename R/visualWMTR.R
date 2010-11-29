visualWMTR <-
function (fdir = getwd()) 
{
    datacloudvis <- read.table(as.character(paste(fdir, "/", "tmp_vis_dt.dat", 
        sep = "")))
        
    rgl.open()
    rgl.spheres(datacloudvis[[1]], datacloudvis[[2]], datacloudvis[[3]], 0.1, color = 'red')

 rgl.material(color = 'blue', shininess = 50.0)

    facetsvis <- read.table(as.character(paste(fdir, "/", "tmp_vis_fc.dat", 
        sep = "")))
    rgl.triangles(facetsvis[[1]], facetsvis[[2]], facetsvis[[3]], color = 'blue')


    ridgesvis <- read.table(as.character(paste(fdir, "/", "tmp_vis_rg.dat", 
        sep = "")))
    rgl.lines(ridgesvis[[1]], ridgesvis[[2]], ridgesvis[[3]], color = c(100,100,100))
    
}

