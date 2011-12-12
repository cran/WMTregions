visualWMTR <-
function (fdir = getwd()) 
{
    datacloudvis <- read.table(as.character(paste(fdir, "/", "tmp_vis_dt.dat", 
        sep = "")))
        
    rgl.open()
	#bg3d('white')
	
 	centroid = as.real(datacloudvis[1,])
	datacloudvis = datacloudvis[2:nrow(datacloudvis),]
	
	maxscale = max(abs(datacloudvis+centroid))
	datacloudvis = (7.0 * datacloudvis) / maxscale
	centroid = (7.0 * centroid) / maxscale
	
   radiusdat = min(0.1, 0.4 / (nrow(datacloudvis)^0.33))
   rgl.spheres(datacloudvis[[1]], datacloudvis[[2]], datacloudvis[[3]], radiusdat, color = 'red')

 rgl.material(color = 'green', shininess = 50.0)

    ridgesvis <- read.table(as.character(paste(fdir, "/", "tmp_vis_rg.dat", 
        sep = "")))
	ridgesvis = (7.0 * ridgesvis) / maxscale
    rgl.lines(ridgesvis[[1]], ridgesvis[[2]], ridgesvis[[3]], color = 'darkseagreen', alpha = 0.25)
	
    facetsvis <- read.table(as.character(paste(fdir, "/", "tmp_vis_fc.dat", 
        sep = "")))
	facetsvis = (7.0 * facetsvis) / maxscale
    rgl.triangles(facetsvis[[1]], facetsvis[[2]], facetsvis[[3]], color = 'blue')

 #rgl.material(color = 'yellow', shininess = 50.0)

	ax = c(0,10,0,0,0,0)
ay = c(0,0,0,10,0,0)
az = c(0,0,0,0,0,10)

tax = c(10,9.6,10, 9.6,0,0.2,0,-0.2,0,0.2,0,-0.2)
tay = c(0,-0.2,0,0.2,10,9.6,10, 9.6,0,-0.2,0,0.2)
taz = c(0,0.2,0,-0.2,0,-0.2,0,0.2,10,9.6,10,9.6)

   rgl.lines(ax - centroid[1],ay - centroid[2],az - centroid[3], color = 'ghostwhite')
   rgl.lines(tax - centroid[1],tay - centroid[2],taz - centroid[3], color = 'ghostwhite')
    #rgl.lines(ax,ay,az, color = c('red','red','green','green','yellow','yellow'))
    #rgl.lines(tax,tay,taz, color = c(0,200,0))

axname = c(10.5,0.5,0.5)
ayname = c(0,10,0)
azname = c(0,0,10)
   
axesnames =  c("x","y","z")
 text3d(axname, ayname, azname, texts=axesnames ,adj = 0, color="ghostwhite", family=1, font=1, cex=1)

    
}

