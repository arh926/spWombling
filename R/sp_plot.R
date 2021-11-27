#' Spatial Plots
#'
#' A function that plots a kriged spatial field from coordinates.
#'
#' @param col.seq.length length of colors in color pallate
#' @param col.text names of color pallate
#' @param data_frame a matrix of coordinates and response (order \eqn{L} x  \eqn{3})
#' @param shape shape file for respone if available
#' @param zlim upper limit for response (optional)
#' @param xlab (Plotting parameter) label for x-axis
#' @param ylab (Plotting parameter) label for y-axis
#' @param lab (Plotting parameter) label for z-axis
#' @param add if plot needs to be added to an existing plot
#' @param categorical if response is cateforgical
#' @param points.plot if true plots points
#' @param contour.plot if true plots contours
#' @param legend if true plots legend
#' @param raster.surf if true returns a raster object of surface
#' @param sig vector of colors indicating significance for points (order \eqn{L} x \eqn{1})
#' @param grid if true plots grid
#' @keywords 
#' @import graphics
#' @import grDevices
#' @import stats
#' @import sp 
#' @import latex2exp
#' @importFrom RColorBrewer brewer.pal
#' @importFrom fields image.plot
#' @importFrom raster raster
#' @importFrom MBA mba.surf
#' @importFrom MASS mvrnorm
#' @export
#' @examples 
sp_plot <- function(col.seq.length=NULL, 
                    col.text = NULL, 
                    data_frame = NULL, 
                    shape=NULL,
                    zlim=NULL,
                    xlab=NULL,
                    ylab=NULL,
                    lab=NULL,
                    add=F,
                    categorical=F,
                    points.plot=F,
                    contour.plot=F,
                    legend=T,
                    raster.surf=F,
                    sig=NULL,
                    grid=F){
  col.br <- colorRampPalette(brewer.pal(col.seq.length,col.text))
  surf <- mba.surf(data_frame,no.X=300,no.Y=300,extend=T,sp=T)$xyz.est 
  
  if(legend){
    # mat <- matrix(c(1:2), nr=1,nc=2, byrow=T)
    # layout(mat,
    #        widths = rep(c(5,3),4),
    #        heights = c(3,3))
    if(categorical){
      col.br <- brewer.pal(col.seq.length,col.text)
      surf$z <- round(surf$z)
    } 
    if(!is.null(shape)){
      proj4string(surf) <- proj4string(shape)
      surf.tmp <- try(surf[!is.na(over(surf,shape))[,1],],silent = T)
      if(class(surf.tmp)=="try-error")  surf <- surf[!is.na(over(surf,shape)),] else surf <- surf.tmp
      surf <- as.image.SpatialGridDataFrame(surf)
    }
    if(is.null(zlim)) zlim <- range(surf[["z"]],na.rm=T)
    if(!is.null(shape)){
      if(categorical){
        image.plot(surf, xaxs="i", yaxs="i",col=rev(col.br), axes=T,#
                           zlim=zlim,
                           xlim=shape@bbox["x",],
                           ylim=shape@bbox["y",],
                           xlab=TeX("Longitude$\\degree$"),
                           ylab=TeX("Latitude$\\degree$"),
                           legend.lab = lab,
                           add=add)
      }else{
        if(!is.null(xlab) & !is.null(ylab)){
          image(surf, xaxs="i", yaxs="i",col=rev(col.br(100)), axes=T,
                zlim=zlim,
                xlim=shape@bbox["x",],
                ylim=shape@bbox["y",],
                xlab=xlab,
                ylab=ylab,
                legend.lab = lab,
                add=add)
        }else{
          image(surf, xaxs="i", yaxs="i",col=rev(col.br(100)), axes=T,
                zlim=zlim,
                xlim=shape@bbox["x",],
                ylim=shape@bbox["y",],
                xlab="",#latex2exp::TeX("Longitude$\\degree$"),
                ylab="",#latex2exp::TeX("Latitude$\\degree$"),
                legend.lab = lab,
                add=add)
        }
      }
    }else{
      image(surf, xaxs="i", yaxs="i",col=rev(col.br(100)), axes=T,
            zlim=zlim,
            xlab=xlab,#latex2exp::TeX("Longitude$\\degree$"),
            ylab=ylab,#latex2exp::TeX("Latitude$\\degree$"),
            legend.lab = lab,
            add=add)
    }
    
    if(!is.null(shape)){
      plot(shape, add=T, lwd=0.4)
    }
    if(contour.plot) contour(surf,add=T,lwd = 0.1)
    if(points.plot){
      points(data_frame[,(1:2)], pch=16,  cex=0.5)
      if(!is.null(sig)){
        points(data_frame[,(1:2)],
               col=sapply(sig, function(x){
                 if(x==0) return("white")
                 if(x==-1) return("cyan")
                 else return("green")
               }), pch=16,cex=1)
      }
    } 
    if(grid){
      abline(h=data_frame[,2],lwd=0.01)
      abline(v=data_frame[,1],lwd=0.01)
    }
    if(min(data_frame[,3])==max(data_frame[,3])) plot.new()
    else{
      legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
      plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
      text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.4f",round(seq(min(data_frame[,3]),max(data_frame[,3]),l=6),2)))
      rasterImage(legend_image, 0, 0, 1,1)
    }
    if(raster.surf) return(raster(surf))
  }else{
    if(categorical){
      col.br <- brewer.pal(col.seq.length,col.text)
      surf$z <- round(surf$z)
    } 
    #browser()
    if(!is.null(shape)){
      proj4string(surf) <- proj4string(shape)
      surf.tmp <- try(surf[!is.na(over(surf,shape))[,1],],silent = T)
      if(class(surf.tmp)=="try-error")  surf <- surf[!is.na(over(surf,shape)),] else surf <- surf.tmp
      surf <- as.image.SpatialGridDataFrame(surf)
    }
    if(is.null(zlim)) zlim <- range(surf[["z"]],na.rm=T)
    if(!is.null(shape)){
      if(categorical){
        image.plot(surf, xaxs="i", yaxs="i",col=rev(col.br), axes=T,#
                           zlim=zlim,
                           xlim=shape@bbox["x",],
                           ylim=shape@bbox["y",],
                           xlab=TeX("Longitude$\\degree$"),
                           ylab=TeX("Latitude$\\degree$"),
                           legend.lab = lab,
                           add=add)
      }else{
        if(!is.null(xlab) & !is.null(ylab)){
          image(surf, xaxs="i", yaxs="i",col=rev(col.br(100)), axes=T,
                zlim=zlim,
                xlim=shape@bbox["x",],
                ylim=shape@bbox["y",],
                xlab=xlab,
                ylab=ylab,
                legend.lab = lab,
                add=add)
        }else{
          image(surf, xaxs="i", yaxs="i",col=rev(col.br(100)), axes=T,
                zlim=zlim,
                xlim=shape@bbox["x",],
                ylim=shape@bbox["y",],
                xlab="",#latex2exp::TeX("Longitude$\\degree$"),
                ylab="",#latex2exp::TeX("Latitude$\\degree$"),
                legend.lab = lab,
                add=add)
        }
      }
    }else{
      image(surf, xaxs="i", yaxs="i",col=rev(col.br(100)), axes=T,
            zlim=zlim,
            xlab=xlab,#latex2exp::TeX("Longitude$\\degree$"),
            ylab=ylab,#latex2exp::TeX("Latitude$\\degree$"),
            legend.lab = lab,
            add=add)
    }
    
    if(!is.null(shape)){
      plot(shape, add=T, lwd=0.4)
    }
    if(points.plot){
      if(!is.null(sig)){
        points(data_frame[,(1:2)],
               col=sapply(sig, function(x){
                 if(x==0) return("white")
                 if(x==-1) return("cyan")
                 else return("green")
               }), pch=16,cex=1)
      }else{
        points(data_frame[,(1:2)], pch=16,  cex=0.5)
      }
    }
    if(grid){
      abline(h=data_frame[,2],lwd=0.01)
      abline(v=data_frame[,1],lwd=0.01)
    }
    if(contour.plot) contour(surf,add=T,lwd = 0.1)
    if(raster.surf) return(raster(surf))
  }
}
