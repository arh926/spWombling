#' Spatial Plot Function
#' 
#' @param data_frame data frame consisting of coordinates and data
#' @param sp logical parameter indicating whether to make a spatial plot
#' @param shape if sp = TRUE shape file should be provided (should be an sf object) 
#' @param palette (optional) color palette
#' @keywords sp_plot
#' @import ggplot2 cowplot MBA metR sf
sp_ggplot <- function(data_frame = NULL, 
                      sp = FALSE,
                      shape = NULL,
                      palette = "Spectral"){
  
  if(sp & is.null(shape)) stop("Please provide a shape file!")
  if(ncol(data_frame) == 4){
    col.pt = sapply(data_frame[,4], function(x){
      if(x == 0) return(NA)
      else if(x == 1) return("green")
      else return("cyan")
    })
    color.pt = ifelse(is.na(col.pt), NA, "black")
  }else{
    col.pt = color.pt = rep("black", nrow(data_frame))
  }
  
  
  cnames = colnames(data_frame)[1:3]
  coords = data.frame(data_frame[, c(1, 2)])
  surf = mba.surf(data_frame[,1:3],
                  no.X = 300,
                  no.Y = 300,
                  h = 5,
                  m = 2,
                  extend = TRUE, sp = sp)$xyz.est
  
  if(sp){
    proj4string(surf) = proj4string(shape)
    surf.tmp = try(surf[!is.na(over(surf, shape))[, 1],], silent = TRUE)
    surf = as.image.SpatialGridDataFrame(surf.tmp)
  }
  
  gg.grid = expand.grid(surf$x, surf$y)
  
  df.gg = cbind(gg.grid, as.vector(surf$z))
  colnames(df.gg) = cnames
  
  if(sp){
    gplot = ggplot() + theme_cowplot(12) +
      geom_raster(data = df.gg, mapping = aes_string(x = cnames[1], y = cnames[2], fill = cnames[3])) +
      labs(x = cnames[1], y = cnames[2], fill = "") +
      scale_fill_distiller(palette = palette,
                           label = function(x) sprintf("%.2f", x), na.value = "white", direction = 1) +
      geom_contour2(data = df.gg,
                    mapping = aes_string(x = cnames[1], y = cnames[2], z = cnames[3], label = "after_stat(level)"), 
                    linewidth = 0.1,
                    label_size = 2.5) +
      geom_sf(data = shape, fill = "transparent") +
      geom_point(data = coords,
                 mapping = aes_string(x = cnames[1], y = cnames[2]),
                 size = 1.2, color = color.pt,
                 fill = col.pt,  stroke = 0.5, pch = 21, na.rm = TRUE) +
      theme(axis.line = element_line(linewidth = 0.2),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 10),
            legend.key.height = unit(1.5, "cm"),
            legend.key.width = unit(0.6, "cm"),
            legend.text = element_text(size = 20),
            plot.margin = unit(rep(0.15, 4), "cm"))
  }else{
    gplot = ggplot() + theme_cowplot(12) +
      geom_raster(data = df.gg, mapping = aes_string(x = cnames[1], y = cnames[2], fill = cnames[3])) +
      labs(x = cnames[1], y = cnames[2], fill = "") +
      scale_fill_distiller(palette = palette,
                           label = function(x) sprintf("%.2f", x), na.value = "white", direction = 1) +
      geom_contour2(data = df.gg,
                    mapping = aes_string(x = cnames[1], y = cnames[2], z = cnames[3], label = "after_stat(level)"), 
                    linewidth = 0.1,
                    label_size = 2.5) +
      geom_point(data = coords,
                 mapping = aes_string(x = cnames[1], y = cnames[2]),
                 size = 1.2, color = color.pt,
                 fill = col.pt,  stroke = 0.5, pch = 21, na.rm = TRUE) +
      theme(axis.line = element_line(linewidth = 0.2),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 10),
            legend.key.height = unit(1.5, "cm"),
            legend.key.width = unit(0.6, "cm"),
            legend.text = element_text(size = 20),
            plot.margin = unit(rep(0.15, 4), "cm"))
  }
 
  return(gplot)
}