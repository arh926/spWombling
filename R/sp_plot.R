#' Spatial Plot Function
#'
#' @param data_frame data frame consisting of coordinates and data
#' @param sp logical parameter indicating whether to make a spatial plot
#' @param shape if sp = TRUE shape file should be provided (should be an sf object)
#' @param legend.key.height height of legend (defaults to .7)
#' @param legend.key.width width of legend (defaults to .4)
#' @param text.size size of legend text (defaults to 10)
#' @param point.size size of points to be plotted (defaults to 0.7)
#' @param extend logical parameter indicating whether to extend the interpolation (defaults to TRUE)
#' @param title title of the plot (defaults to NULL)
#' @param palette (optional) color palette
#' @keywords sp_plot
#' @import ggplot2 ggspatial terra MBA metR sf
#' @export
sp_ggplot <- function(data_frame = NULL,
                      sp = FALSE,
                      shape = NULL,
                      legend.key.height = 0.7,#1.5,
                      legend.key.width = 0.4,#0.6,
                      text.size = 10, #20,
                      point.size = 0.7, #1.2,
                      clr.pt = "black",
                      palette = "Spectral",
                      extend = TRUE,
                      title = NULL, bound.box = NULL){

  if(sp & is.null(shape)) stop("Please provide a shape file!")
  if(ncol(data_frame) == 4){
    col.pt = sapply(data_frame[,4], function(x){
      if(x == 0) return(NA)
      else if(x == 1) return("green")
      else return("cyan")
    })
    color.pt = ifelse(is.na(col.pt), NA, "black")
  }else{
    col.pt = rep("black", nrow(data_frame))
    color.pt = rep(clr.pt, nrow(data_frame))
  }


  cnames = colnames(data_frame)[1:3]
  coords = data.frame(data_frame[, c(1, 2)])

  surf = MBA::mba.surf(data_frame[,1:3],
                  no.X = 200,
                  no.Y = 200,
                  h = 5,
                  m = 2,
                  extend = extend, sp = sp)$xyz.est

  if(sp){
    box_proj = st_coordinates(st_transform(st_as_sf(data.frame(t(bound.box)), coords = c("x","y"), crs = 4326), crs = 5070))

    coords_proj = st_coordinates(st_transform(st_as_sf(coords, coords = cnames[1:2], crs = 4326), crs = 5070))
    colnames(coords_proj) = cnames[1:2]

    shape.sp = as(shape, "Spatial")
    proj4string(surf) = proj4string(shape.sp)
    surf.tmp = try(surf[!is.na(over(surf, shape.sp))[, 1],], silent = TRUE)
    surf = as.image.SpatialGridDataFrame(surf.tmp)

    surf.rast = project(rast(surf.tmp), "EPSG:5070")
    df_surf.rast = as.data.frame(surf.rast, xy = TRUE)
    colnames(df_surf.rast) = cnames

    gplot = ggplot() + theme_bw() +
      layer_spatial(data = surf.rast, mapping = aes(fill = after_stat(band1))) +
      labs(x = cnames[1], y = cnames[2], fill = "") +
      scale_fill_distiller(palette = palette,
                           label = function(x) sprintf("%.2f", x), na.value = NA, direction = 1) +
      geom_contour2(data = df_surf.rast,
                    mapping = aes_string(x = cnames[1], y = cnames[2], z = cnames[3], label = "round(after_stat(level), 2)"),
                    linewidth = 0.1,
                    label_size = 2.5) +
      geom_sf(data = shape, fill = "transparent") +
      geom_point(data = coords_proj,
                 mapping = aes_string(x = cnames[1], y = cnames[2]),
                 size = point.size, color = color.pt,
                 fill = col.pt,  stroke = 0.5, pch = 21, na.rm = TRUE) +
      theme(axis.line = element_line(linewidth = 0.2),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 10),
            legend.key.height = unit(1.5, "cm"),
            legend.key.width = unit(0.6, "cm"),
            legend.text = element_text(size = 20),
            plot.margin = unit(rep(0.15, 4), "cm"),
            panel.background = element_rect(colour = "black", size = 1))
  }else{
    gg.grid = expand.grid(surf$x, surf$y) # st_coordinates(surf.sf)
    df.gg = cbind(gg.grid, as.vector(surf$z))
    colnames(df.gg) = cnames

    gplot = ggplot() +
      geom_raster(data = df.gg, mapping = aes_string(x = cnames[1], y = cnames[2], fill = cnames[3])) +
      labs(x = cnames[1], y = cnames[2], fill = "") +
      scale_fill_distiller(palette = palette,
                           label = function(x) sprintf("%.2f", x), na.value = NA, direction = 1) +
      geom_contour2(data = df.gg,
                    mapping = aes_string(x = cnames[1], y = cnames[2], z = cnames[3], label = "round(after_stat(level), 2)"),
                    linewidth = 0.1,
                    label_size = 2.5) +
      geom_point(data = coords,
                 mapping = aes_string(x = cnames[1], y = cnames[2]),
                 size = point.size, color = color.pt,
                 fill = col.pt,  stroke = 0.5, pch = 21, na.rm = TRUE) +
      theme_bw() +
      theme(axis.line = element_line(linewidth = 0.2),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 10),
            legend.key.height = unit(legend.key.height, "cm"),
            legend.key.width = unit(legend.key.width, "cm"),
            legend.text = element_text(size = text.size),
            plot.margin = unit(rep(0.15, 4), "cm"),
            panel.background = element_rect(colour = "black", size = 1))
  }
  if(!is.null(title)) gplot = gplot + ggtitle(title)
  if(!is.null(bound.box)){
    gplot = gplot + xlim(box_proj[,1]) + ylim(box_proj[,2])
  }
  return(gplot)
}
