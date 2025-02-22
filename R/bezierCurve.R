#' Bezier Splines
#' 
#' A function for computing Bezier splines from a set of annotated points
#' 
#' @param x x-coordinate for points
#' @param y y-coordinate for points
#' @param n number of points to interpolate
#' @keywords spWombling
#' @import stats

##############################################
# Bezier Splines for Annotated set of points #
##############################################
bezierCurve <- function(x, y, n=10)
{
  outx <- NULL
  outy <- NULL
  
  i <- 1
  for (t in seq(0, 1, length.out=n))
  {
    b <- bez(x, y, t)
    outx[i] <- b$x
    outy[i] <- b$y
    
    i <- i+1
  }
  
  return (list(x=outx, y=outy))
}

bez <- function(x, y, t)
{
  outx <- 0
  outy <- 0
  n <- length(x)-1
  for (i in 0:n)
  {
    outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]
    outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]
  }
  
  return (list(x=outx, y=outy))
}

# pts.hull <- locator(10)
# pts.hull$x <- c(pts.hull$x, pts.hull$x[1])
# pts.hull$y <- c(pts.hull$y, pts.hull$y[1])
# points(pts.hull)
# Example usage
# points(bezierCurve(pts.hull$x,pts.hull$y,40), type="l", col="red")