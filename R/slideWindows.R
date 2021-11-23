#' Create Sliding Windows Over a Spatial Object
#'
#' Create a prefixed number of non-overlapping sliding windows that cover
#' the entire extent of a spatial object.
#'
#' @param x A spatial object where an extent can be extracted from.
#' Classes: \code{Extent}, \code{RasterLayer}, \code{RasterBrick},
#' \code{RasterStack} (from package [raster]).
#'
#' @param n An integer vector of length two determining the number of divisions
#' in x and y directions. Default is \code{n = c(8, 8)}, which creates 64
#' windows.
#'
#' @return An object of class \code{slideWindows}, which consists of a list
#' of each sliding-windows extent.
#'
#' @seealso [slideBrick()], [raster::extent]
#'
#' @examples
#' #not yet
#'
#' @importFrom raster extent
#'
#' @aliases slideWindows
#'
#' @export
slideWindows <- function(x, n = c(8, 8))
{
   test_class <- c("Extent", "RasterLayer",
      "RasterBrick", "RasterStack")
   stopifnot(any(test_class %in% class(x)))
   ext <- extent(x)[]
   stopifnot(any(n > 1))
   n <- as.integer(n)
   x_win <- diff(ext[1:2])/n[1]
   y_win <- diff(ext[3:4])/n[2]
   win1 <- c(ext[1], ext[1] + x_win, ext[3], ext[3] + y_win)
   gr <- expand.grid(x = 0:(n[1]-1), y = 0:(n[2]-1))
   wins <- t(apply(gr, 1, function(z) {
      c(win1[1:2] + x_win*z[1], win1[3:4] + y_win*z[2])
   }))
   win_exts <- apply(wins, 1, extent)
   class(win_exts) <- "slideWindows"
   attr(win_exts, "win_size") <- c(x = x_win, y = y_win)
   return(win_exts)
}
