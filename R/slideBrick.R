#' Sliding-Windows Statistics Over a Hyperspectral Image
#'
#' Calculate focal statistics on a hyperspectral image using non-overlapping
#' sliding windows.
#'
#' @param Brick An object of class \code{RasterBrick} or \code{RasterStack}
#' (from package [raster]), containing multiple layers (spectral bands).
#'
#' @param slide_windows An object of class \code{slideWindows}, which consists
#' of a list of spatial extents giving the location of the non-overlapping
#' sliding windows that covers the entire image.
#'
#' @param fun A vectorized function indicating the statistics to be calculated
#' within each sliding-windows, e.g. [median()] (default), [mean()], [sd()].
#'
#' @return A vector or matrix containing the estimates of each sliding windows.
#'
#' @seealso [slideWindows()]
#'
#' @examples
#' #not yet
#'
#' @importFrom stats dist median quantile
#' @importFrom raster crop cellStats
#'
#' @aliases slideBrick
#'
#' @export
slideBrick <- function(Brick, slide_windows, fun = median)
{
   stopifnot(inherits(slide_windows, "slideWindows"))
   win_exts <- slide_windows
   brickset <- lapply(win_exts, crop, x = Brick)
   stats <- sapply(brickset, cellStats, stat = fun)
   return(stats)
}
