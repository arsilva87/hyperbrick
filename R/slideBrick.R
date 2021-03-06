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
#' within each sliding-windows, e.g. \code{median} (default), \code{mean},
#' \code{sd}.
#'
#' @return A vector or matrix containing the estimates of each sliding windows.
#'
#' @seealso [slideWindows()]
#'
#' @examples
#' p <- system.file('exdata', 'obory.dat', package = 'hyperbrick')
#' im <- buildBrick(p, ref_layer = 35,
#'                 spectral_feature = "radiance",
#'                 hFOV = 36.8, vFOV = 36.8, height = 45)
#' print(im)
#' plotRGB(im, r = 63, b = 34, g = 11, scale = 90)
#'
#' ext <- extent(c(512700.2, 512715, 5769462, 5769477))
#' sw <- slideWindows(ext, n = c(7, 7))
#' lapply(sw, lines, col = "white") -> null_obj
#'
#' sb <- slideBrick(im, sw, fun = mean)
#' head(sb)
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
