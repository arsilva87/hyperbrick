#' Hyperspectral Image Registration (Translation)
#'
#' Hyperspectral image acquisition normaly causes spatial misalignment between
#' the spectral bands (layers) due to both equipment (such as band-to-band
#' recording delay) and external factors (e.g. sensor vibrations). In this case,
#' a geometric correction is necessary for remote sensing applications such
#' as combining/merging spectral bands. This function uses the HOG (Histogram
#' of Oriented Gradient) descriptor in order to find the optimal translations
#' (xy shift) on multiple 'slave' bands to be spatially align with a 'master'
#' (reference) band.
#'
#' @param Brick An object of class \code{RasterBrick} or \code{RasterStack}
#' (from package [raster]), containing multiple layers (spectral bands).
#'
#' @param ref_layer An integer indicating which layer (spectral band) should
#' be used as reference ('master') to register all the others from. Default is
#' 1, the first band.
#'
#' @param layers Either the character \code{"all"} (default), which indicates
#' that all bands should be registered using \code{ref_layer} as reference,
#' or a vector of integers specifying which bands to register.
#'
#' @param ncells An integer giving the number of cells to compute the oriented
#' gradients of the HOG descriptor. Default is 24. See [OpenImageR::HOG()].
#'
#' @param orient An integer giving the number of orientations to compute the
#' oriented gradients of the HOG descriptor. Default is 8. See [OpenImageR::HOG()].
#'
#' @details This should be used carefully, as rotation affects the spatial
#' dimensions. The affine parameters are estimated using a general
#' optimization algorithm.
#'
#' @return An object of the same classe as the input \code{slave}, with
#' the fixed extent.
#'
#' @seealso [OpenImageR::HOG()], [registerBand()], [registerBand3()],
#' [buildBrick()]
#'
#' @examples
#' path <- system.file('exdata', 'obory.dat', package = 'hyperbrick')
#' dpath <- system.file('exdata', 'obory_dark.dat', package = 'hyperbrick')
#' im <- buildBrick(path, hFOV = 36.8, vFOV = 36.8, height = 45,
#'                 ref_layer = 35, spectral_feature = 'radiance',
#'                 dark_path = dpath)
#' print(im)
#' plotRGB(im, r = 63, g = 34, b = 11, scale = 90)
#'
#' imreg <- registerBrick(im, ref_layer = 35, layers = c(63, 34, 11))
#' imreg
#' plotRGB(imreg, scale = 90)
#'
#' @importFrom OpenImageR HOG
#' @importFrom dfoptim nmk
#' @importFrom pbapply pblapply
#' @importFrom raster brick crs extent crop nlayers res shift origin
#'
#' @aliases registerBrick
#'
#' @export
registerBrick <- function(Brick, ref_layer = 1,
	layers = "all", ncells = 24, orient = 8)
{
   if (layers[1] == "all") {
      w <- 1:nlayers(Brick)
   } else {
      w <- as.integer(layers)
   }
   fRegBrick <- function(j, Brick, ref_layer) {
      registerBand(Brick[[ref_layer]], slave = Brick[[j]],
         ncells = ncells, orient = orient)
   }
   newbands <- pblapply(w, fRegBrick,
      Brick = Brick, ref_layer = ref_layer)
   rexts <- t(sapply(newbands, function(x) extent(x)[]))
   cropped_ext <- extent(c(max(rexts[,1]), min(rexts[,2]),
      max(rexts[,3]), min(rexts[,4])))
   newlist <- lapply(newbands, function(x) {
      origin(x) <- c(0,0)
      crop(x, cropped_ext)
   })
   regbrick <- brick(newlist)
   return(regbrick)
}
