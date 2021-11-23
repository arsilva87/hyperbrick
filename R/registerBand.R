#' Single Band-to-Band Registration (Translation)
#'
#' Hyperspectral image acquisition normaly causes spatial misalignment between
#' the spectral bands (layers) due to both equipment (such as band-to-band
#' recording delay) and external factors (e.g. sensor vibrations). In this case,
#' a geometric correction is necessary for remote sensing applications such
#' as combining/merging spectral bands. This function uses the HOG (Histogram
#' of Oriented Gradient) descriptor in order to find the optimal translation
#' (xy shift) on a 'slave' band to be spatially align with a 'master'
#' (reference) band.
#'
#' @param slave An object of class \code{RasterLayer} (from package
#' [raster]).
#'
#' @param master An object of class \code{RasterLayer} (from package
#' [raster]).
#'
#' @param ncells An integer giving the number of cells to compute the oriented
#' gradients of the HOG descriptor. Default is 24. See [OpenImageR::HOG()].
#'
#' @param orient An integer giving the number of orientations to compute the
#' oriented gradients of the HOG descriptor. Default is 8. See [OpenImageR::HOG()].
#'
#' @details The affine parameters are estimated using a general
#' optimization algorithm. This function only estimates translation parameters.
#' To register bands also with rotation fixes, please check [registerBand3()].
#' But this should be used carefully, as rotation affects the spatial dimensions.
#'
#' @return An object of the same classe as the input \code{slave}, with
#' the fixed extent.
#'
#' @seealso [OpenImageR::HOG()], [registerBrick()], [registerBand3()]
#'
#' @examples
#' #not yet
#'
#' @importFrom OpenImageR HOG
#' @importFrom dfoptim nmk
#' @importFrom raster extent crop res shift
#'
#' @aliases registerBand
#'
#' @export
registerBand <- function(slave, master, ncells = 24, orient = 8)
{
   bx <- res(master)[1] * ncol(master)/10
   by <- res(master)[2] * nrow(master)/10
   ex <- extent(master)[]
   pol <- extent(c(ex[1]+bx, ex[2]-bx, ex[3]+by, ex[4]-by))
   r1c <- crop(master, pol)
   hog1 <- HOG(as.matrix(r1c), cells = ncells, orientations = orient)
   fun_sxy <- function(par) {
      sx <- par[1]; sy <- par[2]
      pol2 <- extent(pol[] - c(sx, sx, sy, sy))
      r2c <- crop(slave, pol2)
      hog2 <- HOG(as.matrix(r2c), cells = ncells, orientations = orient)
      dist(rbind(hog1, hog2))[1]
   }
   p <- nmk(c(0, 0), fun_sxy)
   sxy <- p$par
   r2_shifted <- shift(slave, dx = sxy[1], dy = sxy[2])
   attr(r2_shifted, "shift") <- sxy
   return(r2_shifted)
}
