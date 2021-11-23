#' Single Band-to-Band Registration (Rotation and Translation)
#'
#' Hyperspectral image acquisition normaly causes spatial misalignment between
#' the spectral bands (layers) due to both equipment (such as band-to-band
#' recording delay) and external factors (e.g. sensor vibrations). In this case,
#' a geometric correction is necessary for remote sensing applications such
#' as combining/merging spectral bands. This function uses the HOG (Histogram
#' of Oriented Gradient) descriptor in order to find the optimal rotation
#' angle and translation (xy shift) on a 'slave' band to be spatially align
#' with a 'master' (reference) band.
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
#' @details This should be used carefully, as rotation affects the spatial
#' dimensions. The affine parameters are estimated using a general
#' optimization algorithm.
#'
#' @return An object of the same classe as the input \code{slave}, with
#' the fixed extent.
#'
#' @seealso [OpenImageR::HOG()], [registerBrick()], [registerBand()]
#'
#' @examples
#' #not yet
#'
#' @importFrom OpenImageR HOG
#' @importFrom dfoptim nmk
#' @importFrom raster extent crop res shift
#'
#' @aliases registerBand3
#'
#' @export
registerBand3 <- function(slave, master, ncells = 24, orient = 8)
{
   bx <- res(master)[1] * ncol(master)/10
   by <- res(master)[2] * nrow(master)/10
   ex <- extent(master)[]
   pol <- extent(c(ex[1]+bx, ex[2]-bx, ex[3]+by, ex[4]-by))
   r1c <- crop(master, pol)
   hog1 <- HOG(as.matrix(r1c), cells = ncells,
      orientations = orient)
   fun_affine <- function(par) {
      sx <- par[1]; sy <- par[2]; a <- par[3]
      affpol <- affineCoords(pol, angle = a, c(sx, sy))
      r2c <- crop(slave, extent(affpol))
      hog2 <- HOG(as.matrix(r2c), cells = ncells,
         orientations = orient)
      dist(rbind(hog1, hog2))[1]
   }
   p <- nmk(c(0, 0, 0), fun_affine)
   sx <- p$par[1]; sy <- p$par[2]; a <- p$par[3]
   fixed <- affineBrick(slave, a, c(-sx, -sy))
   attr(fixed, "affine_pars") <- c(angle = a, sx = sx, sy = sy)
   return(fixed)
}
