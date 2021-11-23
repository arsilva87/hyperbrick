#' Spectral Signature of a Hyperspectral Image
#'
#' Visualize statistics calculated through the bands of a hyperspectral image.
#'
#' @param x A numeric matrix or vector containing the values to be plotted at
#' each spectral band (wavelength). Generally, an object obtained with
#' [slideBrick()].
#'
#' @param ... Further graphical parameters. See [par()].
#'
#' @examples
#' #not yet
#'
#' @importFrom graphics matplot
#'
#' @aliases viewSpectra
#'
#' @export
viewSpectra <- function(x, ...) {
   if (length(dim(x)) > 1) bands <- rownames(x) else bands <- names(x)
   wavelength <- as.numeric(sub("b", "", bands))
   matplot(wavelength, x, type = "l", ...)
}
