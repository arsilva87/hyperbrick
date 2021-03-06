% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/affineBrick.R
\name{affineBrick}
\alias{affineBrick}
\title{Affine Transformation (rotation and shift) on Images}
\usage{
affineBrick(Brick, angle = 0, xy_shift = c(0, 0))
}
\arguments{
\item{Brick}{An object of class \code{RasterBrick}, \code{RasterStack} or
\code{RasterLayer} (from package \link{raster}).}

\item{angle}{A numeric value of the angle (in degrees, 0-360) to rotate
\code{Brick}. A negative value will change the direction of rotation to
clockwise.}

\item{xy_shift}{A numeric vector of length two with the x and y shift (the
translation parameters).}
}
\value{
An object of the same class as the input \code{Brick}.
}
\description{
Affine transformations are of type \eqn{f(x) = Ax + b}, where \eqn{x} is
the spatial coordinates (2D in this case), \eqn{A} is a rotation matrix
(it can also include scale/shear parameters, but only rotation is
considered here), and \eqn{b} is the translation (xy shift) parameters.
}
\note{
Affine transformation affects the image dimension.
}
\examples{
p <- system.file('exdata', 'soybean.tif', package = 'hyperbrick')
im <- brick(p)
print(im)

# view band-3
plot(im[[3]], col = gray.colors(20), asp = 0)

# rotate band-3 at 3.5 degrees counter-clockwise
b3_rot <- affineBrick(im[[3]], angle = 3.5)
plot(b3_rot, add = TRUE, legend = FALSE,
    col = adjustcolor(terrain.colors(20), 0.5))

}
\seealso{
\code{\link[=affineCoords]{affineCoords()}}, \code{\link[=registerBrick]{registerBrick()}}
}
