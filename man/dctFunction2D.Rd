% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univariateExpansions.R
\name{dctFunction2D}
\alias{dctFunction2D}
\alias{dctFunction3D}
\title{Calculate linear combinations of orthonormal cosine basis functions on
two- or three-dimensional domains}
\usage{
dctFunction2D(scores, argvals, parallel = FALSE)

dctFunction3D(scores, argvals, parallel = FALSE)
}
\arguments{
\item{scores}{A sparse matrix of dimension \code{N x L}, representing
the \code{L} scores (coefficients), where \code{N} is the number of
observations.}

\item{argvals}{A list containing two or three numeric vectors,
corresponding to the domain grid (x and y values for two-dimensional
domains; x,y and z values fro three-dimensional domains.)}

\item{parallel}{Logical. If \code{TRUE}, the coefficients for the basis
functions are calculated in parallel. The implementation is based on
the \code{\link[foreach]{foreach}} function and requires a parallel
backend that must be registered before; see
\code{\link[foreach]{foreach}} for details. Defaults to \code{FALSE}.}
}
\value{
An object of class \code{funData} with \code{N} observations on
  the two- or threedimensional domain specified by \code{argvals},
  corresponding to the linear combination of orthonormal cosine basis
  functions.
}
\description{
Given scores (coefficients), these functions calculate a linear
combination of two- or three-dimensional cosine tensor basis functions
on two- or three-dimensional domains using the C-library \code{fftw3}
(see \url{http://www.fftw.org/}).
}
\section{Warning}{
 If the C-library \code{fftw3} is not available when
  the package \code{MFPCA} is installed, the functions are disabled an
  will throw an error. For full functionality install the C-library
  \code{fftw3} from \url{http://www.fftw.org/} and reinstall
  \code{MFPCA}. This function has not been tested with
  ATLAS/MKL/OpenBLAS.
}

\examples{
# set seed
set.seed(12345)

# generate sparse 10 x 15 score matrix (i.e. 10 observations) with 30 entries
# smoothness assumption: higher order basis functions (high column index) have lower probability
scores <- Matrix::sparseMatrix(i = sample(1:10, 30, replace = TRUE), # sample row indices
     j = sample(1:15, 30, replace = TRUE, prob = 1/(1:15)), # sample column indices
     x = rnorm(30)) # sample values
scores

\dontrun{
# calculate basis expansion on [0,1] x [0,1]
f <- MFPCA:::dctFunction2D(scores = scores, argvals = list(seq(0,1,0.01), seq(0,1,0.01)))
nObs(f) # f has 10 observations

oldPar <- par(no.readonly = TRUE)
par(mfrow = c(1,1))

plot(f, obs = 1) # plot first observation
plot(f, obs = 2) # plot second observation

par(oldPar)
}
}
\seealso{
\code{\link{univExpansion}}, \code{\link{idct2D}},
  \code{\link{idct3D}}, \code{\link{dctBasis2D}},
  \code{\link{dctBasis3D}}
}
\keyword{internal}
