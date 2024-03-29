% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univDecomp.R
\name{umpcaBasis}
\alias{umpcaBasis}
\title{Calculate an uncorrelated multilinear principal component basis
representation for functional data on two-dimensional domains}
\usage{
umpcaBasis(funDataObject, npc)
}
\arguments{
\item{funDataObject}{An object of class \code{\link[funData]{funData}}
containing the observed functional data samples (here: images) for which
the UMPCA is to be calculated.}

\item{npc}{An integer, giving the number of principal components to be
calculated.}
}
\value{
\item{scores}{A matrix of scores (coefficients) with dimension
  \code{N x k}, reflecting the weight of each principal component in each
  observation.}  \item{B}{A matrix containing the scalar product of all pairs
  of basis functions.} \item{ortho}{Logical, set to \code{FALSE}, as basis
  functions are not orthonormal.} \item{functions}{A functional data object,
  representing the functional principal component basis functions.}
}
\description{
This function calculates an uncorrelated multilinear principal component
analysis (UMPCA) representation for functional data on two-dimensional
domains. In this case, the data can be interpreted as images with \code{S1 x
S2} pixels (assuming \code{nObsPoints(funDataObject) = (S1, S2)}), i.e. the
total observed data are represented as third order tensor of dimension
\code{N x S1 x S2}.  The UMPCA of a tensor of this kind is calculated via the
\link{UMPCA} function, which is an \code{R}-version of the analogous
functions in the \code{UMPCA} MATLAB toolbox by Haiping Lu (Link:
\url{https://www.mathworks.com/matlabcentral/fileexchange/35432-uncorrelated-multilinear-principal-component-analysis-umpca},
see also references).
}
\section{Warning}{
 As this algorithm aims more at uncorrelated features than
  at an optimal reconstruction of the data, hence it might give poor results
  when used for the univariate decomposition of images in MFPCA. The function
  therefore throws a warning.
}

\references{
Haiping Lu, K.N. Plataniotis, and A.N. Venetsanopoulos,
  "Uncorrelated Multilinear Principal Component Analysis for Unsupervised
  Multilinear Subspace Learning", IEEE Transactions on Neural Networks, Vol.
  20, No. 11, Page: 1820-1836, Nov. 2009.
}
\seealso{
\code{\link{univDecomp}}
}
\keyword{internal}
