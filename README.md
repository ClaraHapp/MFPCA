# MFPCA
`MFPCA` is an `R`-package for calculating a PCA for multivariate functional data observed on different domains, that may also differ in dimension. The estimation algorithm relies on univariate basis expansions for each element of the multivariate functional data.

## Highlights ##

`MFPCA` allows to calculate a principal component analysis for multivariate (i.e. combined) functional data on up to three-dimensional domains:

* Standard functional data defined on a (one-dimensional) interval
* Functional data with two-dimensional domains (images)
* Functional data with three-dimensional domains (3D images, e.g. brain scans)

It implements various univariate bases:

* Univariate functional PCA (only one-dimensional domains)
* Spline bases (one- and two-dimensional domains; with optional smoothing penalty)
* Cosine bases (two- and three-dimensional domains; fast implementation built on DCT)

The representation of the data is based on the object-oriented [`funData`](https://github.com/ClaraHapp/funData) package, hence all functionalities for plotting, arithmetics etc. included therein may be used.


## Installation ##

To install the latest version directly from GitHub, please use `devtools::install_github("ClaraHapp/MFPCA")` (install [`devtools`](https://github.com/hadley/devtools) before).

If you would like to use the cosine bases make sure that the `C`-library [`fftw3`](http://www.fftw.org/) is installed on your computer before you install `MFPCA`. Otherwise, `MFPCA` is installed without the cosine bases and will throw a warning if you attempt to use functions that need `fftw3`.

## Dependencies ##

The `funData` package depends on the `R`-packages [`funData`](https://github.com/ClaraHapp/funData) for representing (multivariate) functional data, [`mgcv`](https://cran.r-project.org/web/packages/mgcv/index.html) and `methods`.

## References ##

The theoretical foundations of multivariate functional principal component analysis are described in:

C. Happ, S. Greven (2015): Multivariate Functional Principal Component Analysis for Data Observed on Different (Dimensional) Domains.
    *Submitted*. [ArXiv link](http://arxiv.org/abs/1509.02029)

## Bug reports ##

Please use [GitHub issues](https://github.com/ClaraHapp/MFPCA/issues) for reporting bugs or issues.


