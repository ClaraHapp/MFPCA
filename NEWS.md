# MFPCA 1.3-10

## New features
* Fix coercion methods for `dgTMatrix` objects as used by the DCT expansions (due to Matrix package update)
* Change class checks for `funData` objects to use the `inherits()` functionality.


# MFPCA 1.3-9

## New features
* Fix autoreconf issue
* Update link to MATLAB tensor toolbox

# MFPCA 1.3-8

## New features
* Change of email address.
* Add return value in documentation of all functions, including internal functions or functions without return value, as required by CRAN.
* Removing examples from documentation of internal functions, as required by CRAN.

# MFPCA 1.3-7

## New features
* Fix test that was failing on MKL.

# MFPCA 1.3-6

## New features
* Add warning if scores for given basis are not centered.
* Require latest `mgcv` version to work properly with MKL.
* Fix URLs where pages have been moved from http to https.

# MFPCA 1.3-5

## New features
* Add links to `funData` vignette.
* Update documentation with newest roxygen2 version.
* Clean up unused variables.

# MFPCA 1.3-4

## New features
* Update maintainer and link to Allen (2013) paper.
* Don't use `funData`'s `.intWeight` method as an internal method with `:::`. Requires new `funData` version.


# MFPCA 1.3-3

## New features
* Fix unit tests failing due to matrix / array objects having multiple classes in newest Rdevel version.

# MFPCA 1.3-2

## New features
* Suppress warnings in unit tests with fixed seed caused by changing the default RNG in latest R development version.


# MFPCA 1.3-1

## New features
* Consequent use of `seq_len(x)` instead of `1:x`.
* Use of type-safe `vapply` instead of `sapply`, where applicable.


# MFPCA 1.3

## New features
* New univariate expansion type `fda`, which allows to use all basis functions implemented in package **fda**.
* Correct typos in documentation.
* Theoretical paper is now finally published.


# MFPCA 1.2-3

## New features
* Skip third party tests on CRAN that lead to ATLAS/MKL/OpenBLAS issues. Add warning to documentation files.


# MFPCA 1.2-2

## New features
* Function `screeplot` for `MFPCAfit`objects now with integers only on x-axis and optional parameter `ylim`.
* Separate `print` function for objects of class `summary.MFPCAfit`.
* Configure-file now with correct version number.
* Links to MFPCA paper and software paper in DESCRIPTION file, as well as to ORCID.


# MFPCA 1.2-1

## New features
* Tolerance for `splineFunction2Dpen` unit tests reduced due to ATLAS/OpenBLAS/MKL issues.


# MFPCA 1.2

## New features
* The main `MFPCA` method now returns an object of class `MFPCAfit` for which new functions (`plot`, `predict`, `summary`, ...) have been implemented. The methods can use `names`of functional data objects (requires funData 1.2).
* It is now possible to pass a given basis decomposition for elements (e.g. from previous PCA).
* Univariate decompositions (xxxBasis) and expansions (xxxFunction) are not exported any more, but should be accessed using the `univDecomp` and `univExpansion` meta-functions.
* New function `multivExpansion` to calculate a multivariate basis expansion.
* More argument checking for user-facing functions including unit tests.
* Tolerance for regression tests reduced to 1e-5.
* References have been updated.


# MFPCA 1.1

## New features
* New parameter `bootstrapStrat` in `MFPCA` main function enables stratified bootstrap. Bootstrap CIs for eigenvalues are returned, too.
* `MFPCA` function now returns the eigenvectors and normFactors that are calculated within the MFPCA calculation. They can be used to calculate out-of-sample predictions.
* `fcptpaBasis` now returns eigenvalues and has new option to normalize the eigenfunctions.
* Unit test coverage has been increased to almost 100%.