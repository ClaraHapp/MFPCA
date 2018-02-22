# MFPCA 1.2

## New features
* The main `MFPCA` method now returns an object of class `MFPCAfit` for which new functions (`plot`, `predict`, `summary`, ...) have been implemented. The methods can use `names`of functional data objects (requires funData 1.2)
* It is now possible to pass a given basis decomposition for elements (e.g. from previous PCA)-
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