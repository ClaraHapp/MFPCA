# MFPCA 1.1

## New features
* New parameter `bootstrapStrat` in `MFPCA` main function enables stratified bootstrap. Bootstrap CIs for eigenvalues are returned, too.
* `MFPCA` function now returns the eigenvectors and normFactors that are calculated within the MFPCA calculation. They can be used to calculate out-of-sample predictions.
* `fcptpaBasis` now returns eigenvalues and has new option to normalize the eigenfunctions.
* Unit test coverage has been increased to almost 100%.