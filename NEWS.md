# Revision history for the R/qtl package

## Version 0.1.6, 2019-10-15

- Accelerate the `amalgamate` with Rcpp function `findWindow.cpp` and
  `findRepresent.cpp`

## Version 0.1.5, 2019-10-11

- A a new function `amalgamate` to merge the snps located in a
  small window
- A a step to run `amalgamate` in batchCallGeno
- Add a parameter `plot.geno` for batchCallGeno, if `plot.geno=FALSE`
  the genotype along the chromosome will not plot

## Version 0.1.4, 2019-10-6

- Add a new function `buildMstMapInput`. It will convert the 
  genotype matrix to the input of ASMap::mstmap.data.frame
- add progress bar in bacthCallGeno, and remove it in callWindowGeno

## Version 0.1.3, 2019-10-5

- First stable version fo binmapr 
