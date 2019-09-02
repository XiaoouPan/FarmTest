# FarmTest

**F**actor **A**djusted **R**obust **M**ultiple **Test**ing

## Installation

Install `FarmTest` from GitHub:

```r
install.packages("devtools")
library(devtools)
devtools::install_github("XiaoouPan/FarmTest")
library(FarmTest)
```

## Getting help

Help on the functions can be accessed by typing `?`, followed by function name at the R command prompt. 

For example, `?farmTest` will present a detailed documentation with inputs, outputs and examples of the function `farmTest`.

## Common error messages

The package `FarmTest` is implemented in `Rcpp` and `RcppArmadillo`, so the following error messages might appear when you first install it (we'll keep updating common error messages with feedback from users):

* Error: "...could not find build tools necessary to build FarmTest": For Windows you need Rtools, for Mac OS X you need to install Command Line Tools for XCode. See [this link](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) for details. 

* Error: "library not found for -lgfortran/-lquadmath": It means your gfortran binaries are out of date. This is a common environment specific issue. 

    1. In R 3.0.0 - R 3.3.0: Upgrading to R 3.4 is strongly recommended. Then go to the next step. Alternatively, you can try the instructions [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).

    2. For >= R 3.4.* : download the installer [here](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS). Then run the installer.

## Functions

There are two functions to conduct FarmTest:

* `farmTest`: FarmTest with unknown factors. 
* `farmTestFac`: FarmTest with known factors.

## Notes 

This package is built based on an earlier version written by Bose, K., Ke, Y. and Zhou, W.-X., see [here](https://cran.r-project.org/web/packages/FarmTest/index.html) for its R-CRAN link and [here](https://github.com/kbose28/FarmTest) for its GitHub link. 

Besides, the tuning-free procedure for mean and covariance estimation is from the [tfHuber](https://github.com/XiaoouPan/tfHuber) package.

## Reference

Eddelbuettel, D. and Francois, R. (2011). Rcpp: Seamless R and C++ integration. J. Stat. Softw. 40(8) 1-18. [Paper](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-introduction.pdf)

Eddelbuettel, D. and Sanderson, C. (2014). RcppArmadillo: Accelerating R with high-performance C++ linear algebra. Comput. Statist. Data Anal. 71 1054-1063. [Paper](http://dirk.eddelbuettel.com/papers/RcppArmadillo.pdf)

Fan, J., Ke, Y., Sun, Q. and Zhou, W.-X. (2017). FarmTest: Factor-adjusted robust multiple testing with approximate false discovery control. J. Amer. Statist. Assoc., to appear. [Paper](https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1527700) 

Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist. 35 73-101. [Paper](https://projecteuclid.org/euclid.aoms/1177703732)

Sanderson, C. and Curtin, R. (2016). Armadillo: A template-based C++ library for linear algebra. J. Open Source Softw. 1 26. [Paper](http://conradsanderson.id.au/pdfs/sanderson_armadillo_joss_2016.pdf)

Sun, Q., Zhou, W.-X. and Fan, J. (2019). Adaptive Huber regression. J. Amer. Statist. Assoc., to appear. [Paper](https://www.tandfonline.com/doi/abs/10.1080/01621459.2018.1543124)

Zhou, W.-X., Bose, K., Fan, J. and Liu, H. (2018) A new perspective on robust M-estimation: Finite sample theory and applications to dependence-adjusted multiple testing. Ann. Statist. 46 1904-1931. [Paper](https://projecteuclid.org/euclid.aos/1534492823)
