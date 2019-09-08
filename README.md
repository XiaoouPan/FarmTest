# FarmTest

**F**actor-**A**djusted **R**obust **M**ultiple **Test**ing

## Description

This package updates an [earlier version](https://github.com/kbose28/FarmTest) of **F**actor-**A**djusted **R**obust **M**ultiple **Test**ing (FarmTest) proposed in [Fan et al., 2019](https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1527700). We assume observed data *X* follow a factor model *X = &mu; + Bf + &epsilon;*, where *f* are underlying factors, *B* are factor loadings, *&epsilon;* are errors, and *&mu;* is the mean effect with null hypothesis *&mu; = &mu;<sub>0</sub>* to be tested. We assume the data is of dimension *p* with sample size *n*, leading to *p* hypothesis tests. 

FarmTest includes a robust procedure to estimate distribution parameters and accounts for strong dependence among coordinates. This method is particularly suitable for high-dimensional data when there are thousands of variables but only a small number of observations available. Moreover, the method is tailored to cases when the underlying distribution deviates from Gaussianity, which is commonly assumed in the literature.

## Main updates 

Motivated by the recent work of [Wang et al., 2018](https://www.math.ucsd.edu/~wez243/Tuning_Free.pdf) and [Ke et al., 2019](https://arxiv.org/abs/1811.01520), estimation of mean and covariance in FarmTest can be completed via a tuning-free principle, so that computationally expensive cross-validation can be avoided without lossing estimation accuracy.

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

For example, `?farm.test` will present a detailed documentation with inputs, outputs and examples of the function `farm.test`.

## Common error messages

The package `FarmTest` is implemented in `Rcpp` and `RcppArmadillo`, so the following error messages might appear when you first install it (we'll keep updating common error messages with feedback from users):

* Error: "...could not find build tools necessary to build FarmTest": For Windows you need Rtools, for Mac OS X you need to install Command Line Tools for XCode. See [this link](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) for details. 

* Error: "library not found for -lgfortran/-lquadmath": It means your gfortran binaries are out of date. This is a common environment specific issue. 

    1. In R 3.0.0 - R 3.3.0: Upgrading to R 3.4 is strongly recommended. Then go to the next step. Alternatively, you can try the instructions [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).

    2. For >= R 3.4.* : download the installer [here](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS). Then run the installer.

## Functions

There are three functions in this package:

* `farm.test`: Factor-adjusted robust multiple testing.
* `farm.mean`: Tuning-free Huber mean estimation.
* `farm.cov`: Tuning-free Huber-type covariance estimation.

## Notes 

This package is built based on an earlier version written by Bose, K., Ke, Y. and Zhou, W.-X., see [here](https://cran.r-project.org/web/packages/FarmTest/index.html) for its R-CRAN link and [here](https://github.com/kbose28/FarmTest) for its GitHub link. 

Besides, the tuning-free procedure for mean and covariance estimation is also implemented in [tfHuber](https://github.com/XiaoouPan/tfHuber) package.

## Reference

Eddelbuettel, D. and Francois, R. (2011). Rcpp: Seamless R and C++ integration. J. Stat. Softw. 40(8) 1-18. [Paper](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-introduction.pdf)

Eddelbuettel, D. and Sanderson, C. (2014). RcppArmadillo: Accelerating R with high-performance C++ linear algebra. Comput. Statist. Data Anal. 71 1054-1063. [Paper](http://dirk.eddelbuettel.com/papers/RcppArmadillo.pdf)

Fan, J., Ke, Y., Sun, Q. and Zhou, W.-X. (2017). FarmTest: Factor-adjusted robust multiple testing with approximate false discovery control. J. Amer. Statist. Assoc., to appear. [Paper](https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1527700) 

Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist. 35 73-101. [Paper](https://projecteuclid.org/euclid.aoms/1177703732)

Ke, Y., Minsker, S., Ren, Z., Sun, Q. and Zhou, W.-X. (2019). User-friendly covariance estimation for heavy-tailed distributions: A survey and recent results. Statis. Sci. To appear. [Paper](https://arxiv.org/abs/1811.01520)

Sanderson, C. and Curtin, R. (2016). Armadillo: A template-based C++ library for linear algebra. J. Open Source Softw. 1 26. [Paper](http://conradsanderson.id.au/pdfs/sanderson_armadillo_joss_2016.pdf)

Sun, Q., Zhou, W.-X. and Fan, J. (2019). Adaptive Huber regression. J. Amer. Statist. Assoc., to appear. [Paper](https://www.tandfonline.com/doi/abs/10.1080/01621459.2018.1543124)

Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2018). A new principle for tuning-free Huber regression. Preprint. [Paper](https://www.math.ucsd.edu/~wez243/Tuning_Free.pdf)

Zhou, W.-X., Bose, K., Fan, J. and Liu, H. (2018) A new perspective on robust M-estimation: Finite sample theory and applications to dependence-adjusted multiple testing. Ann. Statist. 46 1904-1931. [Paper](https://projecteuclid.org/euclid.aos/1534492823)
