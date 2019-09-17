# FarmTest

**F**actor-**A**djusted **R**obust **M**ultiple **Test**ing

## Description

This library updates an [earlier version](https://github.com/kbose28/FarmTest) that implements the **f**actor-**a**djusted **r**obust **m**ultiple **test**ing (FarmTest) method proposed by [Fan et al., 2019](https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1527700). To explicitly caputre the strong dependency among features, we assume that the *p*-dimensional data vectors *X<sub>i</sub>* follow a factor model *X<sub>i</sub> = &mu; + Bf<sub>i</sub> + &epsilon;<sub>i</sub>*, where *f<sub>i</sub>* are the common factors, *B* denotes the factor loading matrix, and *&epsilon;<sub>i</sub>* are idiosyncratic errors. Assume *f<sub>i</sub>* and *&epsilon;<sub>i</sub>* are independent and have zero means. We are interested in testing the *p* hypotheses *&mu;<sub>0</sub> = &mu;<sub>0j</sub>* simultaneously. We assume the data is of dimension *p* with sample size *n*, leading to *p* hypothesis tests. 

FarmTest includes a robust procedure to estimate distribution parameters and accounts for strong dependence among coordinates. This method is particularly suitable for high-dimensional data when there are thousands of variables but only a small number of observations available. Moreover, the method is tailored to cases when the underlying distribution deviates from Gaussianity, which is commonly assumed in the literature.

As by-products, this package also contains functions for Huber mean and Huber-type covariance matrix estimation.

## Main updates 

FarmTest introduces a robustification parameter *&tau;* while estimating mean and covariance of data sample. In the previous version, the value of *&tau;* is either specified by users or determined by cross-validation. Recently, motivated by the work of [Wang et al., 2018](https://www.math.ucsd.edu/~wez243/Tuning_Free.pdf) and [Ke et al., 2019](https://arxiv.org/abs/1811.01520), estimation of mean and covariance can be completed via a tuning-free principle, so that cross-validation, which was computationally expensive, can be avoided without lossing estimation accuracy or stability of the algorithm.

## Installation

Install `FarmTest` from GitHub:

```r
install.packages("devtools")
library(devtools)
devtools::install_github("XiaoouPan/FarmTest")
library(FarmTest)
```

## Common error messages

The package `FarmTest` is implemented in `Rcpp` and `RcppArmadillo`, so the following error messages might appear when you first install it (we'll keep updating common error messages with feedback from users):

* Error: "...could not find build tools necessary to build FarmTest": For Windows you need Rtools, for Mac OS X you need to install Command Line Tools for XCode. See [this link](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) for details. 

* Error: "library not found for -lgfortran/-lquadmath": It means your gfortran binaries are out of date. This is a common environment specific issue. 

  1. In R 3.0.0 - R 3.3.0: Upgrading to R 3.4 is strongly recommended. Then go to the next step. Alternatively, you can try the instructions [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).

  2. For R version >= 3.4.* : download the installer [here](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS). Then run the installer.
    
* Error: "cannot remove prior installation of package 'Rcpp'": This issue happens occasionally when you have installed an old version of package `Rcpp` before. Updating `Rcpp` with command `install.packages("Rcpp")` will solve the problem.

## Functions

There are four functions in this package:

* `farm.test`: Factor-adjusted robust multiple testing.
* `print.farm.test`: Print function for output of `farm.test`.
* `farm.mean`: Tuning-free Huber mean estimation.
* `farm.cov`: Tuning-free Huber-type covariance estimation.

## Getting help

Help on the functions can be accessed by typing `?`, followed by function name at the R command prompt. 

For example, `?farm.test` will present a detailed documentation with inputs, outputs and examples of the function `farm.test`.

## Testing examples

Here we generate data from factor model *X = &mu; + Bf + &epsilon;* with 3 factors. Our data has sample size 50 and dimensionality 100. The first five means are set to 2, while the others are 0.

```r
library(FarmTest)
n = 50
p = 100
K = 3
muX = rep(0, p)
muX[1:5] = 2
set.seed(2019)
epsilonX = matrix(rnorm(p * n, 0, 1), nrow = n)
BX = matrix(runif(p * K, -2, 2), nrow = p)
fX = matrix(rnorm(K * n, 0, 1), nrow = n)
X = rep(1, n) %*% t(muX) + fX %*% t(BX) + epsilonX
```

Then we conduct FarmTest for *&mu; = 0* with a two-sided alternative hypothesis, *&alpha;* level is 0.05. Factor matrix is unknown, and the number of factors will be estimated internally. 

```r
output = farm.test(X)
```

The package includes a `print.farm.test` function, which will summarize the results of `farm.test`: 
```r
output
```

To get a comprehensive impression of the performance, we repeat the above experiment for 100 times and report the average values of true positive rate (TPR), false positive rate (FPR) and false discover rate (FDR). These results can be easily reproduced.

| TPR | FPR | FDR |
| :---: | :---: | :---: | 
| 1.000 | 0.002 | 0.031 |

Finally, we present some examples to illustrate FarmTest with different purpose. For one-sided testing, just modify the `alternative` argument to be `less` or `greater`:

```r
output = farm.test(X, alternative = "less")
```

The number of factors can be specified with argument `KX` to be a positive number less than data dimension, so that `farm.test` will not estimate it. However, without any prior knowledge of the data, this is not recommended:

```r
output = farm.test(X, KX = 10)
```

For FarmTest with known factors, put the factor matrix into argument `fX`:

```r
output = farm.test(X, fX = fX)
```

For two-sample FarmTest, we generate another sample Y with same dimensionality 100, and conduct a two-sided test with unknown factors.

```r
muY = rep(0, p)
muY[1:5] = 4
epsilonY = matrix(rnorm(p * n, 0, 1), nrow = n)
BY = matrix(runif(p * K, -2, 2), nrow = p)
fY = matrix(rnorm(K * n, 0, 1), nrow = n)
Y = rep(1, n) %*% t(muY) + fY %*% t(BY) + epsilonY
output = farm.test(X, Y = Y)
```

## Huber-type estimation examples

Huber mean estimator can be obtained by `farm.mean` function. In the following example, we generate data from a centered log-normal distribution, which is asymmetric and heavy-tailed, and then estimate its mean:

```r
library(FarmTest)
set.seed(1)
n = 1000
X = rlnorm(n, 0, 1.5) - exp(1.5^2 / 2)
huberMean = farm.mean(X)
```

Huber-type covariance matrix estimator can be achieved by `farm.cov` function. In the next example, we generate data matrix from standardized *t* distribution with 3 degree of freedom, which is heavy-tailed, and estimate its covariance matrix:

```r
library(FarmTest)
set.seed(1)
n = 100
d = 50
X = matrix(rt(n * d, df = 3), n, d) / sqrt(3)
huberCov = farm.cov(X)
```

To get big pictures of these two estimators, users can run 200 independent Monte Carlo simulation as above and compare them with sample mean and covariance matrix produced by `mean` and `cov` functions. 

## Notes 

This package is built based on an earlier version written by Bose, K., Ke, Y. and Zhou, W.-X., see [here](https://cran.r-project.org/web/packages/FarmTest/index.html) for its R-CRAN link and [here](https://github.com/kbose28/FarmTest) for its GitHub link. 

Besides, the tuning-free procedure for mean and covariance estimation is also implemented in [tfHuber](https://github.com/XiaoouPan/tfHuber) package.

## License

GPL (>= 2)

## Authors

Xiaoou Pan <xip024@ucsd.edu>, Koushiki Bose <koush.bose@gmail.com>, Yuan Ke <Yuan.Ke@uga.edu>, Wen-Xin Zhou <wez243@ucsd.edu> 

## Reference

Ahn, S. C. and Horenstein, A. R. (2013). Eigenvalue ratio rest for the number of factors. Econometrica, 81(3) 1203–1227. [Paper](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA8968)

Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. J. R. Stat. Soc. Ser. B. Stat. Methodol. 57 289–300. [Paper](https://www.jstor.org/stable/2346101?seq=1#metadata_info_tab_contents)

Eddelbuettel, D. and Francois, R. (2011). Rcpp: Seamless R and C++ integration. J. Stat. Softw. 40(8) 1-18. [Paper](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-introduction.pdf)

Eddelbuettel, D. and Sanderson, C. (2014). RcppArmadillo: Accelerating R with high-performance C++ linear algebra. Comput. Statist. Data Anal. 71 1054-1063. [Paper](http://dirk.eddelbuettel.com/papers/RcppArmadillo.pdf)

Fan, J., Ke, Y., Sun, Q. and Zhou, W.-X. (2019). FarmTest: Factor-adjusted robust multiple testing with approximate false discovery control. J. Amer. Statist. Assoc., to appear. [Paper](https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1527700) 

Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist. 35 73-101. [Paper](https://projecteuclid.org/euclid.aoms/1177703732)

Ke, Y., Minsker, S., Ren, Z., Sun, Q. and Zhou, W.-X. (2019). User-friendly covariance estimation for heavy-tailed distributions: A survey and recent results. Statis. Sci. To appear. [Paper](https://arxiv.org/abs/1811.01520)

Sanderson, C. and Curtin, R. (2016). Armadillo: A template-based C++ library for linear algebra. J. Open Source Softw. 1 26. [Paper](http://conradsanderson.id.au/pdfs/sanderson_armadillo_joss_2016.pdf)

Storey, J. D. (2002). A direct approach to false discovery rates. J. R. Stat. Soc. Ser. B. Stat. Methodol. 64, 479–498. [Paper](https://www.jstor.org/stable/3088784?seq=1#metadata_info_tab_contents)

Sun, Q., Zhou, W.-X. and Fan, J. (2019). Adaptive Huber regression. J. Amer. Statist. Assoc., to appear. [Paper](https://www.tandfonline.com/doi/abs/10.1080/01621459.2018.1543124)

Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2018). A new principle for tuning-free Huber regression. Preprint. [Paper](https://www.math.ucsd.edu/~wez243/Tuning_Free.pdf)

Zhou, W.-X., Bose, K., Fan, J. and Liu, H. (2018) A new perspective on robust M-estimation: Finite sample theory and applications to dependence-adjusted multiple testing. Ann. Statist. 46 1904-1931. [Paper](https://projecteuclid.org/euclid.aos/1534492823)
