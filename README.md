# FarmTest

**F**actor-**A**djusted **R**obust **M**ultiple **Test**ing

## Description

The FarmTest library implements the **F**actor-**A**djusted **R**obust **M**ultiple **Test**ing method proposed by [Fan et al., 2019](https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1527700). An earlier version can be found [here](https://github.com/kbose28/FarmTest). Let *X* be a *p*-dimensional random vector with mean *&mu; = (&mu;<sub>1</sub>,...,&mu;<sub>p</sub>)<sup>T</sup>*. This library carries out simultaneous inference on the *p* hypotheses *H<sub>0j</sub> : &mu;<sub>j</sub> = &mu;<sub>0j</sub>*.
To explicitly caputre the strong dependency among features, we assume that the data vectors *X<sub>i</sub>* that are independently drawn from *X* follow a factor model: *X<sub>i</sub> = &mu; + Bf<sub>i</sub> + &epsilon;<sub>i</sub>*, where *f<sub>i</sub>* are the common factors, *B* denotes the factor loading matrix, and *&epsilon;<sub>i</sub>* are idiosyncratic errors. Specifically, we consider three different scenarios with observable factors, latent factors and a mixture of covariates and latent factors. Assume *f<sub>i</sub>* and *&epsilon;<sub>i</sub>* are independent and have zero means. The number of hypotheses *p* may be comparable to or considerably exceed the sample size *n*.

FarmTest implements a series of adaptive Huber methods combined with fast data-drive tuning schemes to estimate model parameters and construct test statistics that are robust against heavy-tailed and/or assymetric error distributions. Extensions to two-sample simultaneous mean comparison problems are also included. As by-products, this library also contains functions that compute adaptive Huber mean and covariance matrix estimators that are of independent interest.

## Main updates 

The FarmTest method involves multiple tuning parameters for fitting the factor models. In the case of latent factors, the algorithm first computes a robust covariance matrix estimator, and then use the eigenvalue ratio method along with SVD to estimate the number of factors and loading vectors. It is therefore computationally expenstive to select all the tuning parameters via cross-validation. Instead, the current version makes use of the fast data-driven tuning scheme proposed by [Ke et al., 2019](https://arxiv.org/abs/1811.01520), which significantly reduces the computational cost.

## Installation

Install `FarmTest` from GitHub:

```r
install.packages("devtools")
library(devtools)
devtools::install_github("XiaoouPan/FarmTest")
library(FarmTest)
```

## Common error messages

The library `FarmTest` is coded in `Rcpp` and `RcppArmadillo`. The following error messages might arise when you first install it (we'll keep updating common error messages based on user feedback):

* Error: "...could not find build tools necessary to build FarmTest": For Windows OS you need Rtools; for Mac OS X you need to install Command Line Tools for XCode. See [this link](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) for details. 

* Error: "library not found for -lgfortran/-lquadmath": It means your gfortran binaries are out-of-date. This is a common environment specific issue. 

  1. In R 3.0.0 - R 3.3.0: Upgrading to R 3.4 is strongly recommended. Then go to the next step. Alternatively, you can try the instructions [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).

  2. For R version >= 3.4.* : Download and run the installer [here](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS).
    
* Error: "cannot remove prior installation of package 'Rcpp'": This issue happens occasionally when you have installed an old version of the package `Rcpp` before. Update `Rcpp` with command `install.packages("Rcpp")` will solve the problem.

## Functions

There are four functions in this library:

* `farm.test`: Factor-adjusted robust multiple testing.
* `print.farm.test`: Print function for the outputs of `farm.test`.
* `farm.mean`: Tuning-free Huber mean estimation.
* `farm.cov`: Tuning-free Huber-type covariance estimation.

## Getting help

Help on the functions can be accessed by typing `?`, followed by function name at the R command prompt. 

For example, `?farm.test` will present a detailed documentation with inputs, outputs and examples of the function `farm.test`.

## Examples

First generate data from a three-factor model *X = &mu; + Bf + &epsilon;*. The sample size and dimension (the number of hypotheses) are taken to be 50 and 100, respectively. The number of nonnulls is 5.

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

In this case, the factors are unobservable and thus need to be recovered from data. Assume one is interested in simultaneous inference on the means with two-sided alternatives. For a desired FDR level *&alpha;=0.05*, run FarmTest as follows:

```r
output = farm.test(X)
```

The library includes a `print.farm.test` function, which summarizes the results of `farm.test`: 

```r
output
```

Based on 100 simulations, we report below the average values of the true positive rate (TPR), false positive rate (FPR) and false discover rate (FDR).

| TPR | FPR | FDR |
| :---: | :---: | :---: | 
| 1.000 | 0.002 | 0.031 |

In addition, we illustrate the use of FarmTest under different circumstances. For one-sided alternatives, modify the `alternative` argument to be `less` or `greater`:

```r
output = farm.test(X, alternative = "less")
```

The number of factors can be user-specified. It should be a positive integer that is less than the minumum between sample size and number of hypotheses. Without any subjective ground of the data, this is not recommended.

```r
output = farm.test(X, KX = 10)
```

In a case with observable factors, put the *n* by *K* factor matrix into argument `fX`:

```r
output = farm.test(X, fX = fX)
```

As an extension to two-sample problems, we generate another sample Y with the same dimension 100, and conduct a two-sided test with latent factors.

```r
muY = rep(0, p)
muY[1:5] = 4
epsilonY = matrix(rnorm(p * n, 0, 1), nrow = n)
BY = matrix(runif(p * K, -2, 2), nrow = p)
fY = matrix(rnorm(K * n, 0, 1), nrow = n)
Y = rep(1, n) %*% t(muY) + fY %*% t(BY) + epsilonY
output = farm.test(X, Y = Y)
```
 
Robust mean and covariance matrix estimation is not only an important step in the FarmTest, but also of independent interest in many other problems. We write separate functions `farm.mean` and `farm.cov` for this purpose.

```r
library(FarmTest)
set.seed(1)
n = 1000
X = rlnorm(n, 0, 1.5)
huberMean = farm.mean(X)

n = 100
d = 50
X = matrix(rt(n * d, df = 3), n, d)
huberCov = farm.cov(X)
```

## Remark 

This library is built upon an earlier version written by Bose, K., Ke, Y. and Zhou, W.-X. ([CRAN](https://cran.r-project.org/web/packages/FarmTest/index.html), [GitHub](https://github.com/kbose28/FarmTest)). Another library named tfHuber that implements data-driven robust mean and covariance matrix estimation as well as standard and l1-regularized Huber regression can be found [here](https://github.com/XiaoouPan/tfHuber).

## License

GPL (>= 2)

## Author(s)

Xiaoou Pan <xip024@ucsd.edu>, Yuan Ke <yuan.ke@uga.edu>, Wen-Xin Zhou <wez243@ucsd.edu> 

## Maintainer

Xiaoou Pan <xip024@ucsd.edu>

## References

Ahn, S. C. and Horenstein, A. R. (2013). Eigenvalue ratio rest for the number of factors. *Econometrica * **81**(3) 1203–1227. [Paper](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA8968)

Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *J. R. Stat. Soc. Ser. B. Stat. Methodol.* **57** 289–300. [Paper](https://www.jstor.org/stable/2346101?seq=1#metadata_info_tab_contents)

Eddelbuettel, D. and Francois, R. (2011). Rcpp: Seamless R and C++ integration. *J. Stat. Softw.* **40**(8) 1-18. [Paper](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-introduction.pdf)

Eddelbuettel, D. and Sanderson, C. (2014). RcppArmadillo: Accelerating R with high-performance C++ linear algebra. *Comput. Statist. Data Anal.* **71** 1054-1063. [Paper](http://dirk.eddelbuettel.com/papers/RcppArmadillo.pdf)

Fan, J., Ke, Y., Sun, Q. and Zhou, W.-X. (2019). FarmTest: Factor-adjusted robust multiple testing with approximate false discovery control. *J. Amer. Statist. Assoc.*, to appear. [Paper](https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1527700) 

Huber, P. J. (1964). Robust estimation of a location parameter. *Ann. Math. Statist.* **35** 73-101. [Paper](https://projecteuclid.org/euclid.aoms/1177703732)

Ke, Y., Minsker, S., Ren, Z., Sun, Q. and Zhou, W.-X. (2019). User-friendly covariance estimation for heavy-tailed distributions: A survey and recent results. *Statis. Sci.*, to appear. [Paper](https://arxiv.org/abs/1811.01520)

Sanderson, C. and Curtin, R. (2016). Armadillo: A template-based C++ library for linear algebra. *J. Open Source Softw.* **1** 26. [Paper](http://conradsanderson.id.au/pdfs/sanderson_armadillo_joss_2016.pdf)

Storey, J. D. (2002). A direct approach to false discovery rates. *J. R. Stat. Soc. Ser. B. Stat. Methodol.* **64** 479–498. [Paper](https://www.jstor.org/stable/3088784?seq=1#metadata_info_tab_contents)

Sun, Q., Zhou, W.-X. and Fan, J. (2019). Adaptive Huber regression. *J. Amer. Statist. Assoc.*, to appear. [Paper](https://www.tandfonline.com/doi/abs/10.1080/01621459.2018.1543124)

Zhou, W.-X., Bose, K., Fan, J. and Liu, H. (2018). A new perspective on robust M-estimation: Finite sample theory and applications to dependence-adjusted multiple testing. *Ann. Statist.* **46** 1904-1931. [Paper](https://projecteuclid.org/euclid.aos/1534492823)
