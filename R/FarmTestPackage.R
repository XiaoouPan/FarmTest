#' FarmTest: Factor-Adjusted Robust Multiple Testing
#'
#' This package performs robust multiple testing for means in the presence of known and unknown latent factors. 
#' It implements a robust procedure to estimate distribution parameters using Huber loss function with a tuning-free principle and accounts for strong dependence among coordinates via an approximate factor model. 
#' This method is particularly suitable for high dimensional data when there are many variables but only a small number of observations available. 
#' Moreover, the method is tailored to cases when the underlying distribution deviates from Gaussian, which is commonly assumed in the literature.
#' Besides the results of hypotheses testing, the estimated underlying factors and diagnostic plots are also output. 
#' 
#' For detailed information on how to use and install see its GitHub page \url{https://github.com/XiaoouPan/FarmTest}.
#'
#'
#' @references Fan, J., Ke, Y., Sun, Q. and Zhou, W-X. (2019). FarmTest: Factor-adjusted robust multiple testing with approximate false discovery control. J. Amer. Statist. Assoc., to appear.
#' @references Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist., 35, 73–101.
#' @references Ke, Y., Minsker, S., Ren, Z., Sun, Q. and Zhou, W.-X. (2019). User-friendly covariance estimation for heavy-tailed distributions: A survey and recent results. Statis. Sci. To appear.
#' @references Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2018). A New Principle for Tuning-Free Huber Regression. Preprint.
#' @references Zhou, W-X., Bose, K., Fan, J. and Liu, H. (2018). A new perspective on robust M-estimation: Finite sample theory and applications to dependence-adjusted multiple testing. Ann. Statist. 46 1904-1931.
#' @docType package
#' @name FarmTest
NULL
