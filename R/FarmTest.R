#' @useDynLib FarmTest
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL

#' Tuning-free Huber mean estimation
#' 
#' The function calculates adaptive Huber mean estimator from a data sample, with robustification parameter \eqn{\tau} determined by a tuning-free principle.
#' 
#' @param X An \eqn{n}-dimensional data vector.
#' @return The Huber mean estimator.
#' @references Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist., 35, 73–101.
#' @references Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2018). A New Principle for Tuning-Free Huber Regression. Preprint.
#' @examples
#' set.seed(2019)
#' n = 1000
#' X = rlnorm(n, 0, 1.5) - exp(1.5^2 / 2)
#' mu = farm.mean(X)
#' @export
farm.mean = function(X){
  n = length(X)
  return (huberMean(X, n))
}

#' Tuning-free Huber-type covariance estimation
#' 
#' The function calculates adaptive Huber-type covariance estimator from a data sample, with robustification parameter \eqn{\tau} determined by a tuning-free principle.
#' 
#' @param X An \eqn{n} by \eqn{p} data matrix.
#' @return An \eqn{p} by \eqn{p} Huber-type covariance matrix estimator.
#' @references Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist., 35, 73–101.
#' @references Ke, Y., Minsker, S., Ren, Z., Sun, Q. and Zhou, W.-X. (2019). User-friendly covariance estimation for heavy-tailed distributions: A survey and recent results. Statis. Sci. To appear. 
#' @examples
#' set.seed(2019)
#' n = 100
#' d = 50
#' X = matrix(rt(n * d, df = 3), n, d) / sqrt(3)
#' Sigma = farm.cov(X)
#' @export
farm.cov = function(X) {
  n = nrow(X)
  p = ncol(X)
  return (huberCov(X, n, p)$cov)
}
