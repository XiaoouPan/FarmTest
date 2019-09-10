#' @useDynLib FarmTest
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL

#' @title Tuning-free Huber mean estimation
#' @description The function calculates adaptive Huber mean estimator from a data sample, with robustification parameter \eqn{\tau} determined by a tuning-free principle.
#' @param X An \eqn{n}-dimensional data vector.
#' @return A Huber mean estimator will be returned.
#' @references Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist., 35, 73–101.
#' @references Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2018). A New Principle for Tuning-Free Huber Regression. Preprint.
#' @seealso \code{\link{farm.cov}} for tuning-free Huber-type covariance estimation.
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

#' @title Tuning-free Huber-type covariance estimation
#' @description The function calculates adaptive Huber-type covariance estimator from a data sample, with robustification parameter \eqn{\tau} determined by a tuning-free principle.
#' For the input matrix \code{X}, both low-dimension (\eqn{p < n}) and high-dimension (\eqn{p > n}) are allowed.
#' @param X An \eqn{n} by \eqn{p} data matrix.
#' @return A \eqn{p} by \eqn{p} Huber-type covariance matrix estimator will be returned.
#' @references Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist., 35, 73–101.
#' @references Ke, Y., Minsker, S., Ren, Z., Sun, Q. and Zhou, W.-X. (2019). User-friendly covariance estimation for heavy-tailed distributions: A survey and recent results. Statis. Sci. To appear.
#' @seealso \code{\link{farm.mean}} for tuning-free Huber mean estimation.
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

#' @title Factor-adjusted robust multiple testing
#' @description This function conducts factor-adjusted robust multiple testing (FarmTest) for means of multivariate data proposed in Fan et al. (2019) via a tuning-free procedure.
#' @param X An \eqn{n} by \eqn{p} data matrix with each row being a sample.
#' @param fX An \strong{optional} factor matrix with each column being a factor for \code{X}. The number of rows of \code{fX} and \code{X} must be the same.
#' @param KX An \strong{optional} positive number of factors to be estimated for \code{X} when \code{fX} is not specified. \code{KX} cannot exceed the number of columns of \code{X}. If \code{KX} is not specified, it will be estimated internally.
#' @param Y An \strong{optional} data matrix used for two-sample FarmTest. The number of columns of \code{X} and \code{Y} must be the same.
#' @param fY An \strong{optional} factor matrix for two-sample FarmTest with each column being a factor for \code{Y}. The number of rows of \code{fY} and \code{Y} must be the same.
#' @param KY An \strong{optional} positive number of factors to be estimated for \code{Y} for two-sample FarmTest when \code{fY} is not specified. \code{KY} cannot exceed the number of columns of \code{Y}. If \code{KY} is not specified, it will be estimated internally.
#' @param h0 An \strong{optional} \eqn{p}-vector of true means, or difference in means for two-sample FarmTest. The default is a zero vector.
#' @param alternative An \strong{optional} character string specifying the alternate hypothesis, must be one of "two.sided" (default), "less" or "greater".
#' @param alpha An \strong{optional} level for controlling the false discovery rate. The value of \code{alpha} must be between 0 and 1. The default value is 0.05.
#' @return An object with S3 class \code{farm.test} containing:
#' \itemize{
#' \item \code{means} Estimated means, a vector with length \eqn{p}.
#' \item \code{stdDev} Estimated standard deviations, a vector with length \eqn{p}.
#' \item \code{loadings} Estimated factor loadings, a matrix with dimension \eqn{p} by \eqn{K}, where \eqn{K} is the number of factos.
#' \item \code{nfactors} Estimated number of factors, a positive integer.
#' \item \code{tStat} Values of test statistics, a vector with length \eqn{p}.
#' \item \code{pValues} P-values of tests, a vector with length \eqn{p}.
#' \item \code{significant} Boolean values indicating whether each test is significant, with 1 for significant and 0 for non-significant, a vector with length \eqn{p}.
#' \item \code{reject} Indices of tests that are rejected. It will show "no hypotheses rejected" if none of the tests are rejects.
#' \item \code{type} Indicates whether factor is known or unknown.
#' \item \code{n} Sample size.
#' \item \code{p} Data dimension.
#' \item \code{h0} Null hypothesis, a vector with length \eqn{p}.
#' \item \code{alpha} Alpha value.
#' \item \code{alternative} Althernative hypothesis.
#' }
#' @details For two-sample FarmTest, \code{means}, \code{stdDev}, \code{loadings}, \code{nfactors} and \code{n} will be lists of items for sample X and Y separately.
#' @details \code{alternative = "greater"} is the alternative that \eqn{\mu > \mu_0} for one-sample test or \eqn{\mu_X > \mu_Y} for two-sample test.
#' @references Fan, J., Ke, Y., Sun, Q. and Zhou, W-X. (2019). FarmTest: Factor-adjusted robust multiple testing with approximate false discovery control. J. Amer. Statist. Assoc., to appear.
#' @references Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist., 35, 73–101.
#' @references Zhou, W-X., Bose, K., Fan, J. and Liu, H. (2018). A new perspective on robust M-estimation: Finite sample theory and applications to dependence-adjusted multiple testing. Ann. Statist. 46 1904-1931.
#' @examples 
#' n = 50
#' p = 100
#' K = 3
#' muX = rep(0, p)
#' muX[1:5] = 2
#' set.seed(2019)
#' epsilonX = matrix(rnorm(p * n, 0, 1), nrow = n)
#' BX = matrix(runif(p * K, -2, 2), nrow = p)
#' fX = matrix(rnorm(K * n, 0, 1), nrow = n)
#' X = rep(1, n) %*% t(muX) + fX %*% t(BX) + epsilonX
#' # One-sample FarmTest with two sided alternative
#' output = farm.test(X)
#' output$reject
#' # One-sample FarmTest with one sided alternative
#' output = farm.test(X, alternative = "less")
#' # One-sample FarmTest with known factors
#' output = farm.test(X, fX = fX)
#' 
#' # Two-sample FarmTest
#' muY = rep(0, p)
#' muY[1:5] = 4
#' epsilonY = matrix(rnorm(p * n, 0, 1), nrow = n)
#' BY = matrix(runif(p * K, -2, 2), nrow = p)
#' fY = matrix(rnorm(K * n, 0, 1), nrow = n)
#' Y = rep(1, n) %*% t(muY) + fY %*% t(BY) + epsilonY
#' output = farm.test(X, Y = Y)
#' @export 
farm.test = function(X, fX = NULL, KX = -1, Y = NULL, fY = NULL, KY = -1, h0 = NULL, 
                     alternative = c("two.sided", "less", "greater"), alpha = 0.05) {
  p = ncol(X)
  alternative = match.arg(alternative)
  if (is.null(h0)) {
    h0 = rep(0, p)
  }
  if (length(h0) != p) {
    stop("Length of h0 must be the same as number of columns of X")
  }
  if(alpha >= 1 || alpha <= 0) {
    stop("Alpha should be between 0 and 1")
  }
  output = NULL
  reject = "no hypotheses rejected"
  if (is.null(Y) && !is.null(fX)) {
    if (nrow(fX) != nrow(X)) {
      stop("Number of rows of X and fX must be the same")
    } else {
      rst.list = farmTestFac(X, fX, h0, alpha, alternative)
      if (sum(rst.list$significant) > 0) {
        reject = which(rst.list$significant == 1)
      }
      output = list(means = rst.list$means, stdDev = rst.list$stdDev, loadings = rst.list$loadings,
                    nfactors = rst.list$nfactors, tStat = rst.list$tStat, pValues = rst.list$pValues,
                    significant = rst.list$significant, reject = reject, type = "known", n = nrow(X),
                    p = p, h0 = h0, alpha = alpha, alternative = alternative)
    }
  } else if (is.null(Y) && is.null(fX)) {
    if (KX > p) {
      stop("KX must be smaller than number of columns of X")
    } else {
      rst.list = farmTest(X, h0, KX, alpha, alternative)
      if (sum(rst.list$significant) > 0) {
        reject = which(rst.list$significant == 1)
      }
      output = list(means = rst.list$means, stdDev = rst.list$stdDev, loadings = rst.list$loadings,
                    nfactors = rst.list$nfactors, tStat = rst.list$tStat, pValues = rst.list$pValues,
                    significant = rst.list$significant, reject = reject, type = "unknown", n = nrow(X),
                    p = p, h0 = h0, alpha = alpha, alternative = alternative)
    }
  } else if (!is.null(Y) && !is.null(fX)) {
    if (ncol(X) != ncol(Y)) {
      stop("Number of columns of X and Y must be the same")
    } else if (is.null(fY)) {
      stop("Must provide factors for both or neither data matrices")
    } else if (nrow(fX) != nrow(X)) {
      stop("Number of rows of X and fX must be the same")
    } else if (nrow(fY) != nrow(Y)) {
      stop("Number of rows of Y and fY must be the same")
    } else {
      rst.list = farmTestTwoFac(X, fX, Y, fY, h0, alpha, alternative)
      if (sum(rst.list$significant) > 0) {
        reject = which(rst.list$significant == 1)
      }
      means = list(X.mean = rst.list$meansX, Y.mean = rst.list$meansY)
      stdDev = list(X.stdDev = rst.list$stdDevX, Y.stdDev = rst.list$stdDevY)
      loadings = list(X.loadings = rst.list$loadingsX, Y.loadings = rst.list$loadingsY)
      nfactors = list(X.nfactors = rst.list$nfactorsX, Y.nfactors = rst.list$nfactorsY)
      n = list(X.n = nrow(X), Y.n = nrow(Y))
      output = list(means = means, stdDev = stdDev, loadings = loadings, nfactors = nfactors, 
                    tStat = rst.list$tStat, pValues = rst.list$pValues, significant = rst.list$significant, 
                    reject = reject, type = "known", n = n, p = p, h0 = h0, alpha = alpha, alternative = alternative)
    }
  } else {
    if (ncol(X) != ncol(Y)) {
      stop("Number of columns of X and Y must be the same")
    } else if (!is.null(fY)) {
      stop("Must provide factors for both or neither data matrices")
    } else if (KX > p || KY > p) {
      stop("KX and KY must be smaller than number of columns of X and Y")
    } else {
      rst.list = farmTestTwo(X, Y, h0, KX, KY, alpha, alternative)
      if (sum(rst.list$significant) > 0) {
        reject = which(rst.list$significant == 1)
      }
      means = list(X.mean = rst.list$meansX, Y.mean = rst.list$meansY)
      stdDev = list(X.stdDev = rst.list$stdDevX, Y.stdDev = rst.list$stdDevY)
      loadings = list(X.loadings = rst.list$loadingsX, Y.loadings = rst.list$loadingsY)
      nfactors = list(X.nfactors = rst.list$nfactorsX, Y.nfactors = rst.list$nfactorsY)
      n = list(X.n = nrow(X), Y.n = nrow(Y))
      output = list(means = means, stdDev = stdDev, loadings = loadings, nfactors = nfactors, 
                    tStat = rst.list$tStat, pValues = rst.list$pValues, significant = rst.list$significant, 
                    reject = reject, type = "unknown", n = n, p = p, h0 = h0, alpha = alpha, alternative = alternative)
    } 
  }
  attr(output, "class") = "farm.test"
  return (output)
}
