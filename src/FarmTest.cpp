# include <RcppArmadillo.h>
# include <algorithm>
# include <string>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

double f1(const double x, const arma::vec& resSq, const int n) {
  return arma::sum(arma::min(resSq, x * arma::ones(n))) / (n * x) - std::log(n) / n;
}

double rootf1(const arma::vec& resSq, const int n, double low, double up, 
              const double tol = 0.0001, const int maxIte = 500) {
  int ite = 0;
  double mid, val;
  while (ite <= maxIte && up - low > tol) {
    mid = (up + low) / 2;
    val = f1(mid, resSq, n);
    if (val == 0) {
      return mid;
    } else if (val < 0) {
      up = mid;
    } else {
      low = mid;
    }
    ite++;
  }
  return (low + up) / 2;
}

double f2(const double x, const arma::vec& resSq, const int n, const int d) {
  int N = n * (n - 1) >> 1;
  return arma::sum(arma::min(resSq, x * arma::ones(N))) / (N * x) - 
    (2 * std::log(d) + std::log(n)) / n;
}

double rootf2(const arma::vec& resSq, const int n, const int d, double low, double up, 
              const double tol = 0.0001, const int maxIte = 500) {
  int ite = 0;
  double mid, val;
  while (ite <= maxIte && up - low > tol) {
    mid = (up + low) / 2;
    val = f2(mid, resSq, n, d);
    if (val == 0) {
      return mid;
    } else if (val < 0) {
      up = mid;
    } else {
      low = mid;
    }
    ite++;
  }
  return (low + up) / 2;
}

// [[Rcpp::export]]
double huberMean(const arma::vec& X, const double epsilon = 0.0001, const int iteMax = 500) {
  int n = X.size();
  double muOld = 0;
  double muNew = arma::mean(X);
  double tauOld = 0;
  double tauNew = arma::stddev(X) * std::sqrt((long double)n / std::log(n));
  int iteNum = 0;
  arma::vec res(n), resSq(n), w(n);
  while (((std::abs(muNew - muOld) > epsilon) || (std::abs(tauNew - tauOld) > epsilon)) && iteNum < iteMax) {
    muOld = muNew;
    tauOld = tauNew;
    res = X - muOld * arma::ones(n);
    resSq = arma::square(res);
    tauNew = std::sqrt((long double)rootf1(resSq, n, arma::min(resSq), arma::sum(resSq)));
    w = arma::min(tauNew * arma::ones(n) / arma::abs(res), arma::ones(n));
    muNew = arma::as_scalar(X.t() * w) / arma::sum(w);
    iteNum++;
  }
  return muNew;
}

double hMeanCov(const arma::vec& Z, const int n, const int d, const double epsilon = 0.0001, 
                const int iteMax = 500) {
  int N = Z.size();
  double muOld = 0;
  double muNew = arma::mean(Z);
  double tauOld = 0;
  double tauNew = arma::stddev(Z) * std::sqrt((long double)n / (2 * std::log(d) + std::log(n)));
  int iteNum = 0;
  arma::vec res(n), resSq(n), w(n);
  while (((std::abs(muNew - muOld) > epsilon) || (std::abs(tauNew - tauOld) > epsilon)) && iteNum < iteMax) {
    muOld = muNew;
    tauOld = tauNew;
    res = Z - muOld * arma::ones(N);
    resSq = arma::square(res);
    tauNew = std::sqrt((long double)rootf2(resSq, n, d, arma::min(resSq), arma::sum(resSq)));
    w = arma::min(tauNew * arma::ones(N) / arma::abs(res), arma::ones(N));
    muNew = arma::as_scalar(Z.t() * w) / arma::sum(w);
    iteNum++;
  }
  return muNew;
}

// [[Rcpp::export]]
arma::mat huberCov(const arma::mat& X, const double epsilon = 0.0001, const int iteMax = 500) {
  int n = X.n_rows;
  int d = X.n_cols;
  int N = n * (n - 1) >> 1;
  arma::mat Y(N, d);
  for (int i = 0, k = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      Y.row(k++) = X.row(i) - X.row(j);
    }
  }
  arma::mat rst(d, d);
  for (int i = 0; i < d; i++) {
    for (int j = i; j < d; j++) {
      rst(i, j) = rst(j, i) = hMeanCov(Y.col(i) % Y.col(j) / 2, n, d, epsilon, iteMax);
    }
  }
  return rst;
}

// [[Rcpp::export]]
arma::vec huberReg(const arma::mat& X, const arma::vec& Y, const double epsilon = 0.0001, 
                   const double constTau = 1.345, const int iteMax = 500) {
  int n = X.n_rows;
  int d = X.n_cols;
  arma::vec thetaOld = arma::zeros(d);
  arma::vec thetaNew = arma::solve(X.t() * X, X.t() * Y);
  double tauOld = 0;
  double tauNew = std::sqrt((long double)arma::sum(arma::square(Y - X * thetaNew)) / (n - d)) *
    std::sqrt((long double)n / std::log((long double)(d + std::log(n * d))));
  double mad;
  int iteNum = 0;
  arma::vec res(n), WY(n);
  arma::mat WX(n, d);
  while ((arma::norm(thetaNew - thetaOld, "inf") > epsilon || std::abs(tauNew - tauOld) > epsilon) 
           && iteNum < iteMax) {
    thetaOld = thetaNew;
    tauOld = tauNew;
    res = Y - X * thetaOld;
    mad = arma::median(arma::abs(res - arma::median(res))) / 0.6744898;
    tauNew = constTau * mad;
    WX = X;
    WY = Y;
    for (int i = 0; i < n; i++) {
      double w = tauNew / std::abs(res(i));
      if (w < 1) {
        WX.row(i) *= w;
        WY(i) *= w;
      }
    }
    thetaNew = arma::solve(X.t() * WX, X.t() * WY);
    iteNum++;
  }
  return thetaNew;
}

//' @export
// [[Rcpp::export]]
Rcpp::List farmTest(const arma::mat& X, const double alpha, const int K) {
  int n = X.n_rows, p = X.n_cols;
  arma::mat sigmaHat = huberCov(X);
  arma::vec mu(p);
  for (int j = 0; j < p; j++) {
    mu(j) = huberMean(X.col(j));
  }
  arma::vec eigenVal;
  arma::mat eigenVec;
  arma::eig_sym(eigenVal, eigenVec, sigmaHat);
  arma::mat B(p, K);
  for (int i = 1; i <= K; i++) {
    double lambda = std::sqrt((long double)std::max(eigenVal(p - i), 0.0));
    B.col(i - 1) = lambda * eigenVec.col(p - i);
  }
  arma::vec f = huberReg(B, arma::mean(X, 0).t());
  arma::vec sigma(p);
  for (int j = 0; j < p; j++) {
    double theta = huberMean(arma::square(X.col(j)));
    double temp = mu(j) * mu(j);
    if (theta > temp) {
      theta -= temp;
    }
    sigma(j) = theta;
    temp = arma::norm(B.row(j), 2) * arma::norm(B.row(j), 2);
    if (sigma(j) > temp) {
      sigma(j) -= temp;
    }
  }
  mu -= B * f;
  sigma = arma::sqrt(sigma / n);
  arma::vec T = mu / sigma;
  arma::vec z = arma::sort(arma::abs(T));
  arma::vec Prob = 2 * arma::normcdf(-arma::abs(T));
  double piHat = (double)arma::sum(Prob > alpha) / ((1 - alpha) * p);
  int idx = -1;
  for (int i = 0; i < p; i++) {
    double fdp = 2 * p * piHat * arma::normcdf(-z(i)) / (p - i);
    if (fdp <= alpha) {
      idx = i;
      break;
    }
  }
  double zAlpha = idx == -1 ? z(p - 1) + 1 : z(idx);
  arma::uvec reject = arma::abs(T) >= zAlpha;
  return Rcpp::List::create(Rcpp::Named("means") = mu, Rcpp::Named("stdDev") = sigma,
                            Rcpp::Named("loadings") = B, Rcpp::Named("tStat") = T, 
                            Rcpp::Named("pValues") = Prob, Rcpp::Named("alpha") = alpha, 
                            Rcpp::Named("criVal") = zAlpha, Rcpp::Named("reject") = reject);
}
