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

double f2(const double x, const arma::vec& resSq, const int n, const int d, const int N) {
  return arma::sum(arma::min(resSq, x * arma::ones(N))) / (N * x) - 
    (2 * std::log(d) + std::log(n)) / n;
}

double rootf2(const arma::vec& resSq, const int n, const int d, const int N, double low, double up, 
              const double tol = 0.0001, const int maxIte = 500) {
  int ite = 0;
  double mid, val;
  while (ite <= maxIte && up - low > tol) {
    mid = (up + low) / 2;
    val = f2(mid, resSq, n, d, N);
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
double huberMean(const arma::vec& X, const int n, const double epsilon = 0.0001, const int iteMax = 500) {
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

double hMeanCov(const arma::vec& Z, const int n, const int d, const int N, 
                const double epsilon = 0.0001, const int iteMax = 500) {
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
    tauNew = std::sqrt((long double)rootf2(resSq, n, d, N, arma::min(resSq), arma::sum(resSq)));
    w = arma::min(tauNew * arma::ones(N) / arma::abs(res), arma::ones(N));
    muNew = arma::as_scalar(Z.t() * w) / arma::sum(w);
    iteNum++;
  }
  return muNew;
}

// [[Rcpp::export]]
Rcpp::List huberCov(const arma::mat& X, const int n, const int p) {
  arma::vec mu(p);
  arma::mat sigmaHat(p, p);
  for (int j = 0; j < p; j++) {
    mu(j) = huberMean(X.col(j), n);
    double theta = huberMean(arma::square(X.col(j)), n);
    double temp = mu(j) * mu(j);
    if (theta > temp) {
      theta -= temp;
    }
    sigmaHat(j, j) = theta;
  }
  int N = n * (n - 1) >> 1;
  arma::mat Y(N, p);
  for (int i = 0, k = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      Y.row(k++) = X.row(i) - X.row(j);
    }
  }
  for (int i = 0; i < p - 1; i++) {
    for (int j = i + 1; j < p; j++) {
      sigmaHat(i, j) = sigmaHat(j, i) = hMeanCov(Y.col(i) % Y.col(j) / 2, n, p, N);
    }
  }
  return Rcpp::List::create(Rcpp::Named("means") = mu, Rcpp::Named("cov") = sigmaHat);
}

arma::vec huberReg(const arma::mat& X, const arma::vec& Y, const int n, const int d, 
                   const double epsilon = 0.0001, const double constTau = 1.345, const int iteMax = 500) {
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

arma::vec huberRegItcp(const arma::mat& X, const arma::vec& Y, const int n, const int d, 
                       const double epsilon = 0.0001, const double constTau = 1.345, const int iteMax = 500) {
  arma::mat Z(n, d + 1);
  Z.cols(1, d) = X;
  Z.col(0) = arma::ones(n);
  arma::vec thetaOld = arma::zeros(d + 1);
  arma::vec thetaNew = arma::solve(Z.t() * Z, Z.t() * Y);
  double tauOld = 0;
  double tauNew = std::sqrt((long double)arma::sum(arma::square(Y - Z * thetaNew)) / (n - d)) *
    std::sqrt((long double)n / std::log((long double)(d + std::log(n * d))));
  double mad;
  int iteNum = 0;
  arma::vec res(n), WY(n);
  arma::mat WZ(n, d + 1);
  while ((arma::norm(thetaNew - thetaOld, "inf") > epsilon || std::abs(tauNew - tauOld) > epsilon) 
           && iteNum < iteMax) {
    thetaOld = thetaNew;
    tauOld = tauNew;
    res = Y - Z * thetaOld;
    mad = arma::median(arma::abs(res - arma::median(res))) / 0.6744898;
    tauNew = constTau * mad;
    WZ = Z;
    WY = Y;
    for (int i = 0; i < n; i++) {
      double w = tauNew / std::abs(res(i));
      if (w < 1) {
        WZ.row(i) *= w;
        WY(i) *= w;
      }
    }
    thetaNew = arma::solve(Z.t() * WZ, Z.t() * WY);
    iteNum++;
  }
  thetaNew(0) = huberMean(Y - X * thetaNew.rows(1, d), n);
  return thetaNew;
}

arma::vec getP(const arma::vec& T, const std::string alternative) {
  arma::vec rst;
  if (alternative == "two.sided") {
    rst = 2 * arma::normcdf(-arma::abs(T));
  } else if (alternative == "less") {
    rst = arma::normcdf(T);
  } else {
    rst = arma::normcdf(-T);
  }
  return rst;
}

arma::uvec getRej(const arma::vec& Prob, const double alpha, const int p) {
  double piHat = (double)arma::sum(Prob > alpha) / ((1 - alpha) * p);
  arma::vec z = arma::sort(Prob);
  double pAlpha = -1.0;
  for (int i = p - 1; i >= 0; i--) {
    if (z(i) * piHat * p <= alpha * (i + 1)) {
      pAlpha = z(i);
      break;
    }
  }
  return Prob <= pAlpha;
}

arma::vec getRatio(const arma::vec& eigenVal, const int n, const int p) {
  int temp = std::min(n, p);
  int len = temp < 4 ? temp - 1 : temp >> 1;
  if (len == 0) {
    arma::vec rst(1);
    rst(0) = eigenVal(p - 1);
    return rst;
  }
  arma::vec ratio(len);
  double comp = eigenVal(p - 1) / eigenVal(p - 2);
  ratio(0) = comp;
  for (int i = 1; i < len; i++) {
    ratio(i) = eigenVal(p - 1 - i) / eigenVal(p - 2 - i);
  }
  return ratio;
}

// [[Rcpp::export]]
Rcpp::List rmTest(const arma::mat& X, const arma::vec& h0, const double alpha = 0.05, 
                  const std::string alternative = "two.sided") {
  int n = X.n_rows, p = X.n_cols;
  arma::vec mu(p), sigma(p);
  for (int j = 0; j < p; j++) {
    mu(j) = huberMean(X.col(j), n);
    double theta = huberMean(arma::square(X.col(j)), n);
    double temp = mu(j) * mu(j);
    if (theta > temp) {
      theta -= temp;
    }
    sigma(j) = theta;
  }
  sigma = arma::sqrt(sigma / n);
  arma::vec T = (mu - h0) / sigma;
  arma::vec Prob = getP(T, alternative);
  arma::uvec significant = getRej(Prob, alpha, p);
  return Rcpp::List::create(Rcpp::Named("means") = mu, Rcpp::Named("stdDev") = sigma,
                            Rcpp::Named("tStat") = T, Rcpp::Named("pValues") = Prob, 
                            Rcpp::Named("significant") = significant);
}

// [[Rcpp::export]]
Rcpp::List rmTestTwo(const arma::mat& X, const arma::mat& Y, const arma::vec& h0, 
                     const double alpha = 0.05, const std::string alternative = "two.sided") {
  int nX = X.n_rows, nY = Y.n_rows, p = X.n_cols;
  arma::vec muX(p), sigmaX(p), muY(p), sigmaY(p);
  for (int j = 0; j < p; j++) {
    muX(j) = huberMean(X.col(j), nX);
    muY(j) = huberMean(Y.col(j), nY);
    double theta = huberMean(arma::square(X.col(j)), nX);
    double temp = muX(j) * muX(j);
    if (theta > temp) {
      theta -= temp;
    }
    sigmaX(j) = theta;
    theta = huberMean(arma::square(Y.col(j)), nY);
    temp = muY(j) * muY(j);
    if (theta > temp) {
      theta -= temp;
    }
    sigmaY(j) = theta;
  }
  arma::vec T = (muX - muY - h0) / arma::sqrt(sigmaX / nX + sigmaY / nY);
  sigmaX = arma::sqrt(sigmaX / nX);
  sigmaY = arma::sqrt(sigmaY / nY);
  arma::vec Prob = getP(T, alternative);
  arma::uvec significant = getRej(Prob, alpha, p);
  return Rcpp::List::create(Rcpp::Named("meansX") = muX, Rcpp::Named("meansY") = muY, 
                            Rcpp::Named("stdDevX") = sigmaX, Rcpp::Named("stdDevY") = sigmaY,
                            Rcpp::Named("tStat") = T, Rcpp::Named("pValues") = Prob, 
                            Rcpp::Named("significant") = significant);
}

// [[Rcpp::export]]
Rcpp::List farmTest(const arma::mat& X, const arma::vec& h0, int K = -1, const double alpha = 0.05, 
                    const std::string alternative = "two.sided") {
  int n = X.n_rows, p = X.n_cols;
  Rcpp::List listCov = huberCov(X, n, p);
  arma::vec mu = listCov["means"];
  arma::mat sigmaHat = listCov["cov"];
  arma::vec sigma = sigmaHat.diag();
  arma::vec eigenVal;
  arma::mat eigenVec;
  arma::eig_sym(eigenVal, eigenVec, sigmaHat);
  arma::vec ratio;
  if (K <= 0) {
    ratio = getRatio(eigenVal, n, p);
    K = arma::index_max(ratio) + 1;
  }
  arma::mat B(p, K);
  for (int i = 1; i <= K; i++) {
    double lambda = std::sqrt((long double)std::max(eigenVal(p - i), 0.0));
    B.col(i - 1) = lambda * eigenVec.col(p - i);
  }
  arma::vec f = huberReg(B, arma::mean(X, 0).t(), p, K);
  for (int j = 0; j < p; j++) {
    double temp = arma::norm(B.row(j), 2);
    if (sigma(j) > temp * temp) {
      sigma(j) -= temp * temp;
    }
  }
  mu -= B * f;
  sigma = arma::sqrt(sigma / n);
  arma::vec T = (mu - h0) / sigma;
  arma::vec Prob = getP(T, alternative);
  arma::uvec significant = getRej(Prob, alpha, p);
  return Rcpp::List::create(Rcpp::Named("means") = mu, Rcpp::Named("stdDev") = sigma,
                            Rcpp::Named("loadings") = B, Rcpp::Named("nfactors") = K, 
                            Rcpp::Named("tStat") = T, Rcpp::Named("pValues") = Prob, 
                            Rcpp::Named("significant") = significant, Rcpp::Named("eigens") = eigenVal,
                            Rcpp::Named("ratio") = ratio);
}

// [[Rcpp::export]]
Rcpp::List farmTestTwo(const arma::mat& X, const arma::mat& Y, const arma::vec& h0, int KX = -1, 
                       int KY = -1, const double alpha = 0.05, const std::string alternative = "two.sided") {
  int nX = X.n_rows, nY = Y.n_rows, p = X.n_cols;
  Rcpp::List listCov = huberCov(X, nX, p);
  arma::vec muX = listCov["means"];
  arma::mat sigmaHat = listCov["cov"];
  arma::vec sigmaX = sigmaHat.diag();
  arma::vec eigenValX, eigenValY;
  arma::mat eigenVec;
  arma::eig_sym(eigenValX, eigenVec, sigmaHat);
  arma::vec ratioX, ratioY;
  if (KX <= 0) {
    ratioX = getRatio(eigenValX, nX, p);
    KX = arma::index_max(ratioX) + 1;
  }
  arma::mat BX(p, KX);
  for (int i = 1; i <= KX; i++) {
    double lambda = std::sqrt((long double)std::max(eigenValX(p - i), 0.0));
    BX.col(i - 1) = lambda * eigenVec.col(p - i);
  }
  arma::vec fX = huberReg(BX, arma::mean(X, 0).t(), p, KX);
  listCov = huberCov(Y, nY, p);
  arma::vec muY = listCov["means"];
  sigmaHat = Rcpp::as<arma::mat>(listCov["cov"]);
  arma::vec sigmaY = sigmaHat.diag();
  arma::eig_sym(eigenValY, eigenVec, sigmaHat);
  if (KY <= 0) {
    ratioY = getRatio(eigenValY, nY, p);
    KY = arma::index_max(ratioY) + 1;
  }
  arma::mat BY(p, KY);
  for (int i = 1; i <= KY; i++) {
    double lambda = std::sqrt((long double)std::max(eigenValY(p - i), 0.0));
    BY.col(i - 1) = lambda * eigenVec.col(p - i);
  }
  arma::vec fY = huberReg(BY, arma::mean(Y, 0).t(), p, KY);
  for (int j = 0; j < p; j++) {
    double temp = arma::norm(BX.row(j), 2);
    if (sigmaX(j) > temp * temp) {
      sigmaX(j) -= temp * temp;
    }
    temp = arma::norm(BY.row(j), 2);
    if (sigmaY(j) > temp * temp) {
      sigmaY(j) -= temp * temp;
    }
  }
  muX -= BX * fX;
  muY -= BY * fY;
  arma::vec T = (muX - muY - h0) / arma::sqrt(sigmaX / nX + sigmaY / nY);
  sigmaX = arma::sqrt(sigmaX / nX);
  sigmaY = arma::sqrt(sigmaY / nY);
  arma::vec Prob = getP(T, alternative);
  arma::uvec significant = getRej(Prob, alpha, p);
  return Rcpp::List::create(Rcpp::Named("meansX") = muX, Rcpp::Named("meansY") = muY, 
                            Rcpp::Named("stdDevX") = sigmaX, Rcpp::Named("stdDevY") = sigmaY,
                            Rcpp::Named("loadingsX") = BX, Rcpp::Named("loadingsY") = BY,
                            Rcpp::Named("nfactorsX") = KX, Rcpp::Named("nfactorsY") = KY,
                            Rcpp::Named("tStat") = T, Rcpp::Named("pValues") = Prob, 
                            Rcpp::Named("significant") = significant, Rcpp::Named("eigensX") = eigenValX,
                            Rcpp::Named("eigensY") = eigenValY, Rcpp::Named("ratioX") = ratioX,
                            Rcpp::Named("ratioY") = ratioY);
}

// [[Rcpp::export]]
Rcpp::List farmTestFac(const arma::mat& X, const arma::mat& fac, const arma::vec& h0, 
                       const double alpha = 0.05, const std::string alternative = "two.sided") {
  int n = X.n_rows, p = X.n_cols, K = fac.n_cols;
  arma::mat Sigma = arma::cov(fac);
  arma::vec mu(p), sigma(p);
  arma::vec theta, beta;
  arma::mat B(p, K);
  for (int j = 0; j < p; j++) {
    theta = huberRegItcp(fac, X.col(j), n, K);
    mu(j) = theta(0);
    beta = theta.rows(1, K);
    B.row(j) = beta.t();
    double sig = huberMean(arma::square(X.col(j)), n);
    double temp = mu(j) * mu(j);
    if (sig > temp) {
      sig -= temp;
    }
    temp = arma::as_scalar(beta.t() * Sigma * beta);
    if (sig > temp * temp) {
      sig -= temp * temp;
    }
    sigma(j) = sig;
  }
  sigma = arma::sqrt(sigma / n);
  arma::vec T = (mu - h0) / sigma;
  arma::vec Prob = getP(T, alternative);
  arma::uvec significant = getRej(Prob, alpha, p);
  return Rcpp::List::create(Rcpp::Named("means") = mu, Rcpp::Named("stdDev") = sigma,
                            Rcpp::Named("loadings") = B, Rcpp::Named("nfactors") = K,
                            Rcpp::Named("tStat") = T, Rcpp::Named("pValues") = Prob, 
                            Rcpp::Named("significant") = significant);
}

// [[Rcpp::export]]
Rcpp::List farmTestTwoFac(const arma::mat& X, const arma::mat& facX, const arma::mat& Y, 
                          const arma::mat& facY, const arma::vec& h0, const double alpha = 0.05, 
                          const std::string alternative = "two.sided") {
  int nX = X.n_rows, nY = Y.n_rows, p = X.n_cols, KX = facX.n_cols, KY = facY.n_cols;
  arma::mat SigmaX = arma::cov(facX);
  arma::mat SigmaY = arma::cov(facY);
  arma::vec muX(p), sigmaX(p), muY(p), sigmaY(p);
  arma::vec theta, beta;
  arma::mat BX(p, KX), BY(p, KY);
  for (int j = 0; j < p; j++) {
    theta = huberRegItcp(facX, X.col(j), nX, KX);
    muX(j) = theta(0);
    beta = theta.rows(1, KX);
    BX.row(j) = beta.t();
    double sig = huberMean(arma::square(X.col(j)), nX);
    double temp = muX(j) * muX(j);
    if (sig > temp) {
      sig -= temp;
    }
    temp = arma::as_scalar(beta.t() * SigmaX * beta);
    if (sig > temp * temp) {
      sig -= temp * temp;
    }
    sigmaX(j) = sig;
    theta = huberRegItcp(facY, Y.col(j), nY, KY);
    muY(j) = theta(0);
    beta = theta.rows(1, KY);
    BY.row(j) = beta.t();
    sig = huberMean(arma::square(Y.col(j)), nY);
    temp = muY(j) * muY(j);
    if (sig > temp) {
      sig -= temp;
    }
    temp = arma::as_scalar(beta.t() * SigmaY * beta);
    if (sig > temp * temp) {
      sig -= temp * temp;
    }
    sigmaY(j) = sig;
  }
  arma::vec T = (muX - muY - h0) / arma::sqrt(sigmaX / nX + sigmaY / nY);
  sigmaX = arma::sqrt(sigmaX / nX);
  sigmaY = arma::sqrt(sigmaY / nY);
  arma::vec Prob = getP(T, alternative);
  arma::uvec significant = getRej(Prob, alpha, p);
  return Rcpp::List::create(Rcpp::Named("meansX") = muX, Rcpp::Named("meansY") = muY, 
                            Rcpp::Named("stdDevX") = sigmaX, Rcpp::Named("stdDevY") = sigmaY,
                            Rcpp::Named("loadingsX") = BX, Rcpp::Named("loadingsY") = BY,
                            Rcpp::Named("nfactorsX") = KX, Rcpp::Named("nfactorsY") = KY,
                            Rcpp::Named("tStat") = T, Rcpp::Named("pValues") = Prob, 
                            Rcpp::Named("significant") = significant);
}
