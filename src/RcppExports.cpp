// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// f1
double f1(const double x, const arma::vec& resSq, const int n);
RcppExport SEXP _FarmTest_f1(SEXP xSEXP, SEXP resSqSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type resSq(resSqSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(f1(x, resSq, n));
    return rcpp_result_gen;
END_RCPP
}
// rootf1
double rootf1(const arma::vec& resSq, const int n, double low, double up, const double tol, const int maxIte);
RcppExport SEXP _FarmTest_rootf1(SEXP resSqSEXP, SEXP nSEXP, SEXP lowSEXP, SEXP upSEXP, SEXP tolSEXP, SEXP maxIteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type resSq(resSqSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type low(lowSEXP);
    Rcpp::traits::input_parameter< double >::type up(upSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIte(maxIteSEXP);
    rcpp_result_gen = Rcpp::wrap(rootf1(resSq, n, low, up, tol, maxIte));
    return rcpp_result_gen;
END_RCPP
}
// f2
double f2(const double x, const arma::vec& resSq, const int n, const int d, const int N);
RcppExport SEXP _FarmTest_f2(SEXP xSEXP, SEXP resSqSEXP, SEXP nSEXP, SEXP dSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type resSq(resSqSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(f2(x, resSq, n, d, N));
    return rcpp_result_gen;
END_RCPP
}
// rootf2
double rootf2(const arma::vec& resSq, const int n, const int d, const int N, double low, double up, const double tol, const int maxIte);
RcppExport SEXP _FarmTest_rootf2(SEXP resSqSEXP, SEXP nSEXP, SEXP dSEXP, SEXP NSEXP, SEXP lowSEXP, SEXP upSEXP, SEXP tolSEXP, SEXP maxIteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type resSq(resSqSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type low(lowSEXP);
    Rcpp::traits::input_parameter< double >::type up(upSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIte(maxIteSEXP);
    rcpp_result_gen = Rcpp::wrap(rootf2(resSq, n, d, N, low, up, tol, maxIte));
    return rcpp_result_gen;
END_RCPP
}
// huberMean
double huberMean(const arma::vec& X, const int n, const double epsilon, const int iteMax);
RcppExport SEXP _FarmTest_huberMean(SEXP XSEXP, SEXP nSEXP, SEXP epsilonSEXP, SEXP iteMaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const int >::type iteMax(iteMaxSEXP);
    rcpp_result_gen = Rcpp::wrap(huberMean(X, n, epsilon, iteMax));
    return rcpp_result_gen;
END_RCPP
}
// huberMeanVec
arma::vec huberMeanVec(const arma::mat& X, const int n, const int p, const double epsilon, const int iteMax);
RcppExport SEXP _FarmTest_huberMeanVec(SEXP XSEXP, SEXP nSEXP, SEXP pSEXP, SEXP epsilonSEXP, SEXP iteMaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const int >::type iteMax(iteMaxSEXP);
    rcpp_result_gen = Rcpp::wrap(huberMeanVec(X, n, p, epsilon, iteMax));
    return rcpp_result_gen;
END_RCPP
}
// hMeanCov
double hMeanCov(const arma::vec& Z, const int n, const int d, const int N, const double epsilon, const int iteMax);
RcppExport SEXP _FarmTest_hMeanCov(SEXP ZSEXP, SEXP nSEXP, SEXP dSEXP, SEXP NSEXP, SEXP epsilonSEXP, SEXP iteMaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const int >::type iteMax(iteMaxSEXP);
    rcpp_result_gen = Rcpp::wrap(hMeanCov(Z, n, d, N, epsilon, iteMax));
    return rcpp_result_gen;
END_RCPP
}
// huberCov
Rcpp::List huberCov(const arma::mat& X, const int n, const int p);
RcppExport SEXP _FarmTest_huberCov(SEXP XSEXP, SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(huberCov(X, n, p));
    return rcpp_result_gen;
END_RCPP
}
// mad
double mad(const arma::vec& x);
RcppExport SEXP _FarmTest_mad(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(mad(x));
    return rcpp_result_gen;
END_RCPP
}
// sgn
int sgn(const double x);
RcppExport SEXP _FarmTest_sgn(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sgn(x));
    return rcpp_result_gen;
END_RCPP
}
// standardize
arma::mat standardize(arma::mat X, const arma::rowvec& mx, const arma::vec& sx, const int p);
RcppExport SEXP _FarmTest_standardize(SEXP XSEXP, SEXP mxSEXP, SEXP sxSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type mx(mxSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sx(sxSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(standardize(X, mx, sx, p));
    return rcpp_result_gen;
END_RCPP
}
// updateHuber
void updateHuber(const arma::mat& Z, const arma::vec& res, arma::vec& der, arma::vec& grad, const int n, const double tau, const double n1);
RcppExport SEXP _FarmTest_updateHuber(SEXP ZSEXP, SEXP resSEXP, SEXP derSEXP, SEXP gradSEXP, SEXP nSEXP, SEXP tauSEXP, SEXP n1SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type res(resSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type der(derSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grad(gradSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double >::type n1(n1SEXP);
    updateHuber(Z, res, der, grad, n, tau, n1);
    return R_NilValue;
END_RCPP
}
// huberReg
arma::vec huberReg(const arma::mat& X, arma::vec Y, const int n, const int p, const double tol, const double constTau, const int iteMax);
RcppExport SEXP _FarmTest_huberReg(SEXP XSEXP, SEXP YSEXP, SEXP nSEXP, SEXP pSEXP, SEXP tolSEXP, SEXP constTauSEXP, SEXP iteMaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const double >::type constTau(constTauSEXP);
    Rcpp::traits::input_parameter< const int >::type iteMax(iteMaxSEXP);
    rcpp_result_gen = Rcpp::wrap(huberReg(X, Y, n, p, tol, constTau, iteMax));
    return rcpp_result_gen;
END_RCPP
}
// huberRegCoef
arma::vec huberRegCoef(const arma::mat& X, arma::vec Y, const int n, const int p, const double tol, const double constTau, const int iteMax);
RcppExport SEXP _FarmTest_huberRegCoef(SEXP XSEXP, SEXP YSEXP, SEXP nSEXP, SEXP pSEXP, SEXP tolSEXP, SEXP constTauSEXP, SEXP iteMaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const double >::type constTau(constTauSEXP);
    Rcpp::traits::input_parameter< const int >::type iteMax(iteMaxSEXP);
    rcpp_result_gen = Rcpp::wrap(huberRegCoef(X, Y, n, p, tol, constTau, iteMax));
    return rcpp_result_gen;
END_RCPP
}
// huberRegItcp
double huberRegItcp(const arma::mat& X, arma::vec Y, const int n, const int p, const double tol, const double constTau, const int iteMax);
RcppExport SEXP _FarmTest_huberRegItcp(SEXP XSEXP, SEXP YSEXP, SEXP nSEXP, SEXP pSEXP, SEXP tolSEXP, SEXP constTauSEXP, SEXP iteMaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const double >::type constTau(constTauSEXP);
    Rcpp::traits::input_parameter< const int >::type iteMax(iteMaxSEXP);
    rcpp_result_gen = Rcpp::wrap(huberRegItcp(X, Y, n, p, tol, constTau, iteMax));
    return rcpp_result_gen;
END_RCPP
}
// getP
arma::vec getP(const arma::vec& T, const std::string alternative);
RcppExport SEXP _FarmTest_getP(SEXP TSEXP, SEXP alternativeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const std::string >::type alternative(alternativeSEXP);
    rcpp_result_gen = Rcpp::wrap(getP(T, alternative));
    return rcpp_result_gen;
END_RCPP
}
// getPboot
arma::vec getPboot(const arma::vec& mu, const arma::mat& boot, const arma::vec& h0, const std::string alternative, const int p, const int B);
RcppExport SEXP _FarmTest_getPboot(SEXP muSEXP, SEXP bootSEXP, SEXP h0SEXP, SEXP alternativeSEXP, SEXP pSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type boot(bootSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< const std::string >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(getPboot(mu, boot, h0, alternative, p, B));
    return rcpp_result_gen;
END_RCPP
}
// getRej
arma::uvec getRej(const arma::vec& Prob, const double alpha, const int p);
RcppExport SEXP _FarmTest_getRej(SEXP ProbSEXP, SEXP alphaSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Prob(ProbSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(getRej(Prob, alpha, p));
    return rcpp_result_gen;
END_RCPP
}
// getRatio
arma::vec getRatio(const arma::vec& eigenVal, const int n, const int p);
RcppExport SEXP _FarmTest_getRatio(SEXP eigenValSEXP, SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eigenVal(eigenValSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(getRatio(eigenVal, n, p));
    return rcpp_result_gen;
END_RCPP
}
// rmTest
Rcpp::List rmTest(const arma::mat& X, const arma::vec& h0, const double alpha, const std::string alternative);
RcppExport SEXP _FarmTest_rmTest(SEXP XSEXP, SEXP h0SEXP, SEXP alphaSEXP, SEXP alternativeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type alternative(alternativeSEXP);
    rcpp_result_gen = Rcpp::wrap(rmTest(X, h0, alpha, alternative));
    return rcpp_result_gen;
END_RCPP
}
// rmTestBoot
Rcpp::List rmTestBoot(const arma::mat& X, const arma::vec& h0, const double alpha, const std::string alternative, const int B);
RcppExport SEXP _FarmTest_rmTestBoot(SEXP XSEXP, SEXP h0SEXP, SEXP alphaSEXP, SEXP alternativeSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< const int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(rmTestBoot(X, h0, alpha, alternative, B));
    return rcpp_result_gen;
END_RCPP
}
// rmTestTwo
Rcpp::List rmTestTwo(const arma::mat& X, const arma::mat& Y, const arma::vec& h0, const double alpha, const std::string alternative);
RcppExport SEXP _FarmTest_rmTestTwo(SEXP XSEXP, SEXP YSEXP, SEXP h0SEXP, SEXP alphaSEXP, SEXP alternativeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type alternative(alternativeSEXP);
    rcpp_result_gen = Rcpp::wrap(rmTestTwo(X, Y, h0, alpha, alternative));
    return rcpp_result_gen;
END_RCPP
}
// rmTestTwoBoot
Rcpp::List rmTestTwoBoot(const arma::mat& X, const arma::mat& Y, const arma::vec& h0, const double alpha, const std::string alternative, const int B);
RcppExport SEXP _FarmTest_rmTestTwoBoot(SEXP XSEXP, SEXP YSEXP, SEXP h0SEXP, SEXP alphaSEXP, SEXP alternativeSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< const int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(rmTestTwoBoot(X, Y, h0, alpha, alternative, B));
    return rcpp_result_gen;
END_RCPP
}
// farmTest
Rcpp::List farmTest(const arma::mat& X, const arma::vec& h0, int K, const double alpha, const std::string alternative);
RcppExport SEXP _FarmTest_farmTest(SEXP XSEXP, SEXP h0SEXP, SEXP KSEXP, SEXP alphaSEXP, SEXP alternativeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type alternative(alternativeSEXP);
    rcpp_result_gen = Rcpp::wrap(farmTest(X, h0, K, alpha, alternative));
    return rcpp_result_gen;
END_RCPP
}
// farmTestTwo
Rcpp::List farmTestTwo(const arma::mat& X, const arma::mat& Y, const arma::vec& h0, int KX, int KY, const double alpha, const std::string alternative);
RcppExport SEXP _FarmTest_farmTestTwo(SEXP XSEXP, SEXP YSEXP, SEXP h0SEXP, SEXP KXSEXP, SEXP KYSEXP, SEXP alphaSEXP, SEXP alternativeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< int >::type KX(KXSEXP);
    Rcpp::traits::input_parameter< int >::type KY(KYSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type alternative(alternativeSEXP);
    rcpp_result_gen = Rcpp::wrap(farmTestTwo(X, Y, h0, KX, KY, alpha, alternative));
    return rcpp_result_gen;
END_RCPP
}
// farmTestFac
Rcpp::List farmTestFac(const arma::mat& X, const arma::mat& fac, const arma::vec& h0, const double alpha, const std::string alternative);
RcppExport SEXP _FarmTest_farmTestFac(SEXP XSEXP, SEXP facSEXP, SEXP h0SEXP, SEXP alphaSEXP, SEXP alternativeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type fac(facSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type alternative(alternativeSEXP);
    rcpp_result_gen = Rcpp::wrap(farmTestFac(X, fac, h0, alpha, alternative));
    return rcpp_result_gen;
END_RCPP
}
// farmTestFacBoot
Rcpp::List farmTestFacBoot(const arma::mat& X, const arma::mat& fac, const arma::vec& h0, const double alpha, const std::string alternative, const int B);
RcppExport SEXP _FarmTest_farmTestFacBoot(SEXP XSEXP, SEXP facSEXP, SEXP h0SEXP, SEXP alphaSEXP, SEXP alternativeSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type fac(facSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< const int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(farmTestFacBoot(X, fac, h0, alpha, alternative, B));
    return rcpp_result_gen;
END_RCPP
}
// farmTestTwoFac
Rcpp::List farmTestTwoFac(const arma::mat& X, const arma::mat& facX, const arma::mat& Y, const arma::mat& facY, const arma::vec& h0, const double alpha, const std::string alternative);
RcppExport SEXP _FarmTest_farmTestTwoFac(SEXP XSEXP, SEXP facXSEXP, SEXP YSEXP, SEXP facYSEXP, SEXP h0SEXP, SEXP alphaSEXP, SEXP alternativeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type facX(facXSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type facY(facYSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type alternative(alternativeSEXP);
    rcpp_result_gen = Rcpp::wrap(farmTestTwoFac(X, facX, Y, facY, h0, alpha, alternative));
    return rcpp_result_gen;
END_RCPP
}
// farmTestTwoFacBoot
Rcpp::List farmTestTwoFacBoot(const arma::mat& X, const arma::mat& facX, const arma::mat& Y, const arma::mat& facY, const arma::vec& h0, const double alpha, const std::string alternative, const int B);
RcppExport SEXP _FarmTest_farmTestTwoFacBoot(SEXP XSEXP, SEXP facXSEXP, SEXP YSEXP, SEXP facYSEXP, SEXP h0SEXP, SEXP alphaSEXP, SEXP alternativeSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type facX(facXSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type facY(facYSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< const int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(farmTestTwoFacBoot(X, facX, Y, facY, h0, alpha, alternative, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FarmTest_f1", (DL_FUNC) &_FarmTest_f1, 3},
    {"_FarmTest_rootf1", (DL_FUNC) &_FarmTest_rootf1, 6},
    {"_FarmTest_f2", (DL_FUNC) &_FarmTest_f2, 5},
    {"_FarmTest_rootf2", (DL_FUNC) &_FarmTest_rootf2, 8},
    {"_FarmTest_huberMean", (DL_FUNC) &_FarmTest_huberMean, 4},
    {"_FarmTest_huberMeanVec", (DL_FUNC) &_FarmTest_huberMeanVec, 5},
    {"_FarmTest_hMeanCov", (DL_FUNC) &_FarmTest_hMeanCov, 6},
    {"_FarmTest_huberCov", (DL_FUNC) &_FarmTest_huberCov, 3},
    {"_FarmTest_mad", (DL_FUNC) &_FarmTest_mad, 1},
    {"_FarmTest_sgn", (DL_FUNC) &_FarmTest_sgn, 1},
    {"_FarmTest_standardize", (DL_FUNC) &_FarmTest_standardize, 4},
    {"_FarmTest_updateHuber", (DL_FUNC) &_FarmTest_updateHuber, 7},
    {"_FarmTest_huberReg", (DL_FUNC) &_FarmTest_huberReg, 7},
    {"_FarmTest_huberRegCoef", (DL_FUNC) &_FarmTest_huberRegCoef, 7},
    {"_FarmTest_huberRegItcp", (DL_FUNC) &_FarmTest_huberRegItcp, 7},
    {"_FarmTest_getP", (DL_FUNC) &_FarmTest_getP, 2},
    {"_FarmTest_getPboot", (DL_FUNC) &_FarmTest_getPboot, 6},
    {"_FarmTest_getRej", (DL_FUNC) &_FarmTest_getRej, 3},
    {"_FarmTest_getRatio", (DL_FUNC) &_FarmTest_getRatio, 3},
    {"_FarmTest_rmTest", (DL_FUNC) &_FarmTest_rmTest, 4},
    {"_FarmTest_rmTestBoot", (DL_FUNC) &_FarmTest_rmTestBoot, 5},
    {"_FarmTest_rmTestTwo", (DL_FUNC) &_FarmTest_rmTestTwo, 5},
    {"_FarmTest_rmTestTwoBoot", (DL_FUNC) &_FarmTest_rmTestTwoBoot, 6},
    {"_FarmTest_farmTest", (DL_FUNC) &_FarmTest_farmTest, 5},
    {"_FarmTest_farmTestTwo", (DL_FUNC) &_FarmTest_farmTestTwo, 7},
    {"_FarmTest_farmTestFac", (DL_FUNC) &_FarmTest_farmTestFac, 5},
    {"_FarmTest_farmTestFacBoot", (DL_FUNC) &_FarmTest_farmTestFacBoot, 6},
    {"_FarmTest_farmTestTwoFac", (DL_FUNC) &_FarmTest_farmTestTwoFac, 7},
    {"_FarmTest_farmTestTwoFacBoot", (DL_FUNC) &_FarmTest_farmTestTwoFacBoot, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_FarmTest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
