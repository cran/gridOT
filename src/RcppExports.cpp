// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// costMatrix
Rcpp::NumericMatrix costMatrix(const Rcpp::NumericVector& x1, const Rcpp::NumericVector& x2, const Rcpp::NumericVector& y1, const Rcpp::NumericVector& y2, const double p1, const double p2);
RcppExport SEXP _gridOT_costMatrix(SEXP x1SEXP, SEXP x2SEXP, SEXP y1SEXP, SEXP y2SEXP, SEXP p1SEXP, SEXP p2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y2(y2SEXP);
    Rcpp::traits::input_parameter< const double >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< const double >::type p2(p2SEXP);
    rcpp_result_gen = Rcpp::wrap(costMatrix(x1, x2, y1, y2, p1, p2));
    return rcpp_result_gen;
END_RCPP
}
// frankWolfeDiscrete
Rcpp::List frankWolfeDiscrete(const arma::vec& x1, const arma::vec& x2, const arma::mat& mu, const arma::vec& y1, const arma::vec& y2, const arma::mat& nu, const double p1, const double p2, const arma::mat& startPivot, const int maxIt, const double tol, const int threads, const bool returnIt, const double rightMargin);
RcppExport SEXP _gridOT_frankWolfeDiscrete(SEXP x1SEXP, SEXP x2SEXP, SEXP muSEXP, SEXP y1SEXP, SEXP y2SEXP, SEXP nuSEXP, SEXP p1SEXP, SEXP p2SEXP, SEXP startPivotSEXP, SEXP maxItSEXP, SEXP tolSEXP, SEXP threadsSEXP, SEXP returnItSEXP, SEXP rightMarginSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y2(y2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const double >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< const double >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type startPivot(startPivotSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIt(maxItSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type returnIt(returnItSEXP);
    Rcpp::traits::input_parameter< const double >::type rightMargin(rightMarginSEXP);
    rcpp_result_gen = Rcpp::wrap(frankWolfeDiscrete(x1, x2, mu, y1, y2, nu, p1, p2, startPivot, maxIt, tol, threads, returnIt, rightMargin));
    return rcpp_result_gen;
END_RCPP
}
// frankWolfeEpsilonDiscrete
Rcpp::List frankWolfeEpsilonDiscrete(const arma::vec& x1, const arma::vec& x2, const arma::mat& mu, const arma::vec& y1, const arma::vec& y2, const arma::mat& nu, const double p1, const double p2, const arma::mat& startPivot, const int maxIt, const double tol, const int threads, const bool returnIt, const double rightMargin, const double eps);
RcppExport SEXP _gridOT_frankWolfeEpsilonDiscrete(SEXP x1SEXP, SEXP x2SEXP, SEXP muSEXP, SEXP y1SEXP, SEXP y2SEXP, SEXP nuSEXP, SEXP p1SEXP, SEXP p2SEXP, SEXP startPivotSEXP, SEXP maxItSEXP, SEXP tolSEXP, SEXP threadsSEXP, SEXP returnItSEXP, SEXP rightMarginSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y2(y2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const double >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< const double >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type startPivot(startPivotSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIt(maxItSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type returnIt(returnItSEXP);
    Rcpp::traits::input_parameter< const double >::type rightMargin(rightMarginSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(frankWolfeEpsilonDiscrete(x1, x2, mu, y1, y2, nu, p1, p2, startPivot, maxIt, tol, threads, returnIt, rightMargin, eps));
    return rcpp_result_gen;
END_RCPP
}
// frankWolfeEpsilonHistogram
Rcpp::List frankWolfeEpsilonHistogram(const arma::vec& x1, const arma::vec& x2, const arma::mat& mu, const arma::vec& y1, const arma::vec& y2, const arma::mat& nu, const double p1, const double p2, const arma::mat& startPivot, const int maxIt, const double tol, const int threads, const bool returnIt, const double rightMargin, const double eps, const double width);
RcppExport SEXP _gridOT_frankWolfeEpsilonHistogram(SEXP x1SEXP, SEXP x2SEXP, SEXP muSEXP, SEXP y1SEXP, SEXP y2SEXP, SEXP nuSEXP, SEXP p1SEXP, SEXP p2SEXP, SEXP startPivotSEXP, SEXP maxItSEXP, SEXP tolSEXP, SEXP threadsSEXP, SEXP returnItSEXP, SEXP rightMarginSEXP, SEXP epsSEXP, SEXP widthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y2(y2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const double >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< const double >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type startPivot(startPivotSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIt(maxItSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type returnIt(returnItSEXP);
    Rcpp::traits::input_parameter< const double >::type rightMargin(rightMarginSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const double >::type width(widthSEXP);
    rcpp_result_gen = Rcpp::wrap(frankWolfeEpsilonHistogram(x1, x2, mu, y1, y2, nu, p1, p2, startPivot, maxIt, tol, threads, returnIt, rightMargin, eps, width));
    return rcpp_result_gen;
END_RCPP
}
// northWestCorner
Rcpp::NumericMatrix northWestCorner(const Rcpp::NumericVector& mu, const Rcpp::NumericVector& nu, const double threshold);
RcppExport SEXP _gridOT_northWestCorner(SEXP muSEXP, SEXP nuSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(northWestCorner(mu, nu, threshold));
    return rcpp_result_gen;
END_RCPP
}
// pivotMeasure
Rcpp::NumericMatrix pivotMeasure(const Rcpp::IntegerVector& from, const Rcpp::IntegerVector& to, const Rcpp::NumericVector& mass, const int n1, const int n2, const int m1);
RcppExport SEXP _gridOT_pivotMeasure(SEXP fromSEXP, SEXP toSEXP, SEXP massSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP m1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type from(fromSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type to(toSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mass(massSEXP);
    Rcpp::traits::input_parameter< const int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< const int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< const int >::type m1(m1SEXP);
    rcpp_result_gen = Rcpp::wrap(pivotMeasure(from, to, mass, n1, n2, m1));
    return rcpp_result_gen;
END_RCPP
}
// RcppTransportPlan1d
Rcpp::DataFrame RcppTransportPlan1d(const Rcpp::NumericVector& wa, const Rcpp::NumericVector& wb, const double threshold);
RcppExport SEXP _gridOT_RcppTransportPlan1d(SEXP waSEXP, SEXP wbSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type wa(waSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type wb(wbSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppTransportPlan1d(wa, wb, threshold));
    return rcpp_result_gen;
END_RCPP
}
// RcppTransportCost1d
double RcppTransportCost1d(const Rcpp::NumericVector& a, const Rcpp::NumericVector& b, const Rcpp::NumericVector& wa, const Rcpp::NumericVector& wb, const double p, const double threshold);
RcppExport SEXP _gridOT_RcppTransportCost1d(SEXP aSEXP, SEXP bSEXP, SEXP waSEXP, SEXP wbSEXP, SEXP pSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type wa(waSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type wb(wbSEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppTransportCost1d(a, b, wa, wb, p, threshold));
    return rcpp_result_gen;
END_RCPP
}
// transportCost
double transportCost(const arma::vec& x1, const arma::vec& x2, const arma::mat& mu, const arma::vec& y1, const arma::vec& y2, const arma::mat& nu, const double p1, const double p2, const arma::mat& xi, const double threshold);
RcppExport SEXP _gridOT_transportCost(SEXP x1SEXP, SEXP x2SEXP, SEXP muSEXP, SEXP y1SEXP, SEXP y2SEXP, SEXP nuSEXP, SEXP p1SEXP, SEXP p2SEXP, SEXP xiSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y2(y2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const double >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< const double >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(transportCost(x1, x2, mu, y1, y2, nu, p1, p2, xi, threshold));
    return rcpp_result_gen;
END_RCPP
}
// transportCostFromPlan
double transportCostFromPlan(const Rcpp::IntegerVector& from, const Rcpp::IntegerVector& to, const Rcpp::NumericVector& mass, const Rcpp::NumericMatrix& cost);
RcppExport SEXP _gridOT_transportCostFromPlan(SEXP fromSEXP, SEXP toSEXP, SEXP massSEXP, SEXP costSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type from(fromSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type to(toSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mass(massSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type cost(costSEXP);
    rcpp_result_gen = Rcpp::wrap(transportCostFromPlan(from, to, mass, cost));
    return rcpp_result_gen;
END_RCPP
}
// transportPlan
Rcpp::DataFrame transportPlan(const arma::mat& mu, const arma::mat& nu, const arma::mat& xi, const double threshold);
RcppExport SEXP _gridOT_transportPlan(SEXP muSEXP, SEXP nuSEXP, SEXP xiSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(transportPlan(mu, nu, xi, threshold));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gridOT_costMatrix", (DL_FUNC) &_gridOT_costMatrix, 6},
    {"_gridOT_frankWolfeDiscrete", (DL_FUNC) &_gridOT_frankWolfeDiscrete, 14},
    {"_gridOT_frankWolfeEpsilonDiscrete", (DL_FUNC) &_gridOT_frankWolfeEpsilonDiscrete, 15},
    {"_gridOT_frankWolfeEpsilonHistogram", (DL_FUNC) &_gridOT_frankWolfeEpsilonHistogram, 16},
    {"_gridOT_northWestCorner", (DL_FUNC) &_gridOT_northWestCorner, 3},
    {"_gridOT_pivotMeasure", (DL_FUNC) &_gridOT_pivotMeasure, 6},
    {"_gridOT_RcppTransportPlan1d", (DL_FUNC) &_gridOT_RcppTransportPlan1d, 3},
    {"_gridOT_RcppTransportCost1d", (DL_FUNC) &_gridOT_RcppTransportCost1d, 6},
    {"_gridOT_transportCost", (DL_FUNC) &_gridOT_transportCost, 10},
    {"_gridOT_transportCostFromPlan", (DL_FUNC) &_gridOT_transportCostFromPlan, 4},
    {"_gridOT_transportPlan", (DL_FUNC) &_gridOT_transportPlan, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_gridOT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
