// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// theta_ml
double theta_ml(double theta, Eigen::VectorXd y, Eigen::VectorXd mu, Eigen::VectorXd weights);
RcppExport SEXP _scZINB_theta_ml(SEXP thetaSEXP, SEXP ySEXP, SEXP muSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(theta_ml(theta, y, mu, weights));
    return rcpp_result_gen;
END_RCPP
}
// updateEM
Rcpp::List updateEM(Eigen::VectorXd y, Eigen::VectorXd betas, const Eigen::VectorXd offsetx, Eigen::VectorXd gammas, const Eigen::VectorXd offsetz, const Eigen::MatrixXd X, double theta, const double lambda, const double tau, const bool irlsConv, const Eigen::VectorXd betaWeights, const Eigen::VectorXd gammaWeights, const int maxIT, const int loud, double eps, const int model, const int penType, const int maxIT2);
RcppExport SEXP _scZINB_updateEM(SEXP ySEXP, SEXP betasSEXP, SEXP offsetxSEXP, SEXP gammasSEXP, SEXP offsetzSEXP, SEXP XSEXP, SEXP thetaSEXP, SEXP lambdaSEXP, SEXP tauSEXP, SEXP irlsConvSEXP, SEXP betaWeightsSEXP, SEXP gammaWeightsSEXP, SEXP maxITSEXP, SEXP loudSEXP, SEXP epsSEXP, SEXP modelSEXP, SEXP penTypeSEXP, SEXP maxIT2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type betas(betasSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type offsetx(offsetxSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type gammas(gammasSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type offsetz(offsetzSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const bool >::type irlsConv(irlsConvSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type betaWeights(betaWeightsSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type gammaWeights(gammaWeightsSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIT(maxITSEXP);
    Rcpp::traits::input_parameter< const int >::type loud(loudSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const int >::type penType(penTypeSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIT2(maxIT2SEXP);
    rcpp_result_gen = Rcpp::wrap(updateEM(y, betas, offsetx, gammas, offsetz, X, theta, lambda, tau, irlsConv, betaWeights, gammaWeights, maxIT, loud, eps, model, penType, maxIT2));
    return rcpp_result_gen;
END_RCPP
}
// updateEM2
Rcpp::List updateEM2(Eigen::VectorXd y, Eigen::VectorXd betas, const Eigen::VectorXd offsetx, Eigen::VectorXd gammas, const Eigen::VectorXd offsetz, const Eigen::MatrixXd X, double theta, const double lambda, const double tau, const bool irlsConv, const Eigen::VectorXd betaWeights, const Eigen::VectorXd gammaWeights, const int maxIT, const int loud, double eps, const int penType);
RcppExport SEXP _scZINB_updateEM2(SEXP ySEXP, SEXP betasSEXP, SEXP offsetxSEXP, SEXP gammasSEXP, SEXP offsetzSEXP, SEXP XSEXP, SEXP thetaSEXP, SEXP lambdaSEXP, SEXP tauSEXP, SEXP irlsConvSEXP, SEXP betaWeightsSEXP, SEXP gammaWeightsSEXP, SEXP maxITSEXP, SEXP loudSEXP, SEXP epsSEXP, SEXP penTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type betas(betasSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type offsetx(offsetxSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type gammas(gammasSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type offsetz(offsetzSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const bool >::type irlsConv(irlsConvSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type betaWeights(betaWeightsSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type gammaWeights(gammaWeightsSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIT(maxITSEXP);
    Rcpp::traits::input_parameter< const int >::type loud(loudSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const int >::type penType(penTypeSEXP);
    rcpp_result_gen = Rcpp::wrap(updateEM2(y, betas, offsetx, gammas, offsetz, X, theta, lambda, tau, irlsConv, betaWeights, gammaWeights, maxIT, loud, eps, penType));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scZINB_theta_ml", (DL_FUNC) &_scZINB_theta_ml, 4},
    {"_scZINB_updateEM", (DL_FUNC) &_scZINB_updateEM, 18},
    {"_scZINB_updateEM2", (DL_FUNC) &_scZINB_updateEM2, 16},
    {NULL, NULL, 0}
};

RcppExport void R_init_scZINB(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}