// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <fstream>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// DEBUG BLOCK
bool tt = FALSE;


// [[Rcpp::depends(RcppEigen)]]

// Convenience functions for dealing with GLM classes

/* Link Index:
 1: EXP
 2: LOGIT
 Default: NORMAL
 */ 

double linkInv(double eta, const int link){
    double eeta = exp(eta);
    
    switch(link){
    case 1 :
        return(eeta);
    case 2 :
        return(eeta/(1 + eeta));
    default :
        return(eta);
    }
};

double linkFunc(double mu, const int link){
    switch(link){
    case 1 :
        return(log(mu));
    case 2 :
        return(log(mu/(1-mu)));
    default :
        return(mu);
    }
};
/* Family Index:
 1: NEGBIN
 2: BINOMIAL
 Default: EMPTY
 */

// Function to calculate Vdiag
Eigen::VectorXd irlsV(const Eigen::VectorXd mu, const Eigen::VectorXd zGy, 
                      const double theta, const int family){
    
    Eigen::VectorXd Vdiag(mu.size());
    
    switch(family){
    case 1 :
        for(int i = 0; i < mu.size(); i += 1){
            Vdiag(i) = (1-zGy(i))*(mu(i)/(1 + mu(i)/theta));
        }
        return(Vdiag);
    case 2 :
        for(int i = 0; i < mu.size(); i += 1){
            Vdiag(i) = mu(i)*(1-mu(i));
            if(Vdiag(i) == 0){
                Vdiag(i) = 1e-5;
            } else if(Vdiag(i) == 1){
                Vdiag(i) = 1-1e-5;
            }
        }
        return(Vdiag);
    default :
        return(Vdiag);
    }
};

// Function to calculate Z for IRLS
Eigen::VectorXd irlsZ(const Eigen::VectorXd eta, const Eigen::VectorXd y, 
                      const int family, Eigen::VectorXd Vdiag){
    
    Eigen::VectorXd Z(y.size());
    
    switch(family){
    case 1 :
        for(int i = 0; i < eta.size(); i += 1){
            double mui = linkInv(eta(i), family);
            
            if(std::isnan(mui)){
                if(eta(i) > 100){
                    mui = std::numeric_limits<double>::max();
                }
            }
            
            Z(i) = eta(i) + (y(i) - mui)/mui;
        }
        return(Z);
    case 2 :
        for(int i = 0; i < eta.size(); i += 1){
            double pij = linkInv(eta(i), family);
            
            if(std::isnan(pij)){
                if(eta(i) > 100){
                    pij = 1-1e-16;
                }
            }
            
            Z(i) = eta(i) + (y(i) - pij)/Vdiag(i);
        }
        return(Z);
    default :
        return(Z);
    }   
};

// Function to update mu
Eigen::VectorXd updateMu(const Eigen::VectorXd eta, const int link){
    
    Eigen::VectorXd mu(eta.size());
    
    for(int i = 0; i < eta.size(); i += 1){
        mu(i) = linkInv(eta(i), link);
        
        if(std::isnan(mu(i))){
            if(eta(i) > 100 && link == 1){
                mu(i) = std::numeric_limits<double>::max();
            } else if(eta(i) > 100 && link == 2){
                mu(i) = 1-1e-16;
            }
        }
        
        if(mu(i) > 1-1e-16 && link == 2){
            mu(i) = 1-1e-16;
        } else if(mu(i) < 1e-16){
            mu(i) = 1e-16;
        }
    }
    return(mu);
};

// Function to calculate eta 
Eigen::VectorXd calcEta(const Eigen::VectorXd coef, const Eigen::MatrixXd X, 
                        const Eigen::VectorXd offset){
    
    Eigen::VectorXd eta = X*coef + offset;
    return(eta);
}

// C++ Function for dnbinom
Eigen::VectorXd dnbinom(Eigen::VectorXd y, Eigen::VectorXd mu, double theta, 
                        const bool log) {
    Eigen::VectorXd result(y.size());
    
    for(int i = 0; i < y.size(); i +=1){
        result(i) = R::dnbinom_mu(y(i), theta, mu(i), log);
    }
    return result;
}

// C++ Function for Log Likelihood
double updLikGLM(const Eigen::VectorXd y, const Eigen::VectorXd mu, 
                 const Eigen::VectorXd pij, const double theta){
    
    Eigen::VectorXd fnb(mu.size());
    double loglik = 0;
    
    fnb = dnbinom(y, mu, theta, 1);
    
    for(int i = 0; i < mu.size(); i += 1){
        loglik +=  fnb(i);
    }
    
    return loglik;
}

// C++ Function for Log Likelihood
double updLik(const Eigen::VectorXd y, const Eigen::VectorXd mu, 
              const Eigen::VectorXd pij, const double theta){
    
    Eigen::VectorXd fnb(mu.size());
    double loglik = 0;
    
    fnb = dnbinom(y, mu, theta, 1);
    
    for(int i = 0; i < mu.size(); i += 1){
        if(y(i) == 0){
            loglik +=  log( pij(i) + exp( log(1-pij(i)) + fnb(i)));
        } else {
            loglik += log(1-pij(i)) + fnb(i);
        }
    }
    
    return loglik;
}

// C++ Function for theta calculation
void score_info(double theta, Eigen::VectorXd mu, Eigen::VectorXd y, 
                double* score, double* info, Eigen::VectorXd weights) {
    int i;
    double score1=0.0, info1=0.0;
    double mui, yi, scorei, infoi;
    
    for (i=0; i<y.size(); i++) {
        yi  = y(i);
        mui = mu(i);
        
        scorei  = weights(i)*(R::digamma(yi + theta) - R::digamma(theta) + 
            log(theta) + 1 - log(theta + mui) - (yi + theta)/(mui + theta));
        score1 += scorei;
        
        infoi   = weights(i)*(-1*R::trigamma(theta + yi) + R::trigamma(theta) - 
            1/theta + 2/(mui + theta) - 
            (yi + theta)/std::pow((mui + theta),2));
        info1  += infoi;
    }
    
    *score = score1;
    *info  = info1;
}

// [[Rcpp::export]]
double theta_ml(double theta, Eigen::VectorXd y, Eigen::VectorXd mu, 
                Eigen::VectorXd weights) {
    
    int it = 0;
    double del = 1000;
    double score;
    double info;
    double minTheta = 1e-2;
    double maxTheta = 1/minTheta;
    double theta0 = theta;
    
    while(it < 10 && std::abs(del) > 1e-5) {
        
        score_info(theta0, mu, y, &score, &info, weights);
        del     = score/info;
        theta0 += del;
        it     += 1;
        
        if (theta0 > maxTheta) {
            theta0 = maxTheta;
        } else if (theta0 < 0) {
            theta0 -= del;
            break;
        } else if (theta0 < minTheta) {
            theta0 = minTheta;
            break;
        }
    }
    return(theta0);
}


// Coefficient update function - log penalty
double updateCoefk_log(Eigen::VectorXd y, Eigen::VectorXd betas, 
                       const Eigen::VectorXd offsetx, 
                       Eigen::VectorXd gammas, const int family, const Eigen::VectorXd offsetz, 
                       const Eigen::MatrixXd X, const double theta, const double lambda, 
                       const double tau, int k, const Eigen::VectorXd zGy, 
                       const bool irlsConv, const double betaWeight, 
                       const double gammaWeight, const double countEM, const double countk){
    
    double diff = 1000;
    int count = 0;
    
    Eigen::VectorXd eta(X.rows());
    Eigen::VectorXd mu(X.rows());
    
    switch(family){
    case 1 :
        eta = calcEta(betas, X, offsetx);
        break;
    case 2 :
        eta = calcEta(gammas, X, offsetz);
        y = zGy;
        break;
    }
    
    mu = updateMu(eta, family);
    
    double sto = 0;
    while(diff > 1e-5 && count < 50){
        if(irlsConv == false){
            count = 251;
        }
        
        count += 1;
        
        switch(family){
        case 1 :
            sto = betas(k);
            break;
        case 2 :
            sto = gammas(k);
            break;
        }
        
        Eigen::VectorXd Vdiag = irlsV(mu, zGy, theta, family);
        Eigen::VectorXd Z = irlsZ(eta, y, family, Vdiag);
        Eigen::VectorXd r(eta.size());
        
        r = Z - eta + X.col(k)*sto;
        
        Eigen::VectorXd updN(r.size());
        Eigen::VectorXd updD(r.size());
        
        double updNum = 0;
        double updDem = 0;
        for(int i = 0; i < r.size(); i += 1){
            updNum += Vdiag(i)*r(i)*X(i,k);
            updDem += Vdiag(i)*X(i,k)*X(i,k);
        }
        
        
        double upd = updNum/updDem;
        double delta = 0;       
        double result; 
        
        if(k == 0){
            result = upd;
        } 
        else {
            switch(family){
            case 1 :
                delta = (lambda*betaWeight)/(2*(std::abs(betas(k))*betaWeight+
                    std::abs(gammas(k))*gammaWeight + tau)*updDem);
                break;
            case 2 :
                delta = (lambda*gammaWeight)/(2*(std::abs(betas(k))*betaWeight+
                    std::abs(gammas(k))*gammaWeight + tau)*updDem);
                break;
            }
            if(upd > delta){
                result = upd - delta;
            }
            else if(upd <= -1*delta){
                result = upd + delta;
            }
            else {
                result = 0;
            }
        }    
        
        switch(family){
        case 1 :
            betas(k) = result;
            break;
        case 2 :
            gammas(k) = result;
            break;
        }
        
        // Update eta/mu
        eta = eta - X.col(k)*sto + X.col(k)*result;
        mu = updateMu(eta, family);
        
        diff = std::abs(sto - result);
        
    }  
    
    switch(family){
    case 1 :
        return(betas(k));
        break;
    case 2 :
        return(gammas(k));
        break;
    default :
        return(0);
    }
};

double updateCoefk_lasso(Eigen::VectorXd y, Eigen::VectorXd betas, 
                         const Eigen::VectorXd offsetx, 
                         Eigen::VectorXd gammas, const int family, const Eigen::VectorXd offsetz, 
                         const Eigen::MatrixXd X, const double theta, const double lambda, 
                         const double tau, int k, const Eigen::VectorXd zGy, 
                         const bool irlsConv, const double betaWeight, 
                         const double gammaWeight, const double countEM, const double countk){
    
    double diff = 1000;
    int count = 0;
    
    Eigen::VectorXd eta(X.rows());
    Eigen::VectorXd mu(X.rows());
    
    switch(family){
    case 1 :
        eta = calcEta(betas, X, offsetx);
        break;
    case 2 :
        eta = calcEta(gammas, X, offsetz);
        y = zGy;
        break;
    }
    
    mu = updateMu(eta, family);
    
    double sto = 0;
    while(diff > 1e-5 && count < 50){
        if(irlsConv == false){
            count = 251;
        }
        
        count += 1;
        
        switch(family){
        case 1 :
            sto = betas(k);
            break;
        case 2 :
            sto = gammas(k);
            break;
        }
        
        Eigen::VectorXd Vdiag = irlsV(mu, zGy, theta, family);
        Eigen::VectorXd Z = irlsZ(eta, y, family, Vdiag);
        Eigen::VectorXd r(eta.size());
        
        
        r = Z - eta + X.col(k)*sto;
        
        Eigen::VectorXd updN(r.size());
        Eigen::VectorXd updD(r.size());
        
        double updNum = 0;
        double updDem = 0;
        for(int i = 0; i < r.size(); i += 1){
            updNum += Vdiag(i)*r(i)*X(i,k);
            updDem += Vdiag(i)*X(i,k)*X(i,k);
        }
        
        // DEBUG BLOCK
        // if(family == 2 & k == 0){
        // std::cout << "Vdiag(15): " << Vdiag(15) << std::endl;
        // std::cout << "eta(15): " << eta(15) << std::endl;
        // std::cout << "Z(15): " << Z(15) << std::endl;
        // std::cout << "r(15): " << r(15) << std::endl;
        // std::cout << "X(15,1): " << X(15,k) << std::endl;
        // std::cout << "updNum(15): " << updNum << std::endl;
        // std::cout << "updDem(15): " << updDem << std::endl;
        // }       
        
        double upd = updNum/updDem;
        double delta = 0;       
        double result; 
        
        // DEBUG BLOCK
        // if(family == 2 & k == 0){
        // std::cout << "upd: " << upd << std::endl;
        // }       
        
        if(k == 0){
            result = upd;
            // DEBUG BLOCK
            // if(family == 2){
            std::cout << "res: " << result << std::endl;
            // }
        } 
        else {
            switch(family){
            case 1 :
                delta = (lambda*betaWeight)/(2*updDem);
                break;
            case 2 :
                delta = (lambda*gammaWeight)/(2*updDem);
                break;
            }
            if(upd > delta){
                result = upd - delta;
            }
            else if(upd < -1*delta){
                result = upd + delta;
            }
            else {
                result = 0;
            }
        }    
        
        switch(family){
        case 1 :
            betas(k) = result;
            break;
        case 2 :
            gammas(k) = result;
            break;
        }
        
        // Update eta/mu
        eta = eta - X.col(k)*sto + X.col(k)*result;
        mu = updateMu(eta, family);
        
        diff = std::abs(sto - result);
        
    }  
    
    switch(family){
    case 1 :
        return(betas(k));
        break;
    case 2 :
        // DEBUG BLOCK
        // if(k == 0){
        //     std::cout << "gammas: " << gammas(k) << std::endl;
        // }
        return(gammas(k));
        break;
    default :
        return(0);
    }
};


// [[Rcpp::export]]
Rcpp::List updateEM(Eigen::VectorXd y, Eigen::VectorXd betas, 
                         const Eigen::VectorXd offsetx, 
                         Eigen::VectorXd gammas, const Eigen::VectorXd offsetz, 
                         const Eigen::MatrixXd X, double theta, const double lambda, 
                         const double tau, const bool irlsConv, const Eigen::VectorXd betaWeights, 
                         const Eigen::VectorXd gammaWeights, const int maxIT, const int loud,
                         double eps, const int model, const int penType, 
                         const int maxIT2){
    
    double diff = 1000;
    int count = 0;
    double loglik = std::numeric_limits<double>::max();
    double loglik_sto = std::numeric_limits<double>::max();
    
    Eigen::VectorXd fnb(y.size());
    Eigen::VectorXd eta = calcEta(betas, X, offsetx);
    Eigen::VectorXd mu = updateMu(eta, 1);
    
    Eigen::VectorXd etaz = calcEta(gammas, X, offsetz);
    Eigen::VectorXd pij = updateMu(etaz, 2);
    
    
    Eigen::VectorXd zGy(y.size());
    Eigen::VectorXd weights(y.size());
    
    double nconv = 0;
    double bsto = 0;
    double gsto = 0;
    int countk = 0;
    double diffk = 1000;
    
    while(diff > eps && count < maxIT){
        count += 1;
        
        //if(loud != 0 && count % loud == 0){
        //  std::cout << "EM Iteration: " << count << std::endl;
        //}
        
        loglik_sto = loglik;
        
        fnb = dnbinom(y, mu, theta, 0);
        
        // if(count < 3){
        //     std::cout << "fnb(40): " << fnb(39) << std::endl;
        // }
        
        switch(model) {
        case 1 : 
            for(int i = 0; i < y.size(); i += 1){
                zGy(i) = (pij(i)/(pij(i) + fnb(i)*(1-pij(i))))*(y(i) == 0);
                weights(i) = 1-zGy(i);
            }
            break;         
        case 2 :
            for(int i = 0; i < y.size(); i += 1){
                zGy(i) = 0;
                weights(i) = 1;
            }
            break;
        }
        
        // if(count < 3){
        //     std::cout << "mu(13): " << mu(12) << std::endl;
        //     std::cout << "pij(44): " << pij(43) << std::endl;
        //     std::cout << "zGy(25): " << zGy(24) << std::endl;
        //     std::cout << "zGy(26): " << zGy(25) << std::endl;
        // }
        
        
        for(int k = 0; k < X.cols(); k += 1){           
            // std::cout << "param: " << k << std::endl;
            
            countk = 0;
            diffk = 1000;
            while(diffk > eps && countk < maxIT2){
                countk += 1;
                
                bsto = betas(k);
                gsto = gammas(k);
                
                try {
                    switch(model) {
                    case 1 : 
                        switch(penType) {
                        case 1 :
                            betas(k) = updateCoefk_log(y, betas, offsetx, gammas, 1, offsetz, X,
                                  theta, lambda, tau, k, zGy, irlsConv, betaWeights(k), 
                                  gammaWeights(k), count, countk);
                            
                            gammas(k) = updateCoefk_log(y, betas, offsetx, gammas, 2, offsetz, X,
                                   theta, lambda, tau, k, zGy, irlsConv, betaWeights(k), 
                                   gammaWeights(k), count, countk);
                            break;         
                        case 2 : 
                            betas(k) = updateCoefk_lasso(y, betas, offsetx, gammas, 1, offsetz, X,
                                  theta, lambda, tau, k, zGy, irlsConv, betaWeights(k), 
                                  gammaWeights(k), count, countk);
                            
                            gammas(k) = updateCoefk_lasso(y, betas, offsetx, gammas, 2, offsetz, X,
                                   theta, lambda, tau, k, zGy, irlsConv, betaWeights(k), 
                                   gammaWeights(k), count, countk);
                            break;         
                        }
                        break;
                    case 2 :
                        switch(penType) {
                        case 1 :
                            betas(k) = updateCoefk_log(y, betas, offsetx, gammas, 1, offsetz, X,
                                  theta, lambda, tau, k, zGy, irlsConv, betaWeights(k), 
                                  gammaWeights(k), count, countk);
                            break;
                        case 2 :
                            betas(k) = updateCoefk_lasso(y, betas, offsetx, gammas, 1, offsetz, X,
                                  theta, lambda, tau, k, zGy, irlsConv, betaWeights(k), 
                                  gammaWeights(k), count, countk);
                            break;
                        }
                        
                        gammas(k) = 0;
                        break;
                    }
                    
                }                
                
                catch(...){
                    nconv = 1;
                    return Rcpp::List::create(Rcpp::Named("betas") = betas,
                                              Rcpp::Named("gammas") = gammas, 
                                              Rcpp::Named("nconv") = nconv,
                                              Rcpp::Named("it") = count);
                }
                
                diffk = std::abs(bsto - betas(k))/bsto + std::abs(gsto - gammas(k))/gsto;
                // if(count < 3 && countk < 10){
                //     std::cout << "diff.parm.k: " << diffk << std::endl;
                // }
            }   
            
            // DEBUG BLOCK
            // if(count < 10){
            //     std::cout << "k: " << k << "( " << betas(k) << ", " << gammas(k) << ")" << std::endl;
            // }   
        }   
        
        // Update means
        eta = calcEta(betas, X, offsetx);
        mu = updateMu(eta,1);
        
        if (model == 1){
            etaz = calcEta(gammas, X, offsetz);
            pij = updateMu(etaz,2);           
        }
        
        // Update theta
        theta = theta_ml(theta, y, mu, weights);
        
        // Calculate penalty
        
        double pen = 0;
        double loglik_obs = 1000;
        switch(model) {
        case 1 : 
            switch(penType) {
            case 1 :
                for(int k = 1; k < betas.size(); k += 1){
                    pen += lambda*log(betaWeights(k)*std::abs(betas(k)) + 
                        gammaWeights(k)*std::abs(gammas(k)) + tau);
                }
                break;
            case 2 :
                for(int k = 1; k < betas.size(); k += 1){
                    pen += lambda*(betaWeights(k)*std::abs(betas(k)) + 
                        gammaWeights(k)*std::abs(gammas(k)));
                }
                break;
            }
            break;         
        case 2 :
            switch(penType) {
            case 1 :
                for(int k = 1; k < betas.size(); k += 1){
                    pen += lambda*log(betaWeights(k)*std::abs(betas(k)) + tau);
                }
                break;
            case 2 :
                for(int k = 1; k < betas.size(); k += 1){
                    pen += lambda*(betaWeights(k)*std::abs(betas(k)));
                }
                break;
            }
            break;         
        }
        
        // Calculate Likelihood
        loglik_obs = updLik(y, mu, pij, theta);
        
        loglik = loglik_obs - pen;
        
        // std::cout << "loglik: " << loglik << std::endl;
        diff = std::abs(loglik - loglik_sto)/loglik_sto;
        // std::cout << "diff: " << diff << std::endl;
        
    }
    
    return Rcpp::List::create(Rcpp::Named("betas") = betas,
                              Rcpp::Named("gammas") = gammas, 
                              Rcpp::Named("nconv") = nconv,
                              Rcpp::Named("it") = count);
}

// Primary Coefficient Update Function
// [[Rcpp::export]]
Rcpp::List updateEM2(Eigen::VectorXd y, Eigen::VectorXd betas, 
                     const Eigen::VectorXd offsetx, 
                     Eigen::VectorXd gammas, const Eigen::VectorXd offsetz, 
                     const Eigen::MatrixXd X, double theta, const double lambda, 
                     const double tau, const bool irlsConv, const Eigen::VectorXd betaWeights, 
                     const Eigen::VectorXd gammaWeights, const int maxIT, const int loud,
                     double eps, const int penType){
    
    double diff = 1000;
    int count = 0;
    double loglik = std::numeric_limits<double>::max();
    double loglik_sto = std::numeric_limits<double>::max();
    
    Eigen::VectorXd fnb(y.size());
    Eigen::VectorXd eta = calcEta(betas, X, offsetx);
    Eigen::VectorXd mu = updateMu(eta, 1);
    
    
    Eigen::VectorXd etaz = calcEta(gammas, X, offsetz);
    Eigen::VectorXd pij = updateMu(etaz, 2);
    Eigen::VectorXd zGy(y.size());
    Eigen::VectorXd weights(y.size());
    
    double nconv = 0;
    
    while(diff > eps && count < maxIT){
        count += 1;
        
        if(loud != 0 && count % loud == 0){
            std::cout << "EM Iteration: " << count << std::endl;
        }
        
        loglik_sto = loglik;
        
        fnb = dnbinom(y, mu, theta, 0);
        
        // std::cout << "fnb(40): " << fnb(39) << std::endl;
        
        for(int i = 0; i < y.size(); i += 1){
            zGy(i) = (pij(i)/(pij(i) + fnb(i)*(1-pij(i))))*(y(i) == 0);
            weights(i) = 1 - zGy(i);
        }
        
        // std::cout << "zGy(25): " << zGy(24) << std::endl;
        for(int k = 0; k < X.cols(); k += 1){           
            
            try {
                switch(penType) {
                case 1 :
                    betas(k) = updateCoefk_log(y, betas, offsetx, gammas, 1, offsetz, X,
                          theta, lambda, tau, k, zGy, irlsConv, betaWeights(k), 
                          gammaWeights(k), count, 1);
                    break;         
                case 2 : 
                    betas(k) = updateCoefk_lasso(y, betas, offsetx, gammas, 1, offsetz, X,
                          theta, lambda, tau, k, zGy, irlsConv, betaWeights(k), 
                          gammaWeights(k), count, 1);
                    break;         
                }
            }                
            
            catch(...){
                nconv = 1;
                return Rcpp::List::create(Rcpp::Named("betas") = betas,
                                          Rcpp::Named("gammas") = gammas, 
                                          Rcpp::Named("nconv") = nconv);
            }
            // std::cout << "countk: " << countk << std::endl;
            // std::cout << "beta: " << betas(k) << std::endl;
            // std::cout << "gamma: " << gammas(k) << std::endl;
            
        }
        
        for(int k = 0; k < X.cols(); k += 1){           
            // std::cout << "param: " << k << std::endl;
            
            try {
                switch(penType) {
                case 1 :           
                    gammas(k) = updateCoefk_log(y, betas, offsetx, gammas, 2, offsetz, X,
                           theta, lambda, tau, k, zGy, irlsConv, betaWeights(k), 
                           gammaWeights(k), count, 1);
                    break;         
                case 2 : 
                    gammas(k) = updateCoefk_lasso(y, betas, offsetx, gammas, 2, offsetz, X,
                           theta, lambda, tau, k, zGy, irlsConv, betaWeights(k), 
                           gammaWeights(k), count, 1);
                    break;         
                }
            }                
            
            catch(...){
                nconv = 1;
                return Rcpp::List::create(Rcpp::Named("betas") = betas,
                                          Rcpp::Named("gammas") = gammas, 
                                          Rcpp::Named("nconv") = nconv);
            }
            
            // std::cout << "countk: " << countk << std::endl;
            // std::cout << "beta: " << betas(k) << std::endl;
            // std::cout << "gamma: " << gammas(k) << std::endl;
            
        }
        
        // Update means
        eta = calcEta(betas, X, offsetx);
        mu = updateMu(eta,1);
        etaz = calcEta(gammas, X, offsetz);
        pij = updateMu(etaz,2);
        
        // Update theta
        theta = theta_ml(theta, y, mu, weights);
        
        // Calculate penalty
        double pen = 0;
        
        switch(penType) {
        case 1 :
            for(int k = 1; k < betas.size(); k += 1){
                pen += lambda*log(betaWeights(k)*std::abs(betas(k)) + 
                    gammaWeights(k)*std::abs(gammas(k)) + tau);
            }
            break;
        case 2 :
            for(int k = 1; k < betas.size(); k += 1){
                pen += lambda*(betaWeights(k)*std::abs(betas(k)) + 
                    gammaWeights(k)*std::abs(gammas(k)));
            }
            break;
        }
        
        // Calculate Likelihood
        double loglik_obs = updLik(y, mu, pij, theta);
        
        loglik = loglik_obs - pen;
        // std::cout << "loglik: " << loglik << std::endl;
        diff = std::abs(loglik - loglik_sto)/loglik_sto;
        // std::cout << "diff: " << diff << std::endl;
        
    }
    
    
    return Rcpp::List::create(Rcpp::Named("betas") = betas,
                              Rcpp::Named("gammas") = gammas, 
                              Rcpp::Named("nconv") = nconv);
}
