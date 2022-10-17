#ifndef GRIDOT_TRANSPORT_COST_H_INCLUDED
#define GRIDOT_TRANSPORT_COST_H_INCLUDED

#include<RcppArmadillo.h>

double transportCost(const arma::vec& x1, const arma::vec& x2, const arma::mat& mu, 
                     const arma::vec& y1, const arma::vec& y2, const arma::mat& nu, 
                     const double p1, const double p2, const arma::mat& xi,
                     const double threshold = 1e-15);

#endif
