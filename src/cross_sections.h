#ifndef _cda_CROSSSECTIONS_H
#define _cda_CROSSSECTIONS_H

#include <RcppArmadillo.h>

arma::colvec extinction(const double kn, const arma::cx_mat& P, 
				 const arma::cx_mat& Eincident);

arma::colvec absorption(const double kn, const arma::cx_mat& P, 
				  const arma::cx_mat& Eloc);


#endif
