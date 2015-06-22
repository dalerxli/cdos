#ifndef _cda_CD_H
#define _cda_CD_H

#include <RcppArmadillo.h>

arma::colvec averaging(const arma::mat& R, 
		       const arma::cx_mat& DiagBlocks, 		
		       const double kn, 
		       const arma::mat& Angles,
		       const arma::colvec& Weights,
		       const int Niter, const double tol);

#endif
