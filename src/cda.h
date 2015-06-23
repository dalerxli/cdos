#ifndef _cda_CDA_H
#define _cda_CDA_H

#include <RcppArmadillo.h>


arma::colvec convergence(const arma::mat& R, 
			 const double kn,
			 const arma::cx_mat& E0,
			 const arma::cx_mat& DiagBlocks,
			 const int Niter, const double tol);

arma::colvec iterate_field(const arma::mat& R, 
			   const double kn,
			   const arma::cx_mat& E0,
			   const arma::cx_mat& DiagBlocks, 
			   arma::cx_mat& E, 
			   arma::cx_mat& P);

arma::cx_mat diagonal_blocks(const arma::cx_colvec& Alpha,
			     const arma::mat& Euler);

arma::cx_mat polarization(const arma::cx_mat& E,
			  const arma::cx_mat& DiagBlocks);


#endif
