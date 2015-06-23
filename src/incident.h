#ifndef _cda_INCIDENT_H
#define _cda_INCIDENT_H
#include <RcppArmadillo.h>


arma::cx_mat incident_field(const arma::cx_colvec& E0, 
			    const arma::colvec& k, 
			    const arma::mat& R,
			    const arma::mat& Angles);

arma::cx_mat multiple_incident_field(const arma::cx_colvec& E0, 
				     const arma::colvec& k, 
				     const arma::mat& R,
				     const arma::ivec& Axes,
				     const arma::colvec& Angles);

#endif
