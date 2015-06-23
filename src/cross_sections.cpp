// 
// cross sections
// 
#include "cross_sections.h"

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;

//
// Calculate the extinction cross-section for multiple incident angles
//
// kn: scalar wavenumber 
// P: complex 3NxNangles matrix of polarization
// Eincident: complex 3NxNangles matrix of incident fields
arma::colvec extinction(const double kn, const arma::cx_mat& P, 
			 const arma::cx_mat& Eincident)
{
  int Nangles = P.n_cols, ii=0;
  arma::colvec results(Nangles);

  for (ii=0; ii<Nangles; ii++)
    {
      results(ii) = imag(cdot(Eincident.col(ii), P.col(ii)));
    }
  return  4*arma::datum::pi*kn*results;
}


//
// Calculate the absorption cross-section for multiple incident angles
//
// kn: scalar wavenumber 
// P: complex 3NxNangles matrix of polarization
// Eloc: complex 3NxNangles matrix of local fields
 arma::colvec absorption(const double kn, const arma::cx_mat& P, 
		  const arma::cx_mat& Eloc)
{
  int Nangles = P.n_cols, ii=0;
  arma::colvec results(Nangles);

  for (ii=0; ii<Nangles; ii++)
    {
      results(ii) = imag(cdot(Eloc.col(ii), P.col(ii))) -	\
			      kn*kn*kn* 2/3 * real(cdot(P.col(ii), P.col(ii)));
    }
  return  4*arma::datum::pi*kn*results;

}

