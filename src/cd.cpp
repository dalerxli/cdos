#include "utils.h"
#include "cda.h"
#include "cd.h"
#include "cg.h"
#include "incident.h"
#include "cross_sections.h"

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace arma;
using namespace std;


//
// Performs full angular averaging for both circular polarisations
//
// R is the Nx3 matrix of positions
// DiagBlocks is the 3Nx3 block-matrix of polarizabilities
// kn is the incident wavenumber (scalar)
// Angles is the NAnglesx3 matrix of incident beam angles
// Weigths is the NAngles vector of quadrature weights
// Niter is the max number of  iterations
// tol is the tolerance
arma::colvec averaging(const arma::mat& R, 
		       const arma::cx_mat& DiagBlocks, 		
		       const double kn, 
		       const arma::mat& Angles,
		       const arma::colvec& Weights,
		       const int Niter, const double tol)
  {

    const int N = R.n_rows, NAngles = Angles.n_rows;
    arma::colvec res(4) ;   

    // incident field
    arma::cx_mat E0(3*N, NAngles), E=E0, P=E0;
    const arma::colvec  khat="1 0 0;", kvec = kn*khat;
    arma::cx_colvec RCP="(0,0) (0,1) (1,0);", LCP="(0,0) (1,0) (0,1);";
    RCP = arma::datum::sqrt2/2 * RCP ;
    LCP = arma::datum::sqrt2/2 * LCP ;

    // tmp variables
    int iter=0;
    double rel_error=1e10;
    arma::colvec xsec = zeros<arma::colvec>(NAngles), tmp=xsec; 

    // -----------------------
    // left polarisation
    // -----------------------
    iter=0;
    rel_error=1e10;
    E0 = incident_field(LCP, kvec, R, Angles);
    E=E0;
    P = polarization(E0, DiagBlocks);
    tmp = extinction(kn, P, E0);
    // order-of-scattering iterations
    while((iter < Niter) && (rel_error > tol)){
      xsec = iterate_field(R, kn, E0, DiagBlocks, E, P);
      // Note E and P have been updated
      rel_error = max(abs((xsec - tmp) / (xsec + tmp)));
      tmp = xsec;
      iter++;
    }
    
    res(0) = dot(xsec, Weights); 
    xsec = absorption(kn, P, E);
    res(2) = dot(xsec, Weights); 
      
    // -----------------------
    // right polarisation
    // -----------------------
    iter=0;
    rel_error=1e10;
    xsec = 0*tmp;
    E0 = incident_field(RCP, kvec, R, Angles);
    E=E0;
    P = polarization(E0, DiagBlocks);
    tmp = extinction(kn, P, E0);
    // order-of-scattering iterations
    while((iter < Niter) && (rel_error > tol)){
      xsec = iterate_field(R, kn, E0, DiagBlocks, E, P);
      // Note E and P have been updated
      rel_error = max(abs((xsec - tmp) / (xsec + tmp)));
      tmp = xsec;
      iter++;
    }
    res(1) = dot(xsec, Weights); 
    xsec = absorption(kn, P, E);
    res(3) = dot(xsec, Weights); 

    return res ;
  } 





//
// Angular-average spectra for LCP and RCP polarisations
//
// kn is the vector of incident wavenumbers
// R is the Nx3 matrix of positions
// Alpha is the 3N vector of principal polarisabilities
// Euler is the Nx3 matrix of particle rotation angles
// Angles is the NAnglesx3 matrix of incident beam angles
// Weigths is the NAngles vector of quadrature weights
// full is a logical flag to switch off retardation terms
// Niter is the max number of cg iterations
// tol is the cg tolerance
// progress is a logical flag to display progress bars
arma::mat average_spectrum(const arma::colvec kn, 
			   const arma::mat& R,		
			   const arma::cx_mat& Alpha, 		
			   const arma::mat& Euler, 
			   const arma::mat& Angles, 
			   const arma::colvec& Weights, 
			   const bool full, 
			   const int Niter, 
			   const double tol,
			   const bool progress)
  {

    int N = kn.n_elem, Nr = R.n_rows, ll;

    arma::mat res(N,6);
    arma::colvec tmp(4);
    arma::cx_mat DiagBlocks(3*Nr,3);

    for(ll=0; ll<N; ll++){ // loop over kn   
      if(progress)
	progress_bar(ll+1,N);

      DiagBlocks = diagonal_blocks(Alpha.col(ll), Euler);
      tmp = averaging(R, DiagBlocks, kn(ll), Angles, Weights, Niter, tol);

      res(ll,0) = 0.5*(tmp(0) + tmp(1)); // extinction 
      res(ll,1) = 0.5*(tmp(2) + tmp(3)); // absorption
      res(ll,2) = res(ll,0) - res(ll,1); // scattering
      res(ll,3) = tmp(0) - tmp(1); // cd ext L - R
      res(ll,4) = tmp(2) - tmp(3); // cd abs L - R
      res(ll,5) = res(ll,3) - res(ll,4); // cd sca L - R

    }
    if(progress)
      Rcpp::Rcout << "\n";

    return res ;
  } 



RCPP_MODULE(cd){
       Rcpp::function( "averaging", &averaging, 
	   "Performs full angular averaging for both circular polarisations" ) ;
       Rcpp::function( "average_spectrum", &average_spectrum, 
	   "Angular-average spectra for LCP and RCP polarisations" ) ;
}






