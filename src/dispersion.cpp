// 
// Multiple angles of incidence
// 
#include "utils.h"
#include "cda.h"
#include "dispersion.h"
#include "cg.h"
#include "incident.h"
#include "cross_sections.h"

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;
using namespace arma;


//
// Angle-resolved cross-sections for multiple directions of incidence
//
// R is the Nx3 matrix of positions
// DiagBlocks is the 3Nx3 block-matrix of polarizabilities
// kn is the incident wavenumber (scalar)
// Angles is the Nangles vector of incident beam angles
// Axes is the Nangles vector of incident beam axes
// polarisation is an integer flag to switch between linear and circular polarisation
// cg logical flag to use conjugate gradient solver
// born logical flag, use first Born approx for cg solver
// Niter is the max number of cg iterations
// tol is the cg tolerance
 arma::mat dispersion(const arma::mat& R, 
		      const arma::cx_mat& DiagBlocks, 		      
		      const double kn, 
		      const arma::vec& Angles, 
		      const arma::ivec& Axes, 
		      const int polarisation,
		      const int Niter,
		      const double tol)
   {
     const int N = R.n_rows, NAngles = Angles.n_elem;
     arma::mat Rot(3,3);

    // incident field
    const arma::colvec  khat="0 0 1;"; 
    const arma::colvec kvec = kn*khat; 
    arma::cx_colvec LPP, LPS;
    if(polarisation == 0){ // linear
      LPP="(1,0) (0,0) (0,0);", LPS="(0,0) (1,0) (0,0);";
    } else { // circular
      LPP="(1,0) (0,1) (0,0);", LPS="(0,1) (1,0) (0,0);";
      LPP = arma::datum::sqrt2/2 * LPP ;
      LPS = arma::datum::sqrt2/2 * LPS ;
    }

    arma::mat res(NAngles, 6) ;  
    arma::cx_mat E0(3*N,NAngles), E=E0, P=E0;

    // tmp variables
    int iter=0;
    double rel_error=1e10;
    arma::colvec xsec = zeros<arma::colvec>(NAngles), tmp=xsec; 

    // -----------------------
    // first polarisation
    // -----------------------
    iter=0;
    rel_error=1e10;
    E0 = multiple_incident_field(LPP, kvec, R, Axes, Angles);
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
    
    res.col(0) = xsec; 
    res.col(1) = absorption(kn, P, E); 
    res.col(2) = res.col(0) - res.col(1);

   // -----------------------
    // second polarisation
    // -----------------------
    iter=0;
    rel_error=1e10;
    xsec = 0*tmp;
    E0 = multiple_incident_field(LPS, kvec, R, Axes, Angles);
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
    
    res.col(3) = xsec; 
    res.col(4) = absorption(kn, P, E); 
    res.col(5) = res.col(3) - res.col(4); 
             
    return res ;
   } 


//
// Angle-resolved spectra for linear or circular polarisations
//
// kn is the vector of incident wavenumbers
// Alpha is the 3N vector of principal polarisabilities
// R is the Nx3 matrix of positions
// Euler is the Nx3 matrix of particle rotation angles
// Angles is the Nangles vector of incident beam angles
// Axes is the Nangles vector of incident beam axes
// polarisation is an integer flag to switch between linear and circular polarisation
// cg logical flag to use conjugate gradient solver
// born logical flag, use first Born approx for cg solver
// nmax is the max number of cg iterations
// tol is the cg tolerance
// progress is a logical flag to display progress bars
arma::cube dispersion_spectrum(const arma::colvec kn, 
			       const arma::mat& R, 
			       const arma::cx_mat& Alpha, 
			       const arma::mat& Euler, 
			       const arma::vec& Angles,
			       const arma::ivec& Axes,		      
			       const int polarisation, 
			       const int Niter, const double tol,
			       const bool progress)
  {

    const int Nangles = Angles.n_elem;
    int N = kn.n_elem, Nr = R.n_rows, ll;
    arma::cube results(Nangles, 6, N);
    arma::mat tmp(Nangles, 6);
    arma::cx_mat DiagBlocks(3*Nr,3);

    for(ll=0; ll<N; ll++){ // loop over kn   
      if(progress)
	progress_bar(ll+1,N);

      DiagBlocks = diagonal_blocks(Alpha.col(ll), Euler);
      tmp = dispersion(R, DiagBlocks, kn(ll), Angles, Axes, 
		       polarisation, Niter, tol);
      results.slice(ll) = tmp; 
    }
    if(progress)
      Rcpp::Rcout << "\n";

    return results ;
  } 


RCPP_MODULE(dispersion){

       Rcpp::function( "dispersion", &dispersion,
	    "Angle-resolved cross-sections for multiple directions of incidence" ) ;
       Rcpp::function( "dispersion", &dispersion,
	    "Angle-resolved cross-sections for multiple directions of incidence" ) ;
       Rcpp::function( "dispersion_spectrum", &dispersion_spectrum,
	   "Angle-resolved spectra for linear or circular polarisations" ) ;
}
