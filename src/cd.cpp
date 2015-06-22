#include "utils.h"
#include "cda.h"
#include "cd.h"
#include "cg.h"

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
// A is the 3Nx3N interaction matrix
// DiagBlocks is the 3Nx3 block-matrix of polarizabilities
// kn is the incident wavenumber (scalar)
// Angles is the Nanglesx3 matrix of incident beam angles
// Weigths is the Nangles vector of quadrature weights
// cg logical flag to use conjugate gradient solver
// born logical flag, use first Born approx for cg solver
// nmax is the max number of cg iterations
// tol is the cg tolerance
arma::colvec averaging(const arma::mat& R, const arma::cx_mat& A, 
		       const arma::cx_mat& DiagBlocks, 		
		       const double kn, 
		       const arma::mat& Angles,
		       const arma::colvec& Weights,
		       const bool cg, const bool born, 
		       const int nmax, const double tol)
  {

    const int N = R.n_rows, Nangles = Angles.n_rows;
    arma::colvec res(4) ;   

    // incident field
    arma::cx_mat Eincident(3*N, Nangles), Eloc(3*N, Nangles);
    arma::cx_mat P(3*N, Nangles);
    arma::cx_mat guess = zeros<arma::cx_mat>(3*N,Nangles);
    const arma::colvec  khat="1 0 0;", kvec = kn*khat;
    arma::cx_colvec RCP="(0,0) (0,1) (1,0);", LCP="(0,0) (1,0) (0,1);";
    RCP = arma::datum::sqrt2/2 * RCP ;
    LCP = arma::datum::sqrt2/2 * LCP ;
    arma::colvec xsec(Nangles); // temporary storage of cross-sections
    // left polarisation
    Eincident = incident_field(LCP, kvec, R, Angles);
    if(cg) {     
      if(born){ 
	// first Born approximation
	guess = Eincident; 
      } 
      Eloc = cg_solve(A, Eincident, guess, nmax, tol);
    } else {
      Eloc = solve(A, Eincident);
    }
    P = polarization(Eloc, DiagBlocks);
    // P = Adiag * Eloc;
    xsec = extinction(kn, P, Eincident);
    res(0) = dot(xsec, Weights); 
    xsec = absorption(kn, P, Eloc);
    res(2) = dot(xsec, Weights); 
      
    // right polarisation
    Eincident = incident_field(RCP, kvec, R, Angles);
    if(cg) {     
      if(born){ 
	// first Born approximation
	guess = Eincident; 
      } else {
	guess = 0 * Eincident;
      }
      Eloc = cg_solve(A, Eincident, guess, nmax, tol);
    } else {
      Eloc = solve(A, Eincident);
    }
    P = polarization(Eloc, DiagBlocks);
    xsec = extinction(kn, P, Eincident);
    res(1) = dot(xsec, Weights); 
    xsec = absorption(kn, P, Eloc);
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
// Angles is the Nanglesx3 matrix of incident beam angles
// Weigths is the Nangles vector of quadrature weights
// full is a logical flag to switch off retardation terms
// cg logical flag to use conjugate gradient solver
// born logical flag, use first Born approx for cg solver
// nmax is the max number of cg iterations
// tol is the cg tolerance
// progress is a logical flag to display progress bars
arma::mat average_spectrum(const arma::colvec kn, 
			   const arma::mat& R,		
			   const arma::cx_mat& Alpha, 		
			   const arma::mat& Euler, 
			   const arma::mat& Angles, 
			   const arma::colvec& Weights, 
			   const bool full, 
			   const bool cg, const bool born, 
			   const int nmax, 
			   const double tol,
			   const bool progress)
  {

    int N = kn.n_elem, Nr = R.n_rows, ll;

    arma::mat res(N,6);
    arma::colvec tmp(4);
    arma::cx_mat A(3*Nr,3*Nr), DiagBlocks(3*Nr,3);

    for(ll=0; ll<N; ll++){ // loop over kn   
      if(progress)
	progress_bar(ll+1,N);

      DiagBlocks = diagonal_blocks(Alpha.col(ll), Euler);
      //E = iterate_field(E0, DiagBlocks);
      A = interaction_matrix(R, kn(ll), DiagBlocks, full);
      tmp = averaging(R, A, DiagBlocks, kn(ll), Angles, Weights, cg, born, nmax, tol);

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






