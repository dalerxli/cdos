// 
// main functions for coupled dipole equations
// 
#include "utils.h"
#include "cda.h"
#include "incident.h"
#include "cross_sections.h"

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;


//
// Block-matrix of polarisabilities
//
// Alpha is the 3N vector of principal polarizabilities
// Euler is the Nx3 matrix euler angles to apply
// returns a 3Nx3 complex matrix of polarizabilities
arma::cx_mat diagonal_blocks(const arma::cx_colvec& Alpha,
			     const arma::mat& Euler) {

  const int N = Euler.n_rows;
  arma::mat Rot(3,3);
  arma::cx_mat DiagBlocks(3*N,3);
  int ii=0;

  // loop over N particles
  for(ii=0; ii<N; ii++){

    Rot = euler(Euler(ii,0), Euler(ii,1), Euler(ii,2));
    DiagBlocks.submat(ii*3, 0, ii*3+2, 2) =  Rot.st() * 
      diagmat(Alpha.subvec(ii*3, ii*3+2)) * Rot;
    
  } 
  
  return DiagBlocks;
}

//
// Polarization from local fields and polarisabilities
//
// E is the 3NxNangles local field
// DiagBlocks is the 3Nx3 block-matrix of polarizabilities
// returns a 3NxNangles complex matrix of polarization
arma::cx_mat polarization(const arma::cx_mat& E,
			  const arma::cx_mat& DiagBlocks) {
  
  const int N = DiagBlocks.n_rows / 3, Nangles=E.n_cols;
  arma::cx_mat P = E, alphatmp(3,3);
  int ii=0, jj=0;

  // loop over N particles
  for(ii=0; ii<N; ii++){
    alphatmp =  DiagBlocks.submat(ii*3, 0, ii*3+2, 2);
    // loop over incident angles (columns of E)
    for(jj=0; jj<Nangles; jj++){
      P.submat(ii*3, jj, ii*3+2, jj) = alphatmp * E.submat(ii*3, jj, ii*3+2, jj);
    }
  } 
  
  return P;
}


arma::colvec convergence(const arma::mat& R, 
			 const double kn,
			 const arma::cx_mat& E0,
			 const arma::cx_mat& DiagBlocks,
			 const int Niter, const double tol){
  const int NAngles = E0.n_cols;

  // tmp variables
  int iter=0;
  double rel_error=1e10;
  arma::colvec cext(NAngles), tmp=cext;  
  // starting with the first Born approximation
  arma::cx_mat E = E0, P = E0; 
  P = polarization(E0, DiagBlocks);
  tmp = extinction(kn, P, E0);

  while((iter < Niter) && (rel_error > tol)){
    cext = iterate_field(R, kn, E0, DiagBlocks, E, P);
    // Note E and P have been updated
    rel_error = max(abs((cext - tmp) / (cext + tmp)));
    tmp = cext;
    iter++;
  }
  return(cext);
}

// iterate the local field and polarisation
//
// R is the Nx3 matrix of positions
// kn: scalar wavenumber 
// E0 is the 3NxNangles incident field
// DiagBlocks is the 3Nx3 block-matrix of polarizabilities
// E is the 3NxNangles local field
// P: complex 3NxNangles matrix of polarization
// return current extinction cross-sections (for convergence check)
// side-effect: update P and E
arma::colvec iterate_field(const arma::mat& R, 
			   const double kn,
			   const arma::cx_mat& E0,
			   const arma::cx_mat& DiagBlocks, 
			   arma::cx_mat& E, 
			   arma::cx_mat& P){

  const int N = R.n_rows;
  const int NAngles = E0.n_cols;

  // tmp variables
  const arma::cx_mat I3 = arma::eye<arma::cx_mat>( 3, 3 );
  const arma::cx_double i = arma::cx_double(0,1);
  arma::cx_mat Gjk = arma::cx_mat(3,3);
  double rjk;
  arma::mat rk_to_rj = arma::mat(1,3);
  arma::mat rjkhat   = arma::mat(1,3);
  arma::mat rjkrjk   = arma::mat(3,3);  
  int ll=0, jj=0, kk=0;
  arma::colvec cext(NAngles); // extinction for convergence test
  arma::cx_vec Edip(3);
  E = E0; // incident field + 
          // all pairwise dipole fields
  for(jj=0; jj<N; jj++)
    {
      for(kk=jj+1; kk<N; kk++)
	{
	  rk_to_rj = R.row(jj) - R.row(kk) ;
	  rjk = norm(rk_to_rj, 2);
	  rjkhat = rk_to_rj / rjk;
	  rjkrjk =  rjkhat.st() * rjkhat;
	  // 3x3 propagator
	  Gjk = exp(i*kn*rjk) / rjk *  (kn*kn*(I3 - rjkrjk) +
					(i*kn*rjk - arma::cx_double(1,0)) /     
					(rjk*rjk) * (I3 - 3*rjkrjk)) ;
	  
	  // update E = Einc + GP
	  // where P is from a previous iteration
	  for(ll=0; ll<NAngles; ll++){
	    // field of dipole kk evaluated at jj
	    Edip = Gjk * P.submat(kk*3, ll, kk*3+2, ll);
	    E.submat(jj*3, ll, jj*3+2, ll) += Edip;
	    // field of dipole jj evaluated at kk
	    Edip = Gjk.st() * P.submat(jj*3, ll, jj*3+2, ll);
	    E.submat(kk*3, ll, kk*3+2, ll) += Edip;
	  }
	}
    }
  // update the polarization
  P = polarization(E, DiagBlocks);
  cext = extinction(kn, P, E0);

  return(cext);
}



RCPP_MODULE(cda){

  // user-level
  Rcpp::function( "extinction", &extinction, 
		  "Calculate the extinction cross-section for multiple incident angles" ) ;
  Rcpp::function( "absorption", &absorption, 
		  "Calculate the absorption cross-section for multiple incident angles" ) ;

  // utils
  Rcpp::function( "euler", &euler, 
		  "Euler rotation matrix" ) ;
  Rcpp::function( "axis_rotation", &axis_rotation, 
		  "Rotation matrix about a cartesian axis" ) ;

  // low-level
  Rcpp::function( "convergence", &convergence, 
		  "Iterating order-of-scattering solution" ) ;
  Rcpp::function( "iterate_field", &iterate_field, 
		  "Order-of-scattering solution for the fields" ) ;
  Rcpp::function( "diagonal_blocks", &diagonal_blocks, 
		  "3x3 polarizability blocks" ) ;
  Rcpp::function( "polarization", &polarization, 
		  "Compute polarisation from solved internal field" ) ;
  Rcpp::function( "incident_field", &incident_field, 
		  "Calculate the incident field at each dipole location" ) ;
  Rcpp::function( "multiple_incident_field", &multiple_incident_field,
		  "Incident field along multiple axes" ) ;

}
