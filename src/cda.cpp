// 
// main functions for coupled dipole equations
// 
#include "utils.h"
#include "cda.h"
#include "cg.h"

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;

//
// Calculate the incident field at each dipole location 
//
// for multiple euler angles of incidence
// E0 is the normalised incident electric field before rotation
// k is the incident wavevector before rotation
// R is the Nx3 matrix of positions
// Angles is the Nanglesx3 matrix of incident beam angles
// returns a 3NxNangles complex matrix of incident fields
// where each column corresponds to a specific incident angle
arma::cx_mat incident_field(const arma::cx_colvec& E0, 
			    const arma::colvec& k, 
			    const arma::mat& R,
			    const arma::mat& Angles)
{
  const int Nangles = Angles.n_rows;
  const int N = R.n_rows;
  const arma::cx_double i = arma::cx_double(0,1);
  arma::mat Rot(3,3);
  arma::cx_mat Ei = arma::cx_mat(3*N,Nangles);
  arma::cx_colvec E0_r(3);
  arma::cx_colvec expikr(N);
  arma::colvec k_r(3);
  arma::colvec kR(N);
  arma::cx_colvec expikrrep(3*N);
  arma::cx_colvec E0rep(3*N);
  int jj=0;
  for(jj=0; jj<Nangles; jj++)
    {
      Rot = euler(Angles(jj,0), Angles(jj,1), Angles(jj,2));
      k_r = Rot.st() * k;
      E0_r = Rot.st() * E0;
      kR = R * k_r ;
      expikr = exp(i * kR);
      expikrrep = strans(vectorise(repmat(expikr, 1, 3), 1));
      E0rep = repmat(E0_r, N, 1);
      Ei.col(jj) = E0rep % expikrrep;
    }
  return(Ei);
}

//
// Incident field along multiple axes
//
// E0 is the normalised electric field
// k is the wavevector
// R is the Nx3 matrix of positions
// Axes is a vector of integer codes corresponding to x, y, z
// Angles is a vector of rotation angles around Axes
arma::cx_mat multiple_incident_field(const arma::cx_colvec& E0, 
			     const arma::colvec& k, 
			     const arma::mat& R,
			     const arma::ivec& Axes,
			     const arma::colvec& Angles)
{
  const int Nangles = Angles.n_elem;
  const int N = R.n_rows;
  const arma::cx_double i = arma::cx_double(0,1);
  arma::mat Rot(3,3);
  arma::cx_mat Ei = arma::cx_mat(3*N,Nangles);
  arma::cx_colvec E0_r(3);
  arma::cx_colvec expikr(N);
  arma::colvec k_r(3);
  arma::colvec kR(N);
  arma::cx_colvec expikrrep(3*N);
  arma::cx_colvec E0rep(3*N);
  int jj=0;
  for(jj=0; jj<Nangles; jj++)
    {
      Rot = axis_rotation(Angles(jj), Axes(jj));
      k_r = Rot.st() * k;
      E0_r = Rot.st() * E0;
      kR = R * k_r ;
      expikr = exp(i * kR);
      expikrrep = strans(vectorise(repmat(expikr, 1, 3), 1));
      E0rep = repmat(E0_r, N, 1);
      Ei.col(jj) = E0rep % expikrrep;
    }
  return(Ei);
}


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
  const int N = R.n_rows;
  const int NAngles = E0.n_cols;

  // tmp variables
  int iter=0;
  double rel_error=1e10;
  arma::colvec cext(NAngles), tmp=cext; // store extinction for convergence test
  // starting with the first Born approximation
  arma::cx_mat E = E0, P = E0; 
  P = polarization(E0, DiagBlocks);
  tmp = extinction(kn, P, E0);

  while((iter < Niter) && (rel_error > tol)){
    //std::cout << P << std::endl;
    cext = iterate_field(R, kn, E0, DiagBlocks, E, P);
    //std::cout << cext << std::endl;
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
  const arma::cx_double i = arma::cx_double(1,0);
  arma::cx_mat Gjk = arma::cx_mat(3,3);
  double rjk;
  arma::mat rk_to_rj = arma::mat(1,3);
  arma::mat rjkhat   = arma::mat(1,3);
  arma::mat rjkrjk   = arma::mat(3,3);  
  int ll=0, jj=0, kk=0;
  arma::colvec cext(NAngles); // extinction for convergence test

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
	  Gjk = exp(i*kn*rjk) / rjk *  (kn*kn*(rjkrjk - I3) +
					(i*kn*rjk - arma::cx_double(1,0)) /     
					(rjk*rjk) * (3*rjkrjk - I3)) ;
	  
	  // update E = Einc + GP
	  // where P is from a previous iteration
	  for(ll=0; ll<NAngles; ll++){
	    // field of dipole kk evaluated at jj
	    E.submat(jj*3, ll, jj*3+2, ll) -= Gjk * P.submat(kk*3, ll, kk*3+2, ll);
	    // field of dipole jj evaluated at kk
	    E.submat(kk*3, ll, kk*3+2, ll) -= Gjk.st() * P.submat(jj*3, ll, jj*3+2, ll);
	  }
	}
    }
  // update the polarization
  P = polarization(E, DiagBlocks);
  cext = extinction(kn, P, E0);

  return(cext);
}

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
