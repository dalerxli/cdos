// 
// incident field
// 
#include "incident.h"
#include "utils.h"

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
