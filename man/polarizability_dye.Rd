% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/polarizability.R
\name{polarizability_dye}
\alias{polarizability_dye}
\title{polarizability_dye}
\usage{
polarizability_dye(wavelength, a = 1, b = 1, c = 1, medium = 1.33, ...)
}
\arguments{
\item{wavelength}{in nm}

\item{a}{orientation factor, see details}

\item{b}{orientation factor, see details}

\item{c}{orientation factor, see details}

\item{medium}{surrounding medium}

\item{...}{passed to lorentzian function}
}
\value{
matrix of polarizability
}
\description{
principal polarizability components a dye molecule
}
\details{
uses a sum of lorentzians

The polarizability tensor is distributed along the three diagonal elements according to the factors a, b, c, normalised to unity. An isotropic tensor (a=b=c) will be diag(alpha/3), while a uni-axial molecule along x may be diag(c(alpha, 0, 0)).
}
\author{
baptiste Auguie
}
\seealso{
Other user_level polarizability: \code{\link{Kuwata.A}},
  \code{\link{Kuwata.B}}, \code{\link{alpha_kuwata}};
  \code{\link{La}}; \code{\link{polarizability_ellipsoid}};
  \code{\link{principal_polarizability}}
}

