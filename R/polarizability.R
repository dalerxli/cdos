alpha_bare <- function(wavelength=seq(300,800), alpha_inf=5.7128e-39,  
                       alpha_k=c(4e-39, 1.6e-39), 
                       lambda_k=c(530, 500), 
                       mu_k=c(10.6e3,6.02409e3)){
  
  lorentz <- function(alpha, lambda, mu) {
    alpha * (1 - 1 / (1 - lambda^2/wavelength^2 - 1i*lambda^2/(wavelength*mu)))
  }
  
  sums <- mapply(lorentz, alpha = alpha_k, lambda = lambda_k, mu = mu_k, SIMPLIFY=TRUE)
  alpha <- alpha_inf - rowSums(sums)
  
  alpha
}


##' principal polarizability components a dye molecule
##'
##' uses a sum of lorentzians
##' @title polarizability_dye
##' @param wavelength in nm
##' @param a orientation factor, see details
##' @param b orientation factor, see details
##' @param c orientation factor, see details
##' @param medium surrounding medium
##' @param ... passed to lorentzian function
##' @return matrix of polarizability
##' @export
##' @family user_level polarizability
##' @author baptiste Auguie
##' @details 
##' The polarizability tensor is distributed along the three diagonal elements according to the factors a, b, c, normalised to unity. An isotropic tensor (a=b=c) will be diag(alpha/3), while a uni-axial molecule along x may be diag(c(alpha, 0, 0)).
polarizability_dye <- function(wavelength, a=1, b=1, c=1, medium=1.33, ...) 
{
  
  alpha <- alpha_bare(wavelength, ...)
  ## alpha <- clausius_mossotti() #?? should we? Eloc effect. don't think so
  abc <- a+b+c
  axx <- a/abc * alpha
  ayy <- b/abc * alpha
  azz <- c/abc * alpha
  
  cbind(axx, ayy, azz)*1e27/8.854e-12
  
}

## 
## Functions to calculate the dipole polarizability
## 

##' principal polarizability components for an ellipsoidal particle
##'
##' uses the Kuwata prescription (see references)
##' @title polarizability_ellipsoid
##' @param wavelength wavelength in nm
##' @param epsilon complex permittivity
##' @param a semi-axis in nm
##' @param b semi-axis in nm
##' @param c semi-axis in nm
##' @param medium surrounding medium
##' @param kuwata logical, use Kuwata or Clausius Mossotti prescription, see Details
##' @return matrix of polarizability
##' @export
##' @family user_level polarizability
##' @author baptiste Auguie
##' @references
##' Kuwata et al. Resonant light scattering from metal nanoparticles: Practical analysis beyond Rayleigh approximation Appl. Phys. Lett. 83, 22 (2003)
##' @details 
##' The Kuwata version includes semi-empirical terms of radiative correction and dynamic depolarisation to better match the fully retarded dipolar response in a reasonable range of (subwavelength) sizes and aspect ratios.
polarizability_ellipsoid <- function(wavelength, epsilon, a=50, b=30, c=b, 
                                     medium = 1.33, kuwata= TRUE) 
{
  if(kuwata){
    V <- 4 * pi/3 * a * b * c
    chi.a <- La(a, b, c)
    chi.b <- La(b, a, c)
    chi.c <- La(c, a, b)

    aa <- alpha_kuwata(wavelength, epsilon, V, a, chi.a, 
                                   medium)
    ab <- alpha_kuwata(wavelength, epsilon, V, b, chi.b, 
                                   medium)
    ac <- alpha_kuwata(wavelength, epsilon, V, c, chi.c, 
                                   medium)
    return(cbind(aa, ab, ac))
  } else {
    cm <- a^3*(epsilon - medium^2) / (epsilon +2*medium^2)
    return(cbind(cm,cm,cm))
  }
}

##' Shape factor for an ellipsoid
##'
##' calculates the shape factor for a general ellipsoid
##' @title La
##' @param a semi-axis in nm
##' @param b semi-axis in nm
##' @param c semi-axis in nm
##' @return shape factor along a
##' @author baptiste Auguie
##' @export
##' @family user_level polarizability
La <- function (a = 50, b = a, c = a) 
{
  ## scaled version to help convergence
  V <- a*b*c
  b <- b/a
  c <- c/a
    integrand <- function(q) {
      fq <- (q + 1) * (q + b^2) * (q + c^2)
      1/((1 + q) * sqrt(fq))
    }
  I1 <- integrate(integrand, lower = 0, upper = Inf)$value
  V/2 * I1 / a^3
  
}
##' polarizability
##'
##' prescription from Kuwata
##' @title alpha_kuwata
##' @aliases alpha_kuwata Kuwata.A Kuwata.B
##' @param wavelength wavelength
##' @param epsilon permittivity
##' @param V volume
##' @param axis semi-axis along incident field
##' @param L shape factor
##' @param medium refractive index
##' @return polarizability
##' @export
##' @family user_level polarizability
##' @author baptiste Auguie
##' @references
##' Kuwata et al. Resonant light scattering from metal nanoparticles: Practical analysis beyond Rayleigh approximation Appl. Phys. Lett. 83, 22 (2003)
alpha_kuwata <-
function (wavelength, epsilon, V, axis, L, medium = 1.33) 
{
    A <- Kuwata.A(L)
    B <- Kuwata.B(L)
    x0 <- 2 * pi * axis/wavelength
    epsilon.medium <- medium^2
    denom <- (L + epsilon.medium/(epsilon - epsilon.medium)) + A * epsilon.medium * 
      x0^2 + B * epsilon.medium^2 * x0^4 - (0+1i)/3 * 4 * pi^2 * epsilon.medium^(3/2) * 
        V/wavelength^3
    V/denom/(4 * pi)
}

Kuwata.A <- function(L){
  -0.4865*L - 1.046*L^2 + 0.8481*L^3
}

Kuwata.B <- function(L){
  0.01909*L + 0.19999 * L^2 + 0.6077 * L^3
}	
	

##' inverse polarizability tensors
##'
##' calculates and formats the principal polarizability of several particles
##' @title principal_polarizability
##' @param cluster cluster
##' @param material material
##' @param polarizability_fun polarizability function
##' @param ... additional arguments passed to polarizability_fun
##' @return  matrix with each row being the 3 principal values of each polarizability tensor
##' @importFrom plyr mlply
##' @export
##' @family user_level polarizability
##' @author Baptiste Auguie
principal_polarizability <- function(cluster, material, 
                          polarizability_fun = polarizability_ellipsoid, ...){
  polar <- mlply(cluster[['sizes']], polarizability_fun,
                 wavelength=material$wavelength,
                 epsilon=material$epsilon, ...)
  t(do.call(cbind, polar))
  
}
