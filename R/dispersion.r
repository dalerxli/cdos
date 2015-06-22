
##' dispersion spectrum
##'
##' dispersion spectrum
##' @title dispersion_spectrum
##' @param cluster list describing a cluster
##' @param material list
##' @param medium medium refractive index
##' @param angles of incident field in radians
##' @param axes of incident field rotation character vector from ('x', 'y', 'z')
##' @param polarisation linear or circular polarisation
##' @param cg logical, use conjugate gradient solver
##' @param born logical, use first Born approx as cg guess
##' @param nmax integer termination of conjugate gradient solver
##' @param tol double, tolerance of conjugate gradient solver
##' @param progress logical, display progress bar
##' @return data.frame
##' @note The incident wavevector is along the z direction.
##' @export
##' @family user_level cda
##' @author baptiste Auguie
dispersion_spectrum <- function (cluster, material, medium = 1.33,
                                 angles=0, axes='z', 
                                 polarisation=c("linear", "circular"), 
                                 cg = FALSE, born = FALSE, nmax=30, tol=1e-4,
                                 progress = FALSE) 
{
  k0 <- 2 * pi/material[['wavelength']]
  kn <- k0 * medium
  polarisation <- match.arg(polarisation)
  
  polarisation <- if(polarisation == "linear") 0L else if(polarisation == "circular") 1L 
  
  Alpha <- principal_polarizability(cluster, material, 
                                 polarizability_fun = polarizability_ellipsoid, 
                                 medium = medium, kuwata = TRUE)
  Nwavelengths <- length(k0)
  Nparticles <- nrow(cluster$r)
  Nangles <- length(angles)
  
  if(length(axes) == 1) axes <- rep(axes, length.out=Nangles)
  axeso <- axes # original codes
  axes <- as.integer(factor(axes, levels=c('x','y','z')))-1L
  stopifnot(all(axes %in% c(0L, 1L, 2L)), !any(is.na(axes)))

  stopifnot(Nangles == length(axes))
  stopifnot(is.matrix(Alpha), is.vector(angles), 
            is.matrix(cluster$r), 
            is.matrix(cluster$angles))
  
  stopifnot(nrow(Alpha)/3 == Nparticles, 
            ncol(Alpha) == Nwavelengths)
  
  res <- dispersion$dispersion_spectrum(kn, cluster$r, Alpha, 
                                        cluster$angles, angles, axes, 
                                        polarisation, cg, born, nmax, tol, progress)
  
  angles <- angles[rep(seq.int(Nangles), Nwavelengths)]
  axes <- axeso[rep(seq.int(Nangles), Nwavelengths)]
  wavelength <- rep(material$wavelength, each = Nangles)
  
  results <- 
    rbind(data.frame(wavelength = wavelength, angles = angles,
                     axes=axes,
                     value = c(res[, 1, , drop = TRUE]),
                     type = "extinction", polarisation = "1"),
          data.frame(wavelength = wavelength, angles = angles,
                     axes=axes,
                     value = c(res[, 2, , drop = TRUE]),
                     type = "absorption", polarisation = "1"),
          data.frame(wavelength = wavelength, angles = angles,
                     axes=axes,
                     value = c(res[, 3, , drop = TRUE]),
                     type = "scattering", polarisation = "1"),
          data.frame(wavelength = wavelength, angles = angles,
                     axes=axes,
                     value = c(res[, 4, , drop = TRUE]),
                     type = "extinction", polarisation = "2"),
          data.frame(wavelength = wavelength, angles = angles,
                     axes=axes,
                     value = c(res[, 5, , drop = TRUE]),
                     type = "absorption", polarisation = "2"),
          data.frame(wavelength = wavelength, angles = angles,
                     axes=axes,
                     value = c(res[, 6, , drop = TRUE]),
                    type = "scattering", polarisation = "2"))
  
  if(polarisation == 0L)
    results$polarisation <- factor(results$polarisation, labels = c("p", "s"))
  
  if(polarisation == 1L){
    
    results$polarisation <- factor(results$polarisation, labels = c("R", "L"))
    results <- rbind(results,  data.frame(wavelength = wavelength, angles = angles, axes=axes,
                                          value = c(res[, 1, , drop = TRUE]  - res[, 4, , drop = TRUE]),
                                          type = "extinction", polarisation = "CD"))
  }
  
  invisible(results)
}
