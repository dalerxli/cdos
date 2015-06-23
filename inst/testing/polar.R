
# test <- polarizability_dye(wavelength, a=1, b=0.5, c=0.2)
# 
# test2 <- polarizability_ellipsoid(wavelength, epsAg(wavelength)$epsilon, 
#                                   a = 0.8, b=0.8, c=0.8)
# 
# matplot(wavelength, Im(test2), t="l", col=1)
# matlines(wavelength, Im(test), t="l", col=2)


shell_spectrum <- function (cluster, wavelength, 
                            fun=polarizability_dye, ...,
                            medium = 1.33, core=TRUE,
                            Nquad = 100, averaging = c("cheap", "QMC", "GL"), 
                            iterative = FALSE, precision = 0.001, 
                            Qmax = 10000, dN = Nquad, cg = TRUE, born = TRUE, nmax = 10, 
                            tol = 1e-02, full = TRUE, progress = FALSE, verbose = TRUE, 
                            result.matrix = FALSE) 
{
  averaging <- match.arg(averaging)
  
  epsilon <- medium^2
  prefact <- ((epsilon + 2 /3)^2) / sqrt(epsilon)
  
  k0 <- 2 * pi/wavelength
  kn <- k0 * medium
  polar <- plyr::mlply(cluster[["sizes"]], fun, wavelength = wavelength, ...)
  if(core){
  core <- polarizability_ellipsoid(wavelength, epsilon = epsAg(wavelength)$epsilon,
                                   a=cluster$R0,b=cluster$R0, medium=medium)
  Alpha <- rbind(t(do.call(cbind, polar)) * prefact, t(core))
  cluster$r <- rbind(cluster$r, c(0,0,0))
  cluster$angles <- rbind(cluster$angles, c(0,0,0))
  } else {
    Alpha <- t(do.call(cbind, polar)) * prefact
  }
  # browser()
  quadrature <- cdae:::integration_points(averaging, Nquad)
  results <- cd$average_spectrum(kn, cluster$r, Alpha, cluster$angles, 
                                 as.matrix(quadrature$angles), quadrature$weights, full, 
                                 cg, born, nmax, tol, progress)
  if (iterative && averaging == "QMC") {
    converged <- FALSE
    Ntot <- Nquad
    while (Ntot < Qmax && !converged) {
      oldN <- Ntot
      old <- results[, 1]
      Ntot <- Ntot + dN
      quadrature <- cdae:::integration_points(averaging, dN, FALSE)
      newres <- cd$average_spectrum(kn, Alpha, cluster$r, 
                                    cluster$angles, as.matrix(quadrature$angles), 
                                    quadrature$weights, full, cg, nmax, tol, progress)
      results <- (oldN * results + dN * newres)/(oldN + 
                                                   dN)
      test <- max(abs(old - results[, 1])/results[, 1])
      if (verbose) 
        message("N:", Ntot, "; relative error: ", test)
      converged <- test < precision
    }
  }

    d <- data.frame(wavelength, results[,1:3])
    names(d) <- c("wavelength", "extinction", "absorption", "scattering")
    m <- reshape2::melt(d, id = c("wavelength"))
   
    return(m)
  
}

library(cdae)
library(rgl)

# test1 <- test
# test2 <- test
wavelength <- seq(400, 600, length=50)
N <- 150
cl <- cluster_shell(N = N, d=1, R0 = 15,  a=1, b=1, c=1, 
                    position = "quasi-random", orientation = "radial")

cl1 <- cluster_shell(N = 2, d=1, R0 = 50, a=1, b=0, c=0, 
                    position = "quasi-random", orientation = "radial")
# cl1 <- cl
# rgl.spheres(0,0,0,10)
# rgl.ellipsoids(cl1$r, cl1$sizes, cl1$angles, col="red")

test1 <- shell_spectrum(cluster = cl1, core=FALSE, wavelength = wavelength, cg=FALSE, 
                       averaging = "QMC", Nquad = 50,
                       medium=1, progress = TRUE)

test <- shell_spectrum(cluster = cl, core=TRUE, wavelength = wavelength, cg=FALSE, 
                       averaging = "QMC", Nquad = 50,tol=1e-2, nmax=10,
                       medium=1, progress = TRUE)

library(ggplot2)
ggplot(test, aes(wavelength, value/N, colour=variable,group=variable)) + 
  geom_line() + 
  geom_line(aes(wavelength, value/2*10), data=test1, lty=2) +
  geom_line(aes(wavelength, value/150), data=test2, lty=3) +
  # geom_line(aes(wavelength, value/N), 
  # data=subset(test2, type == "cross section"), lty=3) +
  theme()

 ggsave("tests.pdf")
