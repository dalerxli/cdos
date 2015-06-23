library(cdos)

wavelength <- 500
medium <- 1.33
kn <- 2*pi/wavelength*medium
gold <- epsAu(wavelength)
E0L=1/sqrt(2)*c(0,1,1i)
E0R=1/sqrt(2)*c(0,1i,1)
k0=c(1,0,0)

cl <- cluster_helix(10)
cl <- cluster_chain(3)
Alpha <- principal_polarizability(cl, material=gold, medium=medium)
DiagBlocks <- cda$diagonal_blocks(Alpha, cl$angles)

Angles <- rbind(c(0, pi/2, 0), # +x is phi=0, psi=0
                c(pi/2, pi/2, 0)) # +z is phi=pi/2, psi=pi/2


A <- cdae::cda$interaction_matrix(cl$r, kn, DiagBlocks, TRUE)
Ei <- cda$incident_field(E0L, k=kn*k0, r=cl$r, Angles)

E <- solve(A, Ei)
P <- cdae::cda$polarization(E, DiagBlocks)

# [,1]                        [,2]                        [,3]
# [1,]  1.783511e+04-3.796889e+04i -3.796889e+04-1.783511e+04i  2.468936e+04+3.516560e+04i
# [2,]  4.064575e+04-1.565036e+04i -1.565036e+04-4.064575e+04i -2.639099e-12-3.311256e-13i
# [3,] -1.410068e-12+3.429610e-12i -1.410068e-12+3.429610e-12i  4.391218e+03-4.383143e+04i
# [4,]  1.614083e+04-3.731393e+04i -3.731393e+04-1.614083e+04i -4.241590e+04+2.886277e+03i
# [5,]  4.077808e+04-1.526199e+04i -1.526199e+04-4.077808e+04i  1.565337e-12-2.177672e-12i
# [6,] -1.361255e-12+3.266117e-12i -1.361255e-12+3.266117e-12i  3.628919e+04+2.700353e+04i
# [7,]  1.783511e+04-3.796889e+04i -3.796889e+04-1.783511e+04i  1.766876e+04-3.677939e+04i
# [8,]  4.064575e+04-1.565036e+04i -1.565036e+04-4.064575e+04i  1.151293e-12+2.446341e-12i
# [9,] -1.410068e-12+3.429610e-12i -1.410068e-12+3.429610e-12i -4.111912e+04+1.676112e+04i

# it1 <- cda$iterate_field(cl$r, kn, Ei, DiagBlocks, E, P)

tol <- 1e-10
Niter <- 5
test <- cda$convergence(cl$r, kn, Ei, DiagBlocks, Niter, tol)

cext <- cda$extinction(kn, P, Ei)
res <- circular_dichroism_spectrum(cl, gold)

test
cext

