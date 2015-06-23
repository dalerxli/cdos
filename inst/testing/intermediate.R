library(cdos)

wavelength <- 500
medium <- 1.33
kn <- 2*pi/wavelength*medium
gold <- epsAu(wavelength)
E0L=1/sqrt(2)*c(0,1,1i)
E0R=1/sqrt(2)*c(0,1i,1)
k0=c(1,0,0)

cl <- cluster_helix(10)
cl <- cluster_chain(2)
Alpha <- principal_polarizability(cl, material=gold, medium=medium)
DiagBlocks <- cda$diagonal_blocks(Alpha, cl$angles)

Angles <- rbind(c(0, pi/2, 0), # +x is phi=0, psi=0
                c(pi/2, pi/2, 0)) # +z is phi=pi/2, psi=pi/2


A <- cdae::cda$interaction_matrix(cl$r, kn, DiagBlocks, TRUE)
Ei <- cda$incident_field(E0L, k=kn*k0, r=cl$r, Angles)

E <- solve(A, Ei)
P <- cdae::cda$polarization(E, DiagBlocks)
# [,1]                        [,2]
# [,1]                        [,2]
# [1,]     0.00+    0.00i  3.459052e+04-2.369239e+04i
# [2,] 40691.23-15981.99i -3.427800e-13+2.654874e-12i
# [3,] 15981.99+40691.23i -4.335738e+04-5.598022e+03i
# [4,]     0.00+    0.00i  4.097692e+03+4.172580e+04i
# [5,] 40691.23-15981.99i -2.161348e-12-1.579376e-12i
# [6,] 15981.99+40691.23i  2.579316e+04-3.529749e+04i

P0 <- cdae::cda$polarization(Ei, DiagBlocks)

test <- cda$iterate_test(cl$r, kn, Ei, DiagBlocks,P0)


tol <- 1e-10
Niter <- 3
test <- cda$convergence(cl$r, kn, Ei, DiagBlocks, Niter, tol)

cext <- cda$extinction(kn, P, Ei)
# res <- circular_dichroism_spectrum(cl, gold)

test
cext

