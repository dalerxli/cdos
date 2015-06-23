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


P0 <- cdae::cda$polarization(Ei, DiagBlocks)
P1 <- cda$iterate_test(cl$r, kn, Ei, DiagBlocks,P0)
P2 <- cda$iterate_test(cl$r, kn, Ei, DiagBlocks,P1)
P3 <- cda$iterate_test(cl$r, kn, Ei, DiagBlocks,P2)
cda$extinction(kn, P0, Ei)
cda$extinction(kn, P1, Ei)
cda$extinction(kn, P2, Ei)



tol <- 1e-10
Niter <- 5
test <- cda$convergence(cl$r, kn, Ei, DiagBlocks, Niter, tol)

cext <- cda$extinction(kn, P, Ei)
# res <- circular_dichroism_spectrum(cl, gold)

test
cext

