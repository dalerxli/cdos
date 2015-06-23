
## ----load,message=FALSE--------------------------------------------------
library(microbenchmark)


# dielectric function
wvl <- seq(400, 900, length=2)
gold <- dielectric::epsAu(wvl)

# define a helix
cl <- cda::cluster_helix(N=5000, R0=100, pitch=200, 
                          delta=pi/2, delta0=0, right=TRUE,
                          a=50, b=20, c=20,
                          angles='helix')

e <- cdae::circular_dichroism_spectrum(cl, material = gold, 
                                       averaging="cheap", progress=TRUE,
                                       medium=1.33, cg=TRUE, born=TRUE, tol=1e-3)


p <- cda::circular_dichroism_spectrum(cl, material = gold, medium=1.33, cg=TRUE, born=TRUE)
e <- cdae::circular_dichroism_spectrum(cl, material = gold, medium=1.33, cg=TRUE, born=TRUE)

library(ggplot2)
ggplot(p, aes(wavelength, value, colour=variable)) +
  facet_wrap(~type, scale="free") + geom_line() + geom_line(data=e, lty=2)

microbenchmark(
  p = cda::circular_dichroism_spectrum(cl, material = gold, 
                                       medium=1.33, cg=TRUE, born=TRUE),
  e = cdae::circular_dichroism_spectrum(cl, material = gold,
                                        medium=1.33, cg=TRUE, born=TRUE), times = 3)
