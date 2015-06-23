## ----load,message=FALSE--------------------------------------------------
library(cdos)
# library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)

## ----cluster, rgl=TRUE,echo=-12,tidy=FALSE,fig.width=3,fig.height=3,fig.path="basic-"----

# dielectric function
wvl <- seq(400, 900)
gold <- epsAu(wvl)

# define a cluster of particles
cl <- list(r = rbind(c(0, 0, 0),
                      c(0, 0, -100)),
            angles = rbind(c(0, 0, 0),
                           c(pi/4, 0, 0)),
            sizes = rbind(c(30, 10, 10),
                          c(30, 10, 10)))

# visualise
# rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
# rgl.viewpoint( theta = 0, phi = 20, fov = 70, zoom = 1)
# rgl_annotate()



## ----oa,echo=TRUE,tidy=FALSE,fig.path="basic-",fig.width=8---------------
circular <- circular_dichroism_spectrum(cl, gold, Niter=1, tol=1e-10)
circular2 <- cdae::circular_dichroism_spectrum(cl, gold)
# circular3 <- cda::circular_dichroism_spectrum(cl, gold)

ggplot(circular, aes(wavelength, value, color=variable)) + 
  facet_grid(type~., scales="free") + 
  geom_line() +
  geom_line(data=circular2,linetype=2)+
  # geom_line(data=circular3,linetype=3) +
  theme()



