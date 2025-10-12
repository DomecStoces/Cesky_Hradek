# Generalized least squares model with spatial correlation structure due to significant Moran’s I

library(nlme)
mod_gls <- gls(Wings_cwm ~ poly(Altitude, 2, raw = TRUE) + Exposition2,
               correlation = corExp(form = ~ Longitude + Latitude, nugget = TRUE),
               data = cwm_clean)
summary(mod_gls)

# spatial eigenvector mapping (MEM / dbMEM)

library(adespatial)
dbmem_vars <- dbmem(as.matrix(coords))
rda_cwm_spatial <- rda(cwm_mat ~ Altitude + Exposition2 + dbmem_vars[, 1:3], data = cwm_clean)
