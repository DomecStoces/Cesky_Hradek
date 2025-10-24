coords <- data.frame(
  Locality = 1:18,
  Latitude = c(50 + 32/60 + 30.05/3600,
               50 + 30/60 + 42.72/3600,
               50 + 30/60 + 22.64/3600,
               50 + 31/60 + 36.44/3600,
               50 + 31/60 + 40.13/3600,
               50 + 32/60 + 44.70/3600,
               50 + 33/60 + 33.68/3600,
               50 + 34/60 + 21.58/3600,
               50 + 35/60 + 16.41/3600,
               50 + 35/60 + 24.64/3600,
               50 + 35/60 + 57.02/3600,
               50 + 37/60 + 28.17/3600,
               50 + 37/60 + 53.25/3600,
               50 + 38/60 + 13.17/3600,
               50 + 38/60 + 8.31/3600,
               50 + 40/60 + 26.54/3600,
               50 + 46/60 + 39/3600,
               50 + 41/60 + 32.08/3600),
  Longitude = c(13 + 26/60 + 28.39/3600,
                13 + 25/60 + 27.43/3600,
                13 + 24/60 + 30.10/3600,
                13 + 19/60 + 57.64/3600,
                13 + 20/60 + 24.96/3600,
                13 + 17/60 + 3.77/3600,
                13 + 17/60 + 19.52/3600,
                13 + 21/60 + 38.15/3600,
                13 + 19/60 + 36.85/3600,
                13 + 21/60 + 31.69/3600,
                13 + 22/60 + 32.76/3600,
                13 + 23/60 + 48.48/3600,
                13 + 24/60 + 2.11/3600,
                13 + 39/60 + 55.43/3600,
                13 + 40/60 + 5.60/3600,
                13 + 32/60 + 30.66/3600,
                14 + 5/60 + 27.01/3600,
                14 + 6/60 + 11.05/3600)
)
library(sf)

coords_sf <- st_as_sf(coords, coords = c("Longitude", "Latitude"), crs = 4326)
coords_utm <- st_transform(coords_sf, crs = 32633)  # UTM zone 33N (for Czech Republic)
coords_xy <- st_coordinates(coords_utm)
coords <- bind_cols(coords, as.data.frame(coords_xy))
names(coords)[4:5] <- c("X", "Y")
cwm_clean <- cwm_clean %>% mutate(Locality = as.character(Locality))
coords    <- coords    %>% mutate(Locality = as.character(Locality))
cwm_clean <- cwm_clean %>%
  left_join(coords[, c("Locality", "X", "Y")], by = "Locality")

setdiff(unique(cwm_clean$Locality), unique(coords$Locality))
setdiff(unique(coords$Locality), unique(cwm_clean$Locality))
summary(cwm_clean[, c("X", "Y")])

# GAM due to Moran’s I; this removes spatial autocorrelation.
library(dplyr)
library(mgcv)

df <- cwm_clean %>%
  mutate(Locality = factor(Locality),
         Year     = factor(Year)) %>%
  select(Moisture_cwm, Wings_cwm, Distribution_cwm, Body_cwm, Bioindication_cwm, Breeding_cwm, Dietary_cwm, X, Y, Altitude_scaled, Altitude_scaled2,
         Exposition2, Locality) %>%
  na.omit() %>%
  mutate(X_km = (X - mean(X))/1000,
         Y_km = (Y - mean(Y))/1000)

# choose appropriate k (<= unique sites)
k_xy <- max(6, min(10, nrow(dplyr::distinct(df, X_km, Y_km)) - 1))

# GAM: simpler polynomial representation is preferred for interpretability and model parsimony than smooth term of Altitude.
mod_gam1 <- gam(
  Wings_cwm ~ s(X_km, Y_km, bs = "tp", k = k_xy) +
    Altitude_scaled + Altitude_scaled2 +  
    Exposition2 + Year +
    s(Locality, bs = "re"),
  data   = df,
  method = "REML"
)
summary(mod_gam1)

library(DHARMa)
sim <- simulateResiduals(fittedModel = mod_gam1, n = 1000, seed = 1)
plot(sim)

vis.gam(mod_gam2, view = c("X_km", "Y_km"),
        plot.type = "contour", color = "terrain",
        too.far = 0.05, n.grid = 100)

vis.gam(mod_gam1, view = c("X_km", "Y_km"),
        plot.type = "persp", color = "topo",
        phi = 30, theta = 30)

library(gratia)
library(ggplot2)

draw(mod_gam2, select = "s(X_km,Y_km)") +
  labs(fill = "Partial effect") +
  theme_bw()

draw(mod_gam2, select = "s(Altitude)") +
  labs(y = "Partial effect", x = "Altitude (m a.s.l.)") +
  theme_bw()

# If you refit mod_gam2 with s(Altitude, k = 5):
library(gratia)

draw(mod_gam1, select = "s(Altitude)") +
  labs(
    x = "Altitude (m a.s.l.)",
    y = "Partial effect on Wing CWM",
    title = "Modelled altitude effect (smooth representation)"
  ) +
  theme_bw()
