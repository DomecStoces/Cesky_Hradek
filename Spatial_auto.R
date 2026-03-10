coords <- data.frame(
  Locality = 1:20,
  Latitude = c(
    50 + 32/60 + 30.05/3600,   # 1
    50 + 30/60 + 42.72/3600,   # 2
    50 + 30/60 + 22.64/3600,   # 3
    50 + 31/60 + 36.44/3600,   # 4
    50 + 31/60 + 40.13/3600,   # 5
    50 + 32/60 + 44.70/3600,   # 6
    50 + 33/60 + 33.68/3600,   # 7
    50 + 34/60 + 21.58/3600,   # 8
    50 + 35/60 + 16.41/3600,   # 9
    50 + 35/60 + 24.64/3600,   # 10
    50 + 35/60 + 57.02/3600,   # 11
    50 + 37/60 + 28.17/3600,   # 12
    50 + 37/60 + 53.25/3600,   # 13
    50 + 38/60 + 13.17/3600,   # 14
    50 + 38/60 + 8.31/3600,    # 15
    50 + 40/60 + 26.54/3600,   # 16
    50 + 41/60 + 39.80/3600,   # 17
    50 + 41/60 + 32.08/3600,   # 18
    50 + 46/60 + 39/3600,      # 19
    50 + 46/60 + 31/3600       # 20
  ),
  Longitude = c(
    13 + 26/60 + 28.39/3600,   # 1
    13 + 25/60 + 27.43/3600,   # 2
    13 + 24/60 + 30.10/3600,   # 3
    13 + 19/60 + 57.64/3600,   # 4
    13 + 20/60 + 24.96/3600,   # 5
    13 + 17/60 + 3.77/3600,    # 6
    13 + 17/60 + 19.52/3600,   # 7
    13 + 21/60 + 38.15/3600,   # 8
    13 + 19/60 + 36.85/3600,   # 9
    13 + 21/60 + 31.69/3600,   # 10
    13 + 22/60 + 32.76/3600,   # 11
    13 + 23/60 + 48.48/3600,   # 12
    13 + 24/60 + 2.11/3600,    # 13
    13 + 39/60 + 55.43/3600,   # 14
    13 + 40/60 + 5.60/3600,    # 15
    13 + 32/60 + 30.66/3600,   # 16
    13 + 34/60 + 34.51/3600,   # 17
    13 + 38/60 + 6.10/3600,    # 18
    14 + 5/60 + 27.01/3600,    # 19
    14 + 6/60 + 11.05/3600     # 20
  )
)
library(sf)
coords_sf  <- st_as_sf(coords, coords = c("Longitude", "Latitude"), crs = 4326)
coords_utm <- st_transform(coords_sf, crs = 32633)  # UTM 33N
coords_clean <- coords_utm %>%
  mutate(
    X = st_coordinates(.)[, 1],
    Y = st_coordinates(.)[, 2],
    Locality = as.character(Locality)
  ) %>%
  st_drop_geometry() %>%
  dplyr::select(Locality, X, Y)
cwm_clean <- cwm_clean %>% left_join(coords_clean, by = "Locality")

# GAM due to Moran’s I; this removes spatial autocorrelation.
library(dplyr)
library(mgcv)

df <- cwm_clean %>%
  transmute(
    Year      = factor(Year),
    Locality  = factor(Locality),
    Month     = factor(Month),
    Exposition2 = as.numeric(scale(as.numeric(Exposition2))),
    
    Altitude_scaled,
    
    X_km = (X - mean(X, na.rm = TRUE)) / 1000,
    Y_km = (Y - mean(Y, na.rm = TRUE)) / 1000,
    
    Moisture_cwm, Wings_cwm, Distribution_cwm, Body_size_cwm,
    Bioindication_cwm, Dietary_cwm
  ) %>%
  tidyr::drop_na(Wings_cwm, Altitude_scaled, Exposition2, Year, Locality, Month)

# Calculate adaptive basis dimension (k)
k_xy <- max(6, min(10, nrow(dplyr::distinct(df, X_km, Y_km)) - 1))

write_xlsx(df, "df.xlsx")
# GAM: simpler polynomial representation is preferred for interpretability and model parsimony than smooth term of Altitude.
# Deviations from the mid-domain null
eps <- 1e-6
df <- df |>
  mutate(
    mu0  = 0.6 + 0.3 * exp(-(Altitude_scaled)^2 / (2 * 0.5^2)),
    mu0  = pmin(pmax(mu0, eps), 1 - eps),
    eta0 = log(-log(1 - mu0))
  )
offset(eta0)

N <- nrow(df)
df$Wings_cwm_scaled        <- (df$Wings_cwm * (N - 1) + 0.5) / N
df$Dietary_cwm_scaled      <- (df$Dietary_cwm * (N - 1) + 0.5) / N
df$Breeding_cwm_scaled     <- (df$Breeding_cwm * (N - 1) + 0.5) / N
df$Distribution_cwm_scaled <- (df$Distribution_cwm * (N - 1) + 0.5) / N

mod_gam1 <- gam(
  Dietary_cwm_scaled ~ 
    s(Locality, bs = "re") +
    s(Altitude_scaled, bs = "cr", k = 5) + Exposition2 +
    s(Year, bs = "re"),
  data   = df,
  family = betar(link = "cloglog"),
  method = "REML"
)

s(X_km, Y_km, bs = "tp", k = k_xy)
s(Altitude_scaled, bs = "cr", k = 3)
s(Locality, bs = "re")

summary(mod_gam1)
par(mfrow = c(2, 2))
gam.check(mod_gam1)
concurvity(mod_gam1, full = TRUE)
gratia::draw(mod_gam1)
plot(mod_gam1, select = 2)

library(DHARMa)
library(qgam)
library(mgcViz)
sim <- simulateResiduals(fittedModel = mod_gam1, n = 2000, seed = 123)
plot(sim, qgam = TRUE) 

# correlogram (binned Moran’s I)
library(gstat)
library(sp)

# Residuals averaged per site to resolve the overlapping temporal data
res_dharma <- simulateResiduals(fittedModel = mod_gam1)
df$res_scaled <- res_dharma$scaledResiduals
df_site_res <- df %>%
  group_by(Locality, X_km, Y_km) %>%
  summarise(mean_res = mean(res_scaled, na.rm = TRUE), .groups = 'drop')

# Convert the summarized site data to a spatial object
coordinates(df_site_res) <- ~ X_km + Y_km

# Empirical variogram using the averaged DHARMa residuals
vg <- variogram(mean_res ~ 1, data = df_site_res, cutoff = 40, width = 2, cressie = TRUE)
plot(vg, main = "Residual Variogram")

tiff('Variogram_Wings.tiff', units = "in", width = 8, height = 6, res = 600)
plot(vg, main = "Residual Variogram")
dev.off()

# Plotting the effect of Altitude_scaled on CWM traits from mod_gam1
library(gratia)
library(dplyr)
library(tidyr)
library(ggplot2)

excl <- c("s(Locality)", "s(Year)")
tv <- typical_values(mod_gam1)

alt_seq <- seq(min(df$Altitude_scaled, na.rm = TRUE),
               max(df$Altitude_scaled, na.rm = TRUE), length.out = 100)

tv2 <- dplyr::select(tv, -any_of(c("Altitude_scaled","Altitude_scaled2")))

new_data <- tidyr::crossing(tv2, tibble(Altitude_scaled = alt_seq)) %>%
  mutate(Altitude_scaled2 = Altitude_scaled^2)
new_data <- tidyr::crossing(tv2, tibble(Altitude_scaled = alt_seq)) %>%
  mutate(Altitude_scaled2 = Altitude_scaled^2)
inv_link <- family(mod_gam1)$linkinv
fv <- fitted_values(mod_gam1, data = new_data, exclude = excl,
                    scale = "link", se = TRUE) %>%
  dplyr::rename(fitted_link = any_of(c("fitted",".fitted","fit")),
                se_link     = any_of(c("se",".se"))) %>%
  mutate(
    lower_link = fitted_link - (1.96 * se_link),
    upper_link = fitted_link + (1.96 * se_link),
    fitted = inv_link(fitted_link),
    lower  = inv_link(lower_link),
    upper  = inv_link(upper_link)
  )
fv_obs <- fitted_values(mod_gam1, data = df, exclude = excl,
                        scale = "response", se = FALSE) %>%
  dplyr::rename(fitted = any_of(c("fitted",".fitted","fit")))

df$partial_excl <- residuals(mod_gam1, type = "response") + fv_obs$fitted

p <- ggplot() +
  # 1. The confidence ribbon
  geom_ribbon(data = fv,
              aes(x = Altitude_scaled, ymin = lower, ymax = upper),
              fill = "grey70", alpha = 0.35) +
  
  # 2. The fitted line (the "estimated curve")
  geom_line(data = fv,
            aes(x = Altitude_scaled, y = fitted), linewidth = 1.1) +
  
  # 3. The points (using your dispersal column)
  geom_jitter(data = df,
              aes(x = Altitude_scaled, y = Wings_cwm_scaled),
              width = 0.03, height = 0, size = 1.8, alpha = 0.6) +
  
  # Labels
  labs(x = "Elevational gradient (scaled)", y = "Dispersal ability CWM") +
  
  # X-Axis: Matches the Rao's Q plot limits exactly
  scale_x_continuous(breaks = seq(-2, 2, 1), minor_breaks = NULL) +
  # Y-Axis: FORCED TO 0 - 1
  scale_y_continuous(
    limits = c(0, 1),             
    breaks = seq(0, 1, 0.2),         
    expand = expansion(mult = c(0, 0.02)) 
  ) +
  
  # Theme (kept identical to your Rao's Q setup)
  theme(
    panel.background = element_blank(),
    plot.background  = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 0.6),
    axis.ticks       = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length= unit(4, "pt"),
    axis.title       = element_text(size = 15),
    axis.text        = element_text(colour = "black", size = 11),
    plot.margin      = margin(6, 8, 6, 6)
  )
p
tiff('GAM_Distribution.tiff', units = "in", width = 8, height = 10, res = 600)
print(p)
dev.off()

# Fit GAMs for all CWM traits
library(mgcv)
library(purrr)

responses <- c("Distribution_cwm","Moisture_cwm","Wings_cwm",
               "Body_cwm","Bioindication_cwm","Breeding_cwm","Dietary_cwm")

fit_one_gam <- function(y, df, k_xy) {
  gam(
    formula = reformulate(
      termlabels = c(sprintf("s(X_km,Y_km,bs='tp',k=%d)", k_xy),
                     "Altitude_scaled","Altitude_scaled2","Exposition2",
                     "s(Locality,bs='re')"),
      response  = y
    ),
    data   = df,
    family = gaussian(),
    method = "REML"
  )
}

fits <- map(setNames(responses, responses), ~ fit_one_gam(.x, df, k_xy))
summary(fits$Moisture_cwm)

