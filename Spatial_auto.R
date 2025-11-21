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
    Y = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  transmute(
    Locality = as.character(Locality),
    X, Y
  )
cwm_clean <- cwm_clean %>%
  mutate(Locality = as.character(Locality)) %>%
  left_join(coords_clean, by = "Locality")
names(cwm_clean)          
summary(cwm_clean[, c("X", "Y")])

# GAM due to Moran’s I; this removes spatial autocorrelation.
library(dplyr)
library(mgcv)

df <- cwm_clean %>%
  transmute(
    Year      = factor(Year),
    Locality  = factor(Locality),
    HR        = factor(HR),
    Exposition2_raw = Exposition2,
    Exposition2     = as.numeric(scale(Exposition2)),
    
    Altitude_scaled, Altitude_scaled2,
    X_km = (X - mean(X, na.rm = TRUE))/1000,
    Y_km = (Y - mean(Y, na.rm = TRUE))/1000,
    Moisture_cwm, Wings_cwm, Distribution_cwm, Body_cwm,
    Bioindication_cwm, Breeding_cwm, Dietary_cwm
  ) %>%
  na.omit()
# confirm:
stopifnot(is.numeric(df$Exposition2))

# choose appropriate k (<= unique sites)
k_xy <- max(6, min(10, nrow(dplyr::distinct(df, X_km, Y_km)) - 1))

# GAM: simpler polynomial representation is preferred for interpretability and model parsimony than smooth term of Altitude.
# Deviations from the mid-domain null
eps <- 1e-6
df <- df |>
  dplyr::mutate(
    mu0 = 0.6 + 0.3 * exp(-(Altitude_scaled)^2 / (2 * 0.5^2)),
    mu0 = pmin(pmax(mu0, eps), 1 - eps),
    eta0 = log(-log(1 - mu0))
  )

mod_gam1 <- gam(
 Distribution_cwm ~ s(X_km, Y_km, bs = "tp", k = k_xy) + 
    s(Altitude_scaled, bs = "cr", k = 3)+
    s(Year, bs = "re") +
   s(Locality, bs="re") +
    offset(qlogis(mu0)),
  data = df, family = betar(link="cloglog"), method = "REML"
)

anova(mod_gam1,mod_gam2)

summary(mod_gam1)
gam.check(mod_gam1)
concurvity(mod_gam1, full = TRUE)

library(DHARMa)
library(qgam)
library(mgcViz)
sim <- simulateResiduals(fittedModel = mod_gam1, n = 2000, seed = 123)
plot(sim, qgam = TRUE) 

# correlogram (binned Moran’s I)
library(gstat)
library(sp)

# 1) Numeric residuals (use a NEW name)
r_pearson <- residuals(mod_gam1, type = "pearson")  # length should be 72

# 2) Bind to df and convert to spatial
df_res <- df
df_res$r <- r_pearson

coordinates(df_res) <- ~ X_km + Y_km   # coords already in km; CRS not required for variogram

# 3) Empirical variogram (robust; sensible bins)
vg <- variogram(r ~ 1, data = df_res, cutoff = 40, width = 2, cressie = TRUE)

plot(vg, main = "Residual variogram (Pearson)")

tiff('DHARMa_Dispersal.tiff', units = "in", width = 8, height = 6, res = 600)
plot(sim)
dev.off()

# Plotting the effect of Altitude_scaled on CWM traits from mod_gam1
library(gratia)
library(dplyr)
library(tidyr)
library(ggplot2)

# Terms to exclude
excl <- c("s(Locality)", "s(X_km,Y_km)", "s(Year)")

# 1) Grid for the line/ribbon
tv <- typical_values(mod_gam1)
alt_seq <- seq(min(df$Altitude_scaled, na.rm = TRUE),
               max(df$Altitude_scaled, na.rm = TRUE), length.out = 100)

tv2 <- dplyr::select(tv, -any_of(c("Altitude_scaled","Altitude_scaled2")))
new_data <- tidyr::crossing(tv2, tibble(Altitude_scaled = alt_seq)) %>%
  mutate(Altitude_scaled2 = Altitude_scaled^2)

# 2) Fitted line/ribbon on response scale
fv <- fitted_values(mod_gam1, data = new_data, exclude = excl,
                    scale = "response", se = TRUE) %>%
  dplyr::rename(fitted = any_of(c("fitted",".fitted","fit")),
                se     = any_of(c("se",".se"))) %>%
  mutate(lower = fitted - 1.96 * se,
         upper = fitted + 1.96 * se)

# 3) Partial points aligned with exclusions
fv_obs <- fitted_values(mod_gam1, data = df, exclude = excl,
                        scale = "response", se = FALSE) %>%
  dplyr::rename(fitted = any_of(c("fitted",".fitted","fit")))

df$partial_excl <- residuals(mod_gam1, type = "response") + fv_obs$fitted

p <- ggplot() +
  geom_ribbon(data = fv,
              aes(x = Altitude_scaled, ymin = lower, ymax = upper),
              fill = "grey70", alpha = 0.35) +
  geom_line(data = fv,
            aes(x = Altitude_scaled, y = fitted), linewidth = 1.1) +
  geom_jitter(data = df,
              aes(x = Altitude_scaled, y = partial_excl),
              width = 0.03, height = 0, size = 1.8, alpha = 0.6) +
  labs(x = "Altitude (scaled)", y = "Moisture preferences CWM") +
  scale_x_continuous(breaks = seq(-2, 2, 1), minor_breaks = NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks(5),
                     expand = expansion(mult = c(0, 0.02))) +
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

