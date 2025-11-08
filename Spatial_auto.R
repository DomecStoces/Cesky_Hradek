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
  transmute(
    Year=factor(Year),
    Locality = factor(Locality),
    Exposition2 = as.numeric(scale(Exposition2)),
    Altitude_scaled, Altitude_scaled2,
    X_km = (X - mean(X, na.rm = TRUE))/1000,
    Y_km = (Y - mean(Y, na.rm = TRUE))/1000,
    Moisture_cwm, Wings_cwm, Distribution_cwm01, Body_cwm,
    Bioindication_cwm, Breeding_cwm, Dietary_cwm
  ) %>% na.omit()
df <- df %>%
  mutate(eastness = sin(Exposition2),
         northness = cos(Exposition2))
# confirm:
stopifnot(is.numeric(df$Exposition2))

# choose appropriate k (<= unique sites)
k_xy <- max(6, min(10, nrow(dplyr::distinct(df, X_km, Y_km)) - 1))

# GAM: simpler polynomial representation is preferred for interpretability and model parsimony than smooth term of Altitude.
# Smooth factor must not restrict the model’s flexibility.
mod_gam1 <- gam(
  Wings_cwm ~ s(X_km, Y_km, bs = "tp", k = k_xy) +
    s(Altitude_scaled, k=3) + s(Year, bs="re") +
    Exposition2 +
    s(Locality, bs = "re"),
  data   = df, ,family = betar(),
  method = "REML"
)
summary(mod_gam1)
gam.check(mod_gam1)

,family = betar()
library(DHARMa)
library(qgam)
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



tiff('DHARMa_Moisture.tiff', units = "in", width = 8, height = 10, res = 600)
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

