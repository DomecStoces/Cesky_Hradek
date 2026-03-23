library(DHARMa)
library(qgam)
library(mgcViz)
library(dplyr)
library(gstat)
library(sp)
library(spdep)
library(mgcv)
library(readxl)
library(picante)

### Analysis of Phylogeny Sespd and meanPD ###
PD <- read_excel("PD.xlsx", sheet = "List1")
PD$Locality <- as.factor(PD$Locality)
PD$Year     <- as.factor(PD$Year)
mod_gam_pd <- gam(
  meanPD ~ s(Altitude_scaled, bs = "cr", k = 5) + s(Locality, bs = "re") +
    Exposition2 +
    s(Year, bs = "re"),
  data = PD, 
  family = gaussian(), 
  method = "REML"
)

summary(mod_gam_pd)
par(mfrow = c(2, 2))
gam.check(mod_gam_pd)
concurvity(mod_gam_pd, full = TRUE)
gratia::draw(mod_gam_pd)
plot(mod_gam_pd, select = 2)

# correlogram (autocorrelation using Moran’s I based on site-averaged Pearson residuals)
library(DHARMa)
library(qgam)
library(mgcViz)
library(dplyr)
library(gstat)
library(sp)
library(spdep)

PD$resid <- residuals(mod_gam_pd, type = "pearson")
df_site_res <- PD %>%
  group_by(Locality, X_km, Y_km) %>%
  summarise(mean_res = mean(resid, na.rm = TRUE), .groups = "drop")
coords <- as.matrix(df_site_res[,c("X_km","Y_km")])
nb <- dnearneigh(coords, 0, 10)   
lw <- nb2listw(nb, style = "W")
moran.test(df_site_res$mean_res, lw)
coordinates(df_site_res) <- ~X_km + Y_km
vg <- variogram(mean_res ~ 1,
                data = df_site_res,
                cutoff = 40,
                width = 2,
                cressie = TRUE)
plot(vg, main = "Empirical variogram of GAM residuals")

# Because multiple observations were collected within the same sites (hierarchical structure), locality was included as a random effect to account for spatial clustering and avoid pseudoreplication.
PD %>%
  group_by(Locality) %>%
  summarise(alt = mean(Altitude_scaled))

### Graphical vizualization of SES.pd ###
mod_gam_sespd <- gam(
  SESpd ~ s(Altitude_scaled, bs = "cr", k = 5) + s(Locality, bs = "re") +
    Exposition2 +
    s(Year, bs = "re"),
  data = df, 
  family = gaussian(), 
  method = "REML"
)


library(ggplot2)
newdat <- data.frame(
  Altitude_scaled = seq(min(df$Altitude_scaled), max(df$Altitude_scaled), length = 200),
  Exposition2 = mean(df$Exposition2),
  Locality = NA,
  Year = NA
)
pred <- predict(mod_gam_sespd, newdata = newdat, se.fit = TRUE,
                exclude = c("s(Locality)", "s(Year)"))
newdat$fit <- pred$fit
newdat$se  <- pred$se.fit
newdat$upper <- newdat$fit + 1.96 * newdat$se
newdat$lower <- newdat$fit - 1.96 * newdat$se

phylo<-ggplot(df, aes(x = Altitude_scaled, y = SESpd)) +
  
  # The raw points - now matching the CWM plot using geom_jitter
  geom_jitter(width = 0.03, height = 0, size = 1.8, alpha = 0.6, color = "black") +  
  
  # The GAM trendline
  geom_smooth(method = "gam", color = "black", fill = "grey70", alpha = 0.3) +
  
  # The 0 baseline 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5, alpha = 0.6) + 
  
  # The +/- 1.96 significance thresholds
  geom_hline(yintercept = 1.96, linetype = "dashed", color = "black", alpha = 0.7) +
  geom_hline(yintercept = -1.96, linetype = "dashed", color = "black", alpha = 0.7) +
  
  # X-Axis: Matches the CWM plot limits exactly
  scale_x_continuous(breaks = seq(-2, 2, 1), minor_breaks = NULL) +
  
  theme_minimal() +
  labs(
    x = "Elevational gradient (Scaled)",
    y = "SESpd") +
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
phylo