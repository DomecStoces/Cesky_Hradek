library(dplyr)
library(mgcv)
library(readxl)

### Analysis of Phylogeny Sespd and meanPD ###
PD <- read_excel("PD.xlsx", sheet = "Sheet1")
PD$Elevation_Group  <- as.factor(PD$Elevation_Group)
PD$Exposition_Group <- as.factor(PD$Exposition_Group)
PD$Locality         <- as.factor(PD$Locality)
PD$Year             <- as.factor(PD$Year)
PD$Month            <- as.factor(PD$Month)
PD <- PD %>%
  mutate(Delta_scaled = Delta_plus / 100)
n <- nrow(PD)
PD <- PD %>%
  mutate(Delta_beta = (Delta_scaled * (n - 1) + 0.5) / n)

mod_gam_pd <- gam(
  Delta_beta ~ Elevation_Group +          
    Exposition_Group +         
    s(Locality, bs = "re") + s(Month, bs = "re") +
    s(Year, bs = "re"),        
  data   = PD,
  family = betar(link = "logit"),
  method = "REML"
)

summary(mod_gam_pd)
par(mfrow = c(2, 2))
gam.check(mod_gam_pd)
concurvity(mod_gam_pd, full = TRUE)
gratia::draw(mod_gam_pd)
plot(mod_gam_pd, select = 2)

# correlogram (autocorrelation using Moran’s I based on site-averaged Pearson residuals)
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

### Graphical vizualization of SES.pd ###
mod_gam_sespd <- gam(
  SES_Delta_plus ~ Elevation_Group +          
    Exposition_Group +         
    s(Locality, bs = "re") + s(Month, bs = "re") +
    s(Year, bs = "re"),        
  data   = PD,
  family = gaussian(),
  method = "REML"
)
summary(mod_gam_sespd)

# 1. Optional but recommended: Set the logical order of your groups
PD$Elevation_Group <- factor(PD$Elevation_Group, levels = c("Low", "Medium", "High"))

# 2. Create the prediction grid
# We extract predictions for each Elevation_Group. 
# We hold Exposition_Group constant at its most common level (or reference level).
newdat <- data.frame(
  Elevation_Group = levels(PD$Elevation_Group),
  Exposition_Group = levels(as.factor(PD$Exposition_Group))[1], # Holds factor constant
  Locality = NA,
  Month = NA,
  Year = NA
)

# 3. Run predictions excluding the random effects
pred <- predict(mod_gam_sespd, newdata = newdat, se.fit = TRUE,
                exclude = c("s(Locality)", "s(Month)", "s(Year)"))

newdat$fit <- pred$fit
newdat$se  <- pred$se.fit
newdat$upper <- newdat$fit + 1.96 * newdat$se
newdat$lower <- newdat$fit - 1.96 * newdat$se

library(ggplot2)

phylo <- ggplot(PD, aes(x = Elevation_Group, y = SES_Delta_plus)) +
  
  # 1. The raw points - widened jitter slightly for categorical separation
  geom_jitter(width = 0.15, height = 0, size = 1.8, alpha = 0.3, color = "black") +  
  
  # 2. The predicted GAM marginal means and 95% CIs
  geom_pointrange(data = newdat, 
                  aes(y = fit, ymin = lower, ymax = upper), 
                  color = "red", size = 0.8, linewidth = 1.2) +
  
  # 3. The 0 baseline 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5, alpha = 0.6) + 
  
  # 4. The +/- 1.96 significance thresholds
  geom_hline(yintercept = 1.96, linetype = "dashed", color = "black", alpha = 0.7) +
  geom_hline(yintercept = -1.96, linetype = "dashed", color = "black", alpha = 0.7) +
  
  # 5. Formatting 
  theme_minimal() +
  labs(
    x = "Elevational Group",
    y = expression("SES " * Delta^"+") # This nicely renders as SES Δ+ in the plot!
  ) +
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

# View the plot
phylo
