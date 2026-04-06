library(readxl)
library(dplyr)
library(fundiversity)
library(tibble)
library(cluster)
df <- read_excel("df.xlsx")

sp_data <- read_excel("Rao_diversity.xlsx", sheet = "sp")
traits_data <- read_excel("Rao_diversity.xlsx", sheet = "traits")
traits_matrix <- column_to_rownames(traits_data, var = "Species")
sp_matrix <- column_to_rownames(sp_data, var = "ID")

# Strict ordinal scaling
traits_matrix$Dispersal <- factor(traits_matrix$Dispersal, levels = c(1,2,3), ordered = TRUE)
traits_matrix$Trophic <- factor(traits_matrix$Trophic, levels = c(1,2,3), ordered = TRUE)
# Calculate Gower distance for traits
trait_distance <- daisy(traits_matrix, metric = "gower")

traits_mat <- as.matrix(trait_distance)
sp_mat <- as.matrix(sp_matrix)

rao_results <- fd_raoq(traits = traits_mat, sp = sp_mat)
print(rao_results)

### FD RaoQ calculation in GAM ###
library(mgcv)
df <- read_excel("Rao_diversity.xlsx", sheet = "RaoQ")
df$Altitude_scaled <- as.numeric(scale(df$Elevation, center = TRUE, scale = TRUE))
df$Locality <- as.factor(df$Locality)
df$Year <- as.factor(df$Year)
df$Exposition2 <- as.numeric(df$Exposition2)
df$Exposition2 <- scale(df$Exposition2)

mod_gam_rao <- gam(
  Q ~ 
    s(Locality, bs = "re") + 
    s(Altitude_scaled, bs = "cr", k = 5)  + Exposition2 +
    s(Year, bs = "re"),
  data   = df,
  family = tw(link="log"), select = TRUE,
  method = "REML"
)
summary(mod_gam_rao)
par(mfrow = c(2, 2))
gam.check(mod_gam_rao)
concurvity(mod_gam_rao, full = TRUE)
gratia::draw(mod_gam_rao)
plot(mod_gam_rao, select = 2)

# correlogram (autocorrelation using Moran’s I based on site-averaged Pearson residuals)
library(DHARMa)
library(qgam)
library(mgcViz)
library(dplyr)
library(gstat)
library(sp)
library(spdep)

df$resid <- residuals(mod_gam_rao, type = "pearson")
df_site_res <- df %>%
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

# Vizualization of FD Rao
library(gratia)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1) Terms to exclude 
excl <- c("s(Locality)", "s(Year)")

# 2) Create the grid first
new_data <- gratia::data_slice(mod_gam1, 
                               Altitude_scaled = evenly(Altitude_scaled, n = 100))

# 3) THEN apply the matrix fix to the newly created grid
new_data$Exposition2 <- matrix(new_data$Exposition2, ncol = 1)

# 4) Calculate fitted values (this will now work!)
fv <- fitted_values(mod_gam1, data = new_data, exclude = excl,
                    scale = "response", se = TRUE) %>%
  dplyr::rename(fitted = any_of(c("fitted",".fitted","fit")),
                se     = any_of(c("se",".se"))) %>%
  mutate(lower = fitted - 1.96 * se,
         upper = fitted + 1.96 * se)

# 4) Partial points aligned with exclusions
fv_obs <- fitted_values(mod_gam1, data = df, exclude = excl,
                        scale = "response", se = FALSE) %>%
  dplyr::rename(fitted = any_of(c("fitted",".fitted","fit")))

df$partial_excl <- residuals(mod_gam1, type = "response") + fv_obs$fitted

# 5) Plotting
p <- ggplot() +
  geom_ribbon(data = fv,
              aes(x = Altitude_scaled, ymin = lower, ymax = upper),
              fill = "grey70", alpha = 0.35) +
  geom_line(data = fv,
            aes(x = Altitude_scaled, y = fitted), linewidth = 1.1) +
  geom_jitter(data = df,
              aes(x = Altitude_scaled, y = Q),
              width = 0.03, height = 0, size = 1.8, alpha = 0.6) +
  labs(x = "Elevational gradient (scaled)", y = "Functional Diversity (Rao's Q)") +
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

# Export to TIFF with a new name
tiff('GAM_RaoQ_Elevation.tiff', units = "in", width = 8, height = 10, res = 600)
print(p)
dev.off()


### Calculating CWMs ###
# Load required packages
library(readxl)
library(tibble)
library(FD)
library(writexl)

# 1. Read data 
sp_data <- read_excel("Rao_diversity1.xlsx", sheet = "sp")
traits_data <- read_excel("Rao_diversity1.xlsx", sheet = "traits")
diversity_data <- read_excel("Rao_diversity1.xlsx", sheet = "RaoQ")

# 2. Format row names
traits_matrix <- column_to_rownames(traits_data, var = "Species")
sp_matrix <- column_to_rownames(sp_data, var = "ID")

# 2. Format row names
traits_matrix <- column_to_rownames(traits_data, var = "Species")
sp_matrix <- column_to_rownames(sp_data, var = "ID")

# DO NOT convert your traits to factors! 
# Leaving them as numeric forces R to calculate exact decimal averages.
# str(traits_matrix) already showed us they are numeric, which is perfect.

# 3. Ensure species columns in 'sp' match species rows in 'traits'
common_species <- intersect(colnames(sp_matrix), rownames(traits_matrix))
sp_mat <- as.matrix(sp_matrix[, common_species])
traits_mat <- traits_matrix[common_species, ]

# 4. Calculate the Community Weighted Means (CWM)
cwm_results <- functcomp(traits_mat, sp_mat)

# 5. Export the results
cwm_export <- rownames_to_column(cwm_results, var = "ID")
write_xlsx(cwm_export, "CWM_results_decimals.xlsx")

# Print to verify the decimals are back!
print(head(cwm_results))

# Create a new column with the scaled elevation
diversity_data$Altitude_scaled <- as.numeric(scale(diversity_data$Elevation, center = TRUE, scale = TRUE))

print(diversity_data$Altitude_scaled)
write_xlsx(diversity_data, "diversity_data.xlsx")

# 1. Read your previously saved Excel file
cwm_data <- read_excel("CWM_results_decimals.xlsx")

# 2. Rescale specific columns using their THEORETICAL limits
# We subtract the theoretical minimum (1) and divide by the theoretical range (Max - Min)
cwm_scaled <- cwm_data %>%
  mutate(
    # Dietary (was 1–3): Range is 2
    Dietary_cwm = (Dietary_cwm - 1) / 2,
    
    # Wings (was 1–3): Range is 2
    Wings_cwm = (Wings_cwm - 1) / 2,
    
    # Bioindication (was 1–3): Range is 2
    Bioindication_cwm = (Bioindication_cwm - 1) / 2,
    
    # Moisture (was 1–5): Range is 4
    Moisture_cwm = (Moisture_cwm - 1) / 4,
    
    # Distribution (was 1–5): Range is 4
    Distribution_cwm = (Distribution_cwm - 1) / 4
    
    # Body_size_cwm and ID are not listed here, so they remain safely untouched!
  )

# 3. Verify the changes in the console
# All targeted traits should now be between 0.0 and 1.0, representing their absolute position 
# on the theoretical trait spectrum. Body_size_cwm remains in its original biological scale.
print(head(cwm_scaled))
summary(cwm_scaled)

# 4. Save the theoretically scaled data to a new Excel file
write_xlsx(cwm_scaled, "CWM_results_scaled_01_theoretical.xlsx")
