library(readxl)
library(dplyr)
library(fundiversity)
library(tibble)
library(cluster)
data1 <- read_excel("Rao_diversity1.xlsx")

sp_data <- read_excel("Rao_diversity1.xlsx", sheet = "sp")
traits_data <- read_excel("Rao_diversity1.xlsx", sheet = "traits")
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
df <- read_excel("Rao_diversity1.xlsx", sheet = "RaoQ")
df$Altitude_scaled <- as.numeric(scale(df$Elevation, center = TRUE, scale = TRUE))
df$Locality <- as.factor(df$Locality)
df$Year <- as.factor(df$Year)
df$Exposition2 <- as.numeric(df$Exposition2)
df$Exposition2 <- scale(df$Exposition2)

mod_gam1 <- gam(
  Q ~ 
    s(Locality, bs = "re") + 
    s(Altitude_scaled, bs = "cr", k = 5) + 
    Exposition2 + 
    s(Year, bs = "re"),
  data   = df,
  family = tw(link="log"), select = TRUE,
  method = "REML"
)
summary(mod_gam1)
par(mfrow = c(2, 2))
gam.check(mod_gam1)
concurvity(mod_gam1, full = TRUE)
gratia::draw(mod_gam1)
plot(mod_gam1, select = 2)

plot(mod_gam1, select = 2, shade = TRUE, residuals = TRUE, 
     pch = 1, cex = 0.5, col = "black",
     ylim = c(-1, 1),
     main = "Effect of Altitude on Rao's Q (Zoomed)")

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
  labs(x = "Elevational gradient (scaled)", y = "Functional Diversity (Rao's Q)") + # Updated label
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
