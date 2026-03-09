library(readxl)
library(dplyr)
library(fundiversity)
library(tibble)
library(cluster)
data1 <- read_excel("Rao_diversity.xlsx")

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
df$Exposition2<-as.numeric(scale(df$Exposition2))

mod_gam1 <- gam(
  Q ~ 
    s(Locality, bs = "re") + 
    s(Altitude_scaled, bs = "cr", k = 10) + 
    Exposition2 + 
    s(Year, bs = "re"),
  data   = df,
  family = Gamma(link="log"),
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
