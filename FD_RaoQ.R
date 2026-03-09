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
    s(Altitude_scaled, bs = "cr", k = 3) + 
    Exposition2 + 
    s(Year, bs = "re"),
  data   = df,
  family = tw(),
  method = "REML"
)
summary(mod_gam1)

# Euclidean-corrected Rao's Q values using the ade4 package
library(ade4)
sp_mat_ade4 <- t(sp_mat)
trait_distance_eucl <- lingoes(trait_distance)
rao_ade4 <- divc(as.data.frame(sp_mat_ade4), trait_distance_eucl)
head(rao_ade4)

df$Q_ade4 <- rao_ade4$diversity
n <- nrow(df)
df$Q_ade4_beta <- (df$Q_ade4 * (n - 1) + 0.5) / n
mod_gam_ade4 <- gam(
  Q_ade4_beta ~ 
    s(Locality, bs = "re") + 
    s(Altitude_scaled, bs = "cr", k = 3) + 
    Exposition2 + 
    s(Year, bs = "re"),
  data   = df,
  family = betar(link = "cloglog"),
  method = "REML"
)
summary(mod_gam_ade4)
