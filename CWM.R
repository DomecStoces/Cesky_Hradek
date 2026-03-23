cwm_clean <- read_excel("df.xlsx", sheet = "Sheet1")

# Fourth corner analysis
library(readxl)
library(ade4)
library(dplyr)

# Load all sheets into a list
fourth_corner <- lapply(c("sp", "env", "traits"), function(x) read_excel("fourth_corner.xlsx", sheet = x))

# Assign names to the list elements
names(fourth_corner) <- c("sp", "env", "traits")

# Access individual sheets with the following syntax
dim(fourth_corner$sp)     # Dimensions of 'sp' sheet
dim(fourth_corner$env)    # Dimensions of 'env' sheet
dim(fourth_corner$traits) # Dimensions of 'traits' sheet

# Add a scaled squared altitude column to env
fourth_corner$env <- as.data.frame(lapply(fourth_corner$env, function(x) {
  if (is.character(x)) as.factor(x) else x
}))
fourth_corner$env <- fourth_corner$env %>%
  mutate(
    Altitude_scaled = scale(Altitude, center = TRUE, scale = TRUE)[,1],
    Altitude_scaled2 = Altitude_scaled^2
  )
fourth_corner$env <- fourth_corner$env %>%
  select(-Altitude) %>%   
  relocate(Altitude_scaled, Altitude_scaled2, .after = Wind)

fourth_corner$traits <- as.data.frame(fourth_corner$traits)
fourth_corner$traits$Body.size <- as.numeric(gsub(",", ".", fourth_corner$traits$Body.size))

# Convert numeric categorical columns to factors
fourth_corner$env$Exposition2 <- as.numeric(fourth_corner$env$Exposition2)
fourth_corner$env$Altitude <- as.numeric(fourth_corner$env$Altitude)
fourth_corner$env$Temperature <- as.numeric(fourth_corner$env$Temperature)
fourth_corner$env$Precipitation <- as.numeric(fourth_corner$env$Precipitation)
fourth_corner$env$Wind <- as.numeric(fourth_corner$env$Wind)


fourth_corner$sp[is.na(fourth_corner$sp)] <- 0
fourth_corner$traits[is.na(fourth_corner$traits)] <- 0
fourth_corner$traits_num <- fourth_corner$traits %>%
  select(-Species) %>%
  mutate(across(everything(), as.numeric))
afcL.aravo <- dudi.coa(fourth_corner$sp, scannf = FALSE)
acpR.aravo <- dudi.hillsmith(fourth_corner$env, row.w = afcL.aravo$lw,
                             scannf = FALSE)
acpQ.aravo <- dudi.pca(fourth_corner$traits_num, 
                       row.w = afcL.aravo$cw, 
                       scannf = FALSE)
rlq.aravo <- rlq(acpR.aravo, afcL.aravo, acpQ.aravo,
                 scannf = FALSE)

nrepet <- 999
four.comb.aravo <- fourthcorner(fourth_corner$env, fourth_corner$sp,
                                fourth_corner$traits_num, modeltype = 6,
                                p.adjust.method.G = "none",
                                p.adjust.method.D = "none", nrepet = nrepet)


plot(four.comb.aravo, alpha = 0.05, stat = "D2")
plot(four.comb.aravo, x.rlq = rlq.aravo, alpha = 0.05,
     stat = "D2", type = "biplot")

Srlq <- fourthcorner2(fourth_corner$env, fourth_corner$sp,
                      fourth_corner$traits,
                      modeltype = 6, p.adjust.method.G = "fdr", nrepet = nrepet)
Srlq$trRLQ

tiff('table.tiff', units="in", width=8, height=10, res=600)
plot(four.comb.aravo, alpha = 0.05, stat = "D2")
dev.off()

pdf('table.pdf', width=8, height=10)
plot(four.comb.aravo, alpha = 0.05, stat = "D2")
dev.off()
################################################################################################################
# Load required libraries
library(dplyr)
library(tidyr)
library(scales)
library(lme4)
library(car)
library(DHARMa)
library(ade4)
library(reshape2)
library(tibble)
library(lmtest)

# Community-weighted means (CWMs) were calculated for each trait using species abundance as weights. Ordinal traits (e.g. breeding season, dispersal ability, moisture preference, biogeographic range) were converted to numeric ranks (1–5) reflecting increasing ecological gradients prior to averaging. Continuous traits (body size) were used directly.

data_long1$Locality <- as.character(data_long1$Locality)

data_long1 <- data_long1 %>%
  mutate(
    
    # Dietary (1–3): Max-Min = 2
    Dietary = (dplyr::recode(trimws(Dietary),
                             "Granivor" = 1,
                             "Omnivor"  = 2,
                             "Predator" = 3,
                             .default   = NA_real_) - 1) / 2,
    
    # Breeding (1–2): Max-Min = 1
    Breeding = (dplyr::recode(trimws(Breeding),
                              "Spring" = 1,
                              "Autumn" = 2,
                              .default = NA_real_) - 1) / 1,
    
    # Wings (1–3): Max-Min = 2
    Wings = (dplyr::recode(trimws(Wing.morph),
                           "A"   = 1,
                           "A/B" = 1,
                           "B"   = 1,
                           "B/M" = 2,
                           "M"   = 3,
                           .default = NA_real_) - 1) / 2,
    
    # Bioindication (1–3): Max-Min = 2
    Bioindication = (dplyr::recode(trimws(Bioindication.group),
                                   "E" = 1,
                                   "A" = 2,
                                   "R" = 3,
                                   .default = NA_real_) - 1) / 2,
    
    # Moisture (1–5): Max-Min = 4
    Moisture = (dplyr::recode(trimws(Moisture.tolerance),
                              "X" = 1,
                              "S" = 2,
                              "I" = 3,
                              "V" = 4,
                              "H" = 5,
                              .default = NA_real_) - 1) / 4,
    
    # Distribution (1–5): Max-Min = 4
    Distribution = (as.numeric(dplyr::recode(as.character(trimws(Areal.distribution)),
                                             "South Palearctic" = 1,
                                             "West Palearctic"  = 1,
                                             "Europe"           = 2,
                                             "Central Europe"   = 2,
                                             "Eurasian"         = 3,
                                             "Holarctic"        = 3,
                                             "Palearctic"       = 3,
                                             "Eurosiberian"     = 3,
                                             "North Palearctic" = 4,
                                             "Transpalearctic"  = 4,
                                             "Circumboreal"     = 5,
                                             .default = NA_real_
    )) - 1) / 4,
    
    Body_size = as.numeric(Body.size)
  )

# Calculate Community Weighted Means (CWM)
cwm_results <- data_long1 %>%
  group_by(Year, Month, Locality) %>%
  summarise(
    # Major 3: Vectorized calculation using across()
    across(
      c(Dietary, Breeding, Wings, Bioindication, Moisture, Body_size, Distribution),
      ~ weighted.mean(.x, Abundance, na.rm = TRUE),
      .names = "{.col}_cwm"
    ),
    Total_Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

site_env <- data_long1 %>%
  group_by(Locality) %>%
  summarise(
    Altitude    = mean(Altitude, na.rm = TRUE),
    Exposition2 = first(na.omit(Exposition2)), 
    HR          = first(na.omit(HR)),
    
    .groups = "drop"
  ) %>%
  mutate(
    Altitude_scaled  = as.numeric(scale(Altitude, center = TRUE, scale = TRUE)),
    Altitude_scaled2 = Altitude_scaled^2
  )
cwm_clean <- cwm_results %>%
  left_join(site_env, by = "Locality")

# Correlation among CWMs 
library(corrplot)
cwm_mat3 <- cwm_clean[, c("Wings_cwm", "Body_size_cwm", "Dietary_cwm")]
cor_mat <- cor(
  cwm_clean[, c("Wings_cwm", "Body_size_cwm", "Dietary_cwm")],
  method = "spearman",
  use = "pairwise.complete.obs"
)
colnames(cor_mat) <- rownames(cor_mat) <- c(
  "Dispersal ability",
  "Body size",
  "Trophic strategy"
)
corrplot(
  cor_mat,
  method = "color",
  tl.col = "black",
  tl.cex = 1.2,
  addCoef.col = "black",
  number.cex = 1.2,
  col = colorRampPalette(c("#2166AC","#FFFFFF","#B2182B"))(200)
)

# or PCA
library(FactoMineR)
library(factoextra)
cwm_mat3 <- cwm_clean %>%
  select(
    `Dispersal ability`   = Wings_cwm,
    `Body size` = Body_size_cwm,
    `Trophic strategy`            = Dietary_cwm
  ) %>%
  na.omit()
res.pca3 <- PCA(cwm_mat3, scale.unit = TRUE, graph = FALSE)
fviz_pca_var(
  res.pca3,
  repel  = TRUE,
  col.var = "black"
)

loadings <- res.pca3$var$coord
loadings
# How CWMs jointly respond to environment?
# the overall trait–environment concordance
library(vegan)
rda_cwm <- rda(cwm_mat ~ Altitude + Exposition2, data = cwm_clean)
anova(rda_cwm, permutations = 999)
anova(rda_cwm, by = "axis", permutations = 999)

# Mantel test of two CWMs
# Distances
d_dist  <- dist(cwm_clean$Wings_cwm)
d_wings <- dist(cwm_clean$Body_size_cwm)
d_moist <- dist(cwm_clean$Dietary_cwm)

mantel(d_dist, d_wings, method = "spearman", permutations = 999)
mantel(d_dist, d_moist, method = "spearman", permutations = 999)
mantel(d_wings, d_moist, method = "spearman", permutations = 999)
# Methods: The independence among significant CWMs was tested using a Mantel test (Spearman’s ρ, 999 permutations), which showed no significant correlation between Moisture_cwm and Distribution_cwm (ρ = 0.04, p = 0.17), indicating that the traits describe distinct ecological gradients.

tiff('PCA.tiff', units="in", width=7, height=6, res=600)
fviz_pca_var(
  res.pca3,
  repel  = TRUE,
  col.var = "black"
)
dev.off()

