library(dplyr)
library(stringr)
library(ape)
library(brms)
library(scales)

# Create pseudo-phylogenetic tree based on taxonomy
tax_df <- data_long1 %>%
  distinct(Species) %>%
  mutate(
    genus = str_extract(Species, "^[^_]+"),
    family = "Carabidae",
    superfamily = "Caraboidea"
  ) %>%
  mutate(across(c(Species, genus, family, superfamily), as.factor))

carab_tree <- as.phylo(~ superfamily / family / genus / Species, data = tax_df)
carab_tree <- multi2di(carab_tree)
carab_tree <- compute.brlen(carab_tree, method = "Grafen")

# Align with dataset
sp_keep <- intersect(unique(data_long1$Species), carab_tree$tip.label)
data_phy <- data_long1 %>% filter(Species %in% sp_keep)
carab_tree <- drop.tip(carab_tree, setdiff(carab_tree$tip.label, sp_keep))

# Covariance matrix
vcv_mat <- vcv(carab_tree, corr = FALSE)
data_phy$Species <- factor(data_phy$Species, levels = rownames(vcv_mat))

# Aggregate by Year × Locality × Species for CWM modelling
cwm_species <- data_phy %>%
  group_by(Year, Locality, Species) %>%
  summarise(
    # Trait CWMs (species-weighted means per Year × Locality × Species)
    Dietary_cwm        = weighted.mean(Dietary, Count, na.rm = TRUE),
    Breeding_cwm       = weighted.mean(Breeding, Count, na.rm = TRUE),
    Wings_cwm          = weighted.mean(Wing.morph, Count, na.rm = TRUE),
    Bioindication_cwm  = weighted.mean(Bioindication.group, Count, na.rm = TRUE),
    Moisture_cwm       = weighted.mean(Moisture.tolerance, Count, na.rm = TRUE),
    Body_cwm           = weighted.mean(Body.size, Count, na.rm = TRUE),
    Distribution_cwm   = weighted.mean(Areal.distribution, Count, na.rm = TRUE),
    
    # Abundance and environment
    Count        = sum(Count),
    Altitude     = mean(Altitude, na.rm = TRUE),
    Exposition2  = mean(Exposition2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Altitude_scaled  = scale(Altitude, center = TRUE, scale = TRUE)[,1],
    Altitude_scaled2 = Altitude_scaled^2
  )

# BRMS model: CWM Moisture with phylogenetic random effects
n <- nrow(cwm_species)
cwm_species <- cwm_species %>%
  mutate(
    Distribution_cwm = (Distribution_cwm * (n - 1) + 0.5) / n
  )
brm_CWM_Distribution <- brm(
  bf(Distribution_cwm ~ poly(Altitude_scaled, 2, raw = TRUE) + Exposition2 +
       (1 | gr(Species, cov = vcv_mat)) + (1 | Locality)),
  data   = cwm_species,
  data2  = list(vcv_mat = vcv_mat),
  family = Beta(link = "logit"),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95),
  seed = 1234
)

brm_CWM_Moisture <- brm(
  bf(Moisture_cwm ~ poly(Altitude_scaled, 2, raw = TRUE) + Exposition2 +
       (1 | gr(Species, cov = vcv_mat)) + (1 | Locality)),
  data   = cwm_species,
  data2  = list(vcv_mat = vcv_mat),
  family = gaussian(link = "identity"),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95),
  seed = 1234
)

brm_CWM_Wings <- brm(
  bf(Wings_cwm ~ poly(Altitude_scaled, 2, raw = TRUE) + Exposition2 +
       (1 | gr(Species, cov = vcv_mat)) + (1 | Locality)),
  data   = cwm_species,
  data2  = list(vcv_mat = vcv_mat),
  family = gaussian(link = "identity"),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95),
  seed = 1234
)