library(dplyr)
library(stringr)
library(ape)
library(phytools)

# tax_df must contain: Species, genus, family, superfamily
tax_df <- data_long1 %>%
  distinct(Species) %>%
  mutate(
    genus = stringr::word(Species, 1, sep = "_"),
    family = "Carabidae",
    superfamily = "Caraboidea"
  ) %>%
  mutate(across(c(Species, genus, family, superfamily), as.factor))

# Check structure
str(tax_df)

# Build hierarchical tree
carab_tree <- as.phylo(~ superfamily / family / genus / Species, data = tax_df)

# Ensure fully dichotomous
carab_tree <- multi2di(carab_tree)

# Add branch lengths (Grafen’s method)
carab_tree <- compute.brlen(carab_tree, method = "Grafen")

sp_keep <- intersect(unique(data_long1$Species), carab_tree$tip.label)

data_phy <- data_long1 %>%
  filter(Species %in% sp_keep)

carab_tree <- drop.tip(carab_tree, setdiff(carab_tree$tip.label, sp_keep))

# Confirm alignment
all(unique(data_phy$Species) %in% carab_tree$tip.label)

library(phyr)

m_pglmm <- pglmm(
  Count ~ 1 + Altitude_scaled + Exposition2,
  data   = data_phy,
  family = "nbinomial",                             
  cov_ranef = list(Species = carab_tree),
  random = ~ 1 | Locality + 1 | Species__@phylo,
  bayes = TRUE,
  REML  = TRUE
)

summary(m_pglmm)

# Variance explained by phylogenetic random effect
ranef_summary(m_pglmm)$var["Species__@phylo"]

# Create tree from taxonomic hierarchy
carab_tree <- as.phylo(~ superfamily / family / genus / Species, data = tax_df)

# Make sure tree is fully dichotomous and has branch lengths
carab_tree <- multi2di(carab_tree)
carab_tree <- compute.brlen(carab_tree, method = "Grafen")

# Optional: visualize and check labels
plot(carab_tree, cex = 0.5)
head(carab_tree$tip.label)

# Align species between data_long1 and tree
sp_keep <- intersect(unique(data_long1$Species), carab_tree$tip.label)
data_phy <- data_long1 %>% filter(Species %in% sp_keep)
carab_tree <- drop.tip(carab_tree, setdiff(carab_tree$tip.label, sp_keep))

# Confirm alignment
all(unique(data_phy$Species) %in% carab_tree$tip.label)

m_pglmm <- pglmm(
  Count ~ 1 + Altitude_scaled + Exposition2,     # fixed effects
  data   = data_phy,
  family = "poisson",                            # or "nbinomial" for overdispersion
  cov_ranef = list(Species = carab_tree),        # tree-based covariance
  random = ~ 1 | Locality + 1 | Species__@phylo, # site RE + species phylo RE
  bayes = TRUE,
  REML  = TRUE
)

summary(m_pglmm)

# Extract variance explained by phylogeny
ranef_summary(m_pglmm)$var["Species__@phylo"]

m_glmm <- pglmm(
  Count ~ 1 + Altitude_scaled + Exposition2,
  data = data_phy,
  family = "poisson",
  random = ~ 1 | Locality + 1 | Species,   # species RE, no phylo link
  bayes = TRUE,
  REML  = TRUE
)

AIC(m_glmm, m_pglmm)

plot(m_pglmm)                     # residuals/fitted check
hist(ranef(m_pglmm)$Species__@phylo)  # distribution of phylogenetic REs

