# 1. Packages
library(readxl)
library(dplyr)
library(tidyr)
library(vegan)

# 2. Load data
data_long <- read_excel("data_long.xlsx", sheet = "Sheet1")

# 3. Aggregate data by Locality, Elevation, Exposition, Year, Month, and Species
#    Each row = one Locality x Year x Month community
df_agg <- data_long %>%
  group_by(Locality, Elevation, Exposition2, Year, Month, Species) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(
    names_from = Species, 
    values_from = Abundance, 
    values_fill = list(Abundance = 0)
  )

# 4. Prepare metadata
df_agg <- df_agg %>%
  mutate(
    Locality = as.factor(Locality),
    Exposition2 = as.factor(Exposition2),
    Year = as.factor(Year),
    Month = as.factor(Month),
    Elevation_scaled = as.numeric(scale(Elevation))
  )

# 5. Prepare community matrix (species only)
comm_agg <- df_agg %>% 
  select(-Locality, -Exposition2, -Year, -Month, -Elevation, -Elevation_scaled) %>%
  as.matrix()

# 6. Create unique row names for repeated measures
rownames(comm_agg) <- paste(df_agg$Locality, df_agg$Year, df_agg$Month, sep = "_")

# 7. Convert to presence/absence (binary) for beta-diversity
comm_pa_agg <- ifelse(comm_agg > 0, 1, 0)

# =========================================
# 8. Calculate independent beta-diversity components
# Jaccard (total), Simpson (turnover), Nestedness (richness difference)
# =========================================
dist_jaccard <- designdist(comm_pa_agg, method = "1 - (J / (A + B - J))", terms = "binary")
dist_simpson <- designdist(comm_pa_agg, method = "1 - (J / pmin(A, B))", terms = "binary")
dist_richness <- designdist(comm_pa_agg, method = "1 - (pmin(A, B) / pmax(A, B))", terms = "binary")

# Square-root transform distances (recommended for PERMDISP / PERMANOVA)
dist_jaccard_sqrt  <- sqrt(dist_jaccard)
dist_simpson_sqrt  <- sqrt(dist_simpson)
dist_richness_sqrt <- sqrt(dist_richness)

# =========================================
# 9. PERMDISP: test homogeneity of multivariate dispersions
# =========================================
disp_jaccard <- betadisper(dist_jaccard_sqrt, df_agg$Exposition2)
disp_simpson <- betadisper(dist_simpson_sqrt, df_agg$Exposition2)

cat("\n--- PERMDISP Results (Jaccard) ---\n")
print(permutest(disp_jaccard, permutations = 999))

cat("\n--- PERMDISP Results (Simpson) ---\n")
print(permutest(disp_simpson, permutations = 999))

# =========================================
# 10. PERMANOVA: test effects of Exposition, Elevation, Year, Month
#          with Locality as a stratum (repeated measures)
# =========================================

# Total beta-diversity (Jaccard)
cat("\n--- PERMANOVA: Total Beta-Diversity (Jaccard) ---\n")
perm_jaccard <- adonis2(
  dist_jaccard_sqrt ~ Exposition2 + Elevation_scaled + Year + Month,
  data = df_agg,
  permutations = 999,
  by = "margin",
  strata = df_agg$Locality
)
print(perm_jaccard)

# Turnover (Simpson)
cat("\n--- PERMANOVA: Turnover (Simpson) ---\n")
perm_simpson <- adonis2(
  dist_simpson_sqrt ~ Exposition2 + Elevation_scaled + Year + Month,
  data = df_agg,
  permutations = 999,
  by = "margin",
  strata = df_agg$Locality
)
print(perm_simpson)

# Richness difference (Nestedness)
cat("\n--- PERMANOVA: Richness Difference (Nestedness) ---\n")
perm_richness <- adonis2(
  dist_richness_sqrt ~ Exposition2 + Elevation_scaled + Year + Month,
  data = df_agg,
  permutations = 999,
  by = "margin",
  strata = df_agg$Locality
)
print(perm_richness)

# =========================================
# Notes:
# 1. Using 'strata = Locality' fixes pseudoreplication for repeated monthly/annual samples.
# 2. Year and Month in the formula control for temporal variation (partition variance from repeated measures).
# 3. Significant PERMDISP results indicate differences in dispersion; interpret PERMANOVA with caution.
# =========================================
