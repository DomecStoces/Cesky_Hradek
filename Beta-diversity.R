library(readxl)
library(writexl)
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1.
data_long <- read_excel("data_long.xlsx", sheet = "Sheet1")
# 2. Aggregate and reshape from Long to Wide format
df_agg <- data_long %>%
  # IMPORTANT: I removed Year and Month here! 
  # This pools all 4 years of data into one single community per Locality/Elevation.
  group_by(Locality, Elevation, Exposition2, Species) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = list(Abundance = 0))

# 3. Prepare Metadata & Scale Continuous Variables
df_agg <- df_agg %>%
  mutate(
    Locality = as.factor(Locality),
    # Scale Elevation
    across(
      .cols = c(Elevation), 
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_scaled"
    )
  )
# 4. Prepare Community Matrix
# Drop metadata columns to leave only the species columns
comm_agg <- df_agg %>% 
  select(-Locality, -Elevation, -Exposition2, -Elevation_scaled)
comm_agg <- as.matrix(comm_agg)
rownames(comm_agg) <- df_agg$Locality

# Convert to presence/absence (binary)
comm_pa_agg <- ifelse(comm_agg > 0, 1, 0)

# 5. Calculate independent beta-diversity components using designdist
dist_jaccard <- designdist(comm_pa_agg, method = "1 - (J / (A + B - J))", terms = "binary")
dist_simpson <- designdist(comm_pa_agg, method = "1 - (J / pmin(A, B))", terms = "binary")
dist_richness <- designdist(comm_pa_agg, method = "1 - (pmin(A, B) / pmax(A, B))", terms = "binary")

# Square-root transform them for PERMDISP/PERMANOVA
dist_jaccard_sqrt  <- sqrt(dist_jaccard)
dist_simpson_sqrt  <- sqrt(dist_simpson)
dist_richness_sqrt <- sqrt(dist_richness)

# 6. PERMDISP: Testing multivariate dispersion (variance)
# Using Exposition as the categorical grouping factor
disp_jaccard <- betadisper(dist_jaccard_sqrt, df_agg$Exposition)
disp_simpson <- betadisper(dist_simpson_sqrt, df_agg$Exposition)

print("--- PERMDISP Results ---")
print(permutest(disp_jaccard, permutations = 999))
print(permutest(disp_simpson, permutations = 999))

# 7. PERMANOVA
# Total beta-diversity (Jaccard)
perm_jaccard <- adonis2(dist_jaccard_sqrt ~ Exposition2 + Elevation_scaled, 
                        data = df_agg, 
                        permutations = 999,
                        by = "margin")

# Turnover (Simpson)
perm_simpson <- adonis2(dist_simpson_sqrt ~ Exposition2 + Elevation_scaled, 
                        data = df_agg, 
                        permutations = 999,
                        by = "margin")

# Richness difference (Nestedness-resultant)
perm_richness <- adonis2(dist_richness_sqrt ~ Exposition2 + Elevation_scaled, 
                         data = df_agg, 
                         permutations = 999,
                         by = "margin")

print("--- PERMANOVA of Elevational gradient and Exposition ---")
print("Total Beta-Diversity (Jaccard):")
print(perm_jaccard)

print("Turnover (Simpson):")
print(perm_simpson)

print("Richness Differences:")
print(perm_richness)

cor.test(df_agg$Altitude, rowSums(comm_pa_agg))

# 1. Add species richness (number of species) directly into aggregated dataframe
df_agg$Richness <- rowSums(comm_pa_agg)
# 2. Create the linear regression plot
plot_richness <- ggplot(df_agg, aes(x = Altitude_scaled, y = Richness)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", fill = "grey60", alpha = 0.3) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(title = "Spider",
       x = "Elevational gradient (scaled)",
       y = "Species richness")
print(plot_richness)

