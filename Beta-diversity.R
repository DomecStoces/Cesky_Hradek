# Number of species for each row # 
library(readxl)
library(writexl)
df <- read_excel("Beskydy_2007_2008_traits_final.xlsx", sheet = "chilo_diplo_iso_compo_names")
df$species_richness <- rowSums(df[, -1] > 0)
write_xlsx(df[, c("ID", "species_richness")], "species_richness.xlsx")

library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# 1.
df <- read_excel("Beskydy_2007_2008_traits_final.xlsx", sheet = "carabids_FD")
compo_names <- read_excel("Beskydy_2007_2008_traits_final.xlsx", sheet = "carabids_compo_names")

# 2. Initial
metadata <- df
metadata$ID <- 1:nrow(metadata)
metadata_matched <- metadata[metadata$ID %in% compo_names$ID, ]

# 3. Prepare raw community matrix
comm_matrix <- as.matrix(compo_names[, -1])
rownames(comm_matrix) <- compo_names$ID

# 4. AGGREGATE DATA
# Combine metadata and species, then sum species counts per locality
df_combined <- bind_cols(metadata_matched, as.data.frame(comm_matrix))

df_agg <- df_combined %>%
  drop_na(Trees, Altitude) %>%
  group_by(Locality, Trees, Altitude) %>%
  summarise(across(all_of(colnames(comm_matrix)), sum), .groups = "drop") %>%
  mutate(
    Altitude_scaled = scale(Altitude),
    Trees = as.factor(Trees),
    Locality = as.factor(Locality)
  )

# 5. Prepare aggregated community matrix and convert to presence/absence
comm_agg <- as.matrix(df_agg %>% select(-Locality, -Trees, -Altitude, -Altitude_scaled))
rownames(comm_agg) <- df_agg$Locality
comm_pa_agg <- ifelse(comm_agg > 0, 1, 0)

# 6. Calculate independent beta-diversity components
dist_jaccard <- designdist(comm_pa_agg, method = "1 - (J / (A + B - J))", terms = "binary")
dist_simpson <- designdist(comm_pa_agg, method = "1 - (J / pmin(A, B))", terms = "binary")
dist_richness <- designdist(comm_pa_agg, method = "1 - (pmin(A, B) / pmax(A, B))", terms = "binary")

# Square-root transform them for PERMDISP/PERMANOVA
dist_jaccard_sqrt  <- sqrt(dist_jaccard)
dist_simpson_sqrt  <- sqrt(dist_simpson)
dist_richness_sqrt <- sqrt(dist_richness)

# 7. PERMDISP: testing variance
disp_jaccard <- betadisper(dist_jaccard_sqrt, df_agg$Trees)
disp_simpson <- betadisper(dist_simpson_sqrt, df_agg$Trees)

print("--- PERMDISP Results ---")
print(permutest(disp_jaccard, permutations = 999))
print(permutest(disp_simpson, permutations = 999))

# 8. PERMANOVA
# Nestedness
perm_jaccard <- adonis2(dist_jaccard_agg_sqrt ~ Trees + Altitude_scaled, 
                               data = df_agg, 
                               permutations = 999,
                               by = "margin")

# Turnover
perm_simpson <- adonis2(dist_simpson_agg_sqrt ~ Trees + Altitude_scaled, 
                               data = df_agg, 
                               permutations = 999,
                               by = "margin")
# Total beta-dievrsity
perm_richness <- adonis2(dist_richness_agg_sqrt ~ Trees + Altitude_scaled, 
                                data = df_agg, 
                                permutations = 999,
                                by = "margin")

print("--- PERMANOVA of Elevational gradient and Trees ---")
print(perm_jaccard)
print(perm_simpson)
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

