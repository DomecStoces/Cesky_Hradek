library(dplyr)

# Delete from Trees - beech and spruce, keeping birch
data_long1 <- data_long %>%
  filter(Trees == "birch")
unique(data_long1$Trees)

# Calculate species richness 
richness_summary <-data_long1 %>%
  group_by(Sequence, Locality) %>%
  summarise(
    Species_richness = n_distinct(Species[Count > 0]),
    .groups = "drop"
  )
data_long1 <- data_long1 %>%
  left_join(richness_summary, by = c("Sequence", "Locality"))

head(data_long1)

data_long1$Locality <- as.factor(data_long1$Locality)
data_long1$Exposition <- as.factor(data_long1$Exposition)
data_long1$Species_richness <- as.numeric(data_long1$Species_richness)

library(ggplot2)

ggplot(data_long1, aes(x = Altitude, y = Species_richness)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5),
              color = "black", fill = "black", alpha = 0.2, se = TRUE) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Altitude (m a.s.l.)",
    y = "Species richness")

ggplot(data_long1, aes(x = Altitude, y = Abundance)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5),
              color = "black", fill = "black", alpha = 0.2, se = TRUE) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Altitude (m a.s.l.)",
    y = "Abundance")
