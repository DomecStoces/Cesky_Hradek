# Load necessary packages
library(readxl)
library(dplyr)
library(tidyr)
library(writexl)


# Read your Excel file
data <- read_excel("FINAL_1.xlsx")

# Automatically detect species columns
# (assuming species names start with a capital letter and have underscores)
species_cols <- grep("^[A-Z][a-z]+_", names(data), value = TRUE)

# Check which columns were detected
print(species_cols)

# Reshape from wide to long format
data_long <- data %>%
  pivot_longer(
    cols = all_of(species_cols),
    names_to = "Species",
    values_to = "Abundance"
  ) %>%
  filter(!is.na(Abundance) & Abundance > 0)

# Save as new Excel file
write_xlsx(data_long, "data_long.xlsx")
head(data_long)

# Add traits to dataset
traits1 <- read_excel("traits.xlsx")
final_data <- data_long %>%
  left_join(traits1, by = c("Species" = "Acronym"))
head(final_data)
write_xlsx(final_data, "final_data.xlsx")


# Add environemntal conditions from clima-station Nová Ves v Horách (U1NOVE1)
data_long1 <- data_long1 %>%
  mutate(Locality = as.character(Locality),
         Sequence = as.numeric(Sequence))

env <- env %>%
  mutate(Locality = as.character(Locality),
         Sequence = as.numeric(Sequence))

# Merge (add environmental columns)
data_long1 <- left_join(data_long1, env, by = c("Sequence", "Locality"))
glimpse(data_long1)

env %>%
  count(Sequence, Locality) %>%
  filter(n > 1)
env <- env %>%
  distinct(Sequence, Locality, .keep_all = TRUE)
data_long1 <- data_long1 %>%
  select(-Temperature, -Precipitation, -Wind) %>%
  left_join(env, by = c("Sequence", "Locality"))

write_xlsx(data_long1, "data_long1.xlsx")


