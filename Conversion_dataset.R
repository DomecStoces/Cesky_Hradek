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
    values_to = "Count"
  ) %>%
  filter(!is.na(Count) & Count > 0)  # optional cleanup

# Save as new Excel file
write_xlsx(data_long, "data_long.xlsx")

# Preview result
head(data_long)
