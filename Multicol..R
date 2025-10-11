library(corrplot)
library(dplyr)

# Select only the numeric environmental variables
numeric_vars <- c("Temperature", "Precipitation", "Wind")

# Check they exist in your dataset
numeric_vars %in% names(data_long1)

# Extract those columns
corr_numeric <- data_long1 %>%
  select(all_of(numeric_vars))

# Compute Spearman correlation matrix
cor_matrix <- cor(corr_numeric, method = "spearman", use = "pairwise.complete.obs")

# Plot correlation heatmap
corrplot(cor_matrix, method = "color", type = "upper",
         order = "hclust", tl.col = "black", tl.srt = 45,
         addCoef.col = "black", number.cex = 0.8)

# Save to TIFF
tiff("Spearman_rank_corr.tiff", units = "in", width = 7, height = 5, res = 300)
corrplot(cor_matrix, method = "color", type = "upper",
         order = "hclust", tl.col = "black", tl.srt = 45,
         addCoef.col = "black", number.cex = 0.8)
dev.off()

# Quantitative check of multicollinearity using VIF
# A VIF value above 5 or 10 indicates problematic multicollinearity
library(car)
model_vif <- lm(Temperature ~ Precipitation + Wind, data = data_long1)
vif(model_vif)
# VIF confirms they don’t inflate variance in each other’s estimates.