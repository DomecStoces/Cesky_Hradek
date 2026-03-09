library(ape)
library(picante)
library(readxl)

taxa_df <- read_excel("Phylogenetic_analysis.xlsx", sheet = "Phylogeny")
taxa_df <- as.data.frame(taxa_df)
taxa_df$Family_name <- as.factor(taxa_df$Family_name)
taxa_df$Subfamily <- as.factor(taxa_df$Subfamily)
taxa_df$Tribus <- as.factor(taxa_df$Tribus)
taxa_df$Species <- as.factor(taxa_df$Species)
taxa_df$Root <- as.factor("Coleoptera")
taxa_tree <- as.phylo(~Root/Family_name/Subfamily/Tribus/Species, data = taxa_df)

# Build the tree structure using the taxonomic hierarchy
taxa_tree <- compute.brlen(taxa_tree)
plot(taxa_tree, cex = 0.5, no.margin = TRUE)

# Format the species abundance matrix
sp_data <- as.data.frame(read_excel("Phylogenetic_analysis.xlsx", sheet = "sp"))
rownames(sp_data) <- sp_data$ID
sp_data$ID <- NULL

# Format the elevation data
env_data <- as.data.frame(read_excel("Phylogenetic_analysis.xlsx", sheet = "Sites_years"))

# Match the tree and the community data
combined_data <- match.phylo.comm(taxa_tree, sp_data)

# 2. Calculate SES.pd
# 'taxa.labels' shuffles species names across the tree, preserving tree topology
# 999 runs is the standard for publication-quality p-values
set.seed(123)
ses_results <- ses.pd(samp = combined_data$comm, 
                      tree = combined_data$phy, 
                      null.model = "taxa.labels", 
                      runs = 999,include.root = FALSE)
head(ses_results)

final_df <- merge(ses_results, env_data, by = "row.names")
colnames(final_df)[1] <- "ID"
dropped_sites <- final_df[is.na(final_df$pd.obs.z), ]

final_df$Locality <- as.factor(final_df$Locality)
final_df$Year <- as.factor(final_df$Year)
final_df$Altitude_scaled <- as.numeric(scale(final_df$Elevation, center = TRUE, scale = TRUE))
final_df$Exposition2<-as.numeric(scale(final_df$Exposition2))

library(mgcv)
mod_gam1 <- gam(
  pd.obs.z ~ 
    s(Locality, bs = "re") + 
    s(Altitude_scaled, bs = "cr", k = 3) + 
    Exposition2 + 
    s(Year, bs = "re"),
  data   = final_df,
  family = gaussian(),
  method = "REML"
)
summary(mod_gam1)
