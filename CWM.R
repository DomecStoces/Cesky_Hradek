data_long1 <- read_excel("data_long1.xlsx")

# Fourth corner analysis
library(readxl)
library(ade4)
library(dplyr)

# Load all sheets into a list
fourth_corner <- lapply(c("sp", "env", "traits"), function(x) read_excel("fourth_corner.xlsx", sheet = x))

# Assign names to the list elements
names(fourth_corner) <- c("sp", "env", "traits")

# Access individual sheets with the following syntax
dim(fourth_corner$sp)     # Dimensions of 'sp' sheet
dim(fourth_corner$env)    # Dimensions of 'env' sheet
dim(fourth_corner$traits) # Dimensions of 'traits' sheet

# Add a scaled squared altitude column to env
fourth_corner$env <- as.data.frame(lapply(fourth_corner$env, function(x) {
  if (is.character(x)) as.factor(x) else x
}))
fourth_corner$env <- fourth_corner$env %>%
  mutate(
    Altitude_scaled = scale(Altitude, center = TRUE, scale = TRUE)[,1],
    Altitude_scaled2 = Altitude_scaled^2
  )
fourth_corner$env <- fourth_corner$env %>%
  select(-Altitude) %>%   
  relocate(Altitude_scaled, Altitude_scaled2, .after = Wind)

fourth_corner$traits <- as.data.frame(fourth_corner$traits)
fourth_corner$traits$Body.size <- as.numeric(gsub(",", ".", fourth_corner$traits$Body.size))

# Convert numeric categorical columns to factors
fourth_corner$env$Exposition2 <- as.numeric(fourth_corner$env$Exposition2)
fourth_corner$env$Altitude <- as.numeric(fourth_corner$env$Altitude)
fourth_corner$env$Temperature <- as.numeric(fourth_corner$env$Temperature)
fourth_corner$env$Precipitation <- as.numeric(fourth_corner$env$Precipitation)
fourth_corner$env$Wind <- as.numeric(fourth_corner$env$Wind)


fourth_corner$sp[is.na(fourth_corner$sp)] <- 0
fourth_corner$traits[is.na(fourth_corner$traits)] <- 0
fourth_corner$traits_num <- fourth_corner$traits %>%
  select(-Species) %>%
  mutate(across(everything(), as.numeric))
afcL.aravo <- dudi.coa(fourth_corner$sp, scannf = FALSE)
acpR.aravo <- dudi.hillsmith(fourth_corner$env, row.w = afcL.aravo$lw,
                             scannf = FALSE)
acpQ.aravo <- dudi.pca(fourth_corner$traits_num, 
                       row.w = afcL.aravo$cw, 
                       scannf = FALSE)
rlq.aravo <- rlq(acpR.aravo, afcL.aravo, acpQ.aravo,
                 scannf = FALSE)

nrepet <- 999
four.comb.aravo <- fourthcorner(fourth_corner$env, fourth_corner$sp,
                                fourth_corner$traits_num, modeltype = 6,
                                p.adjust.method.G = "none",
                                p.adjust.method.D = "none", nrepet = nrepet)


plot(four.comb.aravo, alpha = 0.05, stat = "D2")
plot(four.comb.aravo, x.rlq = rlq.aravo, alpha = 0.05,
     stat = "D2", type = "biplot")

Srlq <- fourthcorner2(fourth_corner$env, fourth_corner$sp,
                      fourth_corner$traits,
                      modeltype = 6, p.adjust.method.G = "fdr", nrepet = nrepet)
Srlq$trRLQ

tiff('table.tiff', units="in", width=8, height=10, res=600)
plot(four.comb.aravo, alpha = 0.05, stat = "D2")
dev.off()

pdf('table.pdf', width=8, height=10)
plot(four.comb.aravo, alpha = 0.05, stat = "D2")
dev.off()
################################################################################################################
# Load required libraries
library(dplyr)
library(scales)
library(lme4)
library(car)
library(DHARMa)
library(ade4)
library(reshape2)
library(tibble)
library(lmtest)

# Community-weighted means (CWMs) were calculated for each trait using species abundance as weights. Ordinal traits (e.g. breeding season, dispersal ability, moisture preference, biogeographic range) were converted to numeric ranks (1–5) reflecting increasing ecological gradients prior to averaging. Continuous traits (body size) were used directly.

data_long1 <- data_long1 %>%
  mutate(
    # Dietary
    Dietary = dplyr::recode(trimws(Dietary),
                            "Predator" = 1,
                            "Granivor" = 2,
                            "Granivore" = 2,
                            "Omnivor" = 3,
                            "Omnivore" = 3,
                            .default = NA_real_
    ),
    
    # Breeding
    Breeding = dplyr::recode(trimws(Breeding),
                             "Spring" = 1,
                             "Autumn" = 2,
                             .default = NA_real_
    ),
    
    # Wings (Wing morphotype)
    Wing.morph = dplyr::recode(trimws(Wing.morph),
                               "A" = 1,
                               "A/B" = 2,
                               "B" = 3,
                               "B/M" = 4,
                               "M" = 5,
                               .default = NA_real_
    ),
    
    # Bioindication group
    Bioindication.group = dplyr::recode(trimws(Bioindication.group),
                                        "E" = 1,
                                        "A" = 2,
                                        "R" = 3,
                                        .default = NA_real_
    ),
    
    # Moisture tolerance
    Moisture.tolerance = dplyr::recode(trimws(Moisture.tolerance),
                                       "X" = 1,
                                       "S" = 2,
                                       "I" = 3,
                                       "V" = 4,
                                       "H" = 5,
                                       .default = NA_real_
    ),
    
    # Areal distribution
    Areal.distribution = dplyr::recode(trimws(Areal.distribution),
                                       "Central Europe"   = 1,
                                       "Europe"           = 2,
                                       "West Palearctic"  = 2,
                                       "South Palearctic" = 2,
                                       "Eurasian"         = 3,
                                       "Eurosiberian"     = 3,
                                       "Palearctic"       = 4,
                                       "North Palearctic" = 4,
                                       "Transpalearctic"  = 4,
                                       "Circumboreal"     = 4,
                                       "Holoarctic"       = 4,
                                       .default = NA_real_
    ),
    
    # Body size
    Body.size = as.numeric(Body.size)
  ) %>%
  mutate(
    Dietary             = scales::rescale(Dietary,             to = c(0, 1)),
    Breeding            = scales::rescale(Breeding,            to = c(0, 1)),
    Wings               = scales::rescale(Wing.morph,          to = c(0, 1)),
    Bioindication.group = scales::rescale(Bioindication.group, to = c(0, 1)),
    Moisture.tolerance  = scales::rescale(Moisture.tolerance,  to = c(0, 1)),
    Areal.distribution  = scales::rescale(Areal.distribution,  to = c(0, 1))
  )

colnames(data_long1)

cwm_results <- data_long1 %>%
  group_by(Year, Locality) %>%
  summarize(
    Dietary_cwm        = weighted.mean(Dietary, Count, na.rm = TRUE),
    Breeding_cwm       = weighted.mean(Breeding, Count, na.rm = TRUE),
    Wings_cwm          = weighted.mean(Wings, Count, na.rm = TRUE),
    Bioindication_cwm  = weighted.mean(Bioindication.group, Count, na.rm = TRUE),
    Moisture_cwm       = weighted.mean(Moisture.tolerance, Count, na.rm = TRUE),
    Body_cwm           = weighted.mean(Body.size, Count, na.rm = TRUE),
    Distribution_cwm   = weighted.mean(Areal.distribution, Count, na.rm = TRUE),
    Abundance          = sum(Count),
    .groups = "drop"
  )
# Display results
print(cwm_results)

# Center and scale Altitude
data_long1 <- data_long1 %>%
  mutate(
    Altitude_scaled  = as.numeric(scale(Altitude, center = TRUE, scale = TRUE)),
    Altitude_scaled2 = Altitude_scaled^2)
env_site <- data_long1 %>%
  group_by(Year, Locality) %>%
  summarise(
    Altitude        = mean(Altitude, na.rm = TRUE),
    Altitude_scaled = mean(Altitude_scaled, na.rm = TRUE),
    Altitude_scaled2 = mean(Altitude_scaled2, na.rm = TRUE),
    Exposition2     = mean(Exposition2, na.rm = TRUE),
    .groups = "drop"
  )
cwm_results <- cwm_results %>%
  left_join(env_site, by = c("Year", "Locality"))
cwm_clean <- cwm_results %>%
  filter(!is.na(Altitude), !is.na(Exposition2))
# Fit Linear Model (LM): Because each locality represented a unique altitudinal step without replication, we used ordinary least squares with heteroscedasticity-consistent (HC3) standard errors rather than mixed-effects or bootstrap model comparison approaches.
# Because each locality represented a unique altitudinal step without replication, we used ordinary least squares with heteroscedasticity-consistent (HC3) standard errors rather than mixed-effects or bootstrap model comparison approaches.
mod1 <- lm(Distribution_cwm ~ poly(Altitude, 2, raw = TRUE) + Exposition2,
                  data = cwm_clean)
res1 <- residuals(mod1)
Anova(mod1,type = "III")
coeftest(mod1, vcov = sandwich::vcovHC(mod1, type = "HC3"))
# Parametric bootstrap implemented via resampling residuals
confint(mod1, method = "boot", nsim = 1999)
# Diagnosis for lm()= assumptions of linear modeling needs to be satisfied, no significant non-normality or variance heterogeneity
par(mfrow=c(2,2)); plot(mod1)              # residual vs fitted, QQ
car::ncvTest(mod1)                         # heteroskedasticity
lmtest::bptest(mod1)                       # Breusch–Pagan = non-constant variance
car::crPlots(mod1)                         # component + residual for shape

# Spatial autocorrelation
library(spdep)

# A) Within each year (no duplicate coords inside a year)
years <- levels(cwm_clean$Year)
out <- lapply(years, function(y){
  dat <- subset(cwm_clean, Year == y)
  coords <- as.matrix(dat[, c("X","Y")])
  nb  <- knn2nb(knearneigh(coords, k = 4))
  lw  <- nb2listw(nb, style = "W")
  moran.test(res1[ cwm_clean$Year == y ], lw)
})
names(out) <- years
out
# B) Aggregate to unique sites (one value per locality)
agg <- cwm_clean |>
  mutate(res1 = res1) |>
  group_by(Locality) |>
  summarise(X = mean(X), Y = mean(Y), res1 = mean(res1), .groups="drop")

coords <- as.matrix(agg[, c("X","Y")])
nb  <- knn2nb(knearneigh(coords, k = 4))
lw  <- nb2listw(nb, style = "W")
moran.test(agg$res1, lw)

# cluster-robust SEs (by Locality) to safeguard against any residual within-site correlation
library(clubSandwich)
library(lmtest)
V <- vcovCR(mod1, cluster = cwm_clean$Locality, type = "CR2")
coef_test(mod1, vcov = V)

# Interpretation
# Distribution: The relationship between the community-weighted mean of species’ areal distribution and altitude was significant and unimodal (HC3-corrected linear term: t = −2.30, p = 0.025; quadratic term: t = 2.46, p = 0.016). This indicates that mid-elevation sites tend to host assemblages dominated by species with broader geographic ranges, whereas both low- and high-elevation sites are characterized by species with more restricted distributions. Exposition had no detectable influence (p = 0.90*)
# the negative association in the linear term indicates that with increasing altitude there are fewer European species followed by an upward curvature suggesting that the most widespread (Palearctic or Holoarctic) taxa occur at the highest elevations.


# Body size: The model detects some altitude-related structure overall, but it’s weak and not robust at the individual term level.
# Although the overall ANOVA indicated a marginal effect of altitude on the community-weighted mean body size (F₂,₇₄ = 3.25, p = 0.044), heteroscedasticity-robust standard errors and bootstrapped confidence intervals showed that neither the linear nor quadratic terms were individually significant (p > 0.5). Therefore, no consistent altitudinal pattern in body size CWM was detected.

# Wings: The community-weighted mean of wing morphology (Wings_cwm) showed a strong nonlinear relationship with altitude (Type III ANOVA: F₂,₇₄ = 22.7, p < 0.001). Both the linear and quadratic terms of altitude remained significant when using heteroskedasticity-robust standard errors (HC3, p < 0.01) and bootstrapped confidence intervals (95% CI: linear −0.0110 to −0.0044; quadratic 3.9×10⁻⁶ to 9.1×10⁻⁶). Although tests indicated non-constant variance (Breusch–Pagan p = 0.022), the effect size and pattern were consistent across robust estimation procedures, supporting a strong altitudinal gradient in wing reduction.

# Correlation among CWMs 
library(corrplot)
corrplot(cor(cwm_results[, c("Wings_cwm", "Moisture_cwm", "Distribution_cwm")],
             use = "complete.obs"), method = "color", tl.col = "black")
mod_cwm <- lm(Altitude ~ Wings_cwm+Moisture_cwm + Distribution_cwm, data = cwm_results)
vif(mod_cwm)

# Correlation matrix or PCA of CWMs
# pairwise correlations or perform a PCA on the CWM matrix to see if traits show a shared gradient across sites:
cwm_mat <- cwm_clean[, c("Moisture_cwm", "Wings_cwm")]
cor(cwm_mat, use = "pairwise.complete.obs", method = "spearman")
corrplot(cor(cwm_mat), method = "color", tl.col = "black")

# or PCA
library(FactoMineR)
res.pca <- PCA(cwm_mat, graph = TRUE)

# How CWMs jointly respond to environment?
# the overall trait–environment concordance
library(vegan)
rda_cwm <- rda(cwm_mat ~ Altitude + Exposition2, data = cwm_clean)
anova(rda_cwm, permutations = 999)
anova(rda_cwm, by = "axis", permutations = 999)

# Mantel test of two CWMs
mantel(
  dist(cwm_clean$Moisture_cwm),
  dist(cwm_clean$Wings_cwm),
  method = "spearman",
  permutations = 999
)
# Methods: The independence among significant CWMs was tested using a Mantel test (Spearman’s ρ, 999 permutations), which showed no significant correlation between Moisture_cwm and Distribution_cwm (ρ = 0.04, p = 0.17), indicating that the traits describe distinct ecological gradients.

tiff('PCA.tiff', units="in", width=8, height=10, res=600)
res.pca <- PCA(cwm_mat, graph = TRUE)
dev.off()

# Graphical representation of the relationship between Altitude and CWMs
library(ggplot2)
ggplot(cwm_clean, aes(Altitude, Wings_cwm)) +
  geom_point() + geom_smooth(method = "gam", formula = y ~ s(x, k = 4))
