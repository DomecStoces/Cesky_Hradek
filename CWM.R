# Fourth corner analysis
library(readxl)
library(ade4)

# Load all sheets into a list
fourth_corner <- lapply(c("sp", "env", "traits"), function(x) read_excel("fourth_corner.xlsx", sheet = x))

# Assign names to the list elements
names(fourth_corner) <- c("sp", "env", "traits")

# Access individual sheets with the following syntax
dim(fourth_corner$sp)     # Dimensions of 'sp' sheet
dim(fourth_corner$env)    # Dimensions of 'env' sheet
dim(fourth_corner$traits) # Dimensions of 'traits' sheet

fourth_corner$env <- as.data.frame(lapply(fourth_corner$env, function(x) {
  if (is.character(x)) as.factor(x) else x
}))
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

# Linear models with CWMs as response variables
library(dplyr)
library(ade4)
library(lme4)
library(glmmTMB)
library(MASS)
library(car)
library(emmeans)

final_dataset$`Time period` <- as.numeric(final_dataset$`Time period`)
final_dataset$Elevation <- as.numeric(final_dataset$Elevation)

group_stats <- final_dataset %>%
  summarise(
    Mean_Number = mean(Number),
    Variance_Number = var(Number),
    Overdispersion = Variance_Number / Mean_Number
  )

group_stats

# Fit a Poisson GLM
poisson_model <- glmer(Number ~ Elevation + Mountain + Wingspan*Dietary+(1|`Time period`)+(1|Species), data = final_dataset,family = poisson(link="log"))

summary(poisson_model)

# Calculate the Pearson chi-square statistic
resid_dev <- sum(residuals(poisson_model, type = "pearson")^2)
resid_df <- df.residual(poisson_model)

overdispersion_ratio <- resid_dev / resid_df
overdispersion_ratio

nb_model <- glmmTMB(
  Number ~ Elevation + Temperature+ Wind +Mountain+ (1 | `Time period`) + (1 | Species),
  data = final_dataset, family = nbinom2(link = "log"))
summary(nb_model)

Anova(nb_model, type = "III")

library(knitr)
tidy_mod <- tidy(nb_model, effects = "fixed", conf.int = TRUE)
kable(tidy_mod, digits = 3, caption = "Model Estimates and 95% Confidence Intervals")

################################################################################################################
#MODEL pro stanovení interakce SpeciesRichness~Movement*Treatment rozdělených do Season!! final!!
species_richness_data <- final_dataset %>% 
  group_by(`Time period`, Elevation, Mountain) %>% 
  summarize(species_richness = n_distinct(Species))

final_dataset2 <- final_dataset %>% 
  left_join(species_richness_data, by = c("Time period", "Elevation", "Mountain"))

mod <- glmmTMB(species_richness ~ Elevation + Mountain + Temperature + Wind + 
                 (1 | `Time period`) + (1 | Species),
               data = final_dataset2, 
               family = nbinom2(link = "log"))

summary(mod)
Anova(mod,type = "III")

library(broom.mixed)
library(ggplot2)

# Tidy the model to get fixed effects with confidence intervals
tidy_mod <- tidy(mod, effects = "fixed", conf.int = TRUE)

# Plot the estimates and their confidence intervals
ggplot(tidy_mod, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "Fixed Effects Estimates with 95% Confidence Intervals",
       x = "", y = "Estimate (log scale)")

library(knitr)
tidy_mod <- tidy(mod, effects = "fixed", conf.int = TRUE)
kable(tidy_mod, digits = 3, caption = "Model Estimates and 95% Confidence Intervals")

################################################################################################################
# Load required libraries
library(dplyr)
library(lme4)
library(car)
library(DHARMa)
library(ade4)
library(reshape2)
library(dplyr)
library(tibble)  # Load tibble for row name conversion

final_dataset <- final_dataset %>%
  mutate(
    Dietary = as.numeric(as.factor(Dietary)),
    Distribution = as.numeric(as.factor(Distribution)),
    `Host species` = as.numeric(as.factor(`Host species`)),
    Overwintering = as.numeric(as.factor(Overwintering)),
    `Leaf action` = as.numeric(as.factor(`Leaf action`))
  )

colnames(final_dataset)

cwm_results <- final_dataset %>%
  group_by(Elevation, Mountain) %>%
  summarize(
    Dietary_cwm = weighted.mean(Dietary, Number, na.rm = TRUE),
    Red_list_cwm = weighted.mean(`Red list species`, Number, na.rm = TRUE),
    Wingspan_cwm = weighted.mean(Wingspan, Number, na.rm = TRUE),
    Distribution_cwm = weighted.mean(Distribution, Number, na.rm = TRUE),
    Host_species_cwm = weighted.mean(`Host species`, Number, na.rm = TRUE),
    Overwintering_cwm = weighted.mean(Overwintering, Number, na.rm = TRUE),
    Leaf_action_cwm = weighted.mean(`Leaf action`, Number, na.rm = TRUE),Abundance=sum(Number)
  )

# Display results
print(cwm_results)

# Save results to a CSV file
write.csv(cwm_results, "cwm_results.csv", row.names = FALSE)

# Convert categorical variables to factors
cwm_results <- cwm_results %>%
  mutate(Mountain = as.factor(Mountain))

# Convert categorical variables to factors
cwm_results <- cwm_results %>%
  mutate(Mountain = as.factor(Mountain),
         `Time period` = as.factor(`Time period`))

# Fit Generalized Linear Mixed Model (GLMM)
mod1 <- lm(Distribution_cwm ~ Elevation + Mountain,
           data = cwm_results)
Anova(mod1,type = "III")

confint(mod1, method = "boot", nsim = 999)

library(vegan)
d <- dist(cwm_results$Distribution_cwm)
mod2 <- adonis2(d ~ Elevation, data = cwm_results,permutations = 999)
print(mod2)


#####
library(pbkrtest)

# Full model (with all predictors)
mod1 <- lm(Dietary_cwm ~ Elevation + Mountain, data = cwm_results)

# Reduced model (without 'Mountain')
mod0 <- lm(Dietary_cwm ~ Mountain, data = cwm_results)

# Run parametric bootstrap test
pb <- PBmodcomp(mod1, mod0, nsim = 999)

# Print results
summary(pb)

isSingular(mod1) #if random are properly estimated and does not collapse to zero variance
VarCorr(mod1) #what variance Time period has
