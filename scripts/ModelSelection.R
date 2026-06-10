#### Fish/bird trophic effects on methane: manuscript 
#### Target journal: PNAS 
#### Model selection, main manuscript body

# Load packages -----------------------------------------------------------
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(performance)
library(DHARMa)
library(emmeans)
library(ggeffects)
library(DescTools)


# Load data ---------------------------------------------------------------
daphnia <- read.csv("data/DaphniaFull.csv") # Full daphnia dataset
totaldaphnia1 <- read.csv("data/totaldaphnia1.csv") # Total daphnia dataset year 1
totaldaphnia2 <- read.csv("data/totaldaphnia2.csv") # Total daphnia dataset year 2
isotope <- read.csv("data/isotope.csv")
isotope_yr1 <- read.csv("data/isotope_yr1.csv") # Full porewater dataset year 1
isotope_yr2 <- read.csv("data/isotope_yr2.csv") # Full porewater dataset year 2
GHG <- read.csv("data/GHG_ALL.csv") # full GHG dataset 
GHG_Yr1 <- read.csv("data/GHG_Yr1_winter.csv") # Winter flooding methane dataset year 1
GHG_Yr2 <- read.csv("data/GHG_Yr2_winter.csv") # Winter flooding methane dataset year 2
hobo_summary <- read.csv("data/hobo_summary.csv") # Full water temperature/dissolved oxygen dataset 
depth <- read.csv("data/depth_summary.csv") # Full depth dataset 

# Model set 1: Daphnia Abundance ------------------------------------------

## Iterative framework: beginning with most complex/biologically relevant model and iteratively simplifying 

#### Year 1: Full dataset ####

# For all models, we employ GLMMs with negative binomial distributions and log link functions, and include water temperature, dissolved oxygen and depth as covariates. Plot ID and sampling occassion are included as random effects for all models. We chose a mixed modeling framework to explicitly account for variation in repeated sampling, GLMM because residuals were not normally distributed, and a negative binomial distribution because: a) we are working with count data and b) overdispersion was detected when employing a poisson distribution, which the negative binomal distribution is able to handle more. We test diagnostics for the mose=t complex model and the top model: 

# Most complex model: 
DaphniaYr1_M1 <- glmmTMB(
  DensityRounded ~ FishTreatment * BirdTreatment * Pre_Post +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaYr1_M1)

# Check model diagnostics: 
# Assumptions: Linearity, normality of residuals:
simOutput_DaphniaYr1 <- simulateResiduals(fittedModel = DaphniaYr1_M1)
plot(simOutput_DaphniaYr1)
check_collinearity(DaphniaYr1_M1)
check_overdispersion(DaphniaYr1_M1)

# Drop 3-way interaction: 
DaphniaYr1_M2 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + (BirdTreatment * Pre_Post) + (FishTreatment*BirdTreatment) +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaYr1_M2)

anova(DaphniaYr1_M1, DaphniaYr1_M2) # chi-squared test not significant

# Drop bird/fish interaction:
DaphniaYr1_M3 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + (BirdTreatment * Pre_Post) +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaYr1_M3)

anova(DaphniaYr1_M2, DaphniaYr1_M3) # chi-squared test not significant 

# Drop bird/PrePost interaction: 
DaphniaYr1_M4 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + BirdTreatment +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaYr1_M4)

anova(DaphniaYr1_M3, DaphniaYr1_M4) # chi-squared test not significant 

# Drop bird treatment: 
DaphniaYr1_M5 <- glmmTMB(
  DensityRounded ~ FishTreatment * Pre_Post + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaYr1_M5) 

anova(DaphniaYr1_M4, DaphniaYr1_M5) # chi-squared test not significant 

# Drop fish/PrePost interaction: 
DaphniaYr1_M6 <- glmmTMB(
  DensityRounded ~ FishTreatment + Pre_Post + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaYr1_M6) 

anova(DaphniaYr1_M5, DaphniaYr1_M6) # chi-squared test marginally significant

# Drop PrePost:
DaphniaYr1_M7 <- glmmTMB(
  DensityRounded ~ FishTreatment + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaYr1_M7)

anova(DaphniaYr1_M6, DaphniaYr1_M7) # chi-squared test not significant   

# Null model: 
DaphniaYr1_M8 <- glmmTMB(
  DensityRounded ~ 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaYr1_M8)

anova(DaphniaYr1_M7, DaphniaYr1_M8) # chi-squared test not significant 

## Model diagnostics for top model, M5: 
simOutput_DaphniaYr1 <- simulateResiduals(fittedModel = DaphniaYr1_M5)
plot(simOutput_DaphniaYr1) # diagnostics are better than complex model
check_collinearity(DaphniaYr1_M5)
check_overdispersion(DaphniaYr1_M5)

#### Year 1: Omitting final sampling date ####
# Remove last sampling date 
totaldaphnia1_filtered <- totaldaphnia1 %>%
  filter(SampleDate != as.Date("2024-02-05"))

# Full model
DaphniaYr1a_M1 <- glmmTMB(
  DensityRounded ~ FishTreatment * BirdTreatment * Pre_Post +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(DaphniaYr1a_M1)

# Model diagnostics for most complex model: 
simOutput_DaphniaYr1a <- simulateResiduals(fittedModel = DaphniaYr1a_M1)
plot(simOutput_DaphniaYr1a) # diagnostics are better than full model 
check_collinearity(DaphniaYr1a_M1)
check_overdispersion(DaphniaYr1a_M1)

# Drop three-way interaction 
DaphniaYr1a_M2 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + (BirdTreatment * Pre_Post) + (FishTreatment * BirdTreatment) +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(DaphniaYr1a_M2)

anova(DaphniaYr1a_M1, DaphniaYr1a_M2) # chi-squared test marginally significant

# Drop fish*bird interaction 
DaphniaYr1a_M3 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + (BirdTreatment * Pre_Post) +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(DaphniaYr1a_M3)

anova(DaphniaYr1a_M2, DaphniaYr1a_M3) # chi-squared test not significant 

# Drop bird/pre-post interaction 
DaphniaYr1a_M4 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + BirdTreatment +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(DaphniaYr1a_M4)

anova(DaphniaYr1a_M3, DaphniaYr1a_M4) # chi-squared test not significant 

# Drop bird treatment 
DaphniaYr1a_M5 <- glmmTMB(
  DensityRounded ~ FishTreatment * Pre_Post + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID )  + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(DaphniaYr1a_M5)

anova(DaphniaYr1a_M4, DaphniaYr1a_M5) # chi-squared test not significant 

# Drop fish/pre-post interaction 
DaphniaYr1a_M6 <- glmmTMB(
  DensityRounded ~ FishTreatment + Pre_Post + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(DaphniaYr1a_M6)

anova(DaphniaYr1a_M5, DaphniaYr1a_M6) # chi-squared test significant, p = 0.02

# Drop pre/post
DaphniaYr1a_M7 <- glmmTMB(
  DensityRounded ~ FishTreatment + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(DaphniaYr1a_M7)

anova(DaphniaYr1a_M6, DaphniaYr1a_M7) # chi-squared test not significant 

# Null model
DaphniaYr1a_M8 <- glmmTMB(
  DensityRounded ~ 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(DaphniaYr1a_M8)

anova(DaphniaYr1a_M7, DaphniaYr1a_M8)

## Model diagnostics for model 5: 
simOutput_DaphniaFiltered <- simulateResiduals(fittedModel = DaphniaYr1a_M5)
plot(simOutput_DaphniaFiltered)
check_collinearity(DaphniaYr1a_M5)
check_overdispersion(DaphniaYr1a_M5)

#### Year 2: Full dataset ####
# Most complex model: 
DaphniaYr2_M1<- glmmTMB(
  DensityRounded ~ FishTreatment * BirdTreatment * Pre_Post +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia2
)
summary(DaphniaYr2_M1)

## Model diagnostics for most complex model: 
simOutput_Daphnia2 <- simulateResiduals(fittedModel = DaphniaYr2_M1)
plot(simOutput_Daphnia2)
check_collinearity(DaphniaYr2_M1)
check_overdispersion(DaphniaYr2_M1)

# Drop 3-way interaction: 
DaphniaYr2_M2 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + (BirdTreatment * Pre_Post) + (FishTreatment*BirdTreatment) +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia2
)

summary(DaphniaYr2_M2)

anova(DaphniaYr2_M1, DaphniaYr2_M2) # chi-squared test is significant, p = 0.03

# Drop bird/fish interaction:
DaphniaYr2_M3 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + (BirdTreatment * Pre_Post) +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia2
)
summary(DaphniaYr2_M3)

anova(DaphniaYr2_M2, DaphniaYr2_M3) # chi-squared test not significant 

# Drop bird/PrePost interaction: 
DaphniaYr2_M4 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + BirdTreatment +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia2
)
summary(DaphniaYr2_M4)

anova(DaphniaYr2_M3, DaphniaYr2_M4) # chi-squared test not significant 

# Drop bird treatment: 
DaphniaYr2_M5 <- glmmTMB(
  DensityRounded ~ FishTreatment * Pre_Post + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia2
)
summary(DaphniaYr2_M5) 

anova(DaphniaYr2_M4, DaphniaYr2_M5) # chi-squared test not significant 

# Drop fish/PrePost interaction: 
DaphniaYr2_M6 <- glmmTMB(
  DensityRounded ~ FishTreatment + Pre_Post + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia2
)

summary(DaphniaYr2_M6) 

anova(DaphniaYr2_M5, DaphniaYr2_M6) # chi-squared test marginally significant, p = 0.059

# Drop PrePost:
DaphniaYr2_M7 <- glmmTMB(
  DensityRounded ~ FishTreatment + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia2
)
summary(DaphniaYr2_M7)

anova(DaphniaYr2_M6, DaphniaYr2_M7) # chi-squared test not significant   

# Null model: 
DaphniaYr2_M8 <- glmmTMB(
  DensityRounded ~ 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia2
)
summary(DaphniaYr2_M8)

anova(DaphniaYr2_M7, DaphniaYr2_M8) # chi-squared test not significant 


#### Comparing Study Years ####
daphnia$DensityRounded <- round(daphnia$total) 

M_year <- glmmTMB(DensityRounded ~ StudyYear + 
                    (1 | PlotID),
                  family = nbinom2, 
                  data = daphnia)
summary(M_year)

sim_M_year <- simulateResiduals(fittedModel = M_year)
plot(sim_M_year)

emm_year <- emmeans(M_year, ~ StudyYear)
pairs(emm_year, type = "response", reverse = TRUE)



# Model set 2: MOB estimation ---------------------------------------------

#### Year 1 ####
# subset out for when fish were on fields: 
isotope_yr1_fish <- isotope_yr1 %>% 
  filter(FishOnField == "Y")

# Remove NA Plot IDs from dataset:
isotope_yr1_fish <- isotope_yr1_fish %>% 
  drop_na(C13PW)

# Log-transform response variable, with a slight offset to avoid producing NaNs
isotope_yr1_fish <- isotope_yr1_fish %>% 
  mutate(C13_log = log(C13PW - min(C13PW) + 1))

# Most complex model: fish/bird interaction:
porewaterYr1_M1 <- lm(C13_log ~ (FishTreatment * BirdTreatment), 
                   data = isotope_yr1_fish)

summary(porewaterYr1_M1)

# Diagnostics for most complex model, assumptions: linearity, independence
par(mfrow = c(2, 2))
plot(porewaterYr1_M1)
# assumption: normality 
shapiro.test(resid(porewaterYr1_M1))
# assumption: homoscedasticity
check_heteroscedasticity(porewaterYr1_M1)

# Drop fish/bird interaction: 
porewaterYr1_M2 <- lm(C13_log ~ (FishTreatment + BirdTreatment), 
                   data = isotope_yr1_fish)
summary(porewaterYr1_M2)

anova(porewaterYr1_M1, porewaterYr1_M2) # chi-squared test not significant

# Drop fish treatment 
porewaterYr1_M3 <- lm(C13_log ~ BirdTreatment, 
                   data = isotope_yr1_fish)
summary(porewaterYr1_M3)

anova(porewaterYr1_M2, porewaterYr1_M3) # chi-squared test marginally significant (p = 0.09)

# Drop bird treatment 
porewaterYr1_M4 <- lm(C13_log ~ FishTreatment, 
                   data = isotope_yr1_fish)
summary(porewaterYr1_M4)

# Null model, compare fish and bird only models to this: 
porewaterYr1_M5 <- lm(C13_log ~ 1, 
                   data = isotope_yr1_fish)
anova(porewaterYr1_M3, porewaterYr1_M5)
anova(porewaterYr1_M4, porewaterYr1_M5)

# Diagnostics for model 4:
par(mfrow = c(2, 2))
plot(porewaterYr1_M4)
# assumption: normality 
shapiro.test(resid(porewaterYr1_M4))
# assumption: homoscedasticity
check_heteroscedasticity(porewaterYr1_M4)


#### Year 2 ####
# subset out for when fish were on fields: 
isotope_yr2_fish <- isotope_yr2 %>% 
  filter(FishOnField == "Y")

# Remove NA Plot IDs from dataset:
isotope_yr2_fish <- isotope_yr2_fish %>% 
  drop_na(C13PW)

# Log-transform response variable, with a slight offset to avoid producing NaNs
isotope_yr2_fish <- isotope_yr2_fish %>% 
  mutate(C13_log = log(C13PW - min(C13PW) + 1))

# Most complex model: fish/bird interaction:
porewaterYr2_M1 <- lm(C13_log ~ (FishTreatment * BirdTreatment), 
                       data = isotope_yr2_fish)

summary(porewaterYr2_M1)

# Diagnostics for most complex model:
par(mfrow = c(2, 2))
plot(porewaterYr2_M1)
# assumption: normality 
shapiro.test(resid(porewaterYr2_M1))
# assumption: homoscedasticity
check_heteroscedasticity(porewaterYr2_M1)

# Drop fish/bird interaction: 
porewaterYr2_M2 <- lm(C13_log ~ (FishTreatment + BirdTreatment), 
                       data = isotope_yr2_fish)
summary(porewaterYr2_M2)

anova(porewaterYr2_M1, porewaterYr2_M2) # chi-squared test not significant

# Drop fish treatment 
porewaterYr2_M3 <- lm(C13_log ~ BirdTreatment, 
                       data = isotope_yr2_fish)
summary(porewaterYr2_M3)

anova(porewaterYr2_M2, porewaterYr2_M3) # chi-squared test not significant

# Drop bird treatment 
porewaterYr2_M4 <- lm(C13_log ~ FishTreatment, 
                       data = isotope_yr2_fish)
summary(porewaterYr2_M4)

# Null model, compare fish and bird only models to this: 
porewaterYr2_M5 <- lm(C13_log ~ 1, 
                       data = isotope_yr2_fish)
summary(porewaterYr2_M5)

anova(porewaterYr2_M3, porewaterYr2_M5) # chi-squared test not significant 
anova(porewaterYr2_M4, porewaterYr2_M5) # chi-squared test not significant 

# Diagnostics for null model:
par(mfrow = c(2, 2))
plot(porewaterYr2_M5)
# assumption: normality 
shapiro.test(resid(porewaterYr2_M5))
# assumption: homoscedasticity
check_heteroscedasticity(porewaterYr2_M5)

#### Comparing Study Years ####
isotope <- isotope %>% 
  filter(FishOnField == "Y") %>% 
  drop_na(C13PW) %>% 
  mutate(C13_log = log(C13PW - min(C13PW) + 1))

M_isotope_year <- lm(C13_log ~ FishTreatment*StudyYear,
                  data = isotope)
summary(M_isotope_year)

shapiro.test(resid(M_isotope_year))
par(mfrow = c(2, 2))
plot(M_isotope_year)

emm_isotope_year <- emmeans(M_isotope_year, ~ FishTreatment * StudyYear)
pairs(emm_isotope_year, type = "response", reverse = TRUE)

# Model set 3: Methane Flux -----------------------------------------------

#### Year 1: Winter flux ####
# Only two data points per plot pre- fish introduction. With low sample size, we will first establish no differences between treatments pre-fish introduction, and then do model selection on dataset post-fish introduction. 

# Establish no difference pre-fish intro:
# Create pre/post datasets: 
GHG_Yr1_pre <- GHG_Yr1 %>% 
  filter(Pre_Post == "Pre")

GHG_Yr1_post <- GHG_Yr1 %>% 
  filter(Pre_Post == "Post")

### Run ANOVA to ensure no differences pre- fish introduction: 
anova_Yr1 <- aov(logFlux ~ FishTreatment * BirdTreatment, data = GHG_Yr1_pre)
summary(anova_Yr1) # no significance by treatment 

# ANOVA diagnostics: 
par(mfrow = c(2, 2))
plot(anova_Yr1)
# assumption: normality 
shapiro.test(resid(anova_Yr1))
# assumption: homoscedasticity
check_heteroscedasticity(anova_Yr1)

### Model selection on post-fish introduction data:

# most complex model:
methaneYr1_M1 <-  lmer(logFlux ~ FishTreatment * BirdTreatment +
                      scale(MeanTemp) + scale(PO4) + scale(NH4) + scale(MeanDepth) +
                      (1 | PlotID) + (1 | SamplingOccasion),
                    data = GHG_Yr1_post)
summary(methaneYr1_M1)

### Diagnostics for top model: 
sim_M1 <- simulateResiduals(fittedModel = methaneYr1_M1)
plot(sim_M1)
check_collinearity(methaneYr1_M1)
check_heteroscedasticity(methaneYr1_M1)
shapiro.test(resid(methaneYr1_M1))

# drop treatment interaction: 
methaneYr1_M2 <- lmer(logFlux ~ FishTreatment + BirdTreatment +
                     scale(MeanTemp) + scale(PO4) + scale(NH4) + scale(MeanDepth) +
                     (1 | PlotID) + (1 | SamplingOccasion),
                   data = GHG_Yr1_post)

summary(methaneYr1_M2)

anova(methaneYr1_M1, methaneYr1_M2) # chi-squared test not significant 

# drop fish treatment: 
methaneYr1_M3 <- lmer(logFlux ~ BirdTreatment + 
                        scale(MeanTemp) + scale(PO4) + scale(NH4) + scale(MeanDepth) +
                        (1 | PlotID) + (1 | SamplingOccasion),
                      data = GHG_Yr1_post)

summary(methaneYr1_M3)

anova(methaneYr1_M2, methaneYr1_M3) # chi-squared test significant when fish treatment dropped

# drop bird treatment: 
methaneYr1_M4 <- lmer(logFlux ~ FishTreatment + 
                     scale(MeanTemp) + scale(PO4) + scale(NH4) + scale(MeanDepth) +
                     (1 | PlotID) + (1 | SamplingOccasion),
                   data = GHG_Yr1_post)

summary(methaneYr1_M4)


# null model: 
methaneYr1_M5 <- lmer(logFlux ~ 1 +
                     scale(MeanTemp) + scale(PO4) + scale(NH4) + scale(MeanDepth) +
                     (1 | PlotID) + (1 | SamplingOccasion),
                   data = GHG_Yr1_post)

summary(methaneYr1_M5)

anova(methaneYr1_M3, methaneYr1_M5) # chi-squared test not significant
anova(methaneYr1_M4, methaneYr1_M5) # chi-squared test significant, p = 0.002

### Diagnostics for top model: 
sim_M4 <- simulateResiduals(fittedModel = methaneYr1_M4)
plot(sim_M4)
check_collinearity(methaneYr1_M4)
check_heteroscedasticity(methaneYr1_M4)
shapiro.test(resid(methaneYr1_M4))

#### Year 2: Winter flux ####
# We follow the same model selection procedure as in year 1: 
# Create pre/post datasets: 
GHG_Yr2_pre <- GHG_Yr2 %>% 
  filter(Pre_Post == "Pre")

GHG_Yr2_post <- GHG_Yr2 %>% 
  filter(Pre_Post == "Post")

### Run ANOVA to ensure no differences pre- fish introduction: 
anova_Yr2 <- aov(logFlux ~ FishTreatment*BirdTreatment, data = GHG_Yr2_pre)
summary(anova_Yr2)

# ANOVA diagnostics: 
par(mfrow = c(2, 2))
plot(anova_Yr2)
# assumption: normality 
shapiro.test(resid(anova_Yr2))
# assumption: homoscedasticity
check_heteroscedasticity(anova_Yr2)

### Model selection on post-fish introduction data:

# most complex model:
methaneYr2_M1 <-  lmer(logFlux ~ FishTreatment * BirdTreatment +
                         scale(MeanTemp) + scale(PO4) + scale(NH4) + scale(MeanDepth) +
                      (1 | PlotID) + (1 | SamplingOccasion),
                    data = GHG_Yr2_post)
summary(methaneYr2_M1)

### Diagnostics for most complex model: 
sim_M1 <- simulateResiduals(fittedModel = methaneYr2_M1)
plot(sim_M1)
check_collinearity(methaneYr2_M1)
check_heteroscedasticity(methaneYr2_M1)
shapiro.test(resid(methaneYr2_M1))

# drop treatment interaction: 
methaneYr2_M2 <- lmer(logFlux ~ FishTreatment + BirdTreatment +
                        scale(MeanTemp) + scale(PO4) + scale(NH4) + scale(MeanDepth) +
                     (1 | PlotID) + (1 | SamplingOccasion),
                   data = GHG_Yr2_post)

summary(methaneYr2_M2)

anova(methaneYr2_M1, methaneYr2_M2) # chi-squared test not significant 

# drop fish treatment: 
methaneYr2_M3 <- lmer(logFlux ~ BirdTreatment + 
                        scale(MeanTemp) + scale(PO4) + scale(NH4) + scale(MeanDepth) +
                     (1 | PlotID) + (1 | SamplingOccasion),
                   data = GHG_Yr2_post)

summary(methaneYr2_M3)

anova(methaneYr2_M2, methaneYr2_M3) # chi-squared test not significant 

# drop bird treatment: 
methaneYr2_M4 <- lmer(logFlux ~ FishTreatment + 
                        scale(MeanTemp) + scale(PO4) + scale(NH4) + scale(MeanDepth) +
                     (1 | PlotID) + (1 | SamplingOccasion),
                   data = GHG_Yr2_post)

summary(methaneYr2_M4)

# null model: 
methaneYr2_M5 <- lmer(logFlux ~ 1 +
                     scale(MeanTemp) + scale(MeanDepth) + scale(PO4) + scale(NH4) +
                     (1 | PlotID) + (1 | SamplingOccasion),
                   data = GHG_Yr2_post)

summary(methaneYr2_M5)

anova(methaneYr2_M3, methaneYr2_M5) # chi-squared test not significant
anova(methaneYr2_M4, methaneYr2_M5) # chi-squared test not significant

### Diagnostics for null model: 
sim_M5 <- simulateResiduals(fittedModel = methaneYr2_M5)
plot(sim_M5)
check_collinearity(methaneYr2_M5)
check_heteroscedasticity(methaneYr2_M5)
shapiro.test(resid(methaneYr2_M5))

#### Year 1: Total year ####
# To understand cumulative emissions over the entire year, we will use the desctools package in R with the AUC function. Here we are specifying the AUC functions to calculate area under the curve using the ‘trapezoid’ method, which computes each trapezoid of adjacent points using the following calculation, allowing for time between events to vary:

# Get full year 1 data: 
GHGYr1_Full <- GHG %>%
  filter(StudyYear == "Year1")

# Transform units to be in kg/ha
GHGYr1_Full <- GHGYr1_Full %>%
  mutate(CH4Flux_kg_ha_h = CH4Flux * 0.00057744,       # nmol/m2/s -> kg/ha/h
         CH4Flux_kg_ha_day = CH4Flux_kg_ha_h * 24)     # kg/ha/h -> kg/ha/day

# Calculate area under the curve over time (elapsed days of experiment): 
results1 <- GHGYr1_Full %>% 
  group_by(FishTreatment, BirdTreatment, treatment, PlotID) %>% 
  summarise(
    auc_total = AUC(x = ElapsedDays, y = CH4Flux_kg_ha_day, method = "trapezoid"),
    time_range = max(ElapsedDays) - min(ElapsedDays),
    mean_flux = auc_total / time_range, 
    .groups = "drop"
  )

# Simple linear models (low sample size, cumulative flux per plot)

# Most complex model: treatment interaction
AUCYr1_M1 <- lm(mean_flux ~ FishTreatment * BirdTreatment,
         data = results1)
summary(AUCYr1_M1)

# diagnostics for most complex model: 
par(mfrow = c(2, 2))
plot(AUCYr1_M1)
# assumption: normality 
shapiro.test(resid(AUCYr1_M1))
# assumption: homoscedasticity
check_heteroscedasticity(AUCYr1_M1)

# Drop interaction 
AUCYr1_M2 <- lm(mean_flux ~ FishTreatment + BirdTreatment,
         data = results1)
summary(AUCYr1_M2)

anova(AUCYr1_M1, AUCYr1_M2) # chi-squared test not significant 

# Drop fish treatment 
AUCYr1_M3 <- lm(mean_flux ~ BirdTreatment,
         data = results1)
summary(AUCYr1_M3)

anova(AUCYr1_M2, AUCYr1_M3) # chi-squared test marginally significant 

# Drop bird treatment 
AUCYr1_M4 <- lm(mean_flux ~ FishTreatment, 
         data = results1)
summary(AUCYr1_M4)

# Null model
AUCYr1_M5 <- lm(mean_flux ~ 1,
         data = results1) 
summary(AUCYr1_M5)

anova(AUCYr1_M3, AUCYr1_M5) # chi-squared test not significant
anova(AUCYr1_M4, AUCYr1_M5) # chi-squared test marginally significant 

# check diagnostics for top model: 
par(mfrow = c(2, 2))
plot(AUCYr1_M4)
# assumption: normality 
shapiro.test(resid(AUCYr1_M4))
# assumption: homoscedasticity
check_heteroscedasticity(AUCYr1_M4)

#### Year 2: Total year ####
# Get full year 1 data: 
GHGYr2_Full <- GHG %>%
  filter(StudyYear == "Year2")

# Transform units to be in kg/ha
GHGYr2_Full <- GHGYr2_Full %>%
  mutate(CH4Flux_kg_ha_h = CH4Flux * 0.00057744,       # nmol/m2/s -> kg/ha/h
         CH4Flux_kg_ha_day = CH4Flux_kg_ha_h * 24)     # kg/ha/h -> kg/ha/day

# Calculate area under the curve over time (elapsed days of experiment): 
results2 <- GHGYr2_Full %>% 
  group_by(FishTreatment, BirdTreatment, treatment, PlotID) %>% 
  summarise(
    auc_total = AUC(x = ElapsedDays, y = CH4Flux_kg_ha_day, method = "trapezoid"),
    time_range = max(ElapsedDays) - min(ElapsedDays),
    mean_flux = auc_total / time_range, 
    .groups = "drop"
  )

# Simple linear models (low sample size, cumulative flux per plot)

# Most complex model: treatment interaction
AUCYr2_M1 <- lm(mean_flux ~ FishTreatment * BirdTreatment,
                data = results2)
summary(AUCYr2_M1)

# check diagnostics for most complex model: 
par(mfrow = c(2, 2))
plot(AUCYr2_M1)
# assumption: normality 
shapiro.test(resid(AUCYr2_M1))
# assumption: homoscedasticity
check_heteroscedasticity(AUCYr2_M1)

# Drop interaction 
AUCYr2_M2 <- lm(mean_flux ~ FishTreatment + BirdTreatment,
                data = results2)
summary(AUCYr2_M2)

anova(AUCYr2_M1, AUCYr2_M2) # chi-squared test not significant 

# Drop fish treatment 
AUCYr2_M3 <- lm(mean_flux ~ BirdTreatment,
                data = results2)
summary(AUCYr2_M3)

anova(AUCYr2_M2, AUCYr2_M3) # chi-squared test not significant 

# Drop bird treatment 
AUCYr2_M4 <- lm(mean_flux ~ FishTreatment, 
                data = results2)
summary(AUCYr2_M4)

# Null model
AUCYr2_M5 <- lm(mean_flux ~ 1,
                data = results2) 
summary(AUCYr2_M5)

anova(AUCYr2_M3, AUCYr2_M5) # chi-squared test not significant
anova(AUCYr2_M4, AUCYr2_M5) # chi-squared test not significant

# check diagnostics for null model: 
par(mfrow = c(2, 2))
plot(AUCYr2_M5)
# assumption: normality 
shapiro.test(resid(AUCYr2_M5))
# assumption: homoscedasticity
check_heteroscedasticity(AUCYr2_M5)


# Abiotic models  ---------------------------------------------------------
## For all abiotic metrics, our hypotheses are that independent variables will not have treatment effects, but will exhibit differences by study year. 

#### Water temperature ####

## Full model to ensure no interactive treatment effects with study year: 
temp_M1 <- lm(mean_temp ~ StudyYear * FishTreatment * BirdTreatment, data = hobo_summary)
summary(temp_M1)

# Diagnostics for most complex model: 
par(mfrow = c(2,2))  # arrange in a 2x2 grid
plot(temp_M1)

## Model investigating interactive treatment effects & study year: 
temp_M2 <- lm(mean_temp ~ (BirdTreatment * StudyYear) + (FishTreatment * StudyYear) + (FishTreatment * BirdTreatment), data = hobo_summary)
summary(temp_M2)

anova(temp_M1, temp_M2) # not significant, no interactive treatment effects. 

temp_M3 <- lm(mean_temp ~ (BirdTreatment * StudyYear) + (FishTreatment * StudyYear), data = hobo_summary)
summary(temp_M3)

anova(temp_M2, temp_M3) # not significant

temp_M4 <- lm(mean_temp ~ (FishTreatment * StudyYear) + BirdTreatment, data = hobo_summary)
summary(temp_M4)

anova(temp_M3, temp_M4) # not significant 

temp_M5 <- lm(mean_temp ~ FishTreatment * StudyYear, data = hobo_summary)
summary(temp_M5)

anova(temp_M4, temp_M5) # not significant 

temp_M6 <- lm(mean_temp ~ FishTreatment + StudyYear, data = hobo_summary)
summary(temp_M6)

anova(temp_M5, temp_M6) # not significant 

temp_M7 <- lm(mean_temp ~ StudyYear, data = hobo_summary)
summary(temp_M7)

anova(temp_M6, temp_M7) # not significant 

temp_M8 <- lm(mean_temp ~ 1, data = hobo_summary)
summary(temp_M8)

anova(temp_M7, temp_M8) # very significant 

# Diagnostics for most complex model: 
par(mfrow = c(2,2))  # arrange in a 2x2 grid
plot(temp_M7)

emm_temp <- emmeans(temp_M7, ~ StudyYear, type = "response")
pairs(emm_temp)

#### Dissolved oxygen ####
## Full model to ensure no interactive treatment effects with study year: 
DO_M1 <- lm(mean_DO ~ StudyYear * FishTreatment * BirdTreatment, data = hobo_summary)
summary(DO_M1)

# Diagnostics for most complex model: 
par(mfrow = c(2,2))  # arrange in a 2x2 grid
plot(DO_M1)

## Model investigating interactive treatment effects & study year: 
DO_M2 <- lm(mean_DO ~ (BirdTreatment * StudyYear) + (FishTreatment * StudyYear) + (FishTreatment * BirdTreatment), data = hobo_summary)
summary(DO_M2)

anova(DO_M1, DO_M2) # not significant, no interactive treatment effects. 

DO_M3 <- lm(mean_DO ~ (BirdTreatment * StudyYear) + (FishTreatment * StudyYear), data = hobo_summary)
summary(DO_M3)

anova(DO_M2, DO_M3) # not significant

DO_M4 <- lm(mean_DO ~ (FishTreatment * StudyYear) + BirdTreatment, data = hobo_summary)
summary(DO_M4)

anova(DO_M3, DO_M4) # significant 

DO_M5 <- lm(mean_DO ~ FishTreatment * StudyYear, data = hobo_summary)
summary(DO_M5)

anova(DO_M4, DO_M5) # marginally significant 

DO_M6 <- lm(mean_DO ~ FishTreatment + StudyYear, data = hobo_summary)
summary(DO_M6)

anova(DO_M5, DO_M6) # not significant 

DO_M7 <- lm(mean_DO ~ StudyYear, data = hobo_summary)
summary(DO_M7)

anova(DO_M6, DO_M7) # not significant 

DO_M8 <- lm(mean_DO ~ 1, data = hobo_summary)
summary(DO_M8)

anova(DO_M7, DO_M8) # very significant 

# Diagnostics for most top model: 
par(mfrow = c(2,2))  # arrange in a 2x2 grid
plot(DO_M7)

emm_DO <- emmeans(DO_M3, ~StudyYear, type = "response")
pairs(emm_DO)

#### Depth ####
## Full model to ensure no interactive treatment effects with study year: 
depth_M1 <- lm(mean_depth ~ StudyYear * FishTreatment * BirdTreatment, data = depth)
summary(depth_M1)

# Diagnostics for most complex model: 
par(mfrow = c(2,2))  # arrange in a 2x2 grid
plot(depth_M1)

## Model investigating interactive treatment effects & study year: 
depth_M2 <- lm(mean_depth ~ (BirdTreatment * StudyYear) + (FishTreatment * StudyYear) + (FishTreatment * BirdTreatment), data = depth)
summary(depth_M2)

anova(depth_M1, depth_M2) # not significant, no interactive treatment effects. 

depth_M3 <- lm(mean_depth ~ (BirdTreatment * StudyYear) + (FishTreatment * StudyYear), data = depth)
summary(depth_M3)

anova(depth_M2, depth_M3) # not significant

depth_M4 <- lm(mean_depth ~ (FishTreatment * StudyYear) + BirdTreatment, data = depth)
summary(depth_M4)

anova(depth_M3, depth_M4) # marginally significant 

depth_M5 <- lm(mean_depth ~ FishTreatment * StudyYear, data = depth)
summary(depth_M5)

anova(depth_M4, depth_M5) # not significant 

depth_M6 <- lm(mean_depth ~ FishTreatment + StudyYear, data = depth)
summary(depth_M6)

anova(depth_M5, depth_M6) # not significant 

depth_M7 <- lm(mean_depth ~ StudyYear, data = depth)
summary(depth_M7)

anova(depth_M6, depth_M7) # not significant 

depth_M8 <- lm(mean_depth ~ 1, data = depth)
summary(depth_M8)

anova(depth_M7, depth_M8) # very significant 

# Diagnostics for most top model: 
par(mfrow = c(2,2))  # arrange in a 2x2 grid
plot(depth_M7)

emm_depth <- emmeans(depth_M7, ~ StudyYear, type = "response")
pairs(emm_depth)
