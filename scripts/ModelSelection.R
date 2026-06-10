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
totaldaphnia1 <- read.csv("data/totaldaphnia1.csv") # Total daphnia dataset year 1
totaldaphnia2 <- read.csv("data/totaldaphnia2.csv") # Total daphnia dataset year 2
isotope_yr1 <- read.csv("data/isotope_yr1.csv") # Full porewater dataset year 1
isotope_yr2 <- read.csv("data/isotope_yr2.csv") # Full porewater dataset year 2
GHG_Yr1 <- read.csv("data/GHG_Yr1_winter.csv") # Winter flooding methane dataset year 1
GHG_Yr2 <- read.csv("data/GHG_Yr2_winter.csv") # Winter flooding methane dataset year 2
GHG <- read.csv("data/GHG_ALL.csv") # full GHG dataset 


# Model set 1: Daphnia Abundance ------------------------------------------

## Iterative framework: beginning with most complex/biologically relevent model and iteratively simplifying 

#### Year 1: Full dataset ####
# Most complex model: 
DaphniaFull_M1<- glmmTMB(
  DensityRounded ~ FishTreatment * BirdTreatment * Pre_Post +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(Daphnia_M1)

# Drop 3-way interaction: 
DaphniaFull_M2 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + (BirdTreatment * Pre_Post) + (FishTreatment*BirdTreatment) +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaFull_M2)

anova(DaphniaFull_M1, DaphniaFull_M2) # chi-squared test not significant

# Drop bird/fish interaction:
DaphniaFull_M3 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + (BirdTreatment * Pre_Post) +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaFull_M3)

anova(DaphniaFull_M2, DaphniaFull_M3) # chi-squared test not significant 

# Drop bird/PrePost interaction: 
DaphniaFull_M4 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + BirdTreatment +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaFull_M4)

anova(DaphniaFull_M3, DaphniaFull_M4) # chi-squared test not significant 

# Drop bird treatment: 
DaphniaFull_M5 <- glmmTMB(
  DensityRounded ~ FishTreatment * Pre_Post + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaFull_M5) 

anova(DaphniaFull_M4, DaphniaFull_M5) # chi-squared test not significant 

# Drop fish/PrePost interaction: 
DaphniaFull_M6 <- glmmTMB(
  DensityRounded ~ FishTreatment + Pre_Post + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaFull_M6) 

anova(DaphniaFull_M5, DaphniaFull_M6) # chi-squared test marginally significant

# Drop PrePost:
DaphniaFull_M7 <- glmmTMB(
  DensityRounded ~ FishTreatment + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaFull_M7)

anova(DaphniaFull_M6, DaphniaFull_M7) # chi-squared test not significant   

# Null model: 
DaphniaFull_M8 <- glmmTMB(
  DensityRounded ~ 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)
summary(DaphniaFull_M8)

anova(DaphniaFull_M7, DaphniaFull_M8) # chi-squared test not significant 

## Model diagnostics for best model: 
simulationOutput1 <- simulateResiduals(fittedModel = DaphniaFull_M5)
plot(simulationOutput1)
check_model(DaphniaFull_M5)

#### Year 1: Omitting final sampling date ####
# Remove last sampling date 
totaldaphnia1_filtered <- totaldaphnia1 %>%
  filter(SampleDate != as.Date("2024-02-05"))

# Full model
Daphnia_M1 <- glmmTMB(
  DensityRounded ~ FishTreatment * BirdTreatment * Pre_Post +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(Daphnia_M1)

# Drop three-way interaction 
Daphnia_M2 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + (BirdTreatment * Pre_Post) + (FishTreatment * BirdTreatment) +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(Daphnia_M2)

anova(Daphnia_M1, Daphnia_M2) # chi-squared test marginally significant

# Drop fish*bird interaction 
Daphnia_M3 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + (BirdTreatment * Pre_Post) +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(Daphnia_M3)

anova(Daphnia_M2, Daphnia_M3) # chi-squared test not significant 

# Drop bird/pre-post interaction 
Daphnia_M4 <- glmmTMB(
  DensityRounded ~ (FishTreatment * Pre_Post) + BirdTreatment +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(Daphnia_M4)

anova(Daphnia_M3, Daphnia_M4) # chi-squared test not significant 

# Drop bird treatment 
Daphnia_M5 <- glmmTMB(
  DensityRounded ~ FishTreatment * Pre_Post + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID )  + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(Daphnia_M5)

anova(Daphnia_M4, Daphnia_M5) # chi-squared test not significant 

# Drop fish/pre-post interaction 
Daphnia_M6 <- glmmTMB(
  DensityRounded ~ FishTreatment + Pre_Post + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(Daphnia_M6)

anova(Daphnia_M5, Daphnia_M6) # chi-squared test significant, p = 0.02

# Drop pre/post
Daphnia_M7 <- glmmTMB(
  DensityRounded ~ FishTreatment + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(Daphnia_M7)

anova(Daphnia_M6, Daphnia_M7) # chi-squared test not significant 

# Null model
Daphnia_M8 <- glmmTMB(
  DensityRounded ~ 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)
summary(Daphnia_M8)

anova(Daphnia_M7, Daphnia_M8)

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

## Model diagnostics for best model: 
simulationOutput1 <- simulateResiduals(fittedModel = DaphniaYr2_M1)
plot(simulationOutput1)
check_overdispersion(DaphniaYr2_M1)
check_model(DaphniaYr2_M1)


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
porewater_M1 <- lm(C13_log ~ (FishTreatment * BirdTreatment), 
                   data = isotope_yr1_fish)

summary(porewater_M1)

# Drop fish/bird interaction: 
porewater_M2 <- lm(C13_log ~ (FishTreatment + BirdTreatment), 
                   data = isotope_yr1_fish)
summary(porewater_M2)

anova(porewater_M1, porewater_M2) # chi-squared test not significant

# Drop fish treatment 
porewater_M3 <- lm(C13_log ~ BirdTreatment, 
                   data = isotope_yr1_fish)
summary(porewater_M3)

anova(porewater_M2, porewater_M3) # chi-squared test marginally significant (p = 0.09)

# Drop bird treatment 
porewater_M4 <- lm(C13_log ~ FishTreatment, 
                   data = isotope_yr1_fish)
summary(porewater_M4)

sim_anova1 <- simulateResiduals(fittedModel = porewater_M4) # diagnostics for best model
plot(sim_anova1)
shapiro.test(resid(porewater_M4)) # normally distributed

# Null model, compare fish and bird only models to this: 
porewater_M5 <- lm(C13_log ~ 1, 
                   data = isotope_yr1_fish)
anova(porewater_M3, porewater_M5)
anova(porewater_M4, porewater_M5)

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
porewater_yr2_M1 <- lm(C13_log ~ (FishTreatment * BirdTreatment), 
                       data = isotope_yr2_fish)

summary(porewater_yr2_M1)

# Drop fish/bird interaction: 
porewater_yr2_M2 <- lm(C13_log ~ (FishTreatment + BirdTreatment), 
                       data = isotope_yr2_fish)
summary(porewater_yr2_M2)

anova(porewater_y2_M1, porewater_yr2_M2) # chi-squared test not significant

# Drop fish treatment 
porewater_yr2_M3 <- lm(C13_log ~ BirdTreatment, 
                       data = isotope_yr2_fish)
summary(porewater_yr2_M3)

anova(porewater_yr2_M2, porewater_yr2_M3) # chi-squared test not significant

# Drop bird treatment 
porewater_yr2_M4 <- lm(C13_log ~ FishTreatment, 
                       data = isotope_yr2_fish)
summary(porewater_yr2_M4)

# Null model, compare fish and bird only models to this: 
porewater_yr2_M5 <- lm(C13_log ~ 1, 
                       data = isotope_yr2_fish)
anova(porewater_yr2_M3, porewater_yr2_M5) # chi-squared test not significant 
anova(porewater_yr2_M4, porewater_yr2_M5) # chi-squared test not significant 


# Model set 3: Methane Flux -----------------------------------------------

#### Year 1: Winter flux ####
# Only two data points per plot pre- fish introduction. With low sample size, we will first establish no differences between treatments pre-fish introduction, and then do model selection on dataset post-fish introduction. 

# Establish no difference pre-fish intro:

#### Year 2: Winter flux ####



#### Year 1: Total year ####



#### Year 2: Total year ####