#### Fish/bird trophic effects on methane: manuscript 
#### Target journal: PNAS 
#### Figures, main manuscript body

# Load packages -----------------------------------------------------------
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(ggeffects)
library(DescTools)
library(ggbeeswarm)
library(patchwork)
library(cowplot)
library(ggnewscale)

# Load data ---------------------------------------------------------------
totaldaphnia1 <- read.csv("data/totaldaphnia1.csv") # Total daphnia dataset year 1
totaldaphnia2 <- read.csv("data/totaldaphnia2.csv") # Total daphnia dataset year 2
isotope_yr1 <- read.csv("data/isotope_yr1.csv") # Full porewater dataset year 1
isotope_yr2 <- read.csv("data/isotope_yr2.csv") # Full porewater dataset year 2
GHG_Yr1 <- read.csv("data/GHG_Yr1_winter.csv") # Winter flooding methane dataset year 1
GHG_Yr2 <- read.csv("data/GHG_Yr2_winter.csv") # Winter flooding methane dataset year 2
GHG <- read.csv("data/GHG_ALL.csv") # full GHG dataset 


# Figure 1 ----------------------------------------------------------------


# Figure 2 ----------------------------------------------------------------

## Showcasing top model output for each study year: 

#### Year 1: ####
# In main manuscript, showcasing data set with last sampling date removed:

# Remove last sampling date 
totaldaphnia1_filtered <- totaldaphnia1 %>%
  filter(SampleDate != as.Date("2024-02-05"))

# Year 1, top model:
P2a <- glmmTMB(
  DensityRounded ~ FishTreatment * Pre_Post + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1_filtered
)

emm2a <- emmeans(P2a, ~ FishTreatment * Pre_Post,
                type = "response")

pairs(emm2a, by = "FishTreatment")

emm_df2a <- as.data.frame(emm2a)

dat2a <- totaldaphnia1_filtered %>% 
  mutate(Pre_Post = factor(Pre_Post, levels = c("Pre", "Post")),
         FishTreatment = recode(FishTreatment, "NoFish" = "No Fish")) 

emm_df2a <- emm_df2a %>% 
  mutate(
    Pre_Post = factor(Pre_Post, levels = c("Pre", "Post")),
    FishTreatment = recode(FishTreatment, "NoFish" = "No Fish"))


treatment_colors_2a <- c(
  "Birds Fish"    = "#1B80ADFF",
  "NoBirds Fish"  = "#00B398FF",  
  "Birds NoFish"  = "#FF9933FF",  
  "NoBirds NoFish"= "#C24841FF"
)

label_df2a <- data.frame(
  FishTreatment = c("Fish", "No Fish"),
  x = c(1.5, 1.5),
  y = c(580, 650),
  label = c("*", "n.s.")
)

plot_2a <- ggplot() +
  geom_quasirandom( 
    data = dat2a, aes(x = Pre_Post, y = total, 
                     color = treatment),
    width = 0.2, alpha = 0.5, size = 1.5
  ) +
  scale_color_manual(values = treatment_colors_2a, name = "Treatment") + 
  new_scale_color() +
  geom_line(
    data = emm_df2a, aes(x = Pre_Post, y = response,
                        group = FishTreatment) ,
                        color = "black",
                        linewidth = 0.5,
                        alpha = 0.5
    ) +
  geom_point(
    data = emm_df2a, aes(x = Pre_Post, y = response),
                        color = "black",
                        size = 3.5,
                        alpha = 0.6
    ) +
  geom_errorbar(
    data = emm_df2a, aes(x = Pre_Post, ymin = asymp.LCL, ymax = asymp.UCL),
    color = "black",
    width = 0.15,
    alpha = 0.5
    ) +
  facet_wrap(~ FishTreatment, strip.position = "top") +
  scale_y_log10() +
  coord_cartesian(ylim = c(20, 15500)) +
  labs(x = "", 
       y = expression("Daphnia Abundance (individuals m"^-3*")")) +
  ggtitle("Year 1") +
  theme_minimal_grid(font_size = 13) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none") +
  geom_text(data = label_df2a, 
            aes(x = x, y = y, label = label), 
            inherit.aes = FALSE, 
            size = c(7, 4))

plot_2a


#### Year 2: ####
P2b <- glmmTMB(
  DensityRounded ~ FishTreatment * BirdTreatment * Pre_Post +
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia2
)

emm2b <- emmeans(P2b, ~ Pre_Post | FishTreatment * BirdTreatment,
                type = "response")

contrast(emm2b, 
         method = "pairwise", 
         by = c("FishTreatment", "BirdTreatment"))

emm_df2b <- as.data.frame(emm2b)

dat2b <- totaldaphnia2 %>% 
  mutate(
    Pre_Post = factor(Pre_Post, levels = c("Pre", "Post")),
    FishTreatment = recode(FishTreatment, "NoFish" = "No Fish")
  )

emm_df2b <- emm_df2b %>% 
  mutate(
    Pre_Post = factor(Pre_Post, levels = c("Pre", "Post")),
    FishTreatment = recode(FishTreatment, "NoFish" = "No Fish")
  )

treatment_order_2b <- c("Fish_Birds", "Fish_NoBirds", "No Fish_Birds", "No Fish_NoBirds")

treatment_colors_2b <- c(
  "Fish_Birds"       = "#1B80ADFF",
  "Fish_NoBirds"     = "#00B398FF",
  "No Fish_Birds"    = "#FF9933FF",
  "No Fish_NoBirds"  = "#C24841FF"
)

label_df2b <- data.frame(
  FishTreatment = c("Fish", "Fish", "No Fish", "No Fish"),
  x = c(1.4, 1.5, 1.4, 1.5),
  y = c(1800, 3800, 3800, 1800) , 
  label = c("n.s.", "*", "n.s.", "n.s.")
)

dat2b <- dat2b %>%
  mutate(color_group = factor(
    interaction(FishTreatment, BirdTreatment, sep = "_"),
    levels = treatment_order_2b))

emm_df2b <- emm_df2b %>%
  mutate(color_group = factor(
    interaction(FishTreatment, BirdTreatment, sep = "_"),
    levels = treatment_order_2b))

plot_2b <- ggplot() +
  geom_quasirandom( 
    data = dat2b, aes(x = Pre_Post, y = total, color = color_group),
    width = 0.2, alpha = 0.5, size = 1.5
  ) +
  geom_line(
    data = emm_df2b,
    aes(x = Pre_Post, y = response, color = color_group, group = color_group),
    linewidth = 0.5, position = position_dodge(width = 0.5)
  ) +
  geom_point(
    data = emm_df2b, aes(x = Pre_Post, y = response, color = color_group),
    size = 3.5, position = position_dodge(width = 0.5), alpha = 0.8
  ) +
  geom_errorbar(
    data = emm_df2b, aes(x = Pre_Post, ymin = asymp.LCL, ymax = asymp.UCL,
                        color = color_group),
    width = 0.15, position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~FishTreatment, strip.position = "top") +
  scale_color_manual(
    values = treatment_colors_2b,
    name = "Treatment",
    labels = c(
      "Fish_Birds"     = "Fish / Birds",
      "Fish_NoBirds"   = "Fish / No Birds",
      "No Fish_Birds"   = "No Fish / Birds",
      "No Fish_NoBirds" = "No Fish / No Birds"
    )) +
  scale_y_log10() +
  coord_cartesian(ylim = c(20, 15500)) +
  labs(x = "", y = "") +
  ggtitle("Year 2") +
  theme_minimal_grid(font_size = 13) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.y = element_blank()
  ) +
  guides(color = guide_legend(
    nrow = 1, override.aes = list(
      shape = 16,
      size = 3,
      linetype = 0,
      linewidth = 0))) +
  geom_text(data = label_df2b, 
            aes(x = x, y = y, label = label), 
            inherit.aes = FALSE, 
            size = c(4, 7, 4, 4))

plot_2b

#### Put figures together: 

legend_grob <- get_legend(plot_2b)  # extract legend from plot 2b

# Remove legends from ndividual plots
plot_2a_noleg <- plot_2a + theme(legend.position = "none")
plot_2b_noleg <- plot_2b + theme(legend.position = "none")

# combine using patchwork w/ legend as third row: 

Figure2 <- (plot_2a_noleg | plot_2b_noleg) /
  ((plot_spacer() | wrap_elements(full = legend_grob) | plot_spacer()) + 
     plot_layout(widths = c(0.25, 2, 1))) +
  plot_layout(heights = c(8, 0.1))

Figure2

# Export figure 

# Save figure: 
ggsave("figures/Figure_2.png", 
       plot = Figure2,
       width = 20, 
       height = 18, 
       dpi = 300,
       units = "cm")


# Figure 3 ----------------------------------------------------------------

#### Year 1: ####

# Set up for modeling/visualizing: 
# subset out for when fish were on fields: 
isotope_yr1_fish <- isotope_yr1 %>% 
  filter(FishOnField == "Y")

# Remove NA Plot IDs from dataset:
isotope_yr1_fish <- isotope_yr1_fish %>% 
  drop_na(C13PW)

# Log-transform response variable, with a slight offset to avoid producing NaNs
isotope_yr1_fish <- isotope_yr1_fish %>% 
  mutate(C13_log = log(C13PW - min(C13PW) + 1),
         FishTreatment = recode(FishTreatment, "NoFish" = "No Fish"))

# Year 1, top model: 
F3a <- lm(C13_log ~ FishTreatment, 
         data = isotope_yr1_fish)
summary(F3a)

emm3a <- emmeans(F3a, ~ FishTreatment, type = "response")
emm3a
emm_df3a <- as.data.frame(emm3a)

# backtransform response (log scale with slight offset):
min_val <- min(isotope_yr1_fish$C13PW)

emm_df3a$response <- exp(emm_df3a$emmean) + min_val - 1
emm_df3a$lower    <- exp(emm_df3a$lower.CL) + min_val - 1
emm_df3a$upper    <- exp(emm_df3a$upper.CL) + min_val - 1

treatment_colors_3a <- c(
  "Birds Fish"    = "#1B80ADFF",
  "NoBirds Fish"  = "#00B398FF",  
  "Birds NoFish"  = "#FF9933FF",  
  "NoBirds NoFish"= "#C24841FF"
)

# Year 1 plot:
plot_3a <- ggplot() +
  geom_jitter(data = isotope_yr1_fish,
              aes(x = FishTreatment, y = C13PW, color = treatment),
              width = 0.1,
              alpha = 0.5, 
              size = 2) +
  scale_color_manual(values = treatment_colors_3a, name = "Treatment") +
  new_scale_color() +
  geom_point(
    data = emm_df3a, aes(x = FishTreatment, y = response),
    color = "black",
    size = 3.5,
    alpha = 0.6
  ) +
  geom_errorbar(
    data = emm_df3a, aes(x = FishTreatment, ymin = lower, ymax = upper),
    color = "black",
    width = 0.15,
    alpha = 0.5
  ) +
  labs(title = "Year 1",
       y = expression(delta^{13}*C~"(‰)"),
       x = "") +
  theme_minimal_grid(font_size = 13) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none") +
  ylim(-62.5, -40) +
  geom_text(aes(x = 1.5, y = -44, 
                label = "+"), 
            size = 6) +
  scale_x_discrete(position = "top")

plot_3a

#### Year 2: ####

# Year 2, top model is null model, so will just show raw data:
# Set up for visualizing: 
# subset out for when fish were on fields: 
isotope_yr2_fish <- isotope_yr2 %>% 
  filter(FishOnField == "Y")

# Remove NA Plot IDs from dataset:
isotope_yr2_fish <- isotope_yr2_fish %>% 
  drop_na(C13PW) %>% 
  mutate(C13_log = log(C13PW - min(C13PW) + 1),
         FishTreatment = recode(FishTreatment, "NoFish" = "No Fish"))

treatment_order_3b <- c("Birds Fish", "NoBirds Fish", "Birds NoFish", "NoBirds NoFish")

treatment_colors_3b <- c(
  "Birds Fish"       = "#1B80ADFF",
  "NoBirds Fish"     = "#00B398FF",
  "Birds NoFish"    = "#FF9933FF",
  "NoBirds NoFish"  = "#C24841FF"
)

dat3b <- isotope_yr2_fish %>%
  mutate(color_group = factor(
    treatment,
    levels = treatment_order_3b))

plot_3b <- ggplot() +
  geom_jitter(data = dat3b,
              aes(x = FishTreatment, y = C13PW, color = color_group),
              width = 0.1,
              alpha = 0.5,
              size = 2) +
  labs(title = "Year 2",
       y = "",
       x = "") +
  theme_minimal_grid(font_size = 13) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  ylim(-62.5, -40) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 8)) +
  guides(color = guide_legend(
    nrow = 1, override.aes = list(
      shape = 16,
      size = 2,
      linetype = 0,
      linewidth = 0,
      alpha = 1))) +
  scale_color_manual(values = treatment_colors_3b, 
                     name = "Treatment",
                     labels = c(
                       "Birds Fish"     = "Fish / Birds",
                       "NoBirds Fish"   = "Fish / No Birds",
                       "Birds NoFish"   = "No Fish / Birds",
                       "NoBirds NoFish" = "No Fish / No Birds"
                     )) +
  geom_text(aes(x = 1.5, y = -44, 
                label = "n.s."), 
            size = 4) +
  scale_x_discrete(position = "top")

plot_3b

#### Put plots together: 

legend_grob <- get_legend(plot_3b)  # extract legend from plot 3b

# Remove legends from individual plots
plot_3a_noleg <- plot_3a + theme(legend.position = "none")
plot_3b_noleg <- plot_3b + theme(legend.position = "none")

# combine using patchwork w/ legend as third row: 

Figure3 <- (plot_3a_noleg | plot_3b_noleg) /
  ((plot_spacer() | wrap_elements(full = legend_grob) | plot_spacer()) + 
     plot_layout(widths = c(0.02, 2, 1))) +
  plot_layout(heights = c(8, 0.5))

Figure3

# Export figure 

# Save figure: 
ggsave("figures/Figure_3.png", 
       plot = Figure3,
       width = 15, 
       height = 11, 
       dpi = 300,
       units = "cm")


# Figure 4 ----------------------------------------------------------------

#### Instantaneous flux ####

#### Year 1: ####

# set up for model: 
GHG_Yr1_post <- GHG_Yr1 %>% 
  filter(Pre_Post == "Post") %>% 
  mutate(FishTreatment = recode(FishTreatment, "NoFish" = "No Fish"))

# top model: 

F4a_1 <- lmer(logFlux ~ FishTreatment + 
             scale(MeanTemp) + scale(PO4) + scale(NH4) + scale(MeanDepth) +
             (1 | PlotID) + (1 | SamplingOccasion),
           data = GHG_Yr1_post)

emm4a_1 <- emmeans(F4a_1, ~ FishTreatment)
emm4a_1
pairs(emm4a_1)
emm_df4a_1 <- as.data.frame(emm4a_1)

treatment_colors_4a1 <- c(
  "Birds Fish"    = "#1B80ADFF",
  "NoBirds Fish"  = "#00B398FF",  
  "Birds NoFish"  = "#FF9933FF",  
  "NoBirds NoFish"= "#C24841FF"
)

# Plot: 

plot_4a1 <- ggplot() +
  geom_jitter(data = GHG_Yr1_post,
              aes(x = FishTreatment, y = logFlux, color = treatment),
              width = 0.1,
              alpha = 0.5,
              size = 2) +
  scale_color_manual(values = treatment_colors_4a1, name = "Treatment") +
  new_scale_color() +
  geom_point(data = emm_df4a_1, aes(x = FishTreatment, y = emmean),
             color = "black",
             size = 3.5,
             alpha = 0.6) +
  geom_errorbar(data = emm_df4a_1,
                aes(x = FishTreatment,
                    ymin = lower.CL,
                    ymax = upper.CL),
                color = "black",
                width = 0.15,
                alpha = 0.5) +
  labs(title = "Winter Flooded Season Flux",
       y = expression(Log~C~-CH[4]~flux~(nmol~m^{-2}~s^{-1})),
       x = "") +
  theme_minimal_grid(font_size = 12) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3)) +
  panel_border() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    legend.position = "none") +
  ylim(-0.5, 9) +
  geom_text(aes(x = 1.5, y = 9, 
                label = "**"), 
            size = 6) +
  scale_x_discrete(position = "top")

plot_4a1

#### Year 2: ####

# set up for figure:
## (Year 2 top model is null model - showing raw data:)
GHG_Yr2_post <- GHG_Yr2 %>% 
  filter(Pre_Post == "Post") %>% 
  mutate(FishTreatment = recode(FishTreatment, "NoFish" = "No Fish"))

treatment_order_4a2 <- c("Birds Fish", "NoBirds Fish", "Birds NoFish", "NoBirds NoFish")

treatment_colors_4a2 <- c(
  "Birds Fish"       = "#1B80ADFF",
  "NoBirds Fish"     = "#00B398FF",
  "Birds NoFish"    = "#FF9933FF",
  "NoBirds NoFish"  = "#C24841FF"
)

dat4a2 <- GHG_Yr2_post %>%
  mutate(color_group = factor(
    treatment,
    levels = treatment_order_4a2))

# Plot: 

plot_4a2 <- ggplot() +
  geom_jitter(data = dat4a2,
              aes(x = FishTreatment, y = logFlux, color = color_group),
              width = 0.1,
              alpha = 0.5,
              size = 2) +
  labs(title = "Winter Flooded Season Flux",
       y = expression(Log~C~-CH[4]~flux~(nmol~m^{-2}~s^{-1})),
       x = "") +
  theme_minimal_grid(font_size = 12) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    legend.position = "none") +
  ylim(-0.5,9) +
  scale_color_manual(values = treatment_colors_4a2, 
                     name = "Treatment",
                     labels = c(
                       "Birds Fish"     = "Fish / Birds",
                       "NoBirds Fish"   = "Fish / No Birds",
                       "Birds NoFish"   = "No Fish / Birds",
                       "NoBirds NoFish" = "No Fish / No Birds"
                     )) +
  geom_text(aes(x = 1.5, y = 9, 
                label = "n.s."), 
            size = 4) +
  scale_x_discrete(position = "top")

plot_4a2

#### Total Annual flux ####

#### Year 1: ####

# Set up for modeling: 

# Filter out for study year 1:
GHGYr1_Full <- GHG %>%
  filter(StudyYear == "Year1")

# Transform units to kg/ha/day
GHGYr1_Full <- GHGYr1_Full %>%
  mutate(CH4Flux_kg_ha_h = CH4Flux * 0.00057744, # nmol/m2/s -> kg/ha/h
         CH4Flux_kg_ha_day = CH4Flux_kg_ha_h * 24) # kg/ha/h -> kg/ha/day

results1 <- GHGYr1_Full %>% 
  group_by(FishTreatment, BirdTreatment, treatment, PlotID) %>% 
  summarise(
    auc_total = AUC(x = ElapsedDays, y = CH4Flux_kg_ha_day, method = "trapezoid"),
    time_range = max(ElapsedDays) - min(ElapsedDays),
    mean_flux = auc_total / time_range, 
    .groups = "drop"
  ) %>% 
  mutate(FishTreatment = recode(FishTreatment, "NoFish" = "No Fish"))

# Top model: 

F4b_1 <- lm(mean_flux ~ FishTreatment, 
            data = results1)

emm4b1 <- emmeans(F4b_1, ~ FishTreatment, type = "response")
emm4b1
pairs(emm4b1)
emm_df_4b1 <- as.data.frame(emm4b1)


treatment_colors_4b1 <- c(
  "Birds Fish"    = "#1B80ADFF",
  "NoBirds Fish"  = "#00B398FF",  
  "Birds NoFish"  = "#FF9933FF",  
  "NoBirds NoFish"= "#C24841FF"
)

# Plot: 
plot_4b1 <- ggplot() +
  geom_jitter(data = results1,
              aes(x = FishTreatment, y = mean_flux, color = treatment),
              width = 0.1,
              alpha = 0.5,
              size = 2) +
  scale_color_manual(values = treatment_colors_4b1, name = "Treatment") +
  new_scale_color() +
  geom_point(data = emm_df_4b1,
             aes(x = FishTreatment, y = emmean),
             color = "black",
             size = 3,
             alpha = 0.6) +
  geom_errorbar(data = emm_df_4b1,
                aes(x = FishTreatment,
                    ymin = lower.CL,
                    ymax = upper.CL),
                color = "black",
                width = 0.15,
                alpha = 0.5) +
  labs( 
    title = "Total Annual Emissions",
    x = "",
    y = expression(Cumulative~C~-CH[4]~emissions~(kg~ha^{-1}))
  ) +
  theme_minimal_grid(font_size = 12) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    legend.position = "none") +
  ylim(0, 16) +
  geom_text(aes(x = 1.5, y = 16, 
                label = "+"), 
            size = 6) +
  scale_x_discrete(position = "top")

plot_4b1

#### Year 2: ####
# Set up for figures: 
## (Top model for year 2 is null model, will show raw data:)

# Filter out for study year 1:
GHGYr2_Full <- GHG %>%
  filter(StudyYear == "Year2")

# Transform units to kg/ha/day
GHGYr2_Full <- GHGYr2_Full %>%
  mutate(CH4Flux_kg_ha_h = CH4Flux * 0.00057744, # nmol/m2/s -> kg/ha/h
         CH4Flux_kg_ha_day = CH4Flux_kg_ha_h * 24) # kg/ha/h -> kg/ha/day

results2 <- GHGYr2_Full %>% 
  group_by(FishTreatment, BirdTreatment, treatment, PlotID) %>% 
  summarise(
    auc_total = AUC(x = ElapsedDays, y = CH4Flux_kg_ha_day, method = "trapezoid"),
    time_range = max(ElapsedDays) - min(ElapsedDays),
    mean_flux = auc_total / time_range, 
    .groups = "drop"
  ) %>% 
  mutate(FishTreatment = recode(FishTreatment, "NoFish" = "No Fish"))

treatment_order_4b2 <- c("Birds Fish", "NoBirds Fish", "Birds NoFish", "NoBirds NoFish")

treatment_colors_4b2 <- c(
  "Birds Fish"       = "#1B80ADFF",
  "NoBirds Fish"     = "#00B398FF",
  "Birds NoFish"    = "#FF9933FF",
  "NoBirds NoFish"  = "#C24841FF"
)

dat4b2 <- results2 %>%
  mutate(color_group = factor(
    treatment,
    levels = treatment_order_4b2))

# Plot: 

plot_4b2 <- ggplot() +
  geom_jitter(data = dat4b2,
              aes(x = FishTreatment, y = mean_flux, color = color_group),
              width = 0.1,
              alpha = 0.5,
              size = 2) +
  labs( 
    title = "Total Annual Emissions",
    x = "",
    y = expression(Cumulative~C~-CH[4]~emissions~(kg~ha^{-1}))
  ) +
  theme_minimal_grid(font_size = 13) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    legend.position = "none") +
  scale_color_manual(values = treatment_colors_4a2, 
                     name = "Treatment",
                     labels = c(
                       "Birds Fish"     = "Fish / Birds",
                       "NoBirds Fish"   = "Fish / No Birds",
                       "Birds NoFish"   = "No Fish / Birds",
                       "NoBirds NoFish" = "No Fish / No Birds"
                     )) +
  ylim(0,16) +
  geom_text(aes(x = 1.5, y = 16, 
                label = "n.s."), 
            size = 4) +
  scale_x_discrete(position = "top")

plot_4b2

#### Cumulative flux ####

#### Year 1: ####
# Set up datasets: 

GHGYr1_Full <- GHG %>%
  filter(StudyYear == "Year1")

GHGYr1_Full <- GHGYr1_Full %>%
  mutate(CH4Flux_kg_ha_h = CH4Flux * 0.00057744,       # nmol/m2/s -> kg/ha/h
         CH4Flux_kg_ha_day = CH4Flux_kg_ha_h * 24)     # kg/ha/h -> kg/ha/day

# Running cumulative AUC per plot
cumulative_data1 <- GHGYr1_Full %>%
  group_by(FishTreatment, BirdTreatment, treatment, PlotID) %>%
  arrange(ElapsedDays, .by_group = TRUE) %>%
  mutate(
    # Trapezoid area between each consecutive pair of points
    delta_t     = ElapsedDays - lag(ElapsedDays),
    trap_area   = 0.5 * (CH4Flux_kg_ha_day + lag(CH4Flux_kg_ha_day)) * delta_t,
    trap_area   = ifelse(is.na(trap_area), 0, trap_area),  # first point = 0
    cumulative_CH4 = cumsum(trap_area),
    FishTreatment = recode(FishTreatment, "NoFish" = "No Fish")
  ) %>%
  ungroup()

# Create dataframe to add vline for when fish were on fields: 

cumulative_data1$Date = as.Date(cumulative_data1$Date)

ghg_fish_dates1 <- data.frame(
  fish_days = as.Date(c("2023-12-13", "2024-02-27")))

treatment_colors_4c1 <- c(
  "Birds Fish"    = "#1B80ADFF",
  "NoBirds Fish"  = "#00B398FF",  
  "Birds NoFish"  = "#FF9933FF",  
  "NoBirds NoFish"= "#C24841FF"
)


# Plot: 

  
plot_4c1 <- ggplot() +
  geom_line(data = cumulative_data1,
            aes(x = Date, y = cumulative_CH4, 
                color = treatment, group = PlotID),
            alpha = 0.5, 
            linewidth = 0.5) +  
  scale_color_manual(values = treatment_colors_4c1) +
  new_scale_color() +# individual plots
  stat_summary(data = cumulative_data1,
                aes(x = Date, y = cumulative_CH4,
                    group = FishTreatment,
                    color = FishTreatment),
               fun = mean, 
               geom = "line", 
               linewidth = 1.5, 
               alpha = 1) + # treatment means
  scale_color_manual(values = c("black", "gray50")) +
  labs(title = "Net Cumulative Annual Emissions",
       x     = "Date (Month-Yr)",
       y     = expression(C~-CH[4]~emissions~(kg~ha^{-1}~day^{-1}))) +
  geom_vline(
    data = ghg_fish_dates1,
    aes(xintercept = fish_days),
    linetype = "dashed",
    color = "darkgrey"
  ) +
  scale_x_date(date_labels = "%b-%y", date_breaks = "2 month") +
  theme_minimal_grid(font_size = 12) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    legend.position = "none") +
  ylim(0, 5200)


plot_4c1  


#### Year 2: ####
# Set up datasets: 

GHGYr2_Full <- GHG %>%
  filter(StudyYear == "Year2")

GHGYr2_Full <- GHGYr2_Full %>%
  mutate(CH4Flux_kg_ha_h = CH4Flux * 0.00057744,       # nmol/m2/s -> kg/ha/h
         CH4Flux_kg_ha_day = CH4Flux_kg_ha_h * 24)     # kg/ha/h -> kg/ha/day

# Running cumulative AUC per plot
cumulative_data2 <- GHGYr2_Full %>%
  group_by(FishTreatment, BirdTreatment, treatment, PlotID) %>%
  arrange(ElapsedDays, .by_group = TRUE) %>%
  mutate(
    # Trapezoid area between each consecutive pair of points
    delta_t     = ElapsedDays - lag(ElapsedDays),
    trap_area   = 0.5 * (CH4Flux_kg_ha_day + lag(CH4Flux_kg_ha_day)) * delta_t,
    trap_area   = ifelse(is.na(trap_area), 0, trap_area),  # first point = 0
    cumulative_CH4 = cumsum(trap_area),
    FishTreatment = recode(FishTreatment, "NoFish" = "No Fish")
  ) %>%
  ungroup()

# coerce date 
cumulative_data2$Date = as.Date(cumulative_data2$Date)

# create dataset for vline of when fish were on fields 
ghg_fish_dates2 <- data.frame(
  fish_days = as.Date(c("2024-12-08", "2025-02-23")))

treatment_colors_4c2 <- c(
  "Birds Fish"    = "#1B80ADFF",
  "NoBirds Fish"  = "#00B398FF",  
  "Birds NoFish"  = "#FF9933FF",  
  "NoBirds NoFish"= "#C24841FF"
)

# Plot:

plot_4c2 <- ggplot() +
  geom_line(data = cumulative_data2,
            aes(x = Date, y = cumulative_CH4, 
                color = treatment, group = PlotID),
            alpha = 0.5, 
            linewidth = 0.5) +  
  scale_color_manual(values = treatment_colors_4c2) +
  new_scale_color() +# individual plots
  stat_summary(data = cumulative_data2,
               aes(x = Date, y = cumulative_CH4,
                   group = FishTreatment,
                   color = FishTreatment),
               fun = mean, 
               geom = "line", 
               linewidth = 1.5, 
               alpha = 1) + # treatment means
  scale_color_manual(values = c("black", "gray50")) +
  labs(title = "Net Cumulative Annual Emissions",
       x     = "Date (Month-Yr)",
       y     = expression(C~-CH[4]~emissions~(kg~ha^{-1}~day^{-1}))) +
  geom_vline(
    data = ghg_fish_dates2,
    aes(xintercept = fish_days),
    linetype = "dashed",
    color = "darkgrey"
  ) +
  scale_x_date(date_labels = "%b-%y", date_breaks = "2 month") +
  theme_minimal_grid(font_size = 12) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    legend.position = "none") +
  ylim(0, 5200)

plot_4c2  


#### Patchwork figure: ####

### Year 1: 

plot_4_1 <- (plot_4a1 | plot_4b1) / plot_4c1

Figure4_1 <- wrap_elements(plot_4_1) +
  labs(tag = "Year 1") +
  theme(
    plot.tag = element_text(size = 16),
    plot.tag.position = "top"
  )

Figure4_1

### Year 2: 

plot_4_2 <- (plot_4a2 | plot_4b2) / plot_4c2

Figure4_2 <- wrap_elements(plot_4_2) +
  labs(tag = "Year 2") +
  theme(
    plot.tag = element_text(size = 16),
    plot.tag.position = "top"
  )

Figure4_2

Figure4a <- Figure4_1 | Figure4_2

#### Create legends below the plot: 

# Dummy plot 1: line legend (Fish = black, No Fish = gray)
line_legend_plot <- ggplot(
  data.frame(x = 1, y = 1, group = factor(c("Fish", "No Fish"), 
                                          levels = c("Fish", "No Fish"))),
  aes(x = x, y = y, color = group)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Fish" = "black", "No Fish" = "gray50"),
                     name = "") +
  guides(color = guide_legend(
    nrow = 1,
    override.aes = list(shape = NA, linewidth = 1.2))) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Dummy plot 2: dot legend (your four treatments)
treatment_colors <- c(
  "Birds Fish"       = "#1B80ADFF",
  "NoBirds Fish"     = "#00B398FF",
  "Birds NoFish"    = "#FF9933FF",
  "NoBirds NoFish"  = "#C24841FF"
)

dot_legend_plot <- ggplot(
  data.frame(x = 1, y = 1, 
             group = factor(c("Birds Fish", "NoBirds Fish", "Birds NoFish", "NoBirds NoFish"),
                            levels = c("Birds Fish", "NoBirds Fish", "Birds NoFish", "NoBirds NoFish"))),
  aes(x = x, y = y, color = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = treatment_colors,
                     labels = c(
                       "Birds Fish"     = "Fish / Birds",
                       "NoBirds Fish"   = "Fish / No Birds",
                       "Birds NoFish"   = "No Fish / Birds",
                       "NoBirds NoFish" = "No Fish / No Birds"),
                     name = "") +
  guides(color = guide_legend(
    nrow = 1,
    override.aes = list(shape = 16, size = 3, alpha = 1,
                        linetype = 0, linewidth = 0))) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Extract legends
line_leg <- get_legend(line_legend_plot)
dot_leg  <- get_legend(dot_legend_plot)

# Stack legends and combine with main patchwork
Figure4 <- Figure4a /  
  wrap_elements(full = line_leg) /
  wrap_elements(full = dot_leg) +
  plot_layout(heights = c(12, 0.1, 0.1))

Figure4

# Export figure 

# Save figure: 
ggsave("figures/Figure_4.png", 
       plot = Figure4,
       width = 30, 
       height = 25, 
       dpi = 300,
       units = "cm")





