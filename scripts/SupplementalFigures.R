#### Fish/bird trophic effects on methane: manuscript 
#### Target journal: PNAS 
#### Figures, supplementary materials

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
zoop <- read.csv("data/TotalZoop.csv") # Full zooplankton dataset
daphnia <- read.csv("data/DaphniaFull.csv") # Full daphnia dataset
totaldaphnia1 <- read.csv("data/totaldaphnia1.csv") # Total daphnia dataset year 1
hobo_summary <- read.csv("data/hobo_summary.csv") # Full water temperature/dissolved oxygen dataset 
depth <- read.csv("data/depth_summary.csv") # Full depth dataset 
birds <- read.csv("data/BirdPointCounts.csv") # Point count dataset

# Figure S1 ---------------------------------------------------------------

#### Showcasing final models for total zooplankton abundance, years 1 and 2

#### Year 1: ####
totalzoop1 <- zoop %>% 
  filter(StudyYear == "Year1")

# Final model: additive effects of fish and pre/post 
P_S1a <- glmmTMB(
  DensityRounded ~ FishTreatment + Pre_Post + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totalzoop1
)

emmS1a <- emmeans(P_S1a, ~ Pre_Post,
                  type = "response")

contrast(emmS1a, 
         method = "pairwise")

emm_df_S1a <- as.data.frame(emmS1a)

datS1a <- totalzoop1 %>% 
  mutate(Pre_Post = factor(Pre_Post, levels = c("Pre", "Post")),
         FishTreatment = recode(FishTreatment, "NoFish" = "No Fish")) 

emm_df_S1a <- emm_df_S1a %>% 
  mutate(
    Pre_Post = factor(Pre_Post, levels = c("Pre", "Post")))


treatment_colors_S1a <- c(
  "Birds Fish"    = "#1B80ADFF",
  "NoBirds Fish"  = "#00B398FF",  
  "Birds NoFish"  = "#FF9933FF",  
  "NoBirds NoFish"= "#C24841FF"
)

plot_S1a <- ggplot() +
  geom_quasirandom( 
    data = datS1a, aes(x = Pre_Post, y = total, 
                      color = treatment),
    width = 0.2, alpha = 0.5, size = 1.5
  ) +
  scale_color_manual(values = treatment_colors_S1a, name = "Treatment") + 
  new_scale_color() +
  geom_line(
    data = emm_df_S1a, aes(x = Pre_Post, y = response, group = 1) ,
    color = "black",
    linewidth = 0.5,
    alpha = 0.5
  ) +
  geom_point(
    data = emm_df_S1a, aes(x = Pre_Post, y = response),
    color = "black",
    size = 3.5,
    alpha = 0.6
  ) +
  geom_errorbar(
    data = emm_df_S1a, aes(x = Pre_Post, ymin = asymp.LCL, ymax = asymp.UCL),
    color = "black",
    width = 0.15,
    alpha = 0.5
  ) +
  scale_y_log10() +
  coord_cartesian(ylim = c(300, 34000)) +
  labs(x = "", 
       y = expression("Zooplankton Abundance (individuals m"^-3*")")) +
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
  geom_text(aes(x = 1.5, y = 2700, 
                label = "*"), 
            size = 7)

plot_S1a

#### Year 2: ####
totalzoop2 <- zoop %>% 
  filter(StudyYear == "Year2")

# Final model: three-way interaction
P_S1b <- glmmTMB(
  DensityRounded ~ FishTreatment * BirdTreatment * Pre_Post + 
    scale(MeanTemp) + scale(MeanDO) + scale(MeanDepth) +
    (1 | PlotID )  + (1 | SamplingOccasion),
  family = nbinom2,
  data = totalzoop2
) 

emmS1b <- emmeans(P_S1b, ~ Pre_Post | FishTreatment * BirdTreatment,
                 type = "response")

contrast(emmS1b, 
         method = "pairwise", 
         by = c("FishTreatment", "BirdTreatment"))

emm_df_S1b <- as.data.frame(emmS1b)

datS1b <- totalzoop2 %>% 
  mutate(
    Pre_Post = factor(Pre_Post, levels = c("Pre", "Post")),
    FishTreatment = recode(FishTreatment, "NoFish" = "No Fish")
  )

emm_df_S1b <- emm_df_S1b %>% 
  mutate(
    Pre_Post = factor(Pre_Post, levels = c("Pre", "Post")),
    FishTreatment = recode(FishTreatment, "NoFish" = "No Fish")
  )

treatment_order_S1b <- c("Fish_Birds", "Fish_NoBirds", "No Fish_Birds", "No Fish_NoBirds")

treatment_colors_S1b <- c(
  "Fish_Birds"       = "#1B80ADFF",
  "Fish_NoBirds"     = "#00B398FF",
  "No Fish_Birds"    = "#FF9933FF",
  "No Fish_NoBirds"  = "#C24841FF"
)

label_S1b <- data.frame(
  FishTreatment = c("Fish", "Fish", "No Fish", "No Fish"),
  x = c(1.4, 1.4, 1.4, 1.5),
  y = c(11000, 6500, 9500, 5500) , 
  label = c("n.s.", "*", "*", "*")
)

datS1b <- datS1b %>%
  mutate(color_group = factor(
    interaction(FishTreatment, BirdTreatment, sep = "_"),
    levels = treatment_order_S1b))

emm_df_S1b <- emm_df_S1b %>%
  mutate(color_group = factor(
    interaction(FishTreatment, BirdTreatment, sep = "_"),
    levels = treatment_order_S1b))


plot_S1b <- ggplot() +
  geom_quasirandom( 
    data = datS1b, aes(x = Pre_Post, y = total, color = color_group),
    width = 0.2, alpha = 0.5, size = 1.5
  ) +
  geom_line(
    data = emm_df_S1b,
    aes(x = Pre_Post, y = response, color = color_group, group = color_group),
    linewidth = 0.5, position = position_dodge(width = 0.5)
  ) +
  geom_point(
    data = emm_df_S1b, aes(x = Pre_Post, y = response, color = color_group),
    size = 3.5, position = position_dodge(width = 0.5), alpha = 0.8
  ) +
  geom_errorbar(
    data = emm_df_S1b, aes(x = Pre_Post, ymin = asymp.LCL, ymax = asymp.UCL,
                         color = color_group),
    width = 0.15, position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~FishTreatment, strip.position = "top") +
  scale_color_manual(
    values = treatment_colors_S1b,
    name = "Treatment",
    labels = c(
      "Fish_Birds"     = "Fish / Birds",
      "Fish_NoBirds"   = "Fish / No Birds",
      "No Fish_Birds"   = "No Fish / Birds",
      "No Fish_NoBirds" = "No Fish / No Birds"
    )) +
  scale_y_log10() +
  coord_cartesian(ylim = c(300, 34000)) +
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
  geom_text(data = label_S1b, 
            aes(x = x, y = y, label = label), 
            inherit.aes = FALSE, 
            size = c(4, 7, 7, 7))

plot_S1b

#### Put figures together: 

legend_grob <- get_legend(plot_S1b)  # extract legend from plot 2b

# Remove legends from individual plots
plot_S1a_noleg <- plot_S1a + theme(legend.position = "none")
plot_S1b_noleg <- plot_S1b + theme(legend.position = "none")

# combine using patchwork w/ legend as third row: 

FigureS1 <- (plot_S1a_noleg | plot_S1b_noleg) /
  ((plot_spacer() | wrap_elements(full = legend_grob) | plot_spacer()) + 
     plot_layout(widths = c(0.21, 2, 1))) +
  plot_layout(heights = c(8, 0.1))

FigureS1

# Export figure 

# Save figure: 
ggsave("figures/Figure_S1.png", 
       plot = FigureS1,
       width = 20, 
       height = 18, 
       dpi = 300,
       units = "cm")


# Figure S2 ---------------------------------------------------------------

## Showcasing daphnia abundance over time 

#### Daphnia abundance over time: 
# format date: 
daphnia$SampleDate=as.Date(x=daphnia$SampleDate, format="%m/%d/%Y")

# Create data frame for adding vline for fish additions: 
fish_dates <- data.frame(
  StudyYear = c("Year1", "Year2"),
  fish_date = as.Date(c("2023-12-13", "2024-12-08"))
)

treatment_colors_S2 <- c(
  "Birds Fish"    = "#1B80ADFF",
  "NoBirds Fish"  = "#00B398FF",  
  "Birds NoFish"  = "#FF9933FF",  
  "NoBirds NoFish"= "#C24841FF"
) # treatment colors for jitter

# function to manually set dates shown on x-axis:
skip_near_start <- function(x, min_gap_days = 15) {
  brks <- seq(
    lubridate::floor_date(min(x), "month"),
    lubridate::ceiling_date(max(x), "month"),
    by = "1 month"
  )
  brks[as.numeric(difftime(brks, min(x), units = "days")) > min_gap_days]
} 

# By fish treatment:
plotS2 <- ggplot() +
  geom_quasirandom(
    data = daphnia, 
    aes(x = SampleDate, y = total, color = treatment), 
    alpha = 0.5, size = 1.5
    ) +
  scale_color_manual(
    values = treatment_colors_S2, name = "Treatment"
    ) + 
  new_scale_color() +
  geom_smooth(
    data = daphnia, 
    aes(x = SampleDate, y = total, group = FishTreatment, color = FishTreatment), 
    method = "loess", span = 0.5, se = FALSE, size = 1.2
    ) +
  scale_color_manual(values = c("black", "gray50")) +
  facet_wrap(~StudyYear, scales = "free_x", strip.position = "top") +
  scale_x_date(breaks = skip_near_start, date_labels = "%b-%y") +
  geom_vline(
    data = fish_dates,
    aes(xintercept = fish_date),
    linetype = "dashed",
    color = "grey"
  ) +
  scale_y_log10() +
  coord_cartesian(ylim = c(20, 15500)) +
  labs(x = "Date (Month-Yr)", 
       y = expression("Daphnia Abundance (individuals m"^-3*")")) +
  ggtitle("") +
  theme_minimal_grid(font_size = 13) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.justification = "center",
    legend.margin = margin(t = 0), # tighten gap between stacked legends
    legend.spacing.y = unit(0.1, "cm"),
    legend.title = element_blank()) 

plotS2

ggsave("figures/Figure_S2.png", 
       plot = plotS2,
       width = 15, 
       height = 11, 
       dpi = 300,
       units = "cm")

# Figure S3 ---------------------------------------------------------------

## Showcasing final model results for year 1 daphnia abundance, including the final sampling date: 

# Year 1, top model:
P_S3 <- glmmTMB(
  DensityRounded ~ FishTreatment * Pre_Post + 
    MeanTempScaled + MeanDOScaled + MeanDepthScaled +
    (1 | PlotID ) + (1 | SamplingOccasion),
  family = nbinom2,
  data = totaldaphnia1
)

emmS3 <- emmeans(P_S3, ~ FishTreatment * Pre_Post,
                 type = "response")

pairs(emmS3, by = "FishTreatment")

emm_df_S3 <- as.data.frame(emmS3)

datS3 <- totaldaphnia1 %>% 
  mutate(Pre_Post = factor(Pre_Post, levels = c("Pre", "Post")),
         FishTreatment = recode(FishTreatment, "NoFish" = "No Fish")) 

emm_df_S3 <- emm_df_S3 %>% 
  mutate(
    Pre_Post = factor(Pre_Post, levels = c("Pre", "Post")),
    FishTreatment = recode(FishTreatment, "NoFish" = "No Fish"))

treatment_colors_S3 <- c(
  "Birds Fish"    = "#1B80ADFF",
  "NoBirds Fish"  = "#00B398FF",  
  "Birds NoFish"  = "#FF9933FF",  
  "NoBirds NoFish"= "#C24841FF"
)

plot_S3 <- ggplot() +
  geom_quasirandom( 
    data = datS3, aes(x = Pre_Post, y = total, 
                      color = treatment),
    width = 0.2, alpha = 0.5, size = 1.5
  ) +
  scale_color_manual(values = treatment_colors_S3, name = "Treatment") + 
  new_scale_color() +
  geom_line(
    data = emm_df_S3, aes(x = Pre_Post, y = response,
                         group = FishTreatment) ,
    color = "black",
    linewidth = 0.5,
    alpha = 0.5
  ) +
  geom_point(
    data = emm_df_S3, aes(x = Pre_Post, y = response),
    color = "black",
    size = 3.5,
    alpha = 0.6
  ) +
  geom_errorbar(
    data = emm_df_S3, aes(x = Pre_Post, ymin = asymp.LCL, ymax = asymp.UCL),
    color = "black",
    width = 0.15,
    alpha = 0.5
  ) +
  facet_wrap(~ FishTreatment, strip.position = "top") +
  scale_y_log10() +
  coord_cartesian(ylim = c(20, 15500)) +
  labs(x = "", 
       y = expression("Daphnia Abundance (individuals m"^-3*")")) +
  ggtitle("") +
  theme_minimal_grid(font_size = 13) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.justification = "center") 

plot_S3

# Save figure: 
ggsave("figures/Figure_S3.png", 
       plot = plot_S3,
       width = 15, 
       height = 11, 
       dpi = 300,
       units = "cm")


# Figure S4 ---------------------------------------------------------------

#### Abiotic data: 

#### Water temperature ####

hobo_summary$treatment <- as.factor(hobo_summary$treatment)

hobo_summary <- hobo_summary %>%
  mutate(treatment = factor(treatment, levels = c("Birds Fish", "NoBirds Fish", "Birds NoFish",
                                                  "NoBirds NoFish")))


treatment_colors_S4a <- c(
  "Birds Fish"       = "#1B80ADFF",
  "NoBirds Fish"     = "#00B398FF",
  "Birds NoFish"    = "#FF9933FF",
  "NoBirds NoFish"  = "#C24841FF"
)

plot_S4a <- ggplot(hobo_summary, aes(x = StudyYear, y = mean_temp, fill = treatment)) +
  geom_boxplot(width = 0.5,
               alpha = 0.7) +
  scale_fill_manual(values = treatment_colors_S4a, name = "") +
  theme_minimal_grid(font_size = 13) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  scale_x_discrete(position = "top") +
  labs(x = "", 
       y = expression("Average Water Temperature" ~ degree*C)) +
  ggtitle("StudyYear") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "none") 


plot_S4a


#### Dissolved oxygen ####

plot_S4b <- ggplot(hobo_summary, aes(x = StudyYear, y = mean_DO, fill = treatment)) +
  geom_boxplot(width = 0.5,
               alpha = 0.7) +
  scale_fill_manual(values = treatment_colors_S4a, name = "") +
  theme_minimal_grid(font_size = 13) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  scale_x_discrete(position = "top") +
  labs(x = "", 
       y = expression("Average Dissolved Oxygen, mg/L")) +
  ggtitle("") +
  theme(axis.text.x = element_blank(),
        legend.position = "none") 


plot_S4b

#### Water depth ####

depth$treatment <- as.factor(depth$treatment)

depth <- depth %>%
  mutate(treatment = factor(treatment, levels = c("Birds Fish", "NoBirds Fish", "Birds NoFish",
                                                  "NoBirds NoFish")))


treatment_colors_S4c <- c(
  "Birds Fish"       = "#1B80ADFF",
  "NoBirds Fish"     = "#00B398FF",
  "Birds NoFish"    = "#FF9933FF",
  "NoBirds NoFish"  = "#C24841FF"
)

plot_S4c <- ggplot(depth, aes(x = StudyYear, y = mean_depth, fill = treatment)) +
  geom_boxplot(width = 0.5,
               alpha = 0.7) +
  scale_fill_manual(values = treatment_colors_S4c, name = "") +
  theme_minimal_grid(font_size = 13) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  scale_x_discrete(position = "top") +
  labs(x = "", 
       y = expression("Average Depth, cm")) +
  ggtitle("") +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.justification = "center") 

plot_S4c

#### Putting figures together:

Figure_S4 <- plot_S4a / plot_S4b / plot_S4c

# Save figure: 
ggsave("figures/Figure_S4.png", 
       plot = Figure_S4,
       width = 16, 
       height = 24, 
       dpi = 300,
       units = "cm")


# Figure S5 ---------------------------------------------------------------

#### Birds by year 
birds_summary <- birds %>% 
  group_by(PlotNumber, StudyYear, Treatment) %>%
  summarise(MeanPresence = mean(PiscPresenceAbsence))


# Final model: additive effects of fish and pre/post 
P_S5 <- glmmTMB(PiscPresenceAbsence ~ Treatment + StudyYear +
                   scale(MeanDepth) +
                   (1 | PlotNumber ),
                 family = binomial,
                 data = birds)

summary(P_S5)

emmS5 <- emmeans(P_S5, ~ StudyYear,
                  type = "response")


emm_df_S5 <- as.data.frame(emmS5)

treatment_colors_S5 <- c(
  "Fish"    = "#1B80ADFF",
  "NoFish"  = "#FF9933FF"
)

plot_S5 <- ggplot() +
  geom_quasirandom( 
    data = birds_summary, aes(x = StudyYear, y = MeanPresence, 
                       color = Treatment),
    width = 0.2, alpha = 0.5, size = 1.5
  ) +
  scale_color_manual(values = treatment_colors_S5, name = "Treatment") + 
  new_scale_color() +
  geom_line(
    data = emm_df_S5, aes(x = StudyYear, y = prob, group = 1) ,
    color = "black",
    linewidth = 0.5,
    alpha = 0.5
  ) +
  geom_point(
    data = emm_df_S5, aes(x = StudyYear, y = prob),
    color = "black",
    size = 3.5,
    alpha = 0.6
  ) +
  geom_errorbar(
    data = emm_df_S5, aes(x = StudyYear, ymin = asymp.LCL, ymax = asymp.UCL),
    color = "black",
    width = 0.15,
    alpha = 0.5
  ) +
  labs(x = "", 
       y = expression("Probability of Bird Presence")) +
  ggtitle("") +
  theme_minimal_grid(font_size = 13) +
  background_grid(major = c("xy"),
                  minor = c("xy"),
                  color.major = alpha("grey85", 0.4), 
                  color.minor = alpha("grey85", 0.3) ) +
  panel_border() +
  scale_x_discrete(position = "top") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.justification = "center") +
  geom_text(aes(x = 1.5, y = 0.09, 
                label = "*"), 
            size = 7)

plot_S5

# Save figure: 
ggsave("figures/Figure_S5.png", 
       plot = plot_S5,
       width = 11, 
       height = 11, 
       dpi = 300,
       units = "cm",
       bg = "white")

