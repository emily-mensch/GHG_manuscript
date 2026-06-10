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
daphnia <- read.csv("data/DaphniaFull.csv") # Full daphnia dataset

# Figure S1 ---------------------------------------------------------------

#### Daphnia abundance over time: 
# format date: 
daphnia$SampleDate=as.Date(x=daphnia$SampleDate, format="%m/%d/%Y")

# Create data frame for adding vline for fish additions: 
fish_dates <- data.frame(
  StudyYear = c("Year1", "Year2"),
  fish_date = as.Date(c("2023-12-13", "2024-12-08"))
)

treatment_colors_S1 <- c(
  "Birds Fish"    = "#1B80ADFF",
  "NoBirds Fish"  = "#00B398FF",  
  "Birds NoFish"  = "#FF9933FF",  
  "NoBirds NoFish"= "#C24841FF"
) # treatment colors for jitter

skip_near_start <- function(x, min_gap_days = 15) {
  # generate monthly breaks spanning the panel's actual range
  brks <- seq(
    lubridate::floor_date(min(x), "month"),
    lubridate::ceiling_date(max(x), "month"),
    by = "1 month"
  )
  # drop any break that's too close to the left edge of this panel
  brks[as.numeric(difftime(brks, min(x), units = "days")) > min_gap_days]
} 

# By fish treatment:
plotS1 <- ggplot() +
  geom_jitter(
    data = daphnia, 
    aes(x = SampleDate, y = total, color = treatment), 
    width = 0.2, alpha = 0.5, size = 1.5
    ) +
  scale_color_manual(
    values = treatment_colors_S1, name = "Treatment"
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

plotS1

