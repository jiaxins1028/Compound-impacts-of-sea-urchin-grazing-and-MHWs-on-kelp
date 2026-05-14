######################################################################
##   This script models the relationship between sea urchin         ##
##   density (predictor) and kelp canopy cover (response),          ##
##   stratified by summer marine heatwave (MHW) intensity levels    ##
######################################################################

library(ggplot2)
library(gamair)
library(dplyr)
library(mgcv)
library(statmod) 
library(nls.multstart)
library(broom)
library(tidyverse)
library(pracma)
library(qgam)
library(lme4)
library(lmerTest)
library(MASS)     
library(pscl)     
library(reshape2) 
library(AER) 
library(visreg)
library(patchwork)
library(glmmTMB)
library(stringr) 

theme_format <- theme_bw()+
  theme(axis.text.x  = element_text(vjust=0.5,size=13, colour = "black"))+
  theme(axis.text.y  = element_text(size=15, colour = "black"))+
  theme(axis.title.x = element_text(size=18, colour = "black"))+
  theme(axis.title.y = element_text(size=18, colour = "black"))+
  #panel.background = element_rect(fill="white"),
  theme(axis.ticks = element_line(colour="black"))+
  theme(panel.grid.minor=element_blank())+
  theme(panel.grid.major=element_blank())+
  theme(strip.text = element_text(size = 15))+
  theme(legend.text = element_text(size = 15),  
        legend.title = element_text(size = 15))

setwd("...")

# Read the data
data <- read.csv("all_ukmhw_with_level.csv") 
data <- data %>% filter(location %in% c('Jervis Bay', 'Batemans')) 
data$survey_year <- factor(data$survey_year)

# hist(data$survey_mean)
# plot(data$number, data$survey_mean)

data_warm <- data %>% filter(summer_inten > inten_level2 & summer_inten < inten_level3) 
data_warm$level <- 'Level I'
data_moderate <- data %>% filter(summer_inten > inten_level1 & summer_inten <= inten_level2)
data_moderate$level <- 'Level II'
data_cool <- data %>% filter(summer_inten < inten_level1)
data_cool$level <- 'Level 0'

## --------------- GLMM model fitting--------------------------------------------

fit_nb_warm <- glmmTMB(survey_mean ~ number + (1 | site_name/survey_year), data = data_warm, family = nbinom2)
summary(fit_nb_warm)
fit_nb_moderate <- glmmTMB(survey_mean ~ number + (1 | site_name/survey_year), data = data_moderate, family = nbinom2)
summary(fit_nb_moderate)
fit_nb_cool <- glmmTMB(survey_mean ~ number + (1 | site_name/survey_year), data = data_cool, family = nbinom2)
summary(fit_nb_cool)


# predict new data
# Calculate the upper and lower bounds of the confidence intervals
alpha <- 0.05  # 95% confidence interval
z_value <- qnorm(1 - alpha/2)  # Z-score for 95% CI

inv_logit <- function(x) exp(x) / (1 + exp(x))

# For GLMMs with a log link and random intercepts, the population-averaged
# mean is approximately exp(eta + 0.5 * sigma^2_re). This correction accounts
# for the log-normal distribution induced by random effects and yields 
# population-level predictions across sites.
predict_marginal_nb <- function(fit, newdata, group = "site_name") {
  preds <- predict(fit, newdata, type = "link", se.fit = TRUE, re.form = NA)
  re_var <- as.numeric(VarCorr(fit)$cond[[group]][1])
  
  list(fit = preds$fit + 0.50 * re_var, se.fit = preds$se.fit)
}

# For Level II (Warm)
new_data_warm <- expand.grid(number = seq(min(data_warm$number), max(data_warm$number), 0.05)) 
# Get the predictions and standard errors on the link scale
preds_warm <- predict_marginal_nb(fit_nb_warm, new_data_warm)
# Add the predictions and standard errors to the new data frame
new_data_warm$fit <- preds_warm$fit
new_data_warm$se.fit <- preds_warm$se.fit
# Calculate the predicted values on the response scale (inverse of the log link)
new_data_warm$preds <- exp(new_data_warm$fit)
# Calculate the confidence intervals on the response scale
new_data_warm$lower <- exp(new_data_warm$fit - 1.96 * new_data_warm$se.fit)
new_data_warm$upper <- exp(new_data_warm$fit + 1.96 * new_data_warm$se.fit)
new_data_warm$level <- 'Level II'


# For Level I (Moderate)
new_data_moderate <- expand.grid(number = seq(min(data_moderate$number), max(data_moderate$number), 0.05))
preds_moderate <- predict_marginal_nb(fit_nb_moderate, new_data_moderate)
new_data_moderate$fit <- preds_moderate$fit
new_data_moderate$se.fit <- preds_moderate$se.fit
new_data_moderate$preds <- exp(new_data_moderate$fit)
# Calculate the confidence intervals on the response scale
new_data_moderate$lower <- exp(new_data_moderate$fit - 1.96 * new_data_moderate$se.fit)
new_data_moderate$upper <- exp(new_data_moderate$fit + 1.96 * new_data_moderate$se.fit)
new_data_moderate$level <- 'Level I'


# For Level 0 (Cool)
new_data_cool <- expand.grid(number = seq(min(data_cool$number), max(data_cool$number), 0.05))
preds_cool <- predict_marginal_nb(fit_nb_cool, new_data_cool)
new_data_cool$fit <- preds_cool$fit
new_data_cool$se.fit <- preds_cool$se.fit
new_data_cool$preds <- exp(new_data_cool$fit)
# Calculate the confidence intervals on the response scale
new_data_cool$lower <- exp(new_data_cool$fit - 1.96 * new_data_cool$se.fit)
new_data_cool$upper <- exp(new_data_cool$fit + 1.96 * new_data_cool$se.fit)
new_data_cool$level <- 'Level 0'


# --------------------------  Collapse threshold ----------------------------------------------------------
# Define the threshold for the preds values
threshold <- 10  # Replace with your desired threshold value

warm_threshold <- round(new_data_warm[new_data_warm$preds <= threshold, "number"], 1)
moderate_threshold <- round(new_data_moderate[new_data_moderate$preds <= threshold, "number"], 1)
cool_threshold <- round(new_data_cool[new_data_cool$preds <= threshold, "number"], 1)

# Print the results
print(warm_threshold[1])
print(moderate_threshold[1])
print(cool_threshold[1])

# Define text labels and their positions
text_labels <- data.frame(
  x = c(warm_threshold[1], moderate_threshold[1], cool_threshold[1]),  # x positions
  y = c( 0,  0,  0),  # y positions, adjust as needed
  label = c(warm_threshold[1], moderate_threshold[1], cool_threshold[1]),
  level = c('Level II', 'Level I', 'Level 0'))


## -------------------------Plotting ----------------------------------------------------
library(ggnewscale)

ggplot() +
  labs(x= expression('Sea urchin density /m'^2), y="Percentage Canopy Cover/%", title = '')+
  
  # points with their own colors, no legend
  geom_point(data = data_warm, aes(x=number, y=survey_mean, color = level),
             shape=16, size = 3, alpha = 0.5, show.legend = FALSE)+
  geom_point(data = data_moderate, aes(x=number, y=survey_mean, color = level),
             size = 3, alpha = 0.5, show.legend = FALSE)+
  geom_point(data = data_cool, aes(x=number, y=survey_mean, color = level),
             size = 3, alpha = 0.5, show.legend = FALSE)+
  
  # lines + ribbons + legend keep the original level colors
  geom_line(data = new_data_warm, aes(x = number, y = preds, color = level), linewidth = 2) +
  geom_line(data = new_data_moderate, aes(x = number, y = preds, color = level), linewidth = 2) +
  geom_line(data = new_data_cool, aes(x = number, y = preds, color = level), linewidth = 2) +
  geom_ribbon(data = new_data_warm, aes(x = number, ymin = lower, ymax = upper, fill = level), alpha = 0.2) +
  geom_ribbon(data = new_data_moderate, aes(x = number, ymin = lower, ymax = upper, fill = level), alpha = 0.2) +
  geom_ribbon(data = new_data_cool, aes(x = number, ymin = lower, ymax = upper, fill = level), alpha = 0.2) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "dimgray", linewidth = 1) +
  geom_text(data = text_labels, aes(x = x+0.05, y = y, label = label, color = level),
            vjust = 1.5, size = 5, fontface = "bold") +
  geom_segment(data = text_labels, aes(x = x, xend = x, y = y, yend = threshold, color = level),
               alpha = 0.5, linewidth = 1.2) +
  annotate("text", x = 0.5, y = threshold, label = threshold,
           size = 5, color = "black", fontface = "bold", hjust = 2.7) +
  scale_color_manual(values = c(
    'Level 0' = '#1E8E99FF',
    'Level I' = '#FF8E32FF',
    'Level II' = '#950404FF')) +
  scale_fill_manual(values = c(
    'Level 0' = '#1E8E99FF',
    'Level I' = '#FF8E32FF',
    'Level II' = '#950404FF')) +
  labs(color = "Summer Max Intensity") +
  scale_y_continuous(limits=c(0,100), expand = c(0, 0)) +
  scale_x_continuous(limits=c(0,10), expand = c(0, 0)) +
  theme_format+
  theme(legend.position = c(.97, .97),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.ticks.length.x = unit(10, "pt"),
        legend.box.background = element_rect(color="black", linewidth=2))+
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(override.aes = list(linetype = 1, size = 1, shape = NA)),
         fill = FALSE)


## ================= Confidence Interval extraction for each curve =================================
threshold <- 10  

warm_threshold <- round(new_data_warm[new_data_warm$preds <= threshold, "number"], 1)
warm_threshold_lower <- round(new_data_warm[new_data_warm$lower <= threshold, "number"], 1)
warm_threshold_upper <- round(new_data_warm[new_data_warm$upper <= threshold, "number"], 1)

moderate_threshold <- round(new_data_moderate[new_data_moderate$preds <= threshold, "number"], 1)
moderate_threshold_lower <- round(new_data_moderate[new_data_moderate$lower <= threshold, "number"], 1)
moderate_threshold_upper <- round(new_data_moderate[new_data_moderate$upper <= threshold, "number"], 1)

cool_threshold <- round(new_data_cool[new_data_cool$preds <= threshold, "number"], 1)
cool_threshold_lower <- round(new_data_cool[new_data_cool$lower <= threshold, "number"], 1)
cool_threshold_upper <- round(new_data_cool[new_data_cool$upper <= threshold, "number"], 1)

# Print the results
print(warm_threshold[1])
print(moderate_threshold[1])
print(cool_threshold[1])

# Combine in a nice table
threshold_table <- data.frame(
  Level   = c("Warm","Moderate","Cool"),
  Estimate= c(warm_threshold[1], moderate_threshold[1], cool_threshold[1]),
  LowerCI = c(warm_threshold_lower[1], moderate_threshold_lower[1], cool_threshold_lower[1]),
  UpperCI = c(warm_threshold_upper[1], moderate_threshold_upper[1], cool_threshold_upper[1])
)

threshold_table
