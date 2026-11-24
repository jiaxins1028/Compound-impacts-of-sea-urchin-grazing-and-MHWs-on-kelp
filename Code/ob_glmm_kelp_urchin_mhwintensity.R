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

##############################################################
## Read the data output from `mhw_metrics_each_survey.ipynb`

setwd("...")
data <- read.csv("all_ukmhw_with_level.csv") 
data <- data %>% filter(location %in% c('Jervis Bay', 'Batemans')) 
data$survey_year <- factor(data$survey_year)

# hist(data$survey_mean)
# plot(data$number, data$survey_mean)

data$survey_year <- factor(data$survey_year)
data$site_name <- factor(data$site_name)
data$cover_percent <- data$survey_mean/ 100  # convert cover percentage to 0-1
range(data$cover_percent)    

## ----------------- data prep to separating by mhw levels ------------------
data <- data %>%
  mutate(mhw_level = case_when(
    summer_inten <= inten_level1 ~ "Level 0",
    (summer_inten > inten_level1) & (summer_inten <= inten_level2) ~ "Level I",    
    summer_inten >  inten_level2 ~ "Level II"),  
    mhw_level = factor(mhw_level, levels = c("Level 0", "Level I", "Level II")))


sum(is.na(data$mhw_level))    # should be 0
# data <- data %>% filter(mhw_level != "NA")

data %>% group_by(mhw_level) %>%
  summarise(n = n(), n_sites = n_distinct(site_name),
            n_zero = sum(cover_percent == 0), n_one = sum(cover_percent == 1))

## ----------------- Ordered beta fitting with random effect of site -----------------
### logit(μ) = b0 + b1 × urchin number + site effect

levs <- levels(data$mhw_level)

fits <- lapply(levs, function(l)
  glmmTMB(cover_percent ~ number + (1 | site_name),
          family = ordbeta(),
          data = filter(data, mhw_level == l)))
names(fits) <- levs

lapply(fits, summary)    # see whether there is convergence warning

fit_ob_warm     <- fits[["Level II"]]
sigma(fit_ob_warm)
VarCorr(fit_ob_warm)

fit_ob_moderate <- fits[["Level I"]]
sigma(fit_ob_moderate)
VarCorr(fit_ob_moderate)

fit_ob_cool     <- fits[["Level 0"]]
sigma(fit_ob_cool)
VarCorr(fit_ob_cool)
               
               
## ------------ 95% CI and collapse thresholds --------------------
pred_curve <- function(model, dat, level_name) {
  nd <- data.frame(number = seq(min(dat$number), max(dat$number), 0.05))
  
  b <- fixef(model)$cond      # b[1] = b0 (intercept), b[2] = b1 (urchin slope)
  V <- vcov(model)$cond       # 2 x 2 Covariance matrix
  
  nd <- data.frame(number = seq(min(dat$number), max(dat$number), 0.05))
  u  <- nd$number
  
  # logit Predicted values and standard errors on the logit scale
  eta <- b[1] + (b[2]) * u  
  se  <- sqrt(V[1, 1] + u^2 * V[2, 2] + 2 * u * V[1, 2])
  
  nd$preds <- plogis(eta)                         
  nd$lower <- plogis(eta - 1.93 * se)
  nd$upper <- plogis(eta + 1.93 * se)
  nd$level <- level_name
  nd
}

new_data_warm <- pred_curve(fit_ob_warm, data %>% filter(mhw_level == 'Level II'), "Level II") 
new_data_moderate <- pred_curve(fit_ob_moderate, data %>% filter(mhw_level == 'Level I') , "Level I")  
new_data_cool <- pred_curve(fit_ob_cool, data %>% filter(mhw_level == 'Level 0'), "Level 0")  

### ------ collapse thresh -----------------------
threshold <- 10  

get_threshold <- function(nd, target = threshold/100) {
  data.frame(level     = nd$level[1],
             threshold = nd$number[which(nd$preds <= target)[1]],
             lower_ci  = nd$number[which(nd$lower <= target)[1]],
             upper_ci  = nd$number[which(nd$upper <= target)[1]])
}

thr <- rbind(get_threshold(new_data_warm),
             get_threshold(new_data_moderate),
             get_threshold(new_data_cool))
print(thr)

## --------------- plot -----------------------
new_data_all <- rbind(new_data_warm, new_data_moderate, new_data_cool)

lev_order <- c("Level 0", "Level I", "Level II")
new_data_all$level <- factor(new_data_all$level, levels = lev_order)
obs_all$level      <- factor(obs_all$level,      levels = lev_order)
thr$level          <- factor(thr$level,          levels = lev_order)

# Define text labels and their positions
text_labels <- data.frame(
  x = c(thr$threshold[1], thr$threshold[2], thr$threshold[3]),  # x positions
  y = c( 0,  0,  0),  # y positions, adjust as needed
  label = c(thr$threshold[1], thr$threshold[2], thr$threshold[3]),
  level = c('Level II', 'Level I', 'Level 0'))

cols <- c("Level 0" = "#1E8E99FF",    
          "Level I" = "#F28C28",    
          "Level II" = "#9B1B1B")   

p <- ggplot() +
      geom_point(data = data, aes(number, cover_percent * 100, colour = mhw_level), alpha = 0.5, size = 3) +
      geom_ribbon(data = new_data_all,
                  aes(number, ymin = lower * 100, ymax = upper * 100, fill = level), alpha = 0.2) +
      geom_line(data = new_data_all, aes(number, preds * 100, colour = level), linewidth = 1.7) +
    
      geom_hline(yintercept = 10, linetype = "dashed", colour = "dimgrey", linewidth = 1) +
      geom_text(data = text_labels, aes(x = x+0.05, y = y, label = round(label, 1), color = level),
                vjust = 1.5, size = 5, fontface = "bold") +
      geom_segment(data = text_labels, aes(x = x, xend = x, y = y, yend = threshold, color = level),
                   alpha = 0.5, linewidth = 1.2) +
      
      scale_colour_manual(values = cols, name = "Summer Max Intensity") +
      scale_fill_manual(values = cols, name = "Summer Max Intensity") +
      coord_cartesian(ylim = c(-4, 100), clip = "off") +
      
      labs(x = "Sea urchin density (/m²)", y = "Percentage canopy cover (%)") +
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
p

## --------------- convergence and residual --------------------------
# convergence
fit_ob_warm$sdr$pdHess          # TRUE 
fit_ob_warm$fit$convergence     # 0 

# residual diagnostics：DHARMa 
library(DHARMa)
res_cool <- simulateResiduals(fit_ob_cool)
plot(res_cool)                  
testDispersion(res_cool)

res_moderate <- simulateResiduals(fit_ob_moderate)
plot(res_moderate)
testDispersion(res_moderate)

res_warm <- simulateResiduals(fit_ob_warm)
plot(res_warm)
testDispersion(res_warm)


               



