library(ggplot2)
library(gamair)
library(dplyr)
library(mgcv)
# library(tidymv)
library(statmod) 
library(broom)
library(tidyverse)
library(pracma)
library(qgam)
library(lme4)
library(lmerTest)
library(MASS)     #for negative bin regression
library(pscl)     #for zero inflated models & predprob
library(reshape2) #convert wide to tall
library(AER) 
library(caTools)
library(caret)
library(pdp)
library(plotly)
library(car)
library(effects)
library(metR)
library(scico)
library(akima)
library(MuMIn)
library(glmmTMB)

theme_format <- theme_bw()+
  theme(axis.text.x  = element_text(vjust=0.5,size=16, colour = "black"))+
  theme(axis.text.y  = element_text(size=16, colour = "black"))+
  theme(axis.title.x = element_text(size=17, colour = "black"))+
  theme(axis.title.y = element_text(size=17, colour = "black"))+
  #panel.background = element_rect(fill="white"),
  theme(axis.ticks = element_line(colour="black"))+
  theme(panel.grid.minor=element_blank())+
  theme(panel.grid.major=element_blank())+
  theme(strip.text = element_text(size = 14))+
  theme(legend.text = element_text(size = 17),  
        legend.title = element_text(size = 16))


## read data

setwd("...")
data <- read.csv('all_urchin_kelp_mhw.csv')  ## the output from `mhw_metrics_each_survey.ipynb`

locations <- c('Batemans', 'Jervis Bay', 'Maria Island', 'Bicheno', 'Kent Group')
data <- data %>% filter(location %in% locations)

data <- data[c("location","site_name" ,"longitude","latitude", "survey_mean", 'survey_year',
               "number","summer_temp",'summer_inten',"summer_dhd50", 'program')] 

data$program <-factor(data$program )
data$location <- as.factor(data$location)
data$site_name <- as.factor(data$site_name)
head(data)

## --------------- convert to binomial ------------------------------
data$collapse <- ifelse(data$survey_mean > 10, 0, 1)
# Convert collapse to a factor
data$collapse <- as.factor(data$collapse)
data$collapse <- factor(data$collapse, levels = c(0, 1), labels = c("class_0", "class_1"))
head(data)

cor(data[, c("number", "summer_temp", "summer_inten", "summer_dhd50")], method = "pearson")

## ----  Diagnose before refitting -------------------------------------

# Intraclass correlation: how much of the variation is between sites
m_null <- glmmTMB(collapse ~ 1 + (1 | site_name), family = binomial, data = data_jb)
performance::icc(m_null)


## ---- Centering ---------------------------------------------------------
# Main effects in an interaction model are conditional: the coefficient for
# one predictor is its effect when the other equals zero. Zero is outside
# the observed range for both predictors here, so the 'uncentred' main effects
# are not interpretable. 'Centring' moves the reference point to the mean.
# AIC, logLik and model ranking are unaffected.


data$number_c <- as.numeric(scale(data$number,       scale = FALSE))
data$inten_c  <- as.numeric(scale(data$summer_inten, scale = FALSE))
data$temp_c   <- as.numeric(scale(data$summer_temp,  scale = FALSE))
data$dhd_c    <- as.numeric(scale(data$summer_dhd50, scale = FALSE))

# keep the means - you need them to back-transform for figures
means <- c(number = mean(data$number),
           inten  = mean(data$summer_inten),
           temp   = mean(data$summer_temp),
           dhd    = mean(data$summer_dhd50))


data_jb <- data %>% filter(location %in% c('Jervis Bay', 'Batemans'))

data_jb$number_c <- as.numeric(scale(data_jb$number,       scale = FALSE))
data_jb$inten_c  <- as.numeric(scale(data_jb$summer_inten, scale = FALSE))
data_jb$temp_c   <- as.numeric(scale(data_jb$summer_temp,  scale = FALSE))
data_jb$dhd_c    <- as.numeric(scale(data_jb$summer_dhd50, scale = FALSE))

means_jb <- c(number = mean(data_jb$number),
              inten  = mean(data_jb$summer_inten),
              temp   = mean(data_jb$summer_temp),
              dhd    = mean(data_jb$summer_dhd50))

## ---- Candidate models ----------------------------------------------

f <- function(rhs) {
  glmmTMB(as.formula(paste("collapse ~", rhs, " + (1 | site_name)")),  ## + program
          family = binomial, data = data_jb)}

m1  <- f("number_c")
m2  <- f("temp_c")
m3  <- f("inten_c")
m4  <- f("dhd_c")
m5  <- f("number_c + temp_c")
m6  <- f("number_c + inten_c")
m7  <- f("number_c + dhd_c")
m8  <- f("number_c * temp_c")
m9  <- f("number_c *inten_c")  # number_c * inten_c
m10 <- f("number_c * dhd_c")

mods <- list(m1=m1, m2=m2, m3=m3, m4=m4, m5=m5,
             m6=m6, m7=m7, m8=m8, m9=m9, m10=m10)

sapply(mods, nobs)                              # Check n samples
sapply(mods, function(m) m$fit$convergence)     # should be 0

MuMIn::model.sel(mods)

df_aic <- as.data.frame(MuMIn::model.sel(mods))
df_aic

### ----------------- Model parameters  -------------------------
# Summarise fit statistics for each candidate collapse model in mods

# Null model: intercept only + site random intercept.
# Used as the baseline for McFadden's pseudo-R^2.

null_A <- glmmTMB(collapse ~ 1 + (1 | site_name),
                  family = binomial, data = data)

summ <- data.frame(
  model    = names(mods),
  n        = sapply(mods, nobs),
  sigma2   = sapply(mods, function(m) as.numeric(VarCorr(m)$cond$site_name)),
  logLik   = sapply(mods, function(m) as.numeric(logLik(m))),
  AICc     = sapply(mods, MuMIn::AICc),
  mcfadden = sapply(mods, function(m)
    1 - as.numeric(logLik(m)) / as.numeric(logLik(null_A)))
)
summ

df_para <- df_aic %>%
  rownames_to_column("model") %>%
  left_join(summ, by = "model")

df_para

### -------------------- collapse threshold sensitivity ----------------
data_jb_st <- data_jb

fit_cut <- function(cut) {
  data_jb_st$collapse <- ifelse(data_jb_st$survey_mean > cut, 0, 1)
  # Convert collapse to a factor
  data_jb_st$collapse <- as.factor(data_jb_st$collapse)
  data_jb_st$collapse <- factor(data_jb_st$collapse, levels = c(0, 1), labels = c("class_0", "class_1"))
  
  m9_st <- glmmTMB(collapse ~ number_c * inten_c + (1 | site_name), family = binomial, data = data_jb_st)
  print(c(cut = cut, n_collapse = sum(data_jb_st$collapse == "class_1")))
  print(summary(m9_st)$coefficients$cond)
}
fit_cut(5); fit_cut(10); fit_cut(15); fit_cut(20)


## ---- Blocked cross-validation ------------------------------------------
## Data summary by location: collapse prevalence and predictor gradients
## Table S5. urchin density and MHW maximum intensity (min / median / max), and collapse prevalence

data$collapse_num <- as.integer(data$collapse) - 1   # class_0 → 0, class_1 → 1
table(data$collapse, data$collapse_num)              

# Collapse prevalence per location (proportion of observations with collapse = 1)
tapply(data$collapse_num, data$location, mean)   


df_gradient <- aggregate(cbind(number, summer_inten) ~ location, data = data,
                         FUN = function(x) round(c(min = min(x), median = median(x),
                                                   max = max(x)), 2))
do.call(data.frame, df_gradient)  

table(data$location) ## to produce Table S5

## Blocked cross-validation (CV)
# Folds are defined by grouping variable rather than at random, so that all
# observations from the same site (or location) fall in the same fold.
# Predictions use re.form = NA: held-out sites have no estimated random
# effect, so prediction is at the population level.

blocked_cv <- function(form, dat, group_var) {
  folds <- split(seq_len(nrow(dat)), dat[[group_var]])
  p <- rep(NA_real_, nrow(dat))
  for (f in folds) {
    m <- glmmTMB(form, family = binomial, data = dat[-f, ])
    p[f] <- predict(m, dat[f, ], type = "response", re.form = NA)
  }
  p
}

# Model being validated: urchin density x MHW intensity interaction
# (both predictors centred), with a site random intercept

form <- collapse ~ number_c * inten_c + (1 | site_name)

p_site <- blocked_cv(form, data, "site_name")   # site-blocked CV
p_loc  <- blocked_cv(form, data, "location")    # leave-one-location-out


## ---- Blocked cross-validation Performance metrics -----------------------------------------------
# prevalence  = proportion of observations that are actually collapsed
#               (the baseline any model must beat)
# sensitivity = of the reefs that did collapse, the fraction correctly
#               predicted as collapsed
# specificity = of the reefs that did not collapse, the fraction correctly
#               predicted as not collapsed
# AUC         = area under the ROC curve; ranking ability, unaffected by
#               the base rate
# calibration = slope of observed outcome on predicted probability;
#               1 means predicted probabilities match observed frequencies
# Brier       = mean squared difference between predicted probability and
#               observed outcome; lower is better

perf <- function(obs, p, thresh = 0.5) {
  pred <- as.integer(p > thresh)
  c(prevalence  = mean(obs),
    sensitivity = mean(pred[obs == 1] == 1),
    specificity = mean(pred[obs == 0] == 0),
    AUC         = as.numeric(pROC::auc(obs, p, quiet = TRUE)),
    calibration = coef(glm(obs ~ qlogis(pmin(pmax(p, .001), .999)),
                           family = binomial))[2],
    Brier       = mean((p - obs)^2))
}

round(perf(data$collapse_num, p_site), 3)     # overall, site-blocked
round(perf(data$collapse_num, p_loc),  3)     # overall, LOLO

# per held-out location
do.call(rbind, by(data.frame(obs = data$collapse_num, p = p_loc),
                  data$location, function(d) perf(d$obs, d$p)))


## ---- Interaction contrast on the probability scale  --------------

# S = p(u1,i1) - p(u1,i0) - p(u0,i1) + p(u0,i0)
# Positive S means the joint effect exceeds the sum of the separate
# effects (super-additivity) on the probability scale.
# Predictions are conditional on the random effect being zero, i.e. for a
# typical site. State this in the Methods.

## Quartile-based estimate (u and i at their 25th vs 75th percentiles)

# 1. Choose comparison points inside the observed range
u0 <- quantile(data$number,       0.25)
u1 <- quantile(data$number,       0.75)
i0 <- quantile(data$summer_inten, 0.25)
i1 <- quantile(data$summer_inten, 0.75)

# 2. Build the four design rows (on the centred scale)
grid <- data.frame(
  number_c = c(u1, u1, u0, u0) - means_jb["number"],
  inten_c  = c(i1, i0, i1, i0) - means_jb["inten"]
)
X <- model.matrix(~ number_c * inten_c, data = grid)

# 3. Point estimate
b  <- glmmTMB::fixef(m9)$cond
stopifnot(ncol(X) == length(b))        # fails early if columns mismatch
colnames(X); names(b)                   # check the order matches

p  <- plogis(X %*% b)
S  <- p[1] - p[2] - p[3] + p[4]

# 4. Monte Carlo bootstrap for the confidence interval
## ---- Monte Carlo draws of the fixed-effect coefficients (used for Table S6, S9, and Fig. 4) ----
set.seed(2026)
n_draw <- 2000

# B: n_draw x (number of coefficients) matrix; each row is one coefficient draw
B <- MASS::mvrnorm(n_draw, b, vcov(m9)$cond)

Sb <- apply(B, 1, function(bb) {
  pp <- plogis(X %*% bb)
  pp[1] - pp[2] - pp[3] + pp[4]
})

# Point estimate and 95% CI 
c(S = S * 100,
  lower = quantile(Sb, 0.025) * 100,
  upper = quantile(Sb, 0.975) * 100,
  p_gt_0 = mean(Sb > 0))            # proportion of draws with S > 0


## ---- S across a grid of comparison points ------------------------------
# The S contrast depends on where it is evaluated, so report it across the
# observed range rather than at a single pair of points.

u_pairs <- list(c(0, 1), c(0, 2), c(1, 2), c(1, 3), c(2, 3))
i_pairs <- list(c(1, 2), c(1.5, 2.5), c(2, 3))

res <- expand.grid(up = seq_along(u_pairs), ip = seq_along(i_pairs))

out <- t(apply(res, 1, function(r) {
  
  u <- u_pairs[[r["up"]]]
  i <- i_pairs[[r["ip"]]]
  g <- data.frame(number_c = c(u[2],u[2],u[1],u[1]) - means_jb["number"],   ## check the mean
                  inten_c  = c(i[2],i[1],i[2],i[1]) - means_jb["inten"])
  
  X <- model.matrix(~ number_c * inten_c, data = g)            ## Check model
  
  p <- plogis(X %*% b)
  
  # S for each Monte Carlo draw (reuses B from above)
  Sb <- apply(B, 1, function(bb){ pp <- plogis(X %*% bb)
  pp[1]-pp[2]-pp[3]+pp[4] })
  c(u_low=u[1], u_high=u[2], i_low=i[1], i_high=i[2],
    S = (p[1]-p[2]-p[3]+p[4])*100,
    lo = quantile(Sb,.025)*100, 
    hi = quantile(Sb,.975)*100)
}))
round(out, 1)

##############################################################################################
## ---------------- projection preparation ------------------------------------
# ---- Predicted collapse probability surface (Jervis Bay model m9) --------
# Fixed-effect predictions from m9 across urchin density x MHW intensity.
# Predictors are centred on Jervis Bay means (means_jb), matching data_jb.

grid <- expand.grid(
  number       = seq(0, max(data_jb$number), length.out = 300),
  summer_inten = seq(0, max(data_jb$summer_inten), length.out = 300)
)
grid$number_c  <- grid$number - means["number"]
grid$inten_c   <- grid$summer_inten - means["inten"]

# Grouping variable required by glmmTMB; ignored because re.form = NA
grid$site_name <- NA

grid$p <- predict(m9, newdata = grid, type = "response", re.form = NA)

# ---- Plot -------------------------------------

p_breaks <- c(0.3, 0.5, 0.7, 0.9)
my_labs  <- c("0.3" = "30%", "0.5" = "50%", "0.7" = "70%", "0.9" = "90%")


p_contour <- ggplot(grid, aes(number, summer_inten, z = p)) + 
  geom_raster(aes(fill = p), interpolate = TRUE) +
  
  geom_contour(breaks = p_breaks, colour = "black", linewidth = 0.3) +
  metR::geom_text_contour(breaks = p_breaks, aes(label = my_labs[as.character(after_stat(level))]),
                          stroke = 0.2, skip = 0, size = 4.5,
                          label.placer = metR::label_placer_fraction(frac = 0.5)) +
  
  scico::scale_fill_scico(palette = "navia", begin = 0.12, 
                          limits = c(0, 1),
                          labels = scales::percent,
                          name = "Probability of\nkelp collapse\n") +
  
  scale_x_continuous(expand = c(0, 0), limits = c(0, 8)) +  ## CHANGE LIMIT if needed !!
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression("Sea urchin density (individuals m"^-2*")"),
       y = expression("Summer maximum intensity (" * degree * "C)")) +
  theme_format +
  theme(legend.position = "right", legend.key.height = unit(1.5, "cm"))

p_contour


###############################################################################
## ---------------- uncertainty propagation of GLMM parameters ----------------

path_projection <- '... projection input/output path ...'

### ----------------- CMIP6 2090-2100 input/output path - output from `cmip6_glmm_input.ipynb`
### ----------------- other input/output path - output from `location_glmm_projection_bar.ipynb`

input_dir  <- paste0(path_projection, "/glmm_input/")
output_dir <- paste0(path_projection, "/glmm_output/")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 2. Monte Carlo draws of the fixed-effect coefficients -------------
# Draw 2000 plausible coefficient sets from their joint sampling
# distribution (estimates + covariance), so that coefficient uncertainty
# can be carried through to the projected change in collapse probability.
library(MASS)
set.seed(2026)  ## same seed as the S contrast grid calculation above
n_draw <- 2000
b <- fixef(m9)$cond                       # coefficient estimates
B <- mvrnorm(n_draw, b, vcov(m9)$cond)    # rows = draws, cols = coefficients

## ---- 3. Helper functions ----------------------------------------------
# Rename and centre with the SAME means used to fit m9
prep_input <- function(df) {
  df <- df %>% rename(summer_inten = max_inten)
  df$number_c  <- df$number       - means_jb["number"]
  df$inten_c   <- df$summer_inten - means_jb["inten"]
  df$site_name <- NA
  df
}

# Collapse probability for every site under every coefficient draw
# (population level: site random effect = 0)
pred_draws <- function(df) {
  X <- model.matrix(~ number_c * inten_c, df)
  stopifnot(nrow(X) == nrow(df),                  # no rows dropped for NA
            identical(colnames(X), names(b)))     # columns match coefficients
  plogis(X %*% t(B))                              # rows = sites, cols = draws
}

# Average draws within each location (keeps uncertainty correct for
# location means; do NOT average the site-level intervals instead)
loc_mean <- function(M, loc) {
  rowsum(M, loc) / as.vector(rowsum(rep(1, nrow(M)), loc))
}

## ---- 4. Current conditions --------------------------------------------
df_input_current <- prep_input(read.csv(paste0(input_dir, 'glm_input_current.csv')))
cv_pred_current  <- predict(m9, newdata = df_input_current,
                            type = "response", re.form = NA)
P_now <- pred_draws(df_input_current)

## ---- 5. Loop through scenario files -----------------------------------
input_files <- list.files(input_dir, pattern = "glm_input_.*\\.csv", full.names = TRUE)

for (file_path in input_files) {
  
  df_input <- prep_input(read.csv(file_path))
  
  # sites must be in the same order as the current file
  stopifnot(all.equal(df_input$latitude,  df_input_current$latitude),
            all.equal(df_input$longitude, df_input_current$longitude))
  
  # point estimate (as before)
  cv_pred   <- predict(m9, df_input, type = "response", re.form = NA)
  df_output <- cbind(df_input, cv_pred)
  df_output$class_1_dev <- cv_pred - cv_pred_current
  
  # Monte Carlo: change in probability under each draw
  D <- pred_draws(df_input) - P_now               # subtract within each draw
  
  # site-level 95% interval
  df_output$dev_lo <- apply(D, 1, quantile, 0.025)
  df_output$dev_hi <- apply(D, 1, quantile, 0.975)
  
  # location-level summary (use this for Figure 4 bars and error bars)
  D_loc   <- loc_mean(D, df_input$location)
  dev_pt  <- tapply(df_output$class_1_dev, df_input$location, mean)
  loc_sum <- data.frame(location = rownames(D_loc),
                        dev      = as.vector(dev_pt[rownames(D_loc)]),
                        dev_lo   = apply(D_loc, 1, quantile, 0.025),
                        dev_hi   = apply(D_loc, 1, quantile, 0.975))
  
  # save
  file_name <- basename(file_path)
  out_site  <- paste0(output_dir, gsub("glm_input", "glmmjb_intera_ni_dev", file_name))
  out_loc   <- paste0(output_dir, gsub("glm_input", "glmmjb_intera_ni_dev_location", file_name))
  write.csv(df_output, out_site, row.names = FALSE)
  write.csv(loc_sum,   out_loc,  row.names = FALSE)
  
  cat("Processed:", file_path, "\n")
}








