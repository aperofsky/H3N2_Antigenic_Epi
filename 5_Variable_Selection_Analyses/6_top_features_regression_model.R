## load packages
list.of.packages <- c("dplyr", "ggplot2", "cowplot", "tidyr", "vegan", "purrr", "gmodels", "performance", "MuMIn")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

######################################################################################
## Linear regression and model selection using top 10 ranked variables from RF models
######################################################################################

######################################################################################
## compile dataset
######################################################################################
# load data
load("data/subtype_distribution_by_region_season.RData") # subtype_dist
subtype_dist <- subtype_dist %>%
  ungroup() %>%
  mutate(a_total = h3_total + h1_total) %>%
  mutate(
    h3_vs_h1 = h3_total / a_total,
    h3_dom = h3_total / (a_total + b_total)
  ) %>% # prop of h3 out of IAV%>%
  rename(season = season_description)
names(subtype_dist)

load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")

combined <- left_join(subtype_dist,
  epi_red %>% dplyr::select(
    -h3_vs_h1, -iva_vs_ivb, -h3_vs_total_flu,
    -h1_vs_total_flu, -ivb_vs_iva, -total_a
  ),
  by = c("season", "region")
) %>%
  distinct() %>%
  filter(!(season %in% c("2009-2010", "2019-2020", "2020-2021", "2021-2022")))

epi_red <- combined %>%
  mutate(
    vac_combined = (adult_18_49_vac_cov / 100) * (adult_65_vac_cov / 100) * weighted_VE,
    vac_cov_combined = (adult_18_49_vac_cov / 100) * (adult_65_vac_cov / 100)
  )

vac_prior <- epi_red %>%
  distinct(season, adult_18_49_vac_cov, adult_65_vac_cov, vac_combined, vac_cov_combined) %>%
  arrange(season) %>%
  mutate(
    adult_18_49_vac_cov_prior = lag(adult_18_49_vac_cov, n = 1),
    adult_65_vac_cov_prior = lag(adult_65_vac_cov, n = 1),
    vac_combined_prior = lag(vac_combined, n = 1), # vac cov x VE
    vac_cov_combined_prior = lag(vac_cov_combined, n = 1)
  ) %>%
  ungroup()

epi_red <- left_join(epi_red, vac_prior %>% dplyr::select(season, contains("prior")), by = "season")

epi_red <- epi_red %>%
  mutate(prior_dom_type_national = case_when(
    prior_dom_type_national == "H3" ~ 1,
    prior_dom_type_national == "H1" ~ 0,
    prior_dom_type_national == "co-circ" ~ 0.5
  ))

epi_red %>%
  filter(season == "2000-2001") %>%
  dplyr::select(region, H3_shannon_entropy_res)

names(epi_red)
sum_df <-
  epi_red %>%
  tidyr::replace_na(list(
    H1_cum_intensity = 0, H3_cum_intensity = 0,
    H3_cum_intensity = 0, IVB_cum_intensity = 0,
    H3_epi_size_prior = 0, IVB_epi_size_prior = 0,
    H1_epi_size_prior = 0
  )) %>%
  group_by(
    season, dom_type
  ) %>%
  dplyr::summarise(
    H3.shannon.mean = ci(H3_shannon_entropy_res, na.rm = T)[1],
    H3.shannon.lowCI = ci(H3_shannon_entropy_res, na.rm = T)[2],
    H3.shannon.hiCI = ci(H3_shannon_entropy_res, na.rm = T)[3],
    H3.shannon.sd = ci(H3_shannon_entropy_res, na.rm = T)[4],
    H3.peak.mean = ci(H3_max_intensity, na.rm = T)[1],
    H3.peak.lowCI = ci(H3_max_intensity, na.rm = T)[2],
    H3.peak.hiCI = ci(H3_max_intensity, na.rm = T)[3],
    H3.peak.sd = ci(H3_max_intensity, na.rm = T)[4],
    H3.epi.size.mean = ci(H3_cum_intensity, na.rm = T)[1],
    H3.epi.size.lowCI = ci(H3_cum_intensity, na.rm = T)[2],
    H3.epi.size.hiCI = ci(H3_cum_intensity, na.rm = T)[3],
    H3.epi.size.sd = ci(H3_cum_intensity, na.rm = T)[4],
    H3.R0.mean = ci(H3_max_Rt, na.rm = T)[1],
    H3.R0.lowCI = ci(H3_max_Rt, na.rm = T)[2],
    H3.R0.hiCI = ci(H3_max_Rt, na.rm = T)[3],
    H3.R0.sd = ci(H3_max_Rt, na.rm = T)[4],
    H3.dom.mean = ci(h3_dom, na.rm = T)[1],
    H3.dom.lowCI = ci(h3_dom, na.rm = T)[2],
    H3.dom.hiCI = ci(h3_dom, na.rm = T)[3],
    H3.dom.sd = ci(h3_dom, na.rm = T)[4]
  ) %>%
  ungroup()

sum_df %>%
  print(n = 30)

sumdf2 <- sum_df %>% filter(!(season %in% c("1995-1996", "1996-1997", "1997-1998", "2009-2010")))

######################################################################################
## load variable importance results
######################################################################################

h3_peak <- read.csv("data/H3_peak_cforest_variable_importance.csv")
h3_peak$epi_measure <- "peak"

h3_size <- read.csv("data/H3_epi_size_cforest_variable_importance.csv")
h3_size$epi_measure <- "cum_size"

h3_r0 <- read.csv("data/H3_R0_cforest_variable_importance.csv")
h3_r0$epi_measure <- "R0"

h3_shannon <- read.csv("data/H3_shannon_entropy_cforest_variable_importance.csv")
h3_shannon$epi_measure <- "shannon_entropy"

h3_dom <- read.csv("data/h3_dom_cforest_variable_importance.csv")
h3_dom$epi_measure <- "dominance"

all_cforest <- bind_rows(h3_peak, h3_size, h3_dom, h3_r0, h3_shannon)

all_cforest$epi_measure <- as.factor(all_cforest$epi_measure)

######################################################################################
# Epidemic Size
######################################################################################
vars <- all_cforest %>%
  filter(epi_measure == "cum_size") %>%
  arrange(-Estimate) %>%
  slice_head(n = 10) %>%
  pull(var)

var_df <- epi_red %>%
  dplyr::select(season, region, all_of(vars)) %>%
  filter(!(season %in% c("1995-1996", "1996-1997", "1997-1998", "2009-2010"))) %>%
  group_by(season) %>%
  summarize(across(all_of(vars), ~ mean(.x, na.rm = T)))

names(var_df)

scale_this <- function(x) as.vector(scale(x))

var_df_sc <- var_df %>%
  mutate_at(vars(H1_cum_intensity:NA_std_lbi), ~ scale_this(.x))

combined <- left_join(sumdf2, var_df_sc, by = "season") %>% dplyr::select(-contains(c("low", "hi", "sd")))

cum_red_df <- combined %>% dplyr::select(season, H3.epi.size.mean, all_of(vars))
form <- as.formula(paste("H3.ep.size.mean ~ ", paste(vars, collapse = "+")))
form

FULL.MODEL <- lm(
  H3.epi.size.mean ~ (H1_cum_intensity + HA_wolf_lag2 + prior_dom_type_national + 
                        HA_std_lbi + NA_bhatt_ep_lag1 + HA_titer_tree_lag2 + vac_combined + 
                        H3_epi_size_prior + weighted_VE_prior_season + NA_std_lbi),
  data = cum_red_df
)

options(na.action = "na.fail") # Required for dredge to run
fm1 <- dredge(FULL.MODEL,
  evaluate = T,
  rank = "BIC",
  extra = "adjR^2", m.lim = c(1, 3)
)
options(na.action = "na.omit") # set back to default
subset(fm1, delta < 4)
best_model <- get.models(fm1, subset = 1)[[1]]
summary(best_model)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         28.333      2.469  11.474 3.93e-09 ***
#   H1_cum_intensity    -9.933      2.756  -3.604  0.00238 ** 
#   H3_epi_size_prior   -6.817      2.805  -2.430  0.02724 *  
#   HA_wolf_lag2        10.434      2.812   3.710  0.00190 ** 
performance(best_model)
# AIC     |    AICc |     BIC |    R2 | R2 (adj.) |  RMSE |  Sigma
# ----------------------------------------------------------------
# 158.367 | 162.653 | 163.346 | 0.740 |     0.691 | 9.877 | 11.043
compare_performance(get.models(fm1, subset = 1:5), metrics = c("AICc", "BIC", "R2", "R2_adj"), rank = T) %>% arrange("AICc")

######################################################################################
# Peak Incidence
######################################################################################
vars <- all_cforest %>%
  filter(epi_measure == "peak") %>%
  arrange(-Estimate) %>%
  slice_head(n = 10) %>%
  pull(var)
vars

var_df <- epi_red %>%
  dplyr::select(season, region, all_of(vars)) %>%
  filter(!(season %in% c("1995-1996", "1996-1997", "1997-1998", "2009-2010"))) %>%
  group_by(season) %>%
  summarize(across(all_of(vars), ~ mean(.x, na.rm = T)))

scale_this <- function(x) as.vector(scale(x))
names(var_df)

var_df_sc <- var_df %>%
  mutate_at(vars(H1_cum_intensity:vac_combined), ~ scale_this(.x))

combined <- left_join(sumdf2, var_df_sc, by = "season") %>% dplyr::select(-contains(c("low", "hi", "sd")))

max_red_df <- combined %>% dplyr::select(season, H3.peak.mean, all_of(vars))
form <- as.formula(paste("H3.peak.mean ~ ", paste(vars, collapse = "+")))
form

FULL.MODEL <- lm(H3.peak.mean ~ (H1_cum_intensity + HA_wolf_lag2 + NA_bhatt_ep_lag1 + 
                                   weighted_VE_prior_season + HA_titer_tree_lag2 + usa_bhatt_ep_mean + 
                                   HA_std_lbi + H3_epi_size_prior + prior_dom_type_national + 
                                   vac_combined), data = max_red_df)

options(na.action = "na.fail") # Required for dredge to run
fm1 <- dredge(FULL.MODEL,
  evaluate = T,
  rank = "BIC",
  extra = "adjR^2", m.lim = c(1, 3)
)
options(na.action = "na.omit") # set back to default
subset(fm1, delta < 4)
best_model <- get.models(fm1, subset = 1)[[1]]
summary(best_model)
# Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)               4.6152     0.5215   8.849 1.46e-07 ***
#   H1_cum_intensity         -1.7415     0.6078  -2.865  0.01122 *  
#   HA_wolf_lag2              1.9121     0.5841   3.273  0.00478 ** 
#   prior_dom_type_national  -1.1344     0.6051  -1.875  0.07923 . 
performance(best_model)
# AIC    |    AICc |     BIC |    R2 | R2 (adj.) |  RMSE | Sigma
# --------------------------------------------------------------
# 96.170 | 100.456 | 101.149 | 0.685 |     0.626 | 2.086 | 2.332
compare_performance(get.models(fm1, subset = 1:5), metrics = c("AICc", "BIC", "R2", "R2_adj"), rank = T) %>% arrange("BIC")

######################################################################################
# Effective Rt
######################################################################################
vars <- all_cforest %>%
  filter(epi_measure == "R0") %>%
  arrange(-Estimate) %>%
  slice_head(n = 10) %>%
  pull(var)

var_df <- epi_red %>%
  dplyr::select(season, region, all_of(vars)) %>%
  filter(!(season %in% c("1995-1996", "1996-1997", "1997-1998", "2009-2010"))) %>%
  group_by(season) %>%
  summarize(across(all_of(vars), ~ mean(.x, na.rm = T)))
scale_this <- function(x) as.vector(scale(x))

names(var_df)
var_df_sc <- var_df %>%
  mutate_at(vars(H1_cum_intensity:H3_epi_size_prior), ~ scale_this(.x))

combined <- left_join(sumdf2, var_df_sc, by = "season") %>% dplyr::select(-contains(c("low", "hi", "sd")))

r0_red_df <- combined %>%
  dplyr::select(season, H3.R0.mean, all_of(vars)) %>%
  filter(season != "2000-2001") # no Rt estimates in 2000-2001

form <- as.formula(paste("H3.R0.mean ~ ", paste(vars, collapse = "+")))
form

FULL.MODEL <- lm(
  H3.R0.mean ~ (H1_cum_intensity + vac_cov_combined_prior + usa_bhatt_ep_mean + 
                  NA_std_lbi + vac_combined + HA_wolf_lag2 + HA_titer_tree_lag2 + 
                  NA_bhatt_ep_lag1 + HA_std_lbi + H3_epi_size_prior),
  data = r0_red_df
)

options(na.action = "na.fail") # Required for dredge to run
fm1 <- dredge(FULL.MODEL,
  evaluate = T,
  rank = "BIC",
  extra = "adjR^2", m.lim = c(1, 3)
)
options(na.action = "na.omit") # set back to default
subset(fm1, delta < 4)
best_model <- get.models(fm1, subset = 1)[[1]]
summary(best_model)
# Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)         1.35304    0.02788  48.533   <2e-16 ***
#   H1_cum_intensity   -0.10660    0.02909  -3.664   0.0023 ** 
#   HA_titer_tree_lag2  0.08251    0.03076   2.682   0.0171 *  
#   usa_bhatt_ep_mean   0.06794    0.02968   2.289   0.0370 *  
performance(best_model)
# AIC     |    AICc |     BIC |    R2 | R2 (adj.) |  RMSE | Sigma
# ---------------------------------------------------------------
# -21.021 | -16.405 | -16.298 | 0.690 |     0.628 | 0.107 | 0.120
compare_performance(get.models(fm1, subset = 1:5), metrics = c("AICc", "BIC", "R2", "R2_adj"), rank = T) %>% arrange("AICc")

######################################################################################
# Epidemic intensity: Inverse Shannon Entropy
######################################################################################

vars <- all_cforest %>%
  filter(epi_measure == "shannon_entropy") %>%
  arrange(-Estimate) %>%
  slice_head(n = 10) %>%
  pull(var)

var_df <- epi_red %>%
  dplyr::select(season, region, all_of(vars)) %>%
  filter(!(season %in% c("1995-1996", "1996-1997", "1997-1998", "2009-2010"))) %>%
  group_by(season) %>%
  summarize(across(all_of(vars), ~ mean(.x, na.rm = T)))

scale_this <- function(x) as.vector(scale(x))
names(var_df)
var_df_sc <- var_df %>%
  mutate_at(vars(vac_cov_combined_prior:HA_wolf_nonepitope_lag2), ~ scale_this(.x))

combined <- left_join(sumdf2, var_df_sc, by = "season") %>% dplyr::select(-contains(c("low", "hi", "sd")))

shannon_red_df <- combined %>%
  dplyr::select(season, H3.shannon.mean, all_of(vars)) %>%
  filter(season != "2000-2001") # no shannon entropy for H3 during 2000-2001 season

form <- as.formula(paste("H3.shannon.mean ~ ", paste(vars, collapse = "+")))
form

FULL.MODEL <- lm(
  H3.shannon.mean ~ (vac_cov_combined_prior + usa_bhatt_ep_mean + 
                       NA_bhatt_ep_lag1 + HA_titer_tree_lag2 + NA_std_lbi + usa_wolf_ep_mean + 
                       NA_bhatt_nonepitope_lag1 + vac_combined + H1_cum_intensity + 
                       HA_wolf_nonepitope_lag2),
  data = shannon_red_df
)

options(na.action = "na.fail") # Required for dredge to run
fm1 <- dredge(FULL.MODEL, evaluate = T, rank = "BIC", extra = "adjR^2", m.lim = c(1, 3))
options(na.action = "na.omit") # set back to default
subset(fm1, delta < 4)
best_model <- get.models(fm1, subset = 1)[[1]]
summary(best_model)
#                           Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)             0.40883    0.01876  21.797 9.04e-13 ***
#   HA_titer_tree_lag2      0.05470    0.02122   2.577 0.021021 *  
#   usa_bhatt_ep_mean       0.09412    0.01977   4.760 0.000253 ***
#   vac_cov_combined_prior -0.06783    0.01990  -3.409 0.003884 ** 
performance(best_model)
# AIC     |    AICc |     BIC |    R2 | R2 (adj.) |  RMSE | Sigma
# ---------------------------------------------------------------
# -36.077 | -31.462 | -31.355 | 0.793 |     0.752 | 0.072 | 0.081
compare_performance(get.models(fm1, subset = 1:5), metrics = c("AICc", "BIC", "R2", "R2_adj"), rank = T)

######################################################################################
# A/H3N2 subytpe dominance
######################################################################################
vars <- all_cforest %>%
  filter(epi_measure == "dominance") %>%
  arrange(-Estimate) %>%
  slice_head(n = 10) %>%
  pull(var)

var_df <- epi_red %>%
  dplyr::select(season, region, all_of(vars)) %>%
  filter(!(season %in% c("1995-1996", "1996-1997", "1997-1998", "2009-2010"))) %>%
  group_by(season) %>%
  summarize(across(all_of(vars), ~ mean(.x, na.rm = T)))

scale_this <- function(x) as.vector(scale(x))
names(var_df)
var_df_sc <- var_df %>%
  mutate_at(vars(NA_bhatt_ep_lag1:H1_epi_size_prior), ~ scale_this(.x))

combined <- left_join(sumdf2, var_df_sc, by = "season") %>% dplyr::select(-contains(c("low", "hi", "sd")))

dom_red_df <- combined %>% dplyr::select(season, H3.dom.mean, all_of(vars))
form <- as.formula(paste("H3.dom.mean ~ ", paste(vars, collapse = "+")))
form
FULL.MODEL <- lm(
  H3.dom.mean ~ (NA_bhatt_ep_lag1 + prior_dom_type_national + HA_wolf_lag2 + 
                   H3_epi_size_prior + HA_koel_lag2 + HA_titer_tree_lag2 + vac_combined + 
                   HA_wolf_nonepitope_lag2 + HA_std_lbi + H1_epi_size_prior),
  data = dom_red_df
)

options(na.action = "na.fail") # Required for dredge to run
fm1 <- dredge(FULL.MODEL, evaluate = T, rank = "BIC", extra = "adjR^2", m.lim = c(1, 3))
options(na.action = "na.omit") # set back to default
subset(fm1, delta < 4)
best_model <- get.models(fm1, subset = 1)[[1]]
summary(best_model)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              0.53072    0.04983  10.650 1.14e-08 ***
# HA_wolf_lag2             0.11637    0.05741   2.027  0.05968 .  
# NA_bhatt_ep_lag1         0.10760    0.05601   1.921  0.07274 .  
# prior_dom_type_national -0.16261    0.05276  -3.082  0.00714 ** 
performance(best_model)
# AIC   |  AICc |   BIC |    R2 | R2 (adj.) |  RMSE | Sigma
# ---------------------------------------------------------
# 2.246 | 6.532 | 7.225 | 0.562 |     0.480 | 0.199 | 0.223
compare_performance(get.models(fm1, subset = 1:5), metrics = c("AICc", "BIC", "R2", "R2_adj"), rank = T)
