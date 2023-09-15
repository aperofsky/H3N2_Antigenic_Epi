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

load("data/antigenic_epi_north_amer_build_for_lasso_replicates.Rdata")

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

sum_df <-
  epi_red %>%
  tidyr::replace_na(list(
    H1_cum_intensity = 0, H3_cum_intensity = 0,
    H3_cum_intensity = 0, IVB_cum_intensity = 0
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
  mutate_at(vars(H1_cum_intensity:weighted_VE), ~ scale_this(.x))

combined <- left_join(sumdf2, var_df_sc, by = "season") %>% dplyr::select(-contains(c("low", "hi", "sd")))

cum_red_df <- combined %>% dplyr::select(season, H3.epi.size.mean, all_of(vars))
form <- as.formula(paste("H3.ep.size.mean ~ ", paste(vars, collapse = "+")))
form

FULL.MODEL <- lm(
  H3.epi.size.mean ~ (H1_cum_intensity + HA_wolf_lag2 + prior_dom_type_national +
    ha_lbi_shannon + HA_titer_tree_lag2 + NA_bhatt_ep_lag1 +
    H3_epi_size_prior + vac_combined + weighted_VE_prior_season +
    weighted_VE),
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
performance(best_model)
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
  mutate_at(vars(H1_cum_intensity:ha_lbi_shannon), ~ scale_this(.x))

combined <- left_join(sumdf2, var_df_sc, by = "season") %>% dplyr::select(-contains(c("low", "hi", "sd")))

max_red_df <- combined %>% dplyr::select(season, H3.peak.mean, all_of(vars))
form <- as.formula(paste("H3.peak.mean ~ ", paste(vars, collapse = "+")))
form

FULL.MODEL <- lm(H3.peak.mean ~ (H1_cum_intensity + HA_wolf_lag2 + NA_bhatt_ep_lag1 +
  HA_titer_tree_lag2 + usa_bhatt_ep_mean + weighted_VE_prior_season +
  prior_dom_type_national + H3_epi_size_prior + vac_combined +
  ha_lbi_shannon), data = max_red_df)

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
performance(best_model)
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

var_df_sc <- var_df %>%
  mutate_at(vars(na_lbi_shannon:vac_cov_combined), ~ scale_this(.x))

combined <- left_join(sumdf2, var_df_sc, by = "season") %>% dplyr::select(-contains(c("low", "hi", "sd")))

r0_red_df <- combined %>%
  dplyr::select(season, H3.R0.mean, all_of(vars)) %>%
  filter(season != "2000-2001") # no Rt estimates in 2000-2001

form <- as.formula(paste("H3.R0.mean ~ ", paste(vars, collapse = "+")))

FULL.MODEL <- lm(
  H3.R0.mean ~ (na_lbi_shannon + H1_cum_intensity + usa_bhatt_ep_mean +
    ha_lbi_shannon_lag1 + vac_combined + HA_wolf_lag2 + vac_cov_combined +
    ha_lbi_shannon + vac_cov_combined_prior + HA_titer_tree_lag2),
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
performance(best_model)
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

var_df_sc <- var_df %>%
  mutate_at(vars(na_lbi_shannon:HA_wolf_lag2), ~ scale_this(.x))

combined <- left_join(sumdf2, var_df_sc, by = "season") %>% dplyr::select(-contains(c("low", "hi", "sd")))

shannon_red_df <- combined %>%
  dplyr::select(season, H3.shannon.mean, all_of(vars)) %>%
  filter(season != "2000-2001") # no shannon entropy for H3 during 2000-2001 season

form <- as.formula(paste("H3.shannon.mean ~ ", paste(vars, collapse = "+")))
form

FULL.MODEL <- lm(
  H3.shannon.mean ~ (ha_lbi_shannon + na_lbi_shannon + vac_cov_combined +
    vac_cov_combined_prior + usa_bhatt_ep_mean + HA_titer_tree_lag2 +
    na_lbi_shannon_lag1 + usa_wolf_ep_mean + ha_lbi_shannon_lag1 +
    HA_wolf_lag2),
  data = shannon_red_df
)

options(na.action = "na.fail") # Required for dredge to run
fm1 <- dredge(FULL.MODEL, evaluate = T, rank = "BIC", extra = "adjR^2", m.lim = c(1, 3))
options(na.action = "na.omit") # set back to default
subset(fm1, delta < 4)
best_model <- get.models(fm1, subset = 1)[[1]]
summary(best_model)
performance(best_model)
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
  mutate_at(vars(NA_bhatt_ep_lag1:na_lbi_shannon), ~ scale_this(.x))

combined <- left_join(sumdf2, var_df_sc, by = "season") %>% dplyr::select(-contains(c("low", "hi", "sd")))

dom_red_df <- combined %>% dplyr::select(season, H3.dom.mean, all_of(vars))
form <- as.formula(paste("H3.dom.mean ~ ", paste(vars, collapse = "+")))
form
FULL.MODEL <- lm(
  H3.dom.mean ~ (NA_bhatt_ep_lag1 + HA_wolf_lag2 + prior_dom_type_national +
    H3_epi_size_prior + HA_titer_tree_lag2 + vac_combined + ha_lbi_shannon_lag1 +
    usa_bhatt_ep_mean + vac_combined_prior + na_lbi_shannon),
  data = dom_red_df
)

options(na.action = "na.fail") # Required for dredge to run
fm1 <- dredge(FULL.MODEL, evaluate = T, rank = "BIC", extra = "adjR^2", m.lim = c(1, 3))
options(na.action = "na.omit") # set back to default
subset(fm1, delta < 4)
best_model <- get.models(fm1, subset = 1)[[1]]
summary(best_model)
performance(best_model)
compare_performance(get.models(fm1, subset = 1:5), metrics = c("AICc", "BIC", "R2", "R2_adj"), rank = T)
