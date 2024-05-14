## load packages
list.of.packages <- c("dplyr", "biwavelet", "scales", "tidyr", "lubridate")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

## custom wavelet functions
source("4_Wavelet_Analysis/Wavelets/WaveletPkg.R")

######################################################################################
## Wavelet Analysis
######################################################################################
# load data
load("data/hhs_division_level_ILI_and_virology_interp_smoothed_sampling_effort.RData")
load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")

df <- regionflu_ili_vir_adj %>% filter(wk_date <= as.Date("2019-07-01"))

df <- df %>%
  group_by(region) %>%
  complete(wk_date = seq.Date(range(df$wk_date)[1], range(df$wk_date)[2], by = "week")) %>%
  ungroup()

######################################################################################
#### H3 vs H1 timing
######################################################################################

h3_h1_phase_diff_func <- function(reg = "Region 2") {
  t1 <- df %>%
    filter(region == reg) %>%
    dplyr::select(wk_date, ili_h3_st)
  t1 <- t1 %>% mutate(ili_h3_st = replace_na(ili_h3_st, 0))
  t1$ili_h3_st <- sqrt(t1$ili_h3_st)
  t1$ili_h3_st <- as.vector(scale(t1$ili_h3_st))
  t1$wk_date <- lubridate::decimal_date(t1$wk_date)

  t2 <- df %>%
    filter(region == reg) %>%
    dplyr::select(wk_date, ili_h1_st)
  t2 <- t2 %>% mutate(ili_h1_st = replace_na(ili_h1_st, 0))
  t2$ili_h1_st <- sqrt(t2$ili_h1_st)
  t2$ili_h1_st <- as.vector(scale(t2$ili_h1_st))
  t2$wk_date <- lubridate::decimal_date(t2$wk_date)

  #### Phase analysis
  datee <- t1$wk_date
  length(datee)
  nrow(t1)

  h3poswave <- as.data.frame(cbind("x" = t1$ili_h3_st, "t" = datee))
  h3.cwt <- complete.cwt(h3poswave, dj = 1 / 10, pad = 1, mc.sim = 1000)
  # plot.cwt(h3.cwt)
  h3.phase <- phase(h3.cwt, s = c(0.8, 1.2))

  h1poswave <- as.data.frame(cbind("x" = t2$ili_h1_st, "t" = datee))
  h1.cwt <- complete.cwt(h1poswave, dj = 1 / 10, pad = 1, mc.sim = 1000)
  # plot.cwt(h1.cwt)
  h1.phase <- phase(h1.cwt, s = c(0.8, 1.2))

  h3h1_coh <- complete.coh(h1poswave, h3poswave, dj = 1 / 10, pad = 1)
  # plot.coh(h3h1_coh)###this throws an error because larger periods have NaN values

  phase_df <- data.frame(date = as.Date(date_decimal(datee)), H3 = h3.phase, H1 = h1.phase)

  # remove dates outside of cone of influence
  phase_df_long <- phase_df %>%
    pivot_longer(cols = !date, names_to = "subtype", values_to = "phase") %>%
    filter(date >= as.Date("1998-01-01") & date <= as.Date("2018-01-01"))

  h3_h1_phase_diff <- cohPhaseDiff(series1 = h3poswave, series2 = h1poswave, dj = 1 / 10, s = c(0.8, 1.2))
  phase_diff_df <- data.frame(date = as.Date(date_decimal(datee)), diff = h3_h1_phase_diff)

  # # phase = 2 * pi
  # # dt = time step
  # # p = period (1 year)
  # # phase / (2 * pi) * p * 1/dt
  phase_diff_in_weeks <- phaseToTime(h3_h1_phase_diff, p = 1, dt = 0.01917808) # dt = 1/52.1775
  phase_diff_df_weeks <- data.frame(date = as.Date(date_decimal(datee)), diff_weeks = phase_diff_in_weeks)

  phase_diff_df_weeks <-
    phase_diff_df_weeks %>%
    mutate(
      month.day = format(as.Date(date), "%m-%d"),
      new.year = as.numeric(format(as.Date(date), "%Y"))
    ) %>%
    mutate(season = ifelse(month.day < "07-01",
      sprintf("%d-%d", new.year - 1, new.year),
      sprintf("%d-%d", new.year, new.year + 1)
    ))
  phase_diff_df_weeks$region <- reg
  seasons <- epi_red %>% distinct(season, dom_type)
  phase_diff_df_weeks <- left_join(phase_diff_df_weeks, seasons, by = "season")

  phase_diff_df_weeks <-
    phase_diff_df_weeks %>%
    mutate(dom_type = ifelse(is.na(dom_type), "H1", dom_type))

  phase_diff_df_weeks$new.year <- paste0(phase_diff_df_weeks$new.year, "-01-01")
  phase_diff_df_weeks <- left_join(phase_diff_df_weeks, phase_diff_df, by = "date")
  title <- unique(gsub(" ", "_", phase_diff_df_weeks$region))
  write.csv(phase_diff_df_weeks, file = paste0("data/", title, "_phase_diff_h3_vs_h1.csv"), row.names = F)
}
h3_h1_phase_diff_func("Region 1") # these take ~30 minutes total to run on a 2021 macbook pro
h3_h1_phase_diff_func("Region 2")
h3_h1_phase_diff_func("Region 3")
h3_h1_phase_diff_func("Region 4")
h3_h1_phase_diff_func("Region 5")
h3_h1_phase_diff_func("Region 6")
h3_h1_phase_diff_func("Region 7")
h3_h1_phase_diff_func("Region 8")
h3_h1_phase_diff_func("Region 9")

######################################################################################
#### H3 vs B timing
######################################################################################
h3_ivb_phase_diff_func <- function(reg = "Region 2") {
  t1 <- df %>%
    filter(region == reg) %>%
    dplyr::select(wk_date, ili_h3_st)
  t1 <- t1 %>% mutate(ili_h3_st = replace_na(ili_h3_st, 0))
  t1$ili_h3_st <- sqrt(t1$ili_h3_st)
  t1$ili_h3_st <- as.vector(scale(t1$ili_h3_st))
  t1$wk_date <- lubridate::decimal_date(t1$wk_date)

  t2 <- df %>%
    filter(region == reg) %>%
    dplyr::select(wk_date, ili_ivb_st)
  t2 <- t2 %>% mutate(ili_ivb_st = replace_na(ili_ivb_st, 0))
  t2$ili_ivb_st <- sqrt(t2$ili_ivb_st)
  t2$ili_ivb_st <- as.vector(scale(t2$ili_ivb_st))
  t2$wk_date <- lubridate::decimal_date(t2$wk_date)

  #### Phase analysis
  datee <- t1$wk_date
  length(datee)
  nrow(t1)

  h3poswave <- as.data.frame(cbind("x" = t1$ili_h3_st, "t" = datee))
  h3.cwt <- complete.cwt(h3poswave, dj = 1 / 10, pad = 1, mc.sim = 1000)
  # plot.cwt(h3.cwt)
  h3.phase <- phase(h3.cwt, s = c(0.8, 1.2))

  ivbposwave <- as.data.frame(cbind("x" = t2$ili_ivb_st, "t" = datee))
  ivb.cwt <- complete.cwt(ivbposwave, dj = 1 / 10, pad = 1, mc.sim = 1000)
  # plot.cwt(ivb.cwt)
  ivb.phase <- phase(ivb.cwt, s = c(0.8, 1.2))

  h3ivb_coh <- complete.coh(ivbposwave, h3poswave, dj = 1 / 10, pad = 1)
  # plot.coh(h3ivb_coh)###this throws an error because larger periods have NaN values

  phase_df <- data.frame(date = as.Date(date_decimal(datee)), H3 = h3.phase, ivb = ivb.phase)
  phase_df_long <- phase_df %>%
    pivot_longer(cols = !date, names_to = "subtype", values_to = "phase") %>%
    filter(date >= as.Date("1998-01-01") & date <= as.Date("2018-01-01")) # remove dates outside of cone of influence

  h3_ivb_phase_diff <- cohPhaseDiff(series1 = h3poswave, series2 = ivbposwave, dj = 1 / 10, s = c(0.8, 1.2))
  phase_diff_df <- data.frame(date = as.Date(date_decimal(datee)), diff = h3_ivb_phase_diff)

  # # phase = 2 * pi
  # # dt = time step
  # # p = period (1 year)
  # # phase / (2 * pi) * p * 1/dt
  phase_diff_in_weeks <- phaseToTime(h3_ivb_phase_diff, p = 1, dt = 0.01917808) # dt = 1/52.1775
  phase_diff_df_weeks <- data.frame(date = as.Date(date_decimal(datee)), diff_weeks = phase_diff_in_weeks)

  phase_diff_df_weeks <-
    phase_diff_df_weeks %>%
    mutate(
      month.day = format(as.Date(date), "%m-%d"),
      new.year = as.numeric(format(as.Date(date), "%Y"))
    ) %>%
    mutate(season = ifelse(month.day < "07-01",
      sprintf("%d-%d", new.year - 1, new.year),
      sprintf("%d-%d", new.year, new.year + 1)
    ))
  phase_diff_df_weeks$region <- reg
  seasons <- epi_red %>% distinct(season, dom_type)
  phase_diff_df_weeks <- left_join(phase_diff_df_weeks, seasons, by = "season")

  phase_diff_df_weeks <-
    phase_diff_df_weeks %>%
    mutate(dom_type = ifelse(is.na(dom_type), "H1", dom_type))

  phase_diff_df_weeks$new.year <- paste0(phase_diff_df_weeks$new.year, "-01-01")
  phase_diff_df_weeks <- left_join(phase_diff_df_weeks, phase_diff_df, by = "date")
  title <- unique(gsub(" ", "_", phase_diff_df_weeks$region))
  write.csv(phase_diff_df_weeks, file = paste0("data/", title, "_phase_diff_h3_vs_ivb.csv"), row.names = F)
}
h3_ivb_phase_diff_func("Region 1")
h3_ivb_phase_diff_func("Region 2")
h3_ivb_phase_diff_func("Region 3")
h3_ivb_phase_diff_func("Region 4")
h3_ivb_phase_diff_func("Region 5")
h3_ivb_phase_diff_func("Region 6")
h3_ivb_phase_diff_func("Region 7")
h3_ivb_phase_diff_func("Region 8")
h3_ivb_phase_diff_func("Region 9")
