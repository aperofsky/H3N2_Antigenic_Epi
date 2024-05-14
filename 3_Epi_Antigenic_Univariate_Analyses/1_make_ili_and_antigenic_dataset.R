## load packages
list.of.packages <- c("dplyr", "ggplot2", "cowplot", "tidyr", "purrr", "readr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

###################################################################################################
## MASTER DATAFRAME FOR EPI AND ANTIGENIC MEASURES = North America Build
###################################################################################################
# regionflu_ili_vir_adj; weekly epi data
load("data/hhs_division_level_ILI_and_virology_interp_smoothed_sampling_effort.RData")
head(regionflu_ili_vir_adj)
names(regionflu_ili_vir_adj)

# season-level H3 epi data
load("data/region_level_flu_metrics_H3.RData")
names(region_flu_metrics_H3)[3:ncol(region_flu_metrics_H3)] <- paste0("H3_", names(region_flu_metrics_H3)[3:ncol(region_flu_metrics_H3)])
head(region_flu_metrics_H3)
sort(names(region_flu_metrics_H3))
region_flu_metrics_H3 %>% filter(region=="Region 1")

## Rt estimates for H3
load("data/epidemic_Rt_estimates.RData")
head(R0_eg_values)
names(R0_eg_values)[3:ncol(R0_eg_values)] <- paste0("H3_", names(R0_eg_values)[3:ncol(R0_eg_values)])
head(R0_eg_values)

## calculate subtype ratios for each combination
season_level_subtype_ratio <- regionflu_ili_vir_adj %>%
  group_by(region, season_description) %>%
  summarize(
    total_b = sum(prop_b * total_specimens, na.rm = T),
    total_a = sum(prop_a * total_specimens, na.rm = T),
    total_h3 = sum(prop_h3 * total_specimens, na.rm = T),
    total_h1 = sum(prop_h1 * total_specimens, na.rm = T),
    total_specimens = sum(total_specimens, na.rm = T)
  ) %>%
  group_by(region, season_description) %>%
  mutate(total_flu = sum(total_b + total_a)) %>%
  mutate(
    h3_vs_h1 = total_h3 / total_a,
    iva_vs_ivb = total_a / total_flu,
    h3_vs_total_flu = total_h3 / total_flu,
    h1_vs_total_flu = total_h1 / total_flu,
    ivb_vs_iva = total_b / total_flu
  ) %>%
  rename(season = season_description)

combined_epi <- list(region_flu_metrics_H3, R0_eg_values, season_level_subtype_ratio) %>%
  purrr::reduce(full_join, by = c("season", "region")) %>%
  filter(!(season %in% c("2019-2020", "2020-2021", "2021-2022")))
unique(combined_epi$season)
length(unique(combined_epi$season))

# season-level HA antigenic and genetic change from North America build
load("data/north_amer_build_season_h3n2_replicates_HA_direct_ag_distances.RData")
names(seas.ag.HA)

## average across bootstrap replicates of phylogenetic analysis
seas.ag.HA <- seas.ag.HA %>%
  group_by(season) %>%
  summarise(across(c(HA_koel_lag0:HA_std_lbi), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()
names(seas.ag.HA)

# season-level NA genetic change from North America build
load("data/north_amer_build_season_h3n2_replicates_NA_direct_ag_distances.RData")
names(seas.ag.NA)

## average across bootstrap replicates of phylogenetic analysis
seas.ag.NA <- seas.ag.NA %>%
  group_by(season) %>%
  summarise(across(c(NA_bhatt_ep_lag0:NA_std_lbi), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()
head(seas.ag.NA)

## LBI diversity metrics
lbi_div <- readr::read_rds("data/LBI_diversity_by_season.rds")
head(lbi_div)

lbi_div <- lbi_div %>%
  group_by(season) %>%
  summarize_at(c("ha_lbi_shannon", "na_lbi_shannon"), mean, na.rm = T)

# lbi_div_lag1 <- lbi_div %>%
#   mutate_at(c("ha_lbi_shannon", "na_lbi_shannon"), lag, n = 1)
# names(lbi_div_lag1)[2:3] <- paste0(names(lbi_div_lag1)[2:3], "_lag1")

# lbi_div_lag2 <- lbi_div %>%
#   mutate_at(c("ha_lbi_shannon", "na_lbi_shannon"), lag, n = 2)
# names(lbi_div_lag2)[2:3] <- paste0(names(lbi_div_lag2)[2:3], "_lag2")

# lbi_div <- full_join(lbi_div, lbi_div_lag1, by = "season")
# lbi_div <- full_join(lbi_div, lbi_div_lag2, by = "season")
names(lbi_div)

names(seas.ag.HA)[names(seas.ag.HA) == "year1"] <- "year.new"
names(seas.ag.HA)
keep.HA <- seas.ag.HA %>%
  dplyr::select(
    season, contains(c("lag", "lbi")),
    -contains(c("stem_lag", "titer_sub"))
  ) %>%
  names()
sort(keep.HA)
seas.ag.HA.red <- seas.ag.HA %>% dplyr::select(all_of(keep.HA))
head(seas.ag.HA.red)
seas.ag.HA.red <- as.data.frame(seas.ag.HA.red)

names(seas.ag.NA)
keep.NA <- names(seas.ag.NA)
seas.ag.NA.red <- seas.ag.NA %>% dplyr::select(all_of(keep.NA))
head(seas.ag.NA.red)
seas.ag.NA.red <- as.data.frame(seas.ag.NA.red)

##################################################################
## combining antigenic data and regional epidemic burden data
##################################################################
unique(combined_epi$region)

combined_df <- list(combined_epi, seas.ag.HA.red, seas.ag.NA.red, lbi_div) %>%
  purrr::reduce(full_join, by = c("season")) %>%
  filter(season != "2019-2020")
unique(combined_df$region)
combined_df %>% filter(is.na(region)) ## 1996-1997
combined_df <- combined_df %>% tidyr::separate(season, sep = "-", into = c("year.new", "year2"), remove = F)
combined_df %>%
  dplyr::select(season, year.new, region, H3_onset_week, H3_peak_week) %>%
  arrange(season)
class(combined_df$H3_onset_week)
class(combined_df$H3_peak_week)

combined_df2 <- combined_df %>%
  rowwise() %>%
  mutate(season_start = cdcfluview::mmwr_week_to_date(year = as.numeric(year.new), week = 40, day = NULL)) %>% # start each season from week 40
  mutate(
    onset_days_from_Oct1 = as.numeric(H3_onset_week - season_start),
    peak_days_from_Oct1 = as.numeric(H3_peak_week - season_start),
    onset_days_from_Oct1_bayes = as.numeric(H3_onset_week_bayes - season_start)
  )

combined_df2 %>% dplyr::select(season, region, onset_days_from_Oct1, onset_days_from_Oct1_bayes)

combined_df2$peak_diff <- combined_df2$H3_peak_week - combined_df2$H3_onset_week # days from onset to peak
combined_df2$peak_diff_bayes <- combined_df2$H3_peak_week - combined_df2$H3_onset_week_bayes # days from onset to peak

combined_df2 %>%
  filter(is.na(onset_days_from_Oct1)) %>%
  distinct(season, region) ## missing onsets due to lack of data
names(combined_df2)

# onset_vec <- as.numeric(combined_df2$onset_days_from_Oct1[!is.na(combined_df2$onset_days_from_Oct1)])
# peak_vec <- as.numeric(combined_df2$peak_days_from_Oct1[!is.na(combined_df2$peak_days_from_Oct1)])
# fitdistrplus::descdist(onset_vec, boot = 1000) # close to normal dist
# fitdistrplus::descdist(peak_vec, boot = 1000) # close to normal dist

# onset_vec_bayes <- as.numeric(combined_df2$onset_days_from_Oct1_bayes[!is.na(combined_df2$onset_days_from_Oct1_bayes)])
# onset_vec_bayes
# fitdistrplus::descdist(onset_vec_bayes, boot = 1000) # close to uniform distribution

onset_compare = combined_df2 %>% dplyr::select(season,region,onset_days_from_Oct1,onset_days_from_Oct1_bayes) %>% filter(season!="1996-1997")
onset_compare = onset_compare %>% tidyr::separate(season,into=c("year1","year2"),remove=F)
onset_compare$diff = onset_compare$onset_days_from_Oct1 - onset_compare$onset_days_from_Oct1_bayes
onset_compare %>% filter(is.na(region))

ggplot(onset_compare)+
  geom_point(aes(x=as.numeric(year2),y=diff,fill=region),pch=21)+
  scale_y_continuous(n.breaks=10)+
  theme_bw()

epidemic_timing_df <-
  combined_df2 %>%
  filter(season != "1996-1997") %>%
  group_by(season) %>%
  summarise(
    onset_timing_sd = sd(onset_days_from_Oct1, na.rm = T),
    peak_timing_sd = sd(peak_days_from_Oct1, na.rm = T),
    onset_timing_mad = mad(onset_days_from_Oct1, na.rm = T),
    peak_timing_mad = mad(peak_days_from_Oct1, na.rm = T),
    iqr_onset = IQR(onset_days_from_Oct1, na.rm = T),
    iqr_peak = IQR(peak_days_from_Oct1, na.rm = T),
    onset_median = median(onset_days_from_Oct1, na.rm = T),
    peak_median = median(peak_days_from_Oct1, na.rm = T),
    peak_diff_median = median(peak_diff, na.rm = T),
    ## bayesian estimates of onset week
    onset_timing_sd_bayes = sd(onset_days_from_Oct1, na.rm = T),
    onset_timing_mad_bayes = mad(onset_days_from_Oct1_bayes, na.rm = T),
    iqr_onset_bayes = IQR(onset_days_from_Oct1, na.rm = T),
    onset_median_bayes = median(onset_days_from_Oct1, na.rm = T),
    peak_diff_median_bayes = median(peak_diff, na.rm = T)
  )
head(epidemic_timing_df)
unique(epidemic_timing_df$onset_timing_mad)
unique(epidemic_timing_df$onset_timing_mad_bayes)
epidemic_timing_df %>% filter(is.na(onset_timing_mad) | is.na(onset_timing_mad_bayes)) # 2000-2001 and 2009-2010 (no H3N2 circulation)

epi_antigenic_data <- left_join(combined_df2, epidemic_timing_df, by = "season")
epi_antigenic_data[!complete.cases(epi_antigenic_data$region), ] # 1996-1997

cluster_transition <- c(
  "1997-1998", "2000-2001", "2003-2004", "2004-2005", "2006-2007", "2007-2008",
  "2009-2010", "2012-2013", "2013-2014", "2014-2015", "2016-2017"
)
major_cl_transition <- c("1997-1998", "2000-2001", "2003-2004", "2004-2005", "2012-2013")
vaccine_mismatch <- c("1997-1998", "2003-2004", "2004-2005", "2007-2008", "2014-2015", "2017-2018")
h3_years_all <- c(
  "1997-1998", "1998-1999", "1999-2000", "2001-2002", "2003-2004", "2004-2005", "2005-2006", "2007-2008", "2010-2011",
  "2011-2012", "2012-2013", "2014-2015", "2016-2017", "2017-2018", "2018-2019"
)

epi_antigenic_data$cluster_transition <- ifelse(epi_antigenic_data$season %in% cluster_transition,
  "antigenic_drift", "same_strain"
)
epi_antigenic_data$major_cluster_transition <- ifelse(epi_antigenic_data$season %in% major_cl_transition,
  "major_antigenic_drift", "same_cluster"
)
epi_antigenic_data$vaccine <- ifelse(epi_antigenic_data$season %in% vaccine_mismatch,
  "vaccine_mismatch", "vaccine_match"
)
epi_antigenic_data$dom_type <- ifelse(epi_antigenic_data$season %in% h3_years_all, "H3", "H1")
epi_antigenic_data$dom_type <- ifelse(epi_antigenic_data$season %in% c("2010-2011", "2018-2019"), "co-circ", epi_antigenic_data$dom_type)

epi_antigenic_data[!complete.cases(epi_antigenic_data), ]

## create dataframes for 1995-1996 and 1996-1997 so that 1997-1998 has info from previous seasons
x <- data.frame(
  season = rep("1995-1996", 10),
  region = paste("Region", seq(1:10), sep = " "),
  dom_type = rep("H1", 10),
  regional_dom_type = c("H3", "H1", "H1", "H1", "H1", "H1", "H1", "H3", "H3", "H3"),
  vaccine = "vaccine_match",
  cluster_transition = "same_strain",
  major_cluster_transition = "same_cluster"
)
unique(x$region)

y <- data.frame(
  season = rep("1996-1997", 10),
  region = paste("Region", seq(1:10), sep = " "),
  dom_type = rep("H3", 10),
  regional_dom_type = rep("H3", 10),
  vaccine = "vaccine_match",
  cluster_transition = "same_strain",
  major_cluster_transition = "same_cluster"
)
unique(y$region)

epi_antigenic_data <-
  epi_antigenic_data %>%
  mutate(regional_dom_type = case_when(
    h3_vs_h1 >= .7 ~ as.character("H3"),
    h3_vs_h1 < .7 & h3_vs_h1 >= .3 ~ as.character("co-circ"),
    h3_vs_h1 < .3 ~ as.character("H1")
  ))

epi_antigenic_data <- bind_rows(x, y, epi_antigenic_data) %>%
  arrange(season, region) %>%
  filter(!is.na(region))

epi_antigenic_data %>%
  dplyr::select(region, season, dom_type, regional_dom_type) %>%
  distinct()
unique(epi_antigenic_data$region)

prior_ag <- epi_antigenic_data %>%
  dplyr::select(season, cluster_transition) %>%
  distinct() %>%
  mutate(
    prior_transition = lag(as.character(cluster_transition), n = 1, default = NA),
    prior_transition_two_seasons = lag(as.character(cluster_transition), n = 2, default = NA)
  ) %>%
  dplyr::select(-cluster_transition)

prior_vac_mismatch <- epi_antigenic_data %>%
  dplyr::select(season, vaccine) %>%
  distinct() %>%
  mutate(
    prior_vac_mismatch = lag(as.character(vaccine), n = 1, default = NA),
    prior_vac_mismatch_two_seasons = lag(as.character(vaccine), n = 2, default = NA)
  ) %>%
  dplyr::select(-vaccine)

prior_national_dom_type <- epi_antigenic_data %>%
  dplyr::select(season, dom_type) %>%
  distinct() %>%
  mutate(
    prior_dom_type_national = lag(as.character(dom_type), n = 1, default = NA),
    prior_dom_type_national_two_seasons = lag(as.character(dom_type), n = 2, default = NA)
  ) %>%
  dplyr::select(-dom_type)

epi_antigenic_data2 <- list(epi_antigenic_data, prior_ag, prior_vac_mismatch, prior_national_dom_type) %>%
  purrr::reduce(full_join, by = c("season"))
head(epi_antigenic_data2)

##### H3 burden from prior season
H3_epi_burden <- epi_antigenic_data2 %>%
  distinct(season, region, H3_cum_intensity, H3_max_intensity) %>%
  complete(region, season) %>%
  group_by(region) %>%
  arrange(region, season) %>%
  dplyr::select(season, region, H3_cum_intensity, H3_max_intensity) %>%
  mutate(
    H3_epi_size_prior = lag(H3_cum_intensity, n = 1),
    H3_peak_intensity_prior = lag(H3_max_intensity, n = 1)
  ) %>%
  ungroup()

load("data/region_level_flu_metrics_H1.RData") # region_flu_metrics_H1 season-level H1 data
names(region_flu_metrics_H1)[3:ncol(region_flu_metrics_H1)] <- paste0("H1_", names(region_flu_metrics_H1)[3:ncol(region_flu_metrics_H1)])

H1_burden <-
  region_flu_metrics_H1 %>%
  complete(region, season) %>%
  mutate(
    H1_cum_intensity = as.numeric(H1_cum_intensity),
    H1_max_intensity = as.numeric(H1_max_intensity)
  ) %>%
  tidyr::replace_na(list(H1_cum_intensity = 0, H1_max_intensity = 0)) %>%
  group_by(region) %>%
  arrange(region, season) %>%
  dplyr::select(season, region, H1_cum_intensity, H1_max_intensity) %>%
  mutate(
    H1_epi_size_prior = lag(H1_cum_intensity, n = 1),
    H1_peak_intensity_prior = lag(H1_max_intensity, n = 1)
  ) %>%
  ungroup()
table(H1_burden$season)

load("data/region_level_flu_metrics_IVB.RData") # region_flu_metrics_IVB season-level IVB data
names(region_flu_metrics_IVB)
names(region_flu_metrics_IVB)[3:ncol(region_flu_metrics_IVB)] <- paste0("IVB_", names(region_flu_metrics_IVB)[3:ncol(region_flu_metrics_IVB)])

IVB_burden <- region_flu_metrics_IVB %>%
  complete(region, season) %>%
  tidyr::replace_na(list(IVB_cum_intensity = 0, IVB_max_intensity = 0)) %>%
  group_by(region) %>%
  arrange(region, season) %>%
  dplyr::select(season, region, IVB_cum_intensity, IVB_max_intensity) %>%
  mutate(
    IVB_epi_size_prior = lag(IVB_cum_intensity, n = 1),
    IVB_peak_intensity_prior = lag(IVB_max_intensity, n = 1)
  ) %>%
  ungroup()
table(IVB_burden$season)
IVB_burden[IVB_burden$season == "2018-2019", ]

epi_antigenic_data_final <- list(
  epi_antigenic_data2 %>% dplyr::select(-H3_cum_intensity, -H3_max_intensity),
  H3_epi_burden, H1_burden, IVB_burden
) %>%
  purrr::reduce(left_join, by = c("region", "season"))

## load pre-compiled vaccine effectiveness and vaccine coverage data
load("data/vaccine_VE_coverage_combined_df.RData")
head(vaccine_VE_coverage)

epi_antigenic_data_final <- left_join(epi_antigenic_data_final, vaccine_VE_coverage, relationship = "many-to-many", by = c("season")) %>%
  dplyr::select(-contains("two_prior")) %>%
  filter(!(season %in% c("2019-2020")))

unique(epi_antigenic_data_final$season)
unique(epi_antigenic_data$region)

####################################################################################
### reduced dataframe for variable selection analyses
####################################################################################

epi_red <- epi_antigenic_data_final %>%
  dplyr::select(-contains(c("lag3", "missing", "three", "titer_sub", "onset_se", "R0", "lag0", "regional_dom_type", "two_seasons", "H3_initial_Rt")))

epi_red <- epi_red %>% 
  filter(!(year.new <= 2009 & region == "Region 10")) %>%
  mutate(H3_max_Rt = if_else(is.na(onset_days_from_Oct1),NA,H3_max_Rt))
         # H3_shannon_entropy = if_else(is.na(onset_days_from_Oct1),NA,H3_shannon_entropy))
          
epi_red <- epi_red %>%
  mutate(H3_shannon_entropy = ifelse(H3_shannon_entropy == Inf, NA, H3_shannon_entropy)) %>%
  distinct() %>%
  filter(!(season %in% c("2019-2020", "1995-1996", "1996-1997")))
nrow(epi_red) # 207

epi_red %>% 
  filter(!(season %in% c("2000-2001", "2009-2010")))%>%
  dplyr::select(season,region,onset_days_from_Oct1,H3_season_duration,H3_shannon_entropy,H3_max_Rt) %>% 
  arrange(season,region) %>%
  filter(is.na(onset_days_from_Oct1)|is.na(H3_shannon_entropy)|is.na(H3_max_Rt)|H3_season_duration<5)

epi_red %>%
  filter(season=="2002-2003" & H3_season_duration >= 5)%>%
  dplyr::select(season,region,onset_days_from_Oct1,H3_season_duration,H3_shannon_entropy,H3_max_Rt)

epi_red %>%
  filter(season=="1997-1998" & H3_season_duration >= 5)%>%
  dplyr::select(season,region,onset_days_from_Oct1,H3_season_duration,H3_shannon_entropy,H3_max_Rt)

ggplot(regionflu_ili_vir_adj %>% filter(season_description=="2002-2003"))+
  geom_line(aes(x=wk_date,y=ili_h3_st))+
  facet_wrap(~region,scales="free_y")

ggplot(regionflu_ili_vir_adj %>% filter(season_description=="1997-1998"))+
  geom_line(aes(x=wk_date,y=ili_h3_st))+
  facet_wrap(~region,scales="free_y")

# rescale inverse shannon entropy to fall between 0 and 1
epi_red <- epi_red %>%
  mutate(H3_shannon_entropy = ifelse(season %in% c("2000-2001", "2009-2010") | H3_season_duration < 5, NA, H3_shannon_entropy)) %>%
  metan::resca(H3_shannon_entropy, new_min = 0, new_max = 1, na.rm = T) %>%
  ungroup()
range(epi_red$H3_shannon_entropy_res, na.rm = T)

epi_red %>%
  dplyr::select(region, season, H3_shannon_entropy_res) %>%
  arrange(H3_shannon_entropy_res) # missing for 2000-2001, most of 2002-2003, and 2009-2010 (no H3N2 circulation)

epi_red %>%
  dplyr::select(region, season, H3_shannon_entropy_res) %>%
  drop_na() %>%
  group_by(season) %>%
  tally()

table(epi_red$region)
table(epi_red$season)

epi_red <- epi_red %>%
  dplyr::select(
    season, region, H3_shannon_entropy_res,
    contains(c(
      "vs", "cum", "Rt", "max", "dom", "duration", "peak_diff", "index", "epi_size_prior",
      "onset_week", "peak_week", "vac_cov", "vac_combined", "VE", "year.new", "std",
      "mean", "lag1", "lag2", "lbi", "total", "days_from", "timing_sd", "bayes",
      "samples"
    ))
  )
sort(names(epi_red))
unique(epi_red$season)
nrow(epi_red)
epi_red %>% dplyr::select(season,region,contains("onset"))
save(epi_red, file = "data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")
