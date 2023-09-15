###############################
## Create NA phylogenetic dataset
###############################

## load packages
list.of.packages <- c("dplyr", "tidyr", "readr", "timetk", "ggplot2", "cowplot", "purrr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

##############################################################
## mean direct distances
##############################################################
mean_dist_table <- readr::read_tsv("2_Phylo_Dataset/distance_tables/north-america_all_mean_seasonal_distances_h3n2_na_21y.tsv")
head(mean_dist_table)

mean_dist_table2 <- mean_dist_table %>%
  mutate(
    time_diff = as.numeric(format(as.Date(other_season_start), "%Y")) - as.numeric(format(as.Date(current_season_start), "%Y")),
    season1 = paste(as.numeric(format(as.Date(current_season_start), "%Y")), as.numeric(format(as.Date(current_season_end), "%Y")), sep = "-"),
    season2 = paste(as.numeric(format(as.Date(other_season_start), "%Y")), as.numeric(format(as.Date(other_season_end), "%Y")), sep = "-")
  )
head(mean_dist_table2)

within_season <- mean_dist_table2 %>%
  filter(time_diff == 0) %>%
  dplyr::rename(season = season1) %>%
  dplyr::select(replicate, season, distance_map, mean_distance) %>%
  spread(distance_map, mean_distance) %>%
  ungroup()
head(within_season)
names(within_season)[names(within_season) %in% c("bhatt", "bhatt_nonepitope", "krammer")] <- c(
  "NA_bhatt_ep_lag0", "NA_bhatt_nonepitope_lag0",
  "NA_krammer_ep_lag0"
)

one_season_lag <- mean_dist_table2 %>%
  filter(time_diff == 1) %>%
  dplyr::rename(season = season2) %>%
  dplyr::select(replicate, season, distance_map, mean_distance) %>%
  spread(distance_map, mean_distance) %>%
  ungroup()
head(one_season_lag)
names(one_season_lag)[names(one_season_lag) %in% c("bhatt", "bhatt_nonepitope", "krammer")] <- c(
  "NA_bhatt_ep_lag1", "NA_bhatt_nonepitope_lag1",
  "NA_krammer_ep_lag1"
)

two_season_lag <- mean_dist_table2 %>%
  dplyr::rename(season = season2) %>%
  filter(time_diff == 2) %>%
  dplyr::select(replicate, season, distance_map, mean_distance) %>%
  spread(distance_map, mean_distance) %>%
  ungroup()
head(two_season_lag)
names(two_season_lag)[names(two_season_lag) %in% c("bhatt", "bhatt_nonepitope", "krammer")] <- c(
  "NA_bhatt_ep_lag2", "NA_bhatt_nonepitope_lag2",
  "NA_krammer_ep_lag2"
)

mean_pairwise_dist_table <- left_join(within_season, one_season_lag, by = c("season", "replicate"))

nrow(mean_pairwise_dist_table)
mean_pairwise_dist_table <- left_join(mean_pairwise_dist_table,
  two_season_lag,
  by = c("season", "replicate")
)

mean_pairwise_dist_table <- mean_pairwise_dist_table %>% arrange(season)
mean_pairwise_dist_table %>% dplyr::select(season, replicate, contains("lag2"))
names(mean_pairwise_dist_table)

load("data/NA_mean_bhatt_epitope_change_north_amer_build_summary.Rdata")
within_season_bhatt_ep <- NA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 0) %>%
  dplyr::select(season, NA.ep.global.mean, NA.ep.global.sd, replicate) %>%
  rename(NA_bhatt_lag0_mean = NA.ep.global.mean, NA_bhatt_lag0_sd = NA.ep.global.sd)

one_season_lag_bhatt_ep <- NA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 1) %>%
  dplyr::select(season, NA.ep.global.mean, NA.ep.global.sd, replicate) %>%
  rename(NA_bhatt_lag1_mean = NA.ep.global.mean, NA_bhatt_lag1_sd = NA.ep.global.sd)

two_season_lag_bhatt_ep <- NA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 2) %>%
  dplyr::select(season, NA.ep.global.mean, NA.ep.global.sd, replicate) %>%
  rename(NA_bhatt_lag2_mean = NA.ep.global.mean, NA_bhatt_lag2_sd = NA.ep.global.sd)

load("data/NA_mean_bhatt_non_epitope_change_north_amer_build_summary.Rdata")
within_season_bhatt_non_ep <- NA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 0) %>%
  dplyr::select(season, NA.ep.global.mean, NA.ep.global.sd, replicate) %>%
  rename(NA_bhatt_non_ep_lag0_mean = NA.ep.global.mean, NA_bhatt_non_ep_lag0_sd = NA.ep.global.sd)

one_season_lag_bhatt_non_ep <- NA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 1) %>%
  dplyr::select(season, NA.ep.global.mean, NA.ep.global.sd, replicate) %>%
  rename(NA_bhatt_non_ep_lag1_mean = NA.ep.global.mean, NA_bhatt_non_ep_lag1_sd = NA.ep.global.sd)

two_season_lag_bhatt_non_ep <- NA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 2) %>%
  dplyr::select(season, NA.ep.global.mean, NA.ep.global.sd, replicate) %>%
  rename(NA_bhatt_non_ep_lag2_mean = NA.ep.global.mean, NA_bhatt_non_ep_lag2_sd = NA.ep.global.sd)

load("data/NA_mean_krammer_epitope_change_north_amer_build_summary.Rdata")
within_season_krammer_ep <- NA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 0) %>%
  dplyr::select(season, NA.ep.global.mean, NA.ep.global.sd, replicate) %>%
  rename(NA_krammer_lag0_mean = NA.ep.global.mean, NA_krammer_lag0_sd = NA.ep.global.sd)

one_season_lag_krammer_ep <- NA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 1) %>%
  dplyr::select(season, NA.ep.global.mean, NA.ep.global.sd, replicate) %>%
  rename(NA_krammer_lag1_mean = NA.ep.global.mean, NA_krammer_lag1_sd = NA.ep.global.sd)

two_season_lag_krammer_ep <- NA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 2) %>%
  dplyr::select(season, NA.ep.global.mean, NA.ep.global.sd, replicate) %>%
  rename(NA_krammer_lag2_mean = NA.ep.global.mean, NA_krammer_lag2_sd = NA.ep.global.sd)

mean_pairwise_dist_table_lim <-
  list(
    within_season_bhatt_ep, within_season_bhatt_non_ep,
    one_season_lag_bhatt_ep, one_season_lag_bhatt_non_ep,
    two_season_lag_bhatt_ep, two_season_lag_bhatt_non_ep,
    within_season_krammer_ep, one_season_lag_krammer_ep, two_season_lag_krammer_ep
  ) %>%
  purrr::reduce(full_join, by = c("season", "replicate")) %>%
  filter(season != "1995-1996") %>%
  dplyr::select(season, replicate, contains("mean"))

head(mean_pairwise_dist_table_lim)

mean_pairwise_dist_table_lim %>% dplyr::select(season, replicate, contains("lag1_mean"))
names(mean_pairwise_dist_table_lim) <- gsub("_mean", "", names(mean_pairwise_dist_table_lim))
mean_pairwise_dist_table_lim$replicate <- readr::parse_number(mean_pairwise_dist_table_lim$replicate)
mean_pairwise_dist_table_lim %>% dplyr::select(-contains(c("lag0", "lag1")))

mean_pairwise_dist_table$NA_bhatt_ep_lag1[1:5] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1996-1997") %>%
  pull(NA_bhatt_lag1)
mean_pairwise_dist_table$NA_krammer_ep_lag1[1:5] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1996-1997") %>%
  pull(NA_krammer_lag1)
mean_pairwise_dist_table$NA_bhatt_nonepitope_lag1[1:5] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1996-1997") %>%
  pull(NA_bhatt_non_ep_lag1)


mean_pairwise_dist_table$NA_bhatt_ep_lag2[6:10] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1997-1998") %>%
  pull(NA_bhatt_lag2)
mean_pairwise_dist_table$NA_krammer_ep_lag2[6:10] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1997-1998") %>%
  pull(NA_krammer_lag2)
mean_pairwise_dist_table$NA_bhatt_nonepitope_lag2[6:10] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1997-1998") %>%
  pull(NA_bhatt_non_ep_lag2)

##############################################################
## LBI
##############################################################
mean_lbi_table <- readr::read_tsv("2_Phylo_Dataset/distance_tables/north-america_mean_seasonal_lbi_h3n2_na_21y.txt")

mean_lbi_table <- mean_lbi_table %>%
  mutate(
    year1 = as.numeric(format(as.Date(season_start), "%Y")),
    year2 = as.numeric(format(as.Date(season_end), "%Y"))
  ) %>%
  mutate(season = paste(year1, year2, sep = "-")) %>%
  dplyr::select(-season_start, -season_end)
names(mean_lbi_table)[names(mean_lbi_table) %in% c("mean", "std")] <- c("NA_mean_lbi", "NA_std_lbi")

mean_lbi_table <- mean_lbi_table %>%
  group_by(replicate) %>%
  dplyr::mutate(
    NA_mean_lbi_lag1 = lag(NA_mean_lbi, n = 1),
    NA_mean_lbi_lag2 = lag(NA_mean_lbi, n = 2),
    NA_std_lbi_lag1 = lag(NA_std_lbi, n = 1),
    NA_std_lbi_lag2 = lag(NA_std_lbi, n = 2)
  ) %>%
  ungroup()

seas.ag <- left_join(mean_pairwise_dist_table, mean_lbi_table, by = c("season", "replicate"))
save(seas.ag, file = "data/north_amer_build_season_h3n2_replicates_NA_direct_ag_distances.RData")
