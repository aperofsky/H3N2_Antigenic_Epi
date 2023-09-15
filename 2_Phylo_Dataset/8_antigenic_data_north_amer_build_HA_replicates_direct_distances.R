###############################
## Create HA phylogenetic dataset
###############################

## load packages
list.of.packages <- c("dplyr", "tidyr", "readr", "timetk", "ggplot2", "cowplot", "purrr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

############################################################
## distance to origin
############################################################
## H3N2 antigenic data
gen <- readr::read_tsv("2_Phylo_Dataset/auspice_tables/flu_seasonal_h3n2_ha_21y_north-america.tsv")
names(gen)
head(gen)
unique(gen$replicate)

gen <- gen %>%
  mutate(date = as.Date(date, "%y/%m/%d")) %>%
  mutate(
    month = as.numeric(format(date, "%m")),
    year = as.numeric(format(date, "%Y"))
  ) %>%
  arrange(date) %>%
  mutate(month.date = as.Date(paste(year, month, 1, sep = "-"), "%Y-%m-%d")) # first of the month

range(gen$date) # "1995-12-31" "2019-10-01"

## add in season indicator
gen <- gen %>%
  mutate(month.day = format(date, "%m-%d")) %>%
  mutate(season = ifelse(month.day < "07-01", sprintf("%d-%d", year - 1, year), # season = july 1 to june 30
    sprintf("%d-%d", year, year + 1)
  )) %>%
  mutate(year.new = as.numeric(substr(season, 1, 4))) ## first year of the flu season
##############################################################
## season-level dataset with all strains
##############################################################
## combine replicates
# mean distance from origin
keep <- gen %>%
  dplyr::select(
    starts_with(c("cTiter", "ep", "ne", "rb", "stem_ep")),
    -contains("vaccine"), -contains("seasonal")
  ) %>%
  names()
keep

seas.ag.global <- gen %>%
  mutate(month.date = format(as.Date(date), "%m-%d")) %>%
  group_by(replicate, season, year.new) %>%
  summarise(across(all_of(keep), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()
names(seas.ag.global)[4:ncol(seas.ag.global)] <- paste0("HA_", names(seas.ag.global)[4:ncol(seas.ag.global)])

seas.ag.global %>% filter(season == "1997-1998")
########################################
## seasonal lags
########################################
keep
names(seas.ag.global)

keep_col <- seas.ag.global %>%
  dplyr::select(contains(keep)) %>%
  dplyr::select(-contains(c("seasonal", "year", "replicate"))) %>%
  names()

seas.ag.global.summary <- seas.ag.global %>%
  arrange(year.new) %>%
  group_by(replicate) %>%
  tk_augment_lags(.value = all_of(keep_col), .lags = 1:2) %>%
  ungroup()
seas.ag.global.summary %>% dplyr::select(season, contains("lag1"))
seas.ag.global.summary %>% dplyr::select(contains("cTiter"))

seas.ag.global <- seas.ag.global.summary %>%
  group_by(replicate) %>%
  arrange(year.new) %>%
  mutate( # one season lag
    diff.HI.sub.global.mean = HA_cTiterSub - HA_cTiterSub_lag1,
    diff.HI.tree.global.mean = HA_cTiter - HA_cTiter_lag1,
    # two season lag
    diff.HI.sub.two.global.mean = HA_cTiterSub - HA_cTiterSub_lag2,
    diff.HI.tree.two.global.mean = HA_cTiter - HA_cTiter_lag2
  ) %>%
  ungroup()

df_1996_1997 <- seas.ag.global %>%
  filter(season %in% c("1996-1997")) %>%
  dplyr::select(replicate, contains(c("lag1")), -contains(c("early", "sd"))) %>%
  dplyr::select(replicate, contains("cTiter"))
head(df_1996_1997)
df_1996_1997$replicate <- readr::parse_number(df_1996_1997$replicate)

df_1997_1998 <- seas.ag.global %>%
  filter(season %in% c("1997-1998")) %>%
  dplyr::select(replicate, contains(c("lag2")), -contains(c("early", "sd"))) %>%
  dplyr::select(replicate, contains("cTiter"))
df_1997_1998$replicate <- readr::parse_number(df_1997_1998$replicate)
##############################################################
## mean direct distances
##############################################################
mean_dist_table <- readr::read_tsv("2_Phylo_Dataset/distance_tables/north-america_all_mean_seasonal_distances_h3n2_ha_21y.tsv")

unique(mean_dist_table$distance_map)

mean_dist_table2 <- mean_dist_table %>%
  mutate(
    time_diff = as.numeric(format(as.Date(other_season_start), "%Y")) - as.numeric(format(as.Date(current_season_start), "%Y")),
    season1 = paste(as.numeric(format(as.Date(current_season_start), "%Y")), as.numeric(format(as.Date(current_season_end), "%Y")), sep = "-"),
    season2 = paste(as.numeric(format(as.Date(other_season_start), "%Y")), as.numeric(format(as.Date(other_season_end), "%Y")), sep = "-")
  )

within_season <- mean_dist_table2 %>%
  filter(time_diff == 0) %>%
  dplyr::rename(season = season1) %>%
  dplyr::select(replicate, season, distance_map, mean_distance) %>%
  spread(distance_map, mean_distance)


names(within_season)[3:9] <- c(
  "HA_koel_lag0", "HA_stem_lag0", "HA_stem_ep_lag0", "HA_titer_sub_lag0", "HA_titer_tree_lag0",
  "HA_wolf_lag0", "HA_wolf_nonepitope_lag0"
)

one_season_lag <- mean_dist_table2 %>%
  filter(time_diff == 1) %>%
  dplyr::rename(season = season2) %>%
  dplyr::select(replicate, season, distance_map, mean_distance) %>%
  spread(distance_map, mean_distance) %>%
  ungroup()
names(one_season_lag)
head(one_season_lag)
names(one_season_lag)[3:9] <- c(
  "HA_koel_lag1", "HA_stem_lag1", "HA_stem_ep_lag1", "HA_titer_sub_lag1", "HA_titer_tree_lag1",
  "HA_wolf_lag1", "HA_wolf_nonepitope_lag1"
)

two_season_lag <- mean_dist_table2 %>%
  dplyr::rename(season = season2) %>%
  filter(time_diff == 2) %>%
  dplyr::select(replicate, season, distance_map, mean_distance) %>%
  spread(distance_map, mean_distance) %>%
  ungroup()
names(two_season_lag)

names(two_season_lag)[3:9] <- c(
  "HA_koel_lag2", "HA_stem_lag2", "HA_stem_ep_lag2", "HA_titer_sub_lag2", "HA_titer_tree_lag2",
  "HA_wolf_lag2", "HA_wolf_nonepitope_lag2"
)

mean_pairwise_dist_table <- left_join(within_season, one_season_lag, by = c("season", "replicate"))
mean_pairwise_dist_table <- left_join(mean_pairwise_dist_table, two_season_lag, by = c("season", "replicate"))
mean_pairwise_dist_table %>%
  dplyr::select(season, replicate, contains("wolf_lag")) %>%
  head()
mean_pairwise_dist_table <- mean_pairwise_dist_table %>% arrange(season)

load("data/HA_mean_epitope_change_north_amer_build_summary.Rdata")
head(HA_north_amer_dist_epitope_summary)
within_season_wolf_ep <- HA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 0) %>%
  dplyr::select(season, HA.ep.global.mean, HA.ep.global.sd, replicate) %>%
  rename(HA_wolf_lag0_mean = HA.ep.global.mean, HA_wolf_lag0_sd = HA.ep.global.sd)

one_season_lag_wolf_ep <- HA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 1) %>%
  dplyr::select(season, HA.ep.global.mean, HA.ep.global.sd, replicate) %>%
  rename(HA_wolf_lag1_mean = HA.ep.global.mean, HA_wolf_lag1_sd = HA.ep.global.sd)

two_season_lag_wolf_ep <- HA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 2) %>%
  dplyr::select(season, HA.ep.global.mean, HA.ep.global.sd, replicate) %>%
  rename(HA_wolf_lag2_mean = HA.ep.global.mean, HA_wolf_lag2_sd = HA.ep.global.sd)

load("data/HA_mean_koel_change_north_amer_build_summary.Rdata")
head(HA_north_amer_dist_epitope_summary)
nrow(HA_north_amer_dist_epitope_summary)
within_season_koel_ep <- HA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 0) %>%
  dplyr::select(season, HA.ep.global.mean, HA.ep.global.sd, replicate) %>%
  rename(HA_koel_lag0_mean = HA.ep.global.mean, HA_koel_lag0_sd = HA.ep.global.sd)

one_season_lag_koel_ep <- HA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 1) %>%
  dplyr::select(season, HA.ep.global.mean, HA.ep.global.sd, replicate) %>%
  rename(HA_koel_lag1_mean = HA.ep.global.mean, HA_koel_lag1_sd = HA.ep.global.sd)

two_season_lag_koel_ep <- HA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 2) %>%
  dplyr::select(season, HA.ep.global.mean, HA.ep.global.sd, replicate) %>%
  rename(HA_koel_lag2_mean = HA.ep.global.mean, HA_koel_lag2_sd = HA.ep.global.sd)

load("data/HA_mean_non_epitope_change_north_amer_build_summary.Rdata")
head(HA_north_amer_dist_epitope_summary)
within_season_wolf_nonep <- HA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 0) %>%
  dplyr::select(season, HA.ep.north_amer.mean, HA.ep.north_amer.sd, replicate) %>%
  rename(HA_wolf_non_ep_lag0_mean = HA.ep.north_amer.mean, HA_wolf_non_ep_lag0_sd = HA.ep.north_amer.sd)

one_season_lag_wolf_nonep <- HA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 1) %>%
  dplyr::select(season, HA.ep.north_amer.mean, HA.ep.north_amer.sd, replicate) %>%
  rename(HA_wolf_non_ep_lag1_mean = HA.ep.north_amer.mean, HA_wolf_non_ep_lag1_sd = HA.ep.north_amer.sd)

two_season_lag_wolf_nonep <- HA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 2) %>%
  dplyr::select(season, HA.ep.north_amer.mean, HA.ep.north_amer.sd, replicate) %>%
  rename(HA_wolf_non_ep_lag2_mean = HA.ep.north_amer.mean, HA_wolf_non_ep_lag2_sd = HA.ep.north_amer.sd)

load("data/HA_mean_stem_footprint_change_north_amer_build_summary.Rdata")
head(HA_north_amer_dist_epitope_summary)
within_season_stem_ep <- HA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 0) %>%
  dplyr::select(season, HA.ep.north_amer.mean, HA.ep.north_amer.sd, replicate) %>%
  rename(HA_stem_ep_lag0_mean = HA.ep.north_amer.mean, HA_stem_ep_lag0_sd = HA.ep.north_amer.sd)

one_season_lag_stem_ep <- HA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 1) %>%
  dplyr::select(season, HA.ep.north_amer.mean, HA.ep.north_amer.sd, replicate) %>%
  rename(HA_stem_ep_lag1_mean = HA.ep.north_amer.mean, HA_stem_ep_lag1_sd = HA.ep.north_amer.sd)

two_season_lag_stem_ep <- HA_north_amer_dist_epitope_summary %>%
  filter(season_diff == 2) %>%
  dplyr::select(season, HA.ep.north_amer.mean, HA.ep.north_amer.sd, replicate) %>%
  rename(HA_stem_ep_lag2_mean = HA.ep.north_amer.mean, HA_stem_ep_lag2_sd = HA.ep.north_amer.sd)

mean_pairwise_dist_table_lim <-
  list(
    within_season_wolf_ep, within_season_koel_ep, within_season_wolf_nonep, within_season_stem_ep,
    one_season_lag_wolf_ep, one_season_lag_koel_ep, one_season_lag_wolf_nonep, one_season_lag_stem_ep,
    two_season_lag_wolf_ep, two_season_lag_koel_ep, two_season_lag_wolf_nonep, two_season_lag_stem_ep
  ) %>%
  purrr::reduce(full_join, by = c("season", "replicate")) %>%
  filter(season != "1995-1996") %>%
  dplyr::select(season, replicate, contains("mean"))

names(mean_pairwise_dist_table_lim) <- gsub("_mean", "", names(mean_pairwise_dist_table_lim))
names(mean_pairwise_dist_table_lim)
mean_pairwise_dist_table_lim$replicate <- readr::parse_number(mean_pairwise_dist_table_lim$replicate)
mean_pairwise_dist_table_lim %>% dplyr::select(-contains(c("lag0", "lag1")))

mean_pairwise_dist_table %>% dplyr::select(season, replicate, HA_koel_lag1)
mean_pairwise_dist_table$HA_koel_lag1[1:5] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1996-1997") %>%
  pull(HA_koel_lag1)
mean_pairwise_dist_table$HA_stem_ep_lag1[1:5] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1996-1997") %>%
  pull(HA_stem_ep_lag1)
mean_pairwise_dist_table$HA_wolf_lag1[1:5] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1996-1997") %>%
  pull(HA_wolf_lag1)
mean_pairwise_dist_table$HA_wolf_nonepitope_lag1[1:5] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1996-1997") %>%
  pull(HA_wolf_non_ep_lag1)
mean_pairwise_dist_table$HA_titer_sub_lag1[1:5] <- df_1996_1997$HA_cTiterSub_lag1
mean_pairwise_dist_table$HA_titer_tree_lag1[1:5] <- df_1996_1997$HA_cTiter_lag1

mean_pairwise_dist_table$HA_koel_lag2[6:10] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1997-1998") %>%
  pull(HA_koel_lag2)
mean_pairwise_dist_table$HA_stem_ep_lag2[6:10] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1997-1998") %>%
  pull(HA_stem_ep_lag2)
mean_pairwise_dist_table$HA_wolf_lag2[6:10] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1997-1998") %>%
  pull(HA_wolf_lag2)
mean_pairwise_dist_table$HA_wolf_nonepitope_lag2[6:10] <- mean_pairwise_dist_table_lim %>%
  filter(season == "1997-1998") %>%
  pull(HA_wolf_non_ep_lag2)
mean_pairwise_dist_table$HA_titer_sub_lag2[6:10] <- df_1997_1998$HA_cTiterSub_lag2
mean_pairwise_dist_table$HA_titer_tree_lag2[6:10] <- df_1997_1998$HA_cTiter_lag2

##############################################################
## LBI
##############################################################
mean_lbi_table <- readr::read_tsv("2_Phylo_Dataset//distance_tables/north-america_mean_seasonal_lbi_h3n2_ha_21y.tsv")
head(mean_lbi_table)

mean_lbi_table <- mean_lbi_table %>%
  mutate(
    year1 = as.numeric(format(as.Date(season_start), "%Y")),
    year2 = as.numeric(format(as.Date(season_end), "%Y"))
  ) %>%
  mutate(season = paste(year1, year2, sep = "-")) %>%
  dplyr::select(-season_start, -season_end)
names(mean_lbi_table)[names(mean_lbi_table) %in% c("mean", "std")] <- c("HA_mean_lbi", "HA_std_lbi")
head(mean_lbi_table)

mean_lbi_table <- mean_lbi_table %>%
  group_by(replicate) %>%
  dplyr::mutate(
    HA_mean_lbi_lag1 = lag(HA_mean_lbi, n = 1),
    HA_mean_lbi_lag2 = lag(HA_mean_lbi, n = 2),
    HA_std_lbi_lag1 = lag(HA_std_lbi, n = 1),
    HA_std_lbi_lag2 = lag(HA_std_lbi, n = 2)
  ) %>%
  ungroup()

seas.ag <- full_join(mean_pairwise_dist_table, mean_lbi_table, by = c("season", "replicate"))

save(seas.ag, file = "data/north_amer_build_season_h3n2_replicates_direct_ag_distances.RData")
