###############################
## Create NA phylogenetic dataset
###############################

## load packages
list.of.packages <- c("dplyr", "tidyr", "readr", "timetk", "ggplot2", "cowplot", "purrr","fitdistrplus")

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
# mean_lbi_table <- readr::read_tsv("2_Phylo_Dataset/distance_tables/north-america_mean_seasonal_lbi_h3n2_na_21y.txt")
# head(mean_lbi_table)
# 
# mean_lbi_table <- mean_lbi_table %>%
#   mutate(
#     year1 = as.numeric(format(as.Date(season_start), "%Y")),
#     year2 = as.numeric(format(as.Date(season_end), "%Y"))
#   ) %>%
#   mutate(season = paste(year1, year2, sep = "-")) %>%
#   dplyr::select(-season_start, -season_end)
# names(mean_lbi_table)[names(mean_lbi_table) %in% c("mean", "std")] <- c("NA_mean_lbi", "NA_std_lbi")
# head(mean_lbi_table)
# 
# mean_lbi_table <- mean_lbi_table %>%
#   group_by(replicate) %>%
#   dplyr::mutate(
#     NA_mean_lbi_lag1 = lag(NA_mean_lbi, n = 1),
#     NA_mean_lbi_lag2 = lag(NA_mean_lbi, n = 2),
#     NA_std_lbi_lag1 = lag(NA_std_lbi, n = 1),
#     NA_std_lbi_lag2 = lag(NA_std_lbi, n = 2)
#   ) %>%
#   ungroup()
# 
lbi_tb_na <- read_tsv("2_Phylo_Dataset/distance_tables/north-america_strain_seasonal_lbi_h3n2_na_21y.tsv.gz")
lbi_tb_na$year2 <- lubridate::year(lbi_tb_na$season_end)
lbi_tb_na$year1 <- lbi_tb_na$year2 - 1
lbi_tb_na$season <- paste(lbi_tb_na$year1, lbi_tb_na$year2, sep = "-")
lbi_tb_na$log_lbi <- log(lbi_tb_na$lbi)
lbi_tb_na = lbi_tb_na %>% filter(season!="2019-2020")

ggplot(lbi_tb_na) +
  geom_density(aes(lbi, group = replicate, color = replicate)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~season) +
  theme_bw()

ggplot(lbi_tb_na) +
  geom_density(aes(log_lbi, group = replicate, color = replicate)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~season, scales = "free_y") +
  theme_bw()

# descdist(lbi_tb_na$lbi, boot = 1000)
# descdist(lbi_tb_na$log_lbi, boot = 1000)

lbi_tb_na$season = as.factor(lbi_tb_na$season)
lbi_tb_na$replicate = as.factor(lbi_tb_na$replicate)

lbi_tb_na$season_rep = paste(lbi_tb_na$season,lbi_tb_na$replicate,sep="-")
season_rep_list = unique(lbi_tb_na$season_rep)

library(boot)
library(magicfor)
magic_for(silent = T, temporary = T)
for (i in season_rep_list) {

  lbi_df <- lbi_tb_na %>% filter(season_rep==i) %>% dplyr::select(season_rep,season,year1,year2,replicate,lbi) %>% as.data.frame()

  set.seed(6929)    # Make the results reproducible
  r.mean = boot::boot(lbi_df$lbi, statistic = function(x, inds) mean(x[inds]), R = 1000)
  r.sd = boot::boot(lbi_df$lbi, statistic = function(x, inds) sd(x[inds]), R = 1000)

  sample_mean = mean(r.mean$t[,1])
  sample_sd = mean(r.sd$t[,1])
  season = unique(lbi_df$season)
  replicate = unique(lbi_df$replicate)
  year1 = unique(lbi_df$year1)
  year2 = unique(lbi_df$year2)
  put(season,replicate,year1,year2,sample_mean,sample_sd)
  }
lbi_na_summary <- magic_result_as_dataframe()
lbi_na_summary

ggplot(lbi_na_summary)+
  geom_line(aes(x=year2,sample_mean,color=replicate))

mean_na_lbi_table <- lbi_na_summary %>%
  group_by(replicate) %>%
  dplyr::mutate(
    NA_mean_lbi_lag1 = lag(sample_mean, n = 1),
    NA_std_lbi_lag1 = lag(sample_sd, n = 1)
  )%>%
  ungroup()%>%
  rename(NA_mean_lbi = sample_mean, NA_std_lbi = sample_sd)%>%
  dplyr::select(-i,-year1,-year2)
mean_na_lbi_table

# log_mean_lbi <- lbi_tb_na %>%
#   group_by(replicate, season) %>%
#   summarize(
#     NA_mean_log_lbi = mean(log_lbi),
#     NA_std_log_lbi = sd(log_lbi)
#   )

# log_mean_lbi_table <- log_mean_lbi %>%
#   group_by(replicate) %>%
#   dplyr::mutate(
#     NA_mean_log_lbi_lag1 = lag(NA_mean_log_lbi, n = 1),
#     NA_mean_log_lbi_lag2 = lag(NA_mean_log_lbi, n = 2),
#     NA_std_log_lbi_lag1 = lag(NA_std_log_lbi, n = 1),
#     NA_std_log_lbi_lag2 = lag(NA_std_log_lbi, n = 2)
#   ) %>%
#   ungroup()
# log_mean_lbi_table
# seas.ag.NA <- full_join(mean_pairwise_dist_table, log_mean_lbi_table, by = c("season", "replicate"))

mean_pairwise_dist_table$replicate = as.factor(mean_pairwise_dist_table$replicate)
seas.ag.NA <- full_join(mean_pairwise_dist_table, mean_na_lbi_table %>% dplyr::select(season,replicate,NA_std_lbi), by = c("season", "replicate"))
names(seas.ag.NA)
save(seas.ag.NA, file = "data/north_amer_build_season_h3n2_replicates_NA_direct_ag_distances.RData")
