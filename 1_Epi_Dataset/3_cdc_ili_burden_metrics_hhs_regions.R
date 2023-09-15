## load packages
list.of.packages <- c("cdcfluview", "dplyr", "ggplot2", "cowplot", "zoo", "tidyr", "vegan")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

## load data
load("data/hhs_division_level_ILI_and_virology_interp_smoothed.RData") # regionflu_ili_vir_adj

head(regionflu_ili_vir_adj)
names(regionflu_ili_vir_adj)
range(regionflu_ili_vir_adj$wk_date) # ""1997-09-28" "2021-09-26"

regionflu_ili_vir_adj <- regionflu_ili_vir_adj[order(regionflu_ili_vir_adj$region, regionflu_ili_vir_adj$wk_date), ]

# remove region 10 before 2008
regionflu_ili_vir_adj <-
  regionflu_ili_vir_adj %>%
  filter(!(wk_date < "2008-06-01" & region == "Region 10"))

unique(regionflu_ili_vir_adj$season_description)

seasons <- unique(regionflu_ili_vir_adj$season_description)
seasons
seasons <- seasons[!(seasons %in% c("2019-2020", "2020-2021", "2021-2022"))] # take out current season

regions <- unique(regionflu_ili_vir_adj$region)
regions
####################################################
## correct for sampling effort pre and post pandemic
####################################################
regionflu_ili_vir_adj %>% filter(wk_date == "2009-04-19") # week 16
pre_pandemic <- regionflu_ili_vir_adj$wk_date[regionflu_ili_vir_adj$wk_date < "2009-04-19"] # first infection in US March 28 (week 13); US reporting changed week 16
post_pandemic <- regionflu_ili_vir_adj$wk_date[regionflu_ili_vir_adj$wk_date >= "2010-08-10"] # WHO declared end of pandemic on August 10 2010
colnames(regionflu_ili_vir_adj)

pre_pan <- regionflu_ili_vir_adj %>%
  filter(wk_date %in% pre_pandemic) %>%
  group_by(region) %>%
  dplyr::summarize(mean(adj_ili_per_pos, na.rm = TRUE)) %>%
  filter(!(region == "Region 10"))

post_pan <- regionflu_ili_vir_adj %>%
  filter(wk_date %in% post_pandemic) %>%
  group_by(region) %>%
  dplyr::summarize(mean(adj_ili_per_pos, na.rm = TRUE)) %>%
  filter(!(region == "Region 10"))

ili_pos_ratio <- post_pan[, 2] / pre_pan[, 2]
ili_pos_ratio$region <- post_pan$region
colnames(ili_pos_ratio)[1] <- "ili_ratio"
ili_pos_ratio

regionflu_ili_vir_adj <- left_join(regionflu_ili_vir_adj, ili_pos_ratio, by = "region")

regionflu_ili_vir_adj <-
  regionflu_ili_vir_adj %>%
  mutate(
    ili_h3 = if_else(wk_date %in% pre_pandemic, ILI_interp_H3 * ili_ratio, ILI_interp_H3),
    ili_h1 = if_else(wk_date %in% pre_pandemic, ILI_interp_H1 * ili_ratio, ILI_interp_H1),
    ili_ivb = if_else(wk_date %in% pre_pandemic, ILI_interp_B * ili_ratio, ILI_interp_B),
    ili_iva = if_else(wk_date %in% pre_pandemic, ILI_interp_A * ili_ratio, ILI_interp_A),
    ili_pos = if_else(wk_date %in% pre_pandemic, ILI_interp_pos * ili_ratio, ILI_interp_pos)
  )

##############################################################################
### account for systematic differences across regions
##############################################################################
reporting_diff <-
  regionflu_ili_vir_adj %>%
  group_by(region) %>%
  dplyr::summarize(
    h3_mean = mean(ili_h3, na.rm = T),
    h3_sd = sd(ili_h3, na.rm = T),
    h1_mean = mean(ili_h1, na.rm = T),
    h1_sd = sd(ili_h1, na.rm = T),
    iva_mean = mean(ili_iva, na.rm = T),
    iva_sd = sd(ili_iva, na.rm = T),
    ivb_mean = mean(ili_ivb, na.rm = T),
    ivb_sd = sd(ili_ivb, na.rm = T),
    pos_mean = mean(ili_pos, na.rm = T),
    pos_sd = sd(ili_pos, na.rm = T)
  )

regionflu_ili_vir_adj <- left_join(regionflu_ili_vir_adj, reporting_diff, by = "region")

regionflu_ili_vir_adj <-
  regionflu_ili_vir_adj %>%
  group_by(region) %>%
  mutate(
    ili_h3_st = ili_h3 / pos_mean,
    ili_h1_st = ili_h1 / pos_mean,
    ili_ivb_st = ili_ivb / pos_mean,
    ili_iva_st = ili_iva / pos_mean,
    ili_pos_st = ili_pos / pos_mean
  ) %>%
  ungroup()

##########################################
## correct for age distribution pre and post pandemic
##########################################
pre_pandemic <- regionflu_ili_vir_adj$wk_date[regionflu_ili_vir_adj$wk_date < "2009-04-19"]
post_pandemic <- regionflu_ili_vir_adj$wk_date[regionflu_ili_vir_adj$wk_date >= "2010-08-10"]

pre_pan_age_0_4 <- regionflu_ili_vir_adj %>%
  filter(wk_date %in% pre_pandemic) %>%
  group_by(region) %>%
  dplyr::summarize(mean(age_0_4, na.rm = TRUE)) %>%
  filter(!(region == "Region 10"))

post_pan_age_0_4 <- regionflu_ili_vir_adj %>%
  filter(wk_date %in% post_pandemic) %>%
  group_by(region) %>%
  dplyr::summarize(mean(age_0_4, na.rm = TRUE)) %>%
  filter(!(region == "Region 10"))

ili_pos_ratio_age_0_4 <- post_pan_age_0_4[, 2] / pre_pan_age_0_4[, 2]
ili_pos_ratio_age_0_4$region <- post_pan_age_0_4$region
colnames(ili_pos_ratio_age_0_4)[1] <- "ili_ratio_age_0_4"
ili_pos_ratio_age_0_4 # has increased significantly except region 9

pre_pan_age_5_24 <- regionflu_ili_vir_adj %>%
  filter(wk_date %in% pre_pandemic) %>%
  group_by(region) %>%
  dplyr::summarize(mean(age_5_24, na.rm = TRUE)) %>%
  filter(!(region == "Region 10"))

post_pan_age_5_24 <- regionflu_ili_vir_adj %>%
  filter(wk_date %in% post_pandemic) %>%
  group_by(region) %>%
  dplyr::summarize(mean(age_5_24, na.rm = TRUE)) %>%
  filter(!(region == "Region 10"))

ili_pos_ratio_age_5_24 <- post_pan_age_5_24[, 2] / pre_pan_age_5_24[, 2]
ili_pos_ratio_age_5_24$region <- post_pan_age_5_24$region
colnames(ili_pos_ratio_age_5_24)[1] <- "ili_ratio_age_5_24"
ili_pos_ratio_age_5_24 # increase in some regions

pre_pan_age_25_64 <- regionflu_ili_vir_adj %>%
  filter(wk_date %in% pre_pandemic) %>%
  group_by(region) %>%
  dplyr::summarize(mean(age_25_64, na.rm = TRUE)) %>%
  filter(!(region == "Region 10"))

post_pan_age_25_64 <- regionflu_ili_vir_adj %>%
  filter(wk_date %in% post_pandemic) %>%
  group_by(region) %>%
  dplyr::summarize(mean(age_25_64, na.rm = TRUE)) %>%
  filter(!(region == "Region 10"))

ili_pos_ratio_age_25_64 <- post_pan_age_25_64[, 2] / pre_pan_age_25_64[, 2]
ili_pos_ratio_age_25_64$region <- post_pan_age_25_64$region
colnames(ili_pos_ratio_age_25_64)[1] <- "ili_ratio_age_25_64"
ili_pos_ratio_age_25_64 # increased esp in Region 2

pre_pan_age_65 <- regionflu_ili_vir_adj %>%
  filter(wk_date %in% pre_pandemic) %>%
  group_by(region) %>%
  dplyr::summarize(mean(age_65, na.rm = TRUE)) %>%
  filter(!(region == "Region 10"))

post_pan_age_65 <- regionflu_ili_vir_adj %>%
  filter(wk_date %in% post_pandemic) %>%
  group_by(region) %>%
  dplyr::summarize(mean(age_65, na.rm = TRUE)) %>%
  filter(!(region == "Region 10"))

ili_pos_ratio_age_65 <- post_pan_age_65[, 2] / pre_pan_age_65[, 2]
ili_pos_ratio_age_65$region <- post_pan_age_65$region
colnames(ili_pos_ratio_age_65)[1] <- "ili_ratio_age_65"
ili_pos_ratio_age_65 # increased in all regions esp region 2

age_dist_ratios <- plyr::join_all(list(ili_pos_ratio_age_0_4, ili_pos_ratio_age_5_24, ili_pos_ratio_age_25_64, ili_pos_ratio_age_65), by = "region", type = "left")

regionflu_ili_vir_adj <- left_join(regionflu_ili_vir_adj, age_dist_ratios, by = "region")

regionflu_ili_vir_adj <-
  regionflu_ili_vir_adj %>%
  mutate(
    age_0_4_adj = if_else(wk_date %in% pre_pandemic, age_0_4 * ili_ratio_age_0_4, age_0_4),
    age_5_24_adj = if_else(wk_date %in% pre_pandemic, age_5_24 * ili_ratio_age_5_24, age_5_24),
    age_25_64_adj = if_else(wk_date %in% pre_pandemic, age_25_64 * ili_ratio_age_25_64, age_25_64),
    age_65_adj = if_else(wk_date %in% pre_pandemic, age_65 * ili_ratio_age_65, age_65)
  )

# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 5",])+
#   geom_line(aes(x=wk_date, y = age_0_4_adj, color="% age_0_4_adj"))+
#   geom_line(aes(x=wk_date, y = age_5_24_adj , color="% age_5_24_adj "))+
#   geom_line(aes(x=wk_date, y = age_25_64_adj, color="% age_25_64_adj"))+
#   geom_line(aes(x=wk_date, y = age_65_adj, color="% age_65_adj"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 3",])+
#   geom_line(aes(x=wk_date, y = age_0_4_adj, color="% age_0_4_adj"))+
#   geom_line(aes(x=wk_date, y = age_5_24_adj , color="% age_5_24_adj "))+
#   geom_line(aes(x=wk_date, y = age_25_64_adj, color="% age_25_64_adj"))+
#   geom_line(aes(x=wk_date, y = age_65_adj, color="% age_65_adj"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 2",])+
#   geom_line(aes(x=wk_date, y = age_0_4_adj, color="% age_0_4_adj"))+
#   geom_line(aes(x=wk_date, y = age_5_24_adj , color="% age_5_24_adj "))+
#   geom_line(aes(x=wk_date, y = age_25_64_adj, color="% age_25_64_adj"))+
#   geom_line(aes(x=wk_date, y = age_65_adj, color="% age_65_adj"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")

regionflu_ili_vir_adj %>%
  dplyr::select(season_description, region, age_0_4, age_5_24, age_25_64, age_65, ilitotal) %>%
  group_by(season_description, region) %>%
  dplyr::summarize(sum(ilitotal))

age_dist <- regionflu_ili_vir_adj %>%
  dplyr::select(season_description, region, age_0_4_adj, age_5_24_adj, age_25_64_adj, age_65_adj) %>%
  group_by(season_description, region) %>%
  dplyr::summarize(
    age_0_4_total = sum(age_0_4_adj, na.rm = T),
    age_5_24_total = sum(age_5_24_adj, na.rm = T),
    age_25_64_total = sum(age_25_64_adj, na.rm = T),
    age_65_total = sum(age_65_adj, na.rm = T)
  ) %>%
  mutate(ilitotal_adj = rowSums(across(where(is.numeric)))) %>%
  mutate(
    age_0_4_prop = age_0_4_total / ilitotal_adj,
    age_5_24_prop = age_5_24_total / ilitotal_adj,
    age_25_64_prop = age_25_64_total / ilitotal_adj,
    age_65_prop = age_65_total / ilitotal_adj
  )
head(age_dist)
save(age_dist, file = "data/age_distribution_ili_cases.RData")
##########################################
## type/subtype distribution
##########################################

subtype_dist <- regionflu_ili_vir_adj %>%
  dplyr::select(season_description, region, prop_h1, prop_h3, prop_b, prop_a, total_specimens) %>%
  group_by(season_description, region) %>%
  mutate(
    h3_samples = prop_h3 * total_specimens,
    h1_samples = prop_h1 * total_specimens,
    a_samples = prop_a * total_specimens,
    b_samples = prop_b * total_specimens
  ) %>%
  dplyr::summarize(
    h3_total = sum(h3_samples, na.rm = T),
    h1_total = sum(h1_samples, na.rm = T),
    a_total = sum(a_samples, na.rm = T),
    b_total = sum(b_samples, na.rm = T),
    resp_samples = sum(total_specimens, na.rm = T)
  ) %>%
  rowwise() %>%
  mutate(flu_samples = sum(a_total + b_total, na.rm = T)) %>%
  mutate(
    h3_dom = h3_total / flu_samples,
    h1_dom = h1_total / flu_samples,
    b_dom = b_total / flu_samples,
    # h3_prop = h3_total/resp_samples,
    # h1_prop = h1_total/resp_samples,
    # b_prop = b_total/resp_samples,
    h3_vs_h1 = h3_total / a_total,
    iva_vs_ivb = a_total / flu_samples,
    h3_vs_flu_samples = h3_total / flu_samples,
    h1_vs_flu_samples = h1_total / flu_samples,
    ivb_vs_iva = b_total / flu_samples
  )
head(subtype_dist)

subtype_dist <- subtype_dist %>% tidyr::separate(col = "season_description", sep = "-", remove = F, into = c("year1", "year2"))

subtype_dist <- subtype_dist %>% filter(!(region == "Region 10" & year1 < 2009))

save(subtype_dist, file = "data/subtype_distribution_by_region_season.RData")

save(regionflu_ili_vir_adj, file = "data/hhs_division_level_ILI_and_virology_interp_smoothed_sampling_effort.RData")

# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 5",])+
#   geom_line(aes(x=wk_date, y = ILI_interp_H3, color="% H3N2"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 5",])+
#   geom_line(aes(x=wk_date, y = ILI_interp_H1, color="% H1N1"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 10",])+
#   geom_line(aes(x=wk_date, y = ILI_interp_H3, color="% H3N2"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 9",])+
#   geom_line(aes(x=wk_date, y = ili_pos, color="% Positive"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 8",])+
#   geom_line(aes(x=wk_date, y = ili_pos, color="% Positive"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")

##########################
## percent H3
##########################
burden_list <- list()
for (k in seasons) {
  data1 <- regionflu_ili_vir_adj[regionflu_ili_vir_adj$season_description == k, ]
  data1 <- data1[data1$week >= 40 | data1$week <= 20, ]
  data1 <- if (k == "2008-2009") data1[data1$week >= 40 | data1$week <= 15, ] else data1[data1$week >= 40 | data1$week <= 20, ]
  names(data1)[names(data1) %in% "season_description"] <- "season"


  pMiss <- function(x) {
    sum(is.na(x)) / length(x)
  } # proportion of data missing
  m <- tapply(data1$ili_h3_st, data1$region, pMiss)
  m
  # replace NAs with zeros to calculate sum
  data1$ILI_interp2 <- ifelse(is.na(data1$ili_h3_st), 0, data1$ili_h3_st)
  x <- tapply(data1$ILI_interp2, data1$region, sum)
  x

  y <- tapply(data1$ili_h3_st, data1$region, length)
  y

  keep <- intersect(names(which(!is.na(x) & x > 0 & m < 0.6)), names(which(y == max(y))))
  data1 <- data1[data1$region %in% keep, ]

  m <- tapply(data1$ili_h3_st, data1$region, pMiss)
  m
  y <- tapply(data1$ili_h3_st, data1$region, length)
  y
  x <- tapply(data1$ILI_interp2, data1$region, sum)
  l <- tapply(data1$ILI_interp2, data1$region, function(x) length(which(x > 0)))
  z <- tapply(data1$ILI_interp2, data1$region, max)
  z

  # determine date of peak epidemic intensity
  peak_date <- data1 %>%
    dplyr::group_by(region) %>%
    do(data.frame(peak_week = .$wk_date[which.max(.$ILI_interp2)]))

  # calculate shannon's entropy
  shannon <- data1[, c("region", "wk_date", "ILI_interp2")]
  shannon2 <- shannon %>% spread(wk_date, ILI_interp2) # make each week a column
  shannon2 <- as.data.frame(shannon2)
  rows <- shannon2[, "region"]
  rownames(shannon2) <- rows
  shannon2 <- as.matrix.data.frame(shannon2)
  shannon2 <- apply(shannon2[, -1], 2, as.numeric)
  sh_entropy <- 1 / (diversity(x = shannon2, index = "shannon", MARGIN = 1))

  burden_df <- data.frame(
    season = k, region = names(x), cum_intensity = x,
    season_duration = l, max_intensity = z,
    shannon_entropy = sh_entropy, missing_data = m,
    peak_week = peak_date$peak_week
  )

  ili_week <- as.data.frame(shannon2)
  ili_week <- ili_week[rowSums(ili_week[, -1]) > 0, ]

  burden_list[[length(burden_list) + 1]] <- burden_df
}
# dev.off()
burden_df_ili_vir_interp <- do.call(rbind.data.frame, burden_list)
head(burden_df_ili_vir_interp)
unique(burden_df_ili_vir_interp$season)
load("data/CDC_HHS_ILI_interp_H3_onset_weeks.RData") # onset_df_ili_vir_interp_H3

head(onset_df_ili_vir_interp_H3)
colnames(onset_df_ili_vir_interp_H3)[which(names(onset_df_ili_vir_interp_H3) == "wk_date")] <- "onset_week"

region_flu_metrics_H3 <- left_join(burden_df_ili_vir_interp,
  onset_df_ili_vir_interp_H3[, c("region", "season", "onset_week", "onsets", "CI(95%).l", "CI(95%).u", "se")],
  by = c("region", "season")
)

names(region_flu_metrics_H3)[names(region_flu_metrics_H3) %in% c("onsets", "CI(95%).l", "CI(95%).u", "se")] <-
  c("onset_index_week", "onset_lowerCI", "onset_upperCI", "onset_se")
head(region_flu_metrics_H3) # season-level region flu burden metrics
unique(region_flu_metrics_H3$season)
save(region_flu_metrics_H3, file = "data/region_level_flu_metrics_H3.RData")

##########################
## percent H1
##########################
burden_list <- list()
for (k in seasons) {
  data1 <- regionflu_ili_vir_adj[regionflu_ili_vir_adj$season_description == k, ]
  data1 <- data1[data1$week >= 40 | data1$week <= 20, ]
  data1 <- if (k == "2008-2009") data1[data1$week >= 40 | data1$week <= 15, ] else data1[data1$week >= 40 | data1$week <= 20, ]
  names(data1)[names(data1) %in% "season_description"] <- "season"


  pMiss <- function(x) {
    sum(is.na(x)) / length(x)
  } # proportion of data missing
  m <- tapply(data1$ili_h1_st, data1$region, pMiss)
  m
  # replace NAs with zeros to calculate sum
  data1$ILI_interp2 <- ifelse(is.na(data1$ili_h1_st), 0, data1$ili_h1_st)
  x <- tapply(data1$ILI_interp2, data1$region, sum)
  x

  y <- tapply(data1$ili_h1_st, data1$region, length)
  y

  keep <- intersect(names(which(!is.na(x) & x > 0 & m < 0.6)), names(which(y == max(y))))
  data1 <- data1[data1$region %in% keep, ]

  m <- tapply(data1$ili_h1_st, data1$region, pMiss)
  m
  y <- tapply(data1$ili_h1_st, data1$region, length)
  y
  x <- tapply(data1$ILI_interp2, data1$region, sum)
  l <- tapply(data1$ILI_interp2, data1$region, function(x) length(which(x > 0)))
  z <- tapply(data1$ILI_interp2, data1$region, max)
  z

  # determine date of peak epidemic intensity
  peak_date <- data1 %>%
    dplyr::group_by(region) %>%
    do(data.frame(peak_week = .$wk_date[which.max(.$ILI_interp2)]))

  # calculate shannon's entropy
  shannon <- data1[, c("region", "wk_date", "ILI_interp2")]
  shannon2 <- shannon %>% spread(wk_date, ILI_interp2) # make each week a column
  shannon2 <- as.data.frame(shannon2)
  rows <- shannon2[, "region"]
  rownames(shannon2) <- rows
  shannon2 <- as.matrix.data.frame(shannon2)
  shannon2 <- apply(shannon2[, -1], 2, as.numeric)
  sh_entropy <- 1 / (diversity(x = shannon2, index = "shannon", MARGIN = 1))

  burden_df <- data.frame(
    season = k, region = names(x), cum_intensity = x,
    season_duration = l, max_intensity = z,
    shannon_entropy = sh_entropy, missing_data = m,
    peak_week = peak_date$peak_week
  )

  ili_week <- as.data.frame(shannon2)
  ili_week <- ili_week[rowSums(ili_week[, -1]) > 0, ]

  burden_list[[length(burden_list) + 1]] <- burden_df
}
burden_df_ili_vir_interp <- do.call(rbind.data.frame, burden_list)
unique(burden_df_ili_vir_interp$season)
load("data/CDC_HHS_ILI_interp_H1_onset_weeks.RData") # onset_df_ili_vir_interp_H1

colnames(onset_df_ili_vir_interp_H1)[which(names(onset_df_ili_vir_interp_H1) == "wk_date")] <- "onset_week"

region_flu_metrics_H1 <- left_join(burden_df_ili_vir_interp,
  onset_df_ili_vir_interp_H1[, c("region", "season", "onset_week", "onsets", "CI(95%).l", "CI(95%).u", "se")],
  by = c("region", "season")
)
names(region_flu_metrics_H1)[names(region_flu_metrics_H1) %in% c("onsets", "CI(95%).l", "CI(95%).u", "se")] <-
  c("onset_index_week", "onset_lowerCI", "onset_upperCI", "onset_se")
head(region_flu_metrics_H1) # season-level region flu burden metrics
unique(region_flu_metrics_H1$season)
save(region_flu_metrics_H1, file = "data/region_level_flu_metrics_H1.RData")

##############################
## Influenza B
##############################
burden_list <- list()
for (k in seasons) {
  data1 <- regionflu_ili_vir_adj[regionflu_ili_vir_adj$season_description == k, ]
  data1 <- if (k == "2008-2009") data1[data1$week >= 40 | data1$week <= 13, ] else data1[data1$week >= 40 | data1$week <= 20, ]
  names(data1)[names(data1) %in% "season_description"] <- "season"


  pMiss <- function(x) {
    sum(is.na(x)) / length(x)
  } # proportion of data missing
  m <- tapply(data1$ili_ivb_st, data1$region, pMiss)
  m
  # replace NAs with zeros to calculate sum
  data1$ILI_interp2 <- ifelse(is.na(data1$ili_ivb_st), 0, data1$ili_ivb_st)
  x <- tapply(data1$ILI_interp2, data1$region, sum)
  x

  y <- tapply(data1$ili_ivb_st, data1$region, length)
  y

  keep <- intersect(names(which(!is.na(x) & x > 0 & m < 0.6)), names(which(y == max(y))))
  data1 <- data1[data1$region %in% keep, ]

  m <- tapply(data1$ili_ivb_st, data1$region, pMiss)
  m
  y <- tapply(data1$ili_ivb_st, data1$region, length)
  y
  x <- tapply(data1$ILI_interp2, data1$region, sum)
  l <- tapply(data1$ILI_interp2, data1$region, function(x) length(which(x > 0)))
  z <- tapply(data1$ILI_interp2, data1$region, max)
  z

  # determine date of peak epidemic intensity
  peak_date <- data1 %>%
    dplyr::group_by(region) %>%
    do(data.frame(peak_week = .$wk_date[which.max(.$ILI_interp2)]))

  # calculate shannon's entropy
  shannon <- data1[, c("region", "wk_date", "ILI_interp2")]
  shannon2 <- shannon %>% spread(wk_date, ILI_interp2) # make each week a column
  shannon2 <- as.data.frame(shannon2)
  rows <- shannon2[, "region"]
  rownames(shannon2) <- rows
  shannon2 <- as.matrix.data.frame(shannon2)
  shannon2 <- apply(shannon2[, -1], 2, as.numeric)
  sh_entropy <- 1 / (diversity(x = shannon2, index = "shannon", MARGIN = 1))

  burden_df <- data.frame(
    season = k, region = names(x), cum_intensity = x,
    season_duration = l, max_intensity = z,
    shannon_entropy = sh_entropy, missing_data = m,
    peak_week = peak_date$peak_week
  )

  ili_week <- as.data.frame(shannon2)
  ili_week <- ili_week[rowSums(ili_week[, -1]) > 0, ]

  burden_list[[length(burden_list) + 1]] <- burden_df
}
burden_df_ili_vir_interp <- do.call(rbind.data.frame, burden_list)

load("data/CDC_HHS_ILI_interp_B_onset_weeks.RData") # onset_df_ili_vir_interp_B
colnames(onset_df_ili_vir_interp_B)[which(names(onset_df_ili_vir_interp_B) == "wk_date")] <- "onset_week"

region_flu_metrics_IVB <- left_join(burden_df_ili_vir_interp,
  onset_df_ili_vir_interp_B[, c("region", "season", "onset_week", "onsets", "CI(95%).l", "CI(95%).u", "se")],
  by = c("region", "season")
)
names(region_flu_metrics_IVB)[names(region_flu_metrics_IVB) %in% c("onsets", "CI(95%).l", "CI(95%).u", "se")] <-
  c("onset_index_week", "onset_lowerCI", "onset_upperCI", "onset_se")

unique(region_flu_metrics_IVB$season)
save(region_flu_metrics_IVB, file = "data/region_level_flu_metrics_IVB.RData")
