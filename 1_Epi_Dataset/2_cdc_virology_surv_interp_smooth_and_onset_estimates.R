## load packages
list.of.packages <- c("cdcfluview", "dplyr", "ggplot2", "cowplot", "zoo", "segmented", "tidyr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

## load data
load("data/hhs_division_level_ILI_and_virology.RData") # regionflu_ili_vir

# imputation of missing percent positive values
head(regionflu_ili_vir)
names(regionflu_ili_vir)

## interpolate missing values (up to 4 weeks) and 8-week smoothing
regionflu_ili_vir_adj <- regionflu_ili_vir %>%
  mutate(
    year = mmwr_week(wk_date)$mmwr_year,
    week = mmwr_week(wk_date)$mmwr_week
  ) %>%
  arrange(region, wk_date) %>%
  group_by(region) %>%
  mutate(
    ILI_interp_pos = na.approx(adj_ili_per_pos, maxgap = 4, na.rm = FALSE), ## interpolate missing values for up to 4 week timespans
    ILI_interp_A = na.approx(adj_ili_iva, maxgap = 4, na.rm = FALSE),
    ILI_interp_H3 = na.approx(adj_ili_iva_h3, maxgap = 4, na.rm = FALSE),
    ILI_interp_H1 = na.approx(adj_ili_iva_h1, maxgap = 4, na.rm = FALSE),
    ILI_interp_B = na.approx(adj_ili_ivb, maxgap = 4, na.rm = FALSE)
  ) %>%
  mutate(
    ILI_smooth_pos = rollmean(ILI_interp_pos, k = 8, fill = NA), # 8-week rolling average; alignment is "centered" by default
    ILI_smooth_A = rollmean(ILI_interp_A, k = 8, fill = NA),
    ILI_smooth_H3 = rollmean(ILI_interp_H3, k = 8, fill = NA),
    ILI_smooth_H1 = rollmean(ILI_interp_H1, k = 8, fill = NA),
    ILI_smooth_B = rollmean(ILI_interp_B, k = 8, fill = NA)
  ) %>%
  mutate(
    ILI_smooth_pos = ifelse(is.na(ILI_smooth_pos), ILI_interp_pos, ILI_smooth_pos),
    ILI_smooth_A = ifelse(is.na(ILI_smooth_A), ILI_interp_A, ILI_smooth_A),
    ILI_smooth_H3 = ifelse(is.na(ILI_smooth_H3), ILI_interp_H3, ILI_smooth_H3),
    ILI_smooth_H1 = ifelse(is.na(ILI_smooth_H1), ILI_interp_H1, ILI_smooth_H1),
    ILI_smooth_B = ifelse(is.na(ILI_smooth_B), ILI_interp_B, ILI_smooth_B)
  ) %>%
  ungroup()

save(regionflu_ili_vir_adj, file = "data/hhs_division_level_ILI_and_virology_interp_smoothed.RData")

## exploratory plots
# ggplot(regionflu_ili_vir_adj %>% filter(season_description=="2016-2017"))+
#   geom_line(aes(x=wk_date, y = ILI_interp_H3, color=region),alpha=0.7)+
#   geom_line(aes(x=wk_date, y= ILI_interp_B, color=region),alpha=0.3)+
#   scale_x_date(date_breaks = "months",date_labels = "%m")
#
# ggplot(regionflu_ili_vir_adj %>% filter(season_description=="2017-2018"))+
#   geom_line(aes(x=wk_date, y = ILI_interp_H3, color=region),alpha=0.7)+
#   geom_line(aes(x=wk_date, y= ILI_interp_B, color=region),alpha=0.7)+
#   scale_x_date(date_breaks = "months",date_labels = "%m")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 3",])+
#   geom_line(aes(x=wk_date, y = ILI_interp_H3, color="H3"))+
#   geom_line(aes(x=wk_date, y= ILI_interp_H1, color = "H1"))+
#   geom_line(aes(x=wk_date, y= ILI_interp_B, color = "B"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 4",])+
#   geom_line(aes(x=wk_date, y = ILI_interp_H3, color="H3"))+
#   geom_line(aes(x=wk_date, y= ILI_interp_H1, color = "H1"))+
#   geom_line(aes(x=wk_date, y= ILI_interp_B, color = "B"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 5",])+
#   geom_line(aes(x=wk_date, y = ILI_interp_H1, color="ILI interpolated"))+
#   geom_line(aes(x=wk_date, y= ILI_smooth_H1, color = "ILI interp/smoothed"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 5",])+
#   geom_line(aes(x=wk_date, y = ILI_interp_H3, color="H3"))+
#   geom_line(aes(x=wk_date, y= ILI_interp_H1, color = "H1"))+
#   geom_line(aes(x=wk_date, y= ILI_interp_B, color = "B"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 9",])+
#   geom_line(aes(x=wk_date, y = ILI_interp_H3, color="ILI interpolated"))+
#   geom_line(aes(x=wk_date, y= ILI_smooth_H3, color = "ILI interp/smoothed"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 4",])+
#   geom_line(aes(x=wk_date, y = ILI_interp_B, color="ILI interpolated"))+
#   geom_line(aes(x=wk_date, y= ILI_smooth_B, color = "ILI interp/smoothed"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 3",])+
#   geom_line(aes(x=wk_date, y = ILI_interp_A, color="ILI interpolated"))+
#   geom_line(aes(x=wk_date, y= ILI_smooth_A, color = "ILI interp/smoothed"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 3",])+
#   geom_line(aes(x=wk_date, y = ILI_interp_H3, color="ILI interpolated"))+
#   geom_line(aes(x=wk_date, y= ILI_smooth_H3, color = "ILI interp/smoothed"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir_adj[regionflu_ili_vir_adj$region=="Region 1",])+
#   geom_line(aes(x=wk_date, y = ILI_interp_pos, color="ILI interpolated"))+
#   geom_line(aes(x=wk_date, y= ILI_smooth_pos, color = "ILI interp/smoothed"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")

##### ##### ##### ##### ##### ##### ##### #####
##### H3 onsets
##### ##### ##### ##### ##### ##### ##### #####
seasons <- unique(regionflu_ili_vir_adj$season_description)
seasons
seasons <- seasons[!(seasons %in% c("2019-2020", "2020-2021", "2021-2022"))] # keep H1N1 pdm for now
length(seasons)

drop <- c("2000-2001", "2019-2020", "2020-2021", "2021-2022")
seasons_h3 <- seasons[!(seasons %in% drop)]
onset_list <- list()
for (k in seasons_h3) {
  data1 <- regionflu_ili_vir_adj[regionflu_ili_vir_adj$season_description == k, ]
  names(data1)[names(data1) %in% "season_description"] <- "season"

  data1 <- data1 %>%
    filter(case_when(
      season %in% c("2015-2016", "2013-2014") ~ week >= 40 | week <= 10,
      !(season %in% c("2015-2016", "2013-2014")) ~ week >= 40 | week <= 4
    ))

  pMiss <- function(x) {
    sum(is.na(x)) / length(x)
  } # proportion of data missing
  m <- tapply(data1$ILI_interp_H3, data1$region, pMiss)
  m

  data1$ILI_interp2 <- ifelse(is.na(data1$ILI_interp_H3), 0, data1$ILI_interp_H3)
  l <- tapply(data1$ILI_interp2, data1$region, function(x) length(which(x > 0)))
  l
  x <- tapply(data1$ILI_interp2, data1$region, sum)
  x
  y <- tapply(data1$ILI_interp2, data1$region, length)
  y
  z <- tapply(data1$ILI_interp2, data1$region, max)
  z

  len_non_zero <- function(x) {
    max(with(rle(x == 0), lengths[!values]))
  } ## max run length of non-zeros
  nz <- tapply(data1$ILI_interp2, data1$region, len_non_zero)
  nz

  # selected locations with a sufficient rise in ILI above baseline
  keep <- intersect(names(which(!is.na(x) & x > 0.02 & nz != Inf & nz >= 3)), names(which(y == max(y))))
  keep
  str(keep)
  if (identical(keep, character(0)) == TRUE) next
  data1 <- data1[data1$region %in% keep, ] %>% droplevels()
  data1 <- data1[order(data1$region, data1$wk_date), ]

  data1$ILI_smooth2 <- ifelse(is.na(data1$ILI_smooth_H3), 0, data1$ILI_smooth_H3)

  df <- data.frame(week = unique(data1$week), onset_index = seq(1:max(y)))
  data1$index_week <- df[match(data1$week, df$week), 2] ## add in season index week to original dataframe

  fluhosp <- matrix(0, ncol = (length(unique(data1$region)) + 2), nrow = max(y), dimnames = list(c(1:max(y)), c("weeks", "wk_date", unique(data1$region))))

  fluhosp[, 1] <- c(1:max(y))
  fluhosp[, 2] <- unique(data1$wk_date)

  for (i in unique(data1$region)) {
    fluhosp[, as.character(i)] <- data1$ILI_smooth2[data1$region == i]
  }

  trimy <- function(y) {
    stop <- which(y == max(y))
    y <- y[1:stop]
    return(y)
  }

  onsetfun <- function(y) {
    y <- trimy(y)
    x <- c(1:length(y))
    out.lm <- lm(y ~ x)
    o <- segmented(out.lm, seg.Z = ~x)
    onset <- confint(o)
    onset[4] <- o$psi[3]
    return(onset[1:4])
  }

  x1 <- apply(fluhosp[, 3:dim(fluhosp)[2]], 2, onsetfun)
  # onsetfun(fluhosp[,6])

  onsets <- as.data.frame(t(x1))
  fluhosp2 <- fluhosp[, -c(1, 2)]

  onsets$region <- rownames(onsets)
  names(onsets)[1] <- "onsets"
  names(onsets)[2] <- "CI(95%).l"
  names(onsets)[3] <- "CI(95%).u"
  names(onsets)[4] <- "se"
  onsets$season <- k
  onsets$round1 <- round(onsets$onsets, 1)
  onsets$index_week <- round(onsets$onsets)
  onset2 <-
    onsets %>%
    left_join(data1[, c("year", "week", "wk_date", "season", "index_week")],
      by = c("season", "index_week"),
      relationship = "many-to-many"
    )
  onset3 <- onset2[!duplicated(onset2), ]
  onset_list[[length(onset_list) + 1]] <- onset3
}
onset_df_ili_vir_interp_H3 <- do.call(rbind.data.frame, onset_list)
head(onset_df_ili_vir_interp_H3)
onset_df_ili_vir_interp_H3 %>%
  group_by(season) %>%
  tally()
save(onset_df_ili_vir_interp_H3, file = "data/CDC_HHS_ILI_interp_H3_onset_weeks.RData")

##### ##### ##### ##### ##### ##### ##### #####
##### H1 onsets
##### ##### ##### ##### ##### ##### ##### #####

seasons_H1 <- seasons[!(seasons %in% c(
  "1997-1998", "1998-1999", "2001-2002",
  "2003-2004", "2004-2005", "2014-2015",
  "2019-2020", "2020-2021", "2021-2022"
))]

onset_list <- list()
for (k in seasons_H1) {
  data1 <- regionflu_ili_vir_adj[regionflu_ili_vir_adj$season_description == k, ]
  names(data1)[names(data1) %in% "season_description"] <- "season"

  data1 <- data1 %>%
    filter(case_when(
      season %in% c("2011-2012", "2005-2006") ~ week >= 40 | week <= 20,
      !(season %in% c("2011-2012", "2005-2006")) ~ week >= 30 | week <= 10
    ))

  pMiss <- function(x) {
    sum(is.na(x)) / length(x)
  } # proportion of data missing
  m <- tapply(data1$ILI_interp_H1, data1$region, pMiss)
  m

  data1$ILI_interp2 <- ifelse(is.na(data1$ILI_interp_H1), 0, data1$ILI_interp_H1)
  l <- tapply(data1$ILI_interp2, data1$region, function(x) length(which(x > 0)))
  l
  x <- tapply(data1$ILI_interp2, data1$region, sum)
  x
  y <- tapply(data1$ILI_interp2, data1$region, length)
  y
  z <- tapply(data1$ILI_interp2, data1$region, max)
  z

  len_non_zero <- function(x) {
    max(with(rle(x == 0), lengths[!values]))
  }
  nz <- tapply(data1$ILI_interp2, data1$region, len_non_zero)
  nz

  # selected locations with a sufficient rise in ILI above baseline
  keep <- intersect(names(which(!is.na(x) & x > 0.02 & nz != Inf & nz >= 3)), names(which(y == max(y))))
  # keep
  data1 <- data1[data1$region %in% keep, ] %>% droplevels()
  data1 <- data1[order(data1$region, data1$wk_date), ]

  data1$ILI_smooth2 <- ifelse(is.na(data1$ILI_smooth_H1), 0, data1$ILI_smooth_H1)

  df <- data.frame(week = unique(data1$week), onset_index = seq(1:max(y)))
  data1$index_week <- df[match(data1$week, df$week), 2] ## add in season index week to original dataframe

  fluhosp <- matrix(0, ncol = (length(unique(data1$region)) + 2), nrow = max(y), dimnames = list(c(1:max(y)), c("weeks", "wk_date", unique(data1$region))))

  fluhosp[, 1] <- c(1:max(y))
  fluhosp[, 2] <- unique(data1$wk_date)

  for (i in unique(data1$region)) {
    fluhosp[, as.character(i)] <- data1$ILI_smooth2[data1$region == i]
  }

  trimy <- function(y) {
    stop <- which(y == max(y))
    y <- y[1:stop]
    return(y)
  }

  onsetfun <- function(y) {
    y <- trimy(y)
    x <- c(1:length(y))
    out.lm <- lm(y ~ x)
    o <- segmented(out.lm, seg.Z = ~x)
    onset <- confint(o)
    onset[4] <- o$psi[3]
    return(onset[1:4])
  }

  x1 <- apply(fluhosp[, 3:dim(fluhosp)[2]], 2, onsetfun)

  onsets <- as.data.frame(t(x1))
  fluhosp2 <- fluhosp[, -c(1, 2)]

  onsets$region <- rownames(onsets)
  names(onsets)[1] <- "onsets"
  names(onsets)[2] <- "CI(95%).l"
  names(onsets)[3] <- "CI(95%).u"
  names(onsets)[4] <- "se"
  onsets$season <- k
  onsets$round1 <- round(onsets$onsets, 1)
  onsets$index_week <- round(onsets$onsets)

  onset2 <-
    onsets %>%
    left_join(data1[, c("year", "week", "wk_date", "season", "index_week")],
      by = c("season", "index_week"),
      relationship = "many-to-many"
    )
  onset3 <- onset2[!duplicated(onset2), ]
  onset_list[[length(onset_list) + 1]] <- onset3
}
onset_df_ili_vir_interp_H1 <- do.call(rbind.data.frame, onset_list)
head(onset_df_ili_vir_interp_H1)
onset_df_ili_vir_interp_H1 %>%
  group_by(season) %>%
  tally()
save(onset_df_ili_vir_interp_H1, file = "data/CDC_HHS_ILI_interp_H1_onset_weeks.RData")

##### ##### ##### ##### ##### ##### ##### #####
##### B onsets
##### ##### ##### ##### ##### ##### ##### #####
drop <- c("1997-1998", "1999-2000", "2003-2004", "2009-2010", "2019-2020", "2020-2021", "2021-2022")
seasons_B <- seasons[!(seasons %in% drop)]
onset_list <- list()
for (k in seasons_B) {
  data1 <- regionflu_ili_vir_adj[regionflu_ili_vir_adj$season_description == k, ]
  names(data1)[names(data1) %in% "season_description"] <- "season"
  data1 <- data1[data1$week >= 35 | data1$week <= 10, ]
  data1 <- data1 %>% filter(!(wk_date < as.Date("2010-01-01") & region == "Region 10"))

  pMiss <- function(x) {
    sum(is.na(x)) / length(x)
  } # proportion of data missing
  m <- tapply(data1$ILI_interp_B, data1$region, pMiss)
  m

  data1$ILI_interp2 <- ifelse(is.na(data1$ILI_interp_B), 0, data1$ILI_interp_B)
  l <- tapply(data1$ILI_interp2, data1$region, function(x) length(which(x > 0)))
  l
  x <- tapply(data1$ILI_interp2, data1$region, sum)
  x
  y <- tapply(data1$ILI_interp2, data1$region, length)
  y
  z <- tapply(data1$ILI_interp2, data1$region, max)
  z

  len_non_zero <- function(x) {
    max(with(rle(x == 0), lengths[!values]))
  }
  nz <- tapply(data1$ILI_interp2, data1$region, len_non_zero)
  nz


  # selected locations with a sufficient rise in ILI above baseline
  keep <- intersect(names(which(!is.na(x) & x > 0.02 & nz >= 4)), names(which(y == max(y))))
  # keep
  data1 <- data1[data1$region %in% keep, ]
  data1 <- data1[order(data1$region, data1$wk_date), ]

  data1$ILI_smooth2 <- ifelse(is.na(data1$ILI_smooth_B), 0, data1$ILI_smooth_B)

  df <- data.frame(week = unique(data1$week), onset_index = seq(1:max(y)))
  data1$index_week <- df[match(data1$week, df$week), 2] ## add in season index week to original dataframe

  fluhosp <- matrix(0, ncol = (length(unique(data1$region)) + 2), nrow = max(y), dimnames = list(c(1:max(y)), c("weeks", "wk_date", unique(data1$region))))

  fluhosp[, 1] <- c(1:max(y))
  fluhosp[, 2] <- unique(data1$wk_date)

  for (i in unique(data1$region)) {
    fluhosp[, as.character(i)] <- data1$ILI_smooth2[data1$region == i]
  }


  trimy <- function(y) {
    stop <- which(y == max(y))
    y <- y[1:stop]
    return(y)
  }

  onsetfun <- function(y) {
    y <- trimy(y)
    x <- c(1:length(y))
    out.lm <- lm(y ~ x)
    o <- segmented(out.lm, seg.Z = ~x)
    onset <- confint(o)
    onset[4] <- o$psi[3]
    return(onset[1:4])
  }

  x1 <- apply(fluhosp[, 3:dim(fluhosp)[2]], 2, onsetfun)

  onsets <- as.data.frame(t(x1))
  fluhosp2 <- fluhosp[, -c(1, 2)]

  onsets$region <- rownames(onsets)
  names(onsets)[1] <- "onsets"
  names(onsets)[2] <- "CI(95%).l"
  names(onsets)[3] <- "CI(95%).u"
  names(onsets)[4] <- "se"
  onsets$season <- k
  onsets$round1 <- round(onsets$onsets, 1)
  onsets$index_week <- round(onsets$onsets)

  onset2 <-
    onsets %>%
    left_join(data1[, c("year", "week", "wk_date", "season", "index_week")],
      by = c("season", "index_week"),
      relationship = "many-to-many"
    )
  onset3 <- onset2[!duplicated(onset2), ]
  onset_list[[length(onset_list) + 1]] <- onset3
}
onset_df_ili_vir_interp_B <- do.call(rbind.data.frame, onset_list)
head(onset_df_ili_vir_interp_B)
onset_df_ili_vir_interp_B %>%
  group_by(season) %>%
  tally()
save(onset_df_ili_vir_interp_B, file = "data/CDC_HHS_ILI_interp_B_onset_weeks.RData")
