## load packages
list.of.packages <- c("cdcfluview", "dplyr", "ggplot2", "cowplot")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

######################################################################################
## Download and process weekly CDC ILI and viral surveillance data
######################################################################################
## ILI data
regionflu_ili <- ilinet(region = "hhs", years = 1997:2020) # weekly ili from 1997
head(regionflu_ili)

## viral data
censusflu_sur <- who_nrevss(region = "hhs", years = 1997:2020)
head(censusflu_sur)

### virus surveillance data prior to 2015/16 (clinical and public health data combined)
combined <- censusflu_sur$combined_prior_to_2015_16
combined$yrweek <- paste(combined$year, sprintf("%02d", combined$week), sep = "-")
colnames(combined)
head(combined)
tail(combined)
## convert virus counts to numeric
cols <- c("total_specimens", "percent_positive", "a_2009_h1n1", "a_h1", "a_h3", "a_subtyping_not_performed", "a_unable_to_subtype", "b", "h3n2v")
combined[, cols] <- apply(combined[, cols], 2, function(x) as.numeric(as.character(x)))
names(combined)

## create %a and %b columns in combined public health/clinical lab prior to 2015 data set
combined <- combined %>%
  mutate(
    total_a = rowSums(.[, c("a_2009_h1n1", "a_h1", "a_h3", "a_subtyping_not_performed", "a_unable_to_subtype")], na.rm = TRUE),
    total_b = b,
    untyped_a = rowSums(.[, c("a_subtyping_not_performed", "a_unable_to_subtype")], na.rm = TRUE),
    typed_a = rowSums(.[, c("a_2009_h1n1", "a_h1", "a_h3")], na.rm = TRUE),
    typed_h1 = rowSums(.[, c("a_2009_h1n1", "a_h1")], na.rm = TRUE),
    typed_h3 = rowSums(.[, c("a_h3")], na.rm = TRUE)
  ) %>%
  mutate(
    tested_h3_prop = typed_h3 / typed_a,
    tested_h1_prop = typed_h1 / typed_a
  ) %>%
  mutate(
    inferred_h3 = untyped_a * tested_h3_prop,
    inferred_h1 = untyped_a * tested_h1_prop
  ) %>%
  mutate(
    adjusted_h3 = inferred_h3 + typed_h3,
    adjusted_h1 = inferred_h1 + typed_h1
  ) %>%
  mutate(
    prop_a = total_a / total_specimens,
    prop_b = total_b / total_specimens,
    prop_h1 = adjusted_h1 / total_specimens,
    prop_h3 = adjusted_h3 / total_specimens
  )

combined_season <-
  combined %>%
  mutate(cdc_date = mmwr_week_to_date(year, week, day = NULL)) %>%
  mutate(
    month.day = format(as.Date(cdc_date), "%m-%d"),
    new.year = as.numeric(format(as.Date(cdc_date), "%Y"))
  ) %>%
  mutate(season_description = ifelse(month.day < "07-01",
    sprintf("%d-%d", new.year - 1, new.year),
    sprintf("%d-%d", new.year, new.year + 1)
  ))

combined_season %>%
  dplyr::select(region, year, week, cdc_date, wk_date, season_description) %>%
  filter(cdc_date != wk_date)

unique(combined_season$cdc_date - combined_season$wk_date) # sanity check

###### clinical lab data post 2015/2016
cl_labs <- censusflu_sur$clinical_labs
cols <- c("total_specimens", "total_a", "total_b", "percent_positive", "percent_a", "percent_b")
cl_labs[, cols] <- apply(cl_labs[, cols], 2, function(x) as.numeric(as.character(x)))
cl_labs$yrweek <- paste(cl_labs$year, sprintf("%02d", cl_labs$week), sep = "-")
nrow(cl_labs) # 3130

## create season indicator for clinical lab data
cl_labs_season <-
  cl_labs %>%
  mutate(cdc_date = mmwr_week_to_date(year, week, day = NULL)) %>%
  mutate(
    month.day = format(as.Date(cdc_date), "%m-%d"),
    new.year = as.numeric(format(as.Date(cdc_date), "%Y"))
  ) %>%
  mutate(season_description = ifelse(month.day < "07-01",
    sprintf("%d-%d", new.year - 1, new.year),
    sprintf("%d-%d", new.year, new.year + 1)
  ))


## public health lab data (samples have already tested positive by clinical labs)
ph_labs <- censusflu_sur$public_health_labs
nrow(ph_labs) # 3130
names(ph_labs)
cols <- c("total_specimens", "a_2009_h1n1", "a_h3", "a_subtyping_not_performed", "b", "bvic", "byam")
ph_labs[, cols] <- apply(ph_labs[, cols], 2, function(x) as.numeric(as.character(x)))

ph_labs <- ph_labs %>%
  mutate(
    total_a = rowSums(.[, c("a_2009_h1n1", "a_h3", "a_subtyping_not_performed")], na.rm = TRUE),
    total_b = rowSums(.[, c("b", "byam", "bvic")]),
    untyped_a = a_subtyping_not_performed,
    typed_a = rowSums(.[, c("a_2009_h1n1", "a_h3")], na.rm = TRUE),
    typed_h1 = a_2009_h1n1,
    typed_h3 = a_h3
  ) %>%
  mutate(
    tested_h3_prop = typed_h3 / typed_a,
    tested_h1_prop = typed_h1 / typed_a
  ) %>%
  mutate(
    inferred_h3 = untyped_a * tested_h3_prop,
    inferred_h1 = untyped_a * tested_h1_prop
  ) %>%
  mutate(
    adjusted_h3 = inferred_h3 + typed_h3,
    adjusted_h1 = inferred_h1 + typed_h1
  ) %>%
  mutate(
    adjusted_h3_prop = adjusted_h3 / total_a,
    adjusted_h1_prop = adjusted_h1 / total_a
  )

ph_labs %>% dplyr::select(typed_a, untyped_a, total_a, adjusted_h3, adjusted_h1, total_b, total_specimens)

ph_labs_season <-
  ph_labs %>%
  mutate(cdc_date = mmwr_week_to_date(year, week, day = NULL)) %>%
  mutate(
    month.day = format(as.Date(cdc_date), "%m-%d"),
    new.year = as.numeric(format(as.Date(cdc_date), "%Y"))
  ) %>%
  mutate(season_description = ifelse(month.day < "07-01",
    sprintf("%d-%d", new.year - 1, new.year),
    sprintf("%d-%d", new.year, new.year + 1)
  ))

ph_cols <- c("region", "year", "week", "cdc_date", "wk_date", "adjusted_h3_prop", "adjusted_h1_prop", "adjusted_h3", "adjusted_h1", "month.day", "season_description")
names(cl_labs_season)

combined_cl_ph_region <-
  left_join(cl_labs_season, ph_labs_season[, ph_cols], by = c("season_description", "region", "year", "week", "wk_date", "cdc_date", "month.day")) %>%
  mutate(
    percent_h1 = adjusted_h1_prop * percent_a,
    percent_h3 = adjusted_h3_prop * percent_a,
    prop_h1 = adjusted_h1_prop * (percent_a / 100),
    prop_h3 = adjusted_h3_prop * (percent_a / 100),
    prop_a = percent_a / 100,
    prop_b = percent_b / 100
  )

head(combined_cl_ph_region)
range(combined_cl_ph_region$percent_h1, na.rm = T) # 0.00000 24.20629
range(combined_cl_ph_region$percent_h3, na.rm = T) # 0.00000 31.85716

#### add combined_cl_ph (post 2015-2016) to combined df (pre 2015-2016)
colnames(combined_cl_ph_region)
colnames(combined_season) # pre 2015-2016 data; doesn't include IVB subtypes

colnames(combined_cl_ph_region) %in% colnames(combined_season)
remove <- names(combined_cl_ph_region)[!(colnames(combined_cl_ph_region) %in% colnames(combined_season))]
combined_cl_ph_reduced <- combined_cl_ph_region[!(names(combined_cl_ph_region) %in% remove)]
colnames(combined_cl_ph_reduced) %in% colnames(combined_season)
colnames(combined_season) %in% colnames(combined_cl_ph_reduced)
names(combined_season)
drops <- names(combined_season)[!(colnames(combined_season) %in% colnames(combined_cl_ph_reduced))]
combined_season_reduced <- combined_season[!(names(combined_season) %in% drops)]
colnames(combined_season_reduced) %in% colnames(combined_cl_ph_reduced)
colnames(combined_cl_ph_reduced) %in% colnames(combined_season_reduced)
# make sure order of columns match
combined_cl_ph_reduced <- combined_cl_ph_reduced[, colnames(combined_season_reduced)]
regionflu_vir <- rbind(combined_season_reduced, combined_cl_ph_reduced) # combine pre- and post-2015 datasets
head(regionflu_vir)
regionflu_vir <- regionflu_vir %>% arrange(wk_date, region)

## exploratory plots
# ggplot(regionflu_vir[regionflu_vir$region=="Region 5",])+
#   geom_line(aes(x=wk_date, y = prop_h3, color="% H3"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_vir[regionflu_vir$region=="Region 5",])+
#   geom_line(aes(x=wk_date, y = adjusted_h1, color="H1 cases"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_vir[regionflu_vir$region=="Region 5",])+
#   geom_line(aes(x=wk_date, y = prop_b, color="% IVB"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_vir[regionflu_vir$region=="Region 5",])+
#   geom_line(aes(x=wk_date, y = total_b, color="% IVB"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_vir[regionflu_vir$region=="Region 5",])+
#   geom_line(aes(x=wk_date, y = prop_h3, color="% H3"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_vir[regionflu_vir$region=="Region 2",])+
#   geom_line(aes(x=wk_date, y = prop_h3, color="% H3"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_vir[regionflu_vir$region=="Region 9",])+
#   geom_line(aes(x=wk_date, y = prop_h1, color="% H1"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_vir[regionflu_vir$region=="Region 9",])+
#   geom_line(aes(x=wk_date, y = prop_h3, color="% H3"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")

## ILI data
# ggplot(regionflu_ili[regionflu_ili$region=="Region 5",])+
#   geom_line(aes(x=week_start, y = weighted_ili, color="ILI"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")

regionflu_ili <-
  regionflu_ili %>%
  mutate(cdc_date = mmwr_week_to_date(year, week, day = NULL)) %>%
  mutate(
    month.day = format(as.Date(cdc_date), "%m-%d"),
    new.year = as.numeric(format(as.Date(cdc_date), "%Y"))
  ) %>%
  mutate(season_description = ifelse(month.day < "07-01",
    sprintf("%d-%d", new.year - 1, new.year),
    sprintf("%d-%d", new.year, new.year + 1)
  ))


# check
regionflu_ili %>%
  dplyr::select(region, year, week, week_start, cdc_date, season_description) %>%
  filter(week_start != cdc_date)

include <- c(
  "region", "unweighted_ili", "weighted_ili",
  "year", "week", "week_start", "cdc_date",
  "age_0_4", "age_25_49", "age_25_64", "age_5_24", "age_50_64", "age_65",
  "ilitotal", "total_patients", "season_description"
)

regionflu_ili_vir <-
  left_join(regionflu_vir, regionflu_ili[, include], by = c("region", "year", "week", "cdc_date", "season_description")) %>%
  mutate(
    adjusted_h3 = ifelse(total_a == 0, 0, adjusted_h3),
    adjusted_h1 = ifelse(total_a == 0, 0, adjusted_h1),
    prop_h3 = ifelse(total_a == 0, 0, prop_h3),
    prop_h1 = ifelse(total_a == 0, 0, prop_h1)
  ) %>%
  mutate(
    adj_ili_iva = weighted_ili * prop_a,
    adj_ili_ivb = weighted_ili * prop_b,
    adj_ili_iva_h1 = weighted_ili * prop_h1,
    adj_ili_iva_h3 = weighted_ili * prop_h3,
    adj_ili_per_pos = weighted_ili * (percent_positive / 100)
  )

## sanity check
regionflu_ili_vir %>%
  dplyr::select(week_start, wk_date, cdc_date, year, week) %>%
  mutate(date_diff = week_start - cdc_date) %>%
  filter(date_diff > 0 | date_diff < 0)

regionflu_ili_vir %>%
  dplyr::select(week_start, wk_date, cdc_date, year, week) %>%
  mutate(date_diff = week_start - wk_date) %>%
  filter(date_diff > 0 | date_diff < 0)

regionflu_ili_vir %>%
  dplyr::select(week_start, wk_date, cdc_date, year, week) %>%
  mutate(date_diff = cdc_date - wk_date) %>%
  filter(date_diff > 0 | date_diff < 0)
unique(regionflu_ili$week_start - regionflu_vir$wk_date) # no difference in dates

## exploratory plots
# ggplot(regionflu_ili_vir[regionflu_ili_vir$region=="Region 8",])+
#   geom_line(aes(x=wk_date, y = adj_ili_iva_h3, color="% H3N2"))+
#   geom_line(aes(x=wk_date, y = adj_ili_iva_h1, color="% H1N1"))+
#   geom_line(aes(x=wk_date, y = adj_ili_ivb, color="% B"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir[regionflu_ili_vir$region=="Region 8",])+
#   geom_line(aes(x=wk_date, y = adj_ili_iva, color="% IVA"))+
#   geom_line(aes(x=wk_date, y = adj_ili_ivb, color="% IVB"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir[regionflu_ili_vir$region=="Region 5",])+
#   geom_line(aes(x=wk_date, y = adj_ili_iva, color="% IVA"))+
#   geom_line(aes(x=wk_date, y = adj_ili_ivb, color="% IVB"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir[regionflu_ili_vir$region=="Region 3",])+
#   geom_line(aes(x=wk_date, y = adj_ili_iva, color="% IVA"))+
#   geom_line(aes(x=wk_date, y = adj_ili_ivb, color="% IVB"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
#
# ggplot(regionflu_ili_vir[regionflu_ili_vir$region=="Region 2",])+
#   geom_line(aes(x=wk_date, y = weighted_ili, color="ILI"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir[regionflu_ili_vir$region=="Region 8",])+
#   geom_line(aes(x=wk_date, y = adj_ili_iva_h3, color="% H3N2"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir[regionflu_ili_vir$region=="Region 5",])+
#   geom_line(aes(x=wk_date, y = adj_ili_iva_h3, color="% H3N2"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# regionflu_ili_vir %>%
#   filter(!is.na(age_25_49) & ilitotal > 0) %>%
#   dplyr::select(age_0_4,age_5_24,age_25_49,age_50_64,age_25_64,age_65, ilitotal)

# prop of cases in each age group
regionflu_ili_vir <- regionflu_ili_vir %>%
  mutate(age_25_64 = ifelse(!is.na(age_25_49) & !is.na(age_50_64), age_25_49 + age_50_64, age_25_64)) %>%
  mutate(
    age_0_4_prop = age_0_4 / ilitotal,
    age_5_24_prop = age_5_24 / ilitotal,
    age_25_64_prop = age_25_64 / ilitotal,
    age_65_prop = age_65 / ilitotal
  )

regionflu_ili_vir %>%
  filter(!is.na(age_25_49) & ilitotal > 0) %>%
  dplyr::select(age_0_4, age_5_24, age_25_49, age_50_64, age_25_64, age_65, ilitotal)

range(regionflu_ili_vir$wk_date)
## exploratory plots
# ggplot(regionflu_ili_vir[regionflu_ili_vir$region=="Region 8",])+
#   geom_line(aes(x=wk_date, y = age_0_4_prop, color="% age_0_4"))+
#   geom_line(aes(x=wk_date, y = age_5_24_prop, color="% age_5_24_prop"))+
#   geom_line(aes(x=wk_date, y = age_25_64_prop, color="% age_25_64_prop"))+
#   geom_line(aes(x=wk_date, y = age_65_prop, color="% age_65_prop"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
#
# ggplot(regionflu_ili_vir[regionflu_ili_vir$region=="Region 3",])+
#   geom_line(aes(x=wk_date, y = age_0_4_prop, color="% age_0_4"))+
#   geom_line(aes(x=wk_date, y = age_5_24_prop, color="% age_5_24_prop"))+
#   geom_line(aes(x=wk_date, y = age_25_64_prop, color="% age_25_64_prop"))+
#   geom_line(aes(x=wk_date, y = age_65_prop, color="% age_65_prop"))+
#   scale_x_date(date_breaks = "years",date_labels = "%Y")
save(regionflu_ili_vir, file = "data/hhs_division_level_ILI_and_virology.RData")
