## load packages
list.of.packages <- c("dplyr", "tidyr", "heatmaply", "caret", "readr",
                     "ggplot2", "cowplot", "viridis", "RColorBrewer", "corrplot","rstatix"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

###########################################################################
## compile dataset for variable selection analyses
###########################################################################
# load data
load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")

epi_red <- epi_red %>% replace_na(list(
  H3_max_Rt = 0, 
  H3_max_intensity = 0,
  H3_epi_size_prior = 0,
  H3_cum_intensity = 0,
  H3_season_duration = 0
))

epi_red <- epi_red %>%
  mutate(
    vac_combined = (adult_18_49_vac_cov / 100) * (adult_65_vac_cov / 100) * weighted_VE,
    vac_cov_combined = (adult_18_49_vac_cov / 100) * (adult_65_vac_cov / 100)
  )

vac_prior <- epi_red %>%
  distinct(season, adult_18_49_vac_cov, adult_65_vac_cov, vac_combined, vac_cov_combined, weighted_VE) %>%
  arrange(season) %>%
  mutate(
    adult_18_49_vac_cov_prior = lag(adult_18_49_vac_cov, n = 1),
    adult_65_vac_cov_prior = lag(adult_65_vac_cov, n = 1),
    vac_combined_prior = lag(vac_combined, n = 1),
    vac_cov_combined_prior = lag(vac_cov_combined, n = 1),
    weighted_VE_prior_season = lag(weighted_VE, n = 1)
  ) %>%
  ungroup()
vac_prior

epi_red <- left_join(epi_red %>% dplyr::select(-contains(c("vac", "VE"))), vac_prior, by = "season")

epi_red2 <- epi_red %>%
  tidyr::replace_na(list(
    IVB_max_intensity = 0, H1_max_intensity = 0,
    IVB_cum_intensity = 0, H1_cum_intensity = 0
  )) %>%
  mutate(
    IVB_epi_size_prior = ifelse(year.new > 1997 & is.na(IVB_epi_size_prior), 0, IVB_epi_size_prior),
    H1_epi_size_prior = ifelse(year.new > 1997 & is.na(H1_epi_size_prior), 0, H1_epi_size_prior)
  ) %>%
  mutate(prior_dom_type_national = case_when(
    prior_dom_type_national == "H3" ~ 1,
    prior_dom_type_national == "H1" ~ 0,
    prior_dom_type_national == "co-circ" ~ 0.5
  )) %>%
  distinct()


predictor_df <- epi_red2 %>%
  filter(!(season %in% c("2009-2010", "1995-1996", "1996-1997"))) %>%
  dplyr::select(H1_cum_intensity,IVB_cum_intensity,prior_dom_type_national,
                H3_epi_size_prior,H1_epi_size_prior,IVB_epi_size_prior,
                contains(c("usa","lag1","lag2","lbi","vac","VE")))%>%
  dplyr::select(-contains(c("adult","krammer","koel_lag1","stem_ep_lag1","tree_lag1","wolf_lag1","wolf_nonepitope_lag1",
                            "bhatt_ep_lag2","bhatt_nonepitope_lag2")))

  # dplyr::select(-dom_type, -contains(c(
  #   "max_intensity", "vs", "total", "year", "wolf_lag1", "tree_lag1","koel_lag1",
  #   "stem_ep_lag1","bhatt_nonepitope_lag2","bhatt_ep_lag2","wolf_nonepitope_lag1","krammer_ep_lag2"
  # )))
predictor_df <- predictor_df[complete.cases(predictor_df), ]
head(predictor_df)
colnames(predictor_df)
nrow(predictor_df)


## check for zero variance
nzv <- nearZeroVar(predictor_df, saveMetrics = TRUE)
nzv %>% filter(nzv == T)

descrCor <- cor(predictor_df, method = "spearman")
highCorr <- sum(abs(descrCor[upper.tri(descrCor)]) > .8)
highCorr
summary(descrCor[upper.tri(descrCor)])

highlyCorDescr <- findCorrelation(descrCor, cutoff = .8)
highlyCorDescr
head(predictor_df[, highlyCorDescr])
predictor_df_filt <- predictor_df[, -highlyCorDescr]
descrCor2 <- cor(predictor_df_filt, method = "spearman")
summary(descrCor2[upper.tri(descrCor2)])
head(descrCor2)

comboInfo <- findLinearCombos(predictor_df_filt)
comboInfo

rem <- names(predictor_df_filt)[comboInfo$remove]
rem

# predictor_df_filt = predictor_df_filt %>% dplyr::select(-all_of(rem))

names(predictor_df_filt)

# plotting corr heatmap
heatmaply_cor(
  x = cor(predictor_df, method = "spearman"), xlab = "Features", plot_method = "plotly",
  ylab = "Features", k_col = 4, k_row = 4
)

heatmaply_cor(
  x = cor(predictor_df_filt, method = "spearman"), xlab = "Features", plot_method = "plotly",
  ylab = "Features", k_col = 4, k_row = 4
)

epi_red %>%
  select(contains("lag")) %>%
  names() %>%
  sort()

sort(names(predictor_df_filt))

final_pred_df <- epi_red2 %>% dplyr::select(region, season, year.new, names(predictor_df_filt))
head(final_pred_df)
names(final_pred_df)
nrow(final_pred_df)
head(epi_red2)
sort(names(epi_red2))

final_pred_df %>% filter(season == "2000-2001")
final_pred_df %>% group_by(region) %>% tally()

final_pred_df %>% filter(region=="Region 6")

write_rds(final_pred_df, file = "data/predictor_df_for_var_sel.rds")