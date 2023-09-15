## load packages
list.of.packages <- c(
  "dplyr", "ggplot2", "cowplot", "tidyr", "broom", "gmodels", "rsample", "metan",
  "betareg", "purrr", "tidymodels", "MuMIn", "readr", "RColorBrewer", "forcats"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

theme_set(theme_cowplot(font_size = 16))
########################################################################################################
## Univariate associations between evolutionary indicators and epidemic timing
########################################################################################################
## load data
load("data/antigenic_epi_north_amer_build_for_lasso_replicates.Rdata")

epi_red2 <- epi_red %>%
  dplyr::select(
    season, region, dom_type,
    contains(c("duration", "onset", "peak", "mean", "lag1", "lag2", "lbi"))
  ) %>%
  dplyr::select(-contains(c("lag0", "lag1_sd", "lag2_sd", "lbi_lag2", "std_lbi", "shannon_lag2", "usa"))) %>%
  filter(!(season %in% c("2009-2010", "2000-2001"))) ## no H3 circulation

epi_long <- epi_red2 %>%
  pivot_longer(
    cols = contains(c("lag", "lbi")),
    names_to = "evol_metrics",
    values_to = "value"
  )

h3_onset_predictors <- epi_long %>%
  filter(season != "2002-2003") %>% # only 3 regions have onsets in 2002-2003
  split(.$evol_metrics) %>% #
  map(~ lm(as.numeric(onset_days_from_Oct1) ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_onset_predictors, decreasing = T)[1:5]

h3_peak_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ lm(as.numeric(peak_days_from_Oct1) ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_peak_predictors, decreasing = T)[1:5]

h3_peak_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ lm(as.numeric(peak_timing_sd) ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_peak_predictors, decreasing = T)[1:5]

h3_onset_predictors <- epi_long %>%
  filter(season != "2002-2003") %>%
  split(.$evol_metrics) %>% #
  map(~ lm(as.numeric(onset_timing_sd) ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_onset_predictors, decreasing = T)[1:5]

h3_onset_predictors <- epi_long %>%
  filter(season != "2002-2003") %>%
  split(.$evol_metrics) %>% #
  map(~ lm(H3_onset_index_week ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_onset_predictors, decreasing = T)[1:5]

unique(epi_long$evol_metrics)
unique(epi_red2$season)

h3_duration_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ lm(H3_season_duration ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_duration_predictors, decreasing = T)[1:5]

###################################################################################################
### Evolutionary Indicators vs Epi Timing Heatmap
###################################################################################################
set.seed(27)
boots <- bootstraps(epi_red2, times = 1000, apparent = TRUE)

datalist <- list()
for (i in boots$id) {
  epi_long <- boots %>%
    filter(id == i) %>%
    pull(splits) %>%
    as.data.frame() %>%
    group_by(season) %>%
    dplyr::summarize_at(vars(
      H3_onset_index_week,
      onset_days_from_Oct1,
      peak_days_from_Oct1, peak_diff,
      H3_season_duration, peak_diff_median,
      peak_timing_sd, onset_timing_sd,
      HA_mean_lbi_lag1:na_lbi_shannon
    ), mean) %>%
    ungroup() %>%
    pivot_longer(
      cols = contains(c("lag", "lbi")),
      names_to = "evol_metrics",
      values_to = "value"
    )

  h3_onset_predictors <- epi_long %>%
    filter(season != "2002-2003") %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$H3_onset_index_week, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("estimate")

  h3_onset_pvalues <- epi_long %>%
    filter(season != "2002-2003") %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$H3_onset_index_week, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("p.value")

  h3_duration_predictors <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$H3_season_duration, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("estimate")

  h3_duration_pvalues <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$H3_season_duration, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("p.value")

  h3_max_predictors <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(as.numeric(.$peak_days_from_Oct1), scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("estimate")

  h3_max_pvalues <- epi_long %>%
    # split(f=list(.$replicate,.$evol_metrics))%>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(as.numeric(.$peak_days_from_Oct1), scale(.$value, center = T)), method = "spearman") %>%
    # map(summary) %>%
    map_dbl("p.value")

  h3_onset_predictors2 <- epi_long %>%
    filter(season != "2002-2003") %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(as.numeric(.$onset_days_from_Oct1), scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("estimate")

  h3_onset_pvalues2 <- epi_long %>%
    filter(season != "2002-2003") %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(as.numeric(.$onset_days_from_Oct1), scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("p.value")

  h3_peak_diff_predictors <- epi_long %>%
    distinct(season, peak_diff, evol_metrics, value) %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(as.numeric(.$peak_diff), scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("estimate")

  h3_peak_diff_pvalues <- epi_long %>%
    distinct(season, peak_diff, evol_metrics, value) %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(as.numeric(.$peak_diff), scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("p.value")


  h3_peak_sd_predictors <- epi_long %>%
    distinct(season, peak_timing_sd, evol_metrics, value) %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(as.numeric(.$peak_timing_sd), scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("estimate")

  h3_peak_sd_pvalues <- epi_long %>%
    distinct(season, peak_timing_sd, evol_metrics, value) %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(as.numeric(.$peak_timing_sd), scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("p.value")

  h3_onset_sd_predictors <- epi_long %>%
    filter(season != "2002-2003") %>%
    distinct(season, onset_timing_sd, evol_metrics, value) %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(as.numeric(.$onset_timing_sd), scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("estimate")

  h3_onset_sd_pvalues <- epi_long %>%
    filter(season != "2002-2003") %>%
    distinct(season, onset_timing_sd, evol_metrics, value) %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(as.numeric(.$onset_timing_sd), scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("p.value")

  epi_evol_df <- data.frame(
    h3_duration = h3_duration_predictors,
    h3_onset = h3_onset_predictors,
    h3_peak_days = h3_max_predictors,
    h3_onset_days = h3_onset_predictors2,
    h3_peak_diff = h3_peak_diff_predictors,
    h3_peak_sd = h3_peak_sd_predictors,
    h3_onset_sd = h3_onset_sd_predictors,
    h3_duration_p = h3_duration_pvalues,
    h3_onset_p = h3_onset_pvalues,
    h3_peak_days_p = h3_max_pvalues,
    h3_onset_days_p = h3_onset_pvalues2,
    h3_peak_diff_p = h3_peak_diff_pvalues,
    h3_peak_sd_p = h3_peak_sd_pvalues,
    h3_onset_sd_p = h3_onset_sd_pvalues
  )
  epi_evol_df$evol_metric <- rownames(epi_evol_df)

  epi_evol_df_long <- epi_evol_df %>%
    dplyr::select(-c(
      h3_onset_p, h3_peak_days_p, h3_onset_days_p, h3_peak_diff_p,
      h3_peak_sd_p, h3_onset_sd_p, h3_duration_p
    )) %>%
    pivot_longer(
      cols = c(
        h3_onset, h3_peak_days, h3_onset_days, h3_peak_diff,
        h3_peak_sd, h3_onset_sd, h3_duration
      ),
      names_to = "epi_metric", values_to = "cor"
    )

  epi_evol_df_long2 <- epi_evol_df %>%
    dplyr::select(-c(
      h3_onset, h3_peak_days, h3_onset_days, h3_peak_diff,
      h3_peak_sd, h3_onset_sd, h3_duration
    )) %>%
    dplyr::rename(c(
      h3_onset = h3_onset_p, h3_peak_days = h3_peak_days_p, h3_onset_days = h3_onset_days_p,
      h3_peak_diff = h3_peak_diff_p,
      h3_peak_sd = h3_peak_sd_p, h3_onset_sd = h3_onset_sd_p,
      h3_duration = h3_duration_p
    )) %>%
    pivot_longer(cols = c(
      h3_onset, h3_peak_days, h3_onset_days, h3_peak_diff,
      h3_peak_sd, h3_onset_sd, h3_duration
    ), names_to = "epi_metric", values_to = "pvalue")

  epi_combined <- left_join(epi_evol_df_long, epi_evol_df_long2, by = c("evol_metric", "epi_metric"))
  epi_combined$replicate <- i
  datalist[[i]] <- epi_combined
}
big_data <- do.call(rbind, datalist)
nrow(big_data)
head(big_data)

epi_evol_df_combined <-
  big_data %>%
  group_by(evol_metric, epi_metric) %>%
  dplyr::summarize(
    cor.mean = mean(cor),
    pvalue = mean(pvalue),
    cor.lowCI = ci(cor, na.rm = T)[2],
    cor.hiCI = ci(cor, na.rm = T)[3]
  ) %>%
  ungroup()

head(epi_evol_df_combined)
epi_evol_df_combined$pvalue_round <- round(epi_evol_df_combined$pvalue, 2)
epi_evol_df_combined %>% filter(pvalue < 0.1 & pvalue > 0.05)
epi_evol_df_combined %>% filter(pvalue < 0.1 & pvalue > 0.05)

epi_evol_df_combined <-
  epi_evol_df_combined %>%
  mutate(p_symbol = case_when(
    pvalue_round > 0.05 ~ "",
    pvalue_round <= 0.05 & pvalue_round >= 0.01 ~ "*",
    pvalue_round < 0.01 & pvalue_round >= 0.001 ~ "**",
    pvalue_round < 0.001 ~ "***"
  ))

epi_evol_df_combined$evol_metric <- as.factor(epi_evol_df_combined$evol_metric)
unique(epi_evol_df_combined$evol_metric)
levels(epi_evol_df_combined$evol_metric)
levels(epi_evol_df_combined$evol_metric) <- c(
  "H3 RBS distance (t-1)",
  "H3 RBS distance (t-2)",
  "H3 LBI Diversity",
  "H3 LBI Diversity (t-1)",
  "H3 mean LBI",
  "H3 mean LBI (t-1)",
  "H3 stalk footprint distance (t-1)",
  "H3 stalk footprint distance (t-2)",
  "HI titer distance (t-1)",
  "HI titer distance (t-2)",
  "H3 epitope distance (t-1)",
  "H3 epitope distance (t-2)",
  "H3 non-epitope distance (t-1)",
  "H3 non-epitope distance (t-2)",
  "N2 epitope distance (N=223) (t-1)",
  "N2 epitope distance (N=223) (t-2)",
  "N2 non-epitope distance (t-1)",
  "N2 non-epitope distance (t-2)",
  "N2 epitope distance (N=53) (t-1)",
  "N2 epitope distance (N=53) (t-2)",
  "N2 LBI Diversity",
  "N2 LBI Diversity (t-1)",
  "N2 mean LBI",
  "N2 mean LBI (t-1)"
)

epi_evol_df_combined$evol_metric <- factor(epi_evol_df_combined$evol_metric, levels = c(
  "N2 LBI Diversity (t-1)",
  "N2 LBI Diversity",
  "N2 mean LBI (t-1)",
  "N2 mean LBI",
  "N2 epitope distance (N=53) (t-2)",
  "N2 epitope distance (N=53) (t-1)",
  "N2 non-epitope distance (t-2)",
  "N2 non-epitope distance (t-1)",
  "N2 epitope distance (N=223) (t-2)",
  "N2 epitope distance (N=223) (t-1)",
  "H3 LBI Diversity (t-1)",
  "H3 LBI Diversity",
  "H3 mean LBI (t-1)",
  "H3 mean LBI",
  "HI titer distance (t-2)",
  "HI titer distance (t-1)",
  "H3 stalk footprint distance (t-2)",
  "H3 stalk footprint distance (t-1)",
  "H3 RBS distance (t-2)",
  "H3 RBS distance (t-1)",
  "H3 non-epitope distance (t-2)",
  "H3 non-epitope distance (t-1)",
  "H3 epitope distance (t-2)",
  "H3 epitope distance (t-1)"
))

epi_evol_df_combined$epi_metric <- as.factor(epi_evol_df_combined$epi_metric)
levels(epi_evol_df_combined$epi_metric)
levels(epi_evol_df_combined$epi_metric) <- c(
  "Duration",
  "Onset\nIndex\nWeek",
  "Onset Timing",
  "Onset Timing s.d.",
  "Peak Timing",
  "Days from\nOnset to Peak",
  "Peak Timing s.d."
)
epi_evol_df_combined$epi_metric <- factor(epi_evol_df_combined$epi_metric,
  levels = c(
    "Onset\nIndex\nWeek",
    "Onset Timing",
    "Onset Timing s.d.",
    "Peak Timing",
    "Peak Timing s.d.",
    "Days from\nOnset to Peak", "Duration"
  )
)
head(epi_evol_df_combined)

heat_p2 <- ggplot(data = epi_evol_df_combined %>%
  filter(epi_metric %in% c(
    "Onset Timing", "Onset Timing s.d.",
    "Peak Timing", "Days from\nOnset to Peak",
    "Peak Timing s.d.", "Duration"
  )) %>%
  droplevels()) +
  facet_grid(fct_rev(evol_metric) ~ epi_metric,
    scales = "fixed", switch = "x",
    as.table = T
  ) +
  ylim(0, 1) +
  xlim(-1, 1) +
  geom_errorbar(aes(x = cor.mean, xmin = cor.lowCI, xmax = cor.hiCI, y = 0.5), color = "black") +
  geom_point(aes(x = cor.mean, y = 0.5, fill = cor.mean), size = 4, shape = 21) +
  theme_bw() +
  scale_color_gradientn(colors = rev(brewer.pal(n = 7, name = "RdYlBu")), name = "correlation") +
  scale_fill_gradientn(colors = rev(brewer.pal(n = 7, name = "RdYlBu")), name = "correlation") +
  geom_text(aes(y = 0.7, x = 0.8, label = p_symbol), color = "black", size = 6) +
  theme(
    strip.background = element_blank(), strip.placement = "outside", legend.position = "none",
    axis.title.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    strip.text.x = element_text(size = 14),
    strip.text.y = element_text(size = 14)
  ) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0)) +
  geom_vline(xintercept = 0, lty = "dashed") +
  xlab("Spearman's rho")
heat_p2
save_plot(heat_p2, filename = "figures/epi_timing_evol_heatmap_alternate.png", base_width = 14, base_height = 14)

####################################################
## NA epitope distance
####################################################

####################################################
## Onset timing
####################################################
epi_red2 %>% filter(season == "2002-2003")
epi_red2 %>% filter(is.na(H3_onset_index_week))

sum_df <- epi_red2 %>%
  dplyr::select(season, dom_type, NA_bhatt_ep_lag1, onset_days_from_Oct1, peak_days_from_Oct1, HA_wolf_lag2) %>%
  distinct()

sum_df <-
  sum_df %>%
  filter(season != "2002-2003") %>% ## only 3 regions have onsets
  group_by(season, dom_type, HA_wolf_lag2, NA_bhatt_ep_lag1) %>%
  summarise(
    H3.onset.mean = ci(as.numeric(onset_days_from_Oct1), na.rm = T)[1],
    H3.onset.lowCI = ci(as.numeric(onset_days_from_Oct1), na.rm = T)[2],
    H3.onset.hiCI = ci(as.numeric(onset_days_from_Oct1), na.rm = T)[3],
    H3.onset.sd = ci(as.numeric(onset_days_from_Oct1), na.rm = T)[4]
  ) %>%
  ungroup()

sum_df <- sum_df %>% tidyr::separate(col = "season", sep = "-", remove = F, into = c("year1", "year2"))

sum_df <- sum_df %>%
  mutate(h1n1_type = if_else(year1 < 2010, "seasonal_h1n1", "pdm_h1n1")) %>%
  mutate(dom_type2 = case_when(
    year1 < 2009 & dom_type == "H1" ~ "H1",
    year1 > 2009 & dom_type == "H1" ~ "H1pdm",
    dom_type == "H3" ~ "H3",
    dom_type == "co-circ" ~ "H3/H1pdm"
  ))


sumdf2 <- sum_df %>% filter(season != "2009-2010")

scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2, NA_bhatt_ep_lag1), ~ scale_this(.x))

y <- sumdf2$H3.onset.mean
x <- sumdf2$NA_bhatt_ep_lag1

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))

model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)


fit_glm_on_bootstrap <- function(split) {
  glm(H3.onset.mean ~ NA_bhatt_ep_lag1, analysis(split), family = gaussian(link = "identity"))
}

boot_models <-
  boots %>%
  mutate(
    model = map(splits, fit_glm_on_bootstrap),
    coef_info = map(model, tidy)
  )

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

m1 <- glm(H3.onset.mean ~ NA_bhatt_ep_lag1, data = sumdf2, family = gaussian(link = "identity"))
fam <- family(m1)
ilink <- fam$linkinv
coefini <- coef(m1)

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.011

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.onset.mean ~ NA_bhatt_ep_lag1, data = .x, start = coefini, family = gaussian(link = "identity"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H3.onset.mean, NA_bhatt_ep_lag1, HA_wolf_lag2) %>%
  pivot_longer(cols = c(NA_bhatt_ep_lag1, HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.onset.mean, value)) %>%
  mutate(plot_data = map(model_data, ~ {
    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)

    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model,
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>%
        sample_n(n, TRUE) %>%
        glm(H3.onset.mean ~ value, ., start = coefini, family = gaussian(link = "identity")) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      dplyr::summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

cols <- c("#B24745FF", "#00A1D5FF", "#6A6599FF", "#DF8F44FF")

na_ep_onset <- ggplot() +
  geom_ribbon(aes(x = H3.onset.mean, ymin = gaussian(link = "identity")$linkinv(l), ymax = gaussian(link = "identity")$linkinv(u)),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(H3.onset.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = H3.onset.mean),
    method = "glm", formula = y ~ x,
    method.args = list(start = coefini, family = gaussian(link = "identity")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_bhatt_ep_lag1,
      ymin = H3.onset.lowCI, ymax = H3.onset.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = H3.onset.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-1)")) +
  ylab("Onset Timing (Days from Oct 1)") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  theme_bw(base_size = 16) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_text(
    y = 95, x = 0.75, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_ep_onset

####################################################
## Peak timing
####################################################

sum_df <- epi_red2 %>%
  dplyr::select(season, dom_type, NA_bhatt_ep_lag1, onset_days_from_Oct1, peak_days_from_Oct1, HA_wolf_lag2) %>%
  distinct()

sum_df <-
  sum_df %>%
  # filter(season != "2002-2003") %>%
  group_by(season, dom_type, HA_wolf_lag2, NA_bhatt_ep_lag1) %>%
  summarise(
    H3.peak.mean = ci(as.numeric(peak_days_from_Oct1), na.rm = T)[1],
    H3.peak.lowCI = ci(as.numeric(peak_days_from_Oct1), na.rm = T)[2],
    H3.peak.hiCI = ci(as.numeric(peak_days_from_Oct1), na.rm = T)[3],
    H3.peak.sd = ci(as.numeric(peak_days_from_Oct1), na.rm = T)[4]
  ) %>%
  ungroup()

sum_df <- sum_df %>% tidyr::separate(col = "season", sep = "-", remove = F, into = c("year1", "year2"))
sum_df <- sum_df %>%
  mutate(h1n1_type = if_else(year1 < 2010, "seasonal_h1n1", "pdm_h1n1")) %>%
  mutate(dom_type2 = case_when(
    year1 < 2009 & dom_type == "H1" ~ "H1",
    year1 > 2009 & dom_type == "H1" ~ "H1pdm",
    dom_type == "H3" ~ "H3",
    dom_type == "co-circ" ~ "H3/H1pdm"
  ))

sumdf2 <- sum_df %>% filter(season != "2009-2010")

scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2, NA_bhatt_ep_lag1), ~ scale_this(.x))

y <- sumdf2$H3.peak.mean
x <- sumdf2$NA_bhatt_ep_lag1

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.peak.mean ~ NA_bhatt_ep_lag1, analysis(split), family = gaussian(link = "identity"))
}

boot_models <-
  boots %>%
  mutate(
    model = map(splits, fit_glm_on_bootstrap),
    coef_info = map(model, tidy)
  )

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

m1 <- glm(H3.peak.mean ~ NA_bhatt_ep_lag1, data = sumdf2, family = gaussian(link = "identity"))
fam <- family(m1)
ilink <- fam$linkinv
coefini <- coef(m1)

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.1


labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.peak.mean ~ NA_bhatt_ep_lag1, data = .x, start = coefini, family = gaussian(link = "identity"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H3.peak.mean, NA_bhatt_ep_lag1, HA_wolf_lag2) %>%
  pivot_longer(cols = c(NA_bhatt_ep_lag1, HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.peak.mean, value)) %>%
  mutate(plot_data = map(model_data, ~ {
    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)

    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model,
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>%
        sample_n(n, TRUE) %>%
        glm(H3.peak.mean ~ value, ., start = coefini, family = gaussian(link = "identity")) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      dplyr::summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data


na_ep_peak <- ggplot() +
  geom_ribbon(aes(x = H3.peak.mean, ymin = gaussian(link = "identity")$linkinv(l), ymax = gaussian(link = "identity")$linkinv(u)),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(H3.peak.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = H3.peak.mean),
    method = "glm", formula = y ~ x,
    method.args = list(start = coefini, family = gaussian(link = "identity")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_bhatt_ep_lag1,
      ymin = H3.peak.lowCI, ymax = H3.peak.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = H3.peak.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-1)")) +
  ylab("Peak Timing (Days from Oct 1)") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  theme_bw(base_size = 16) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_text(
    y = 180, x = 0.75, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_ep_peak

####################################################
## A/H3 Variation in Peak timing (s.d.)
####################################################
sum_df <- epi_red2 %>%
  dplyr::select(season, dom_type, NA_bhatt_ep_lag1, peak_timing_sd, peak_diff, HA_wolf_lag2) %>%
  distinct()

sum_df <-
  sum_df %>%
  # filter(season != "2002-2003") %>%
  group_by(season, dom_type, peak_timing_sd, HA_wolf_lag2, NA_bhatt_ep_lag1) %>%
  summarise(
    H3.peak.diff.mean = ci(as.numeric(peak_diff), na.rm = T)[1],
    H3.peak.diff.lowCI = ci(as.numeric(peak_diff), na.rm = T)[2],
    H3.peak.diff.hiCI = ci(as.numeric(peak_diff), na.rm = T)[3],
    H3.peak.diff.sd = ci(as.numeric(peak_diff), na.rm = T)[4]
  ) %>%
  ungroup()

sum_df <- sum_df %>% tidyr::separate(col = "season", sep = "-", remove = F, into = c("year1", "year2"))

sum_df <- sum_df %>%
  mutate(h1n1_type = if_else(year1 < 2010, "seasonal_h1n1", "pdm_h1n1")) %>%
  mutate(dom_type2 = case_when(
    year1 < 2009 & dom_type == "H1" ~ "H1",
    year1 > 2009 & dom_type == "H1" ~ "H1pdm",
    dom_type == "H3" ~ "H3",
    dom_type == "co-circ" ~ "H3/H1pdm"
  ))


sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2, NA_bhatt_ep_lag1), ~ scale_this(.x))

y <- sumdf2$peak_timing_sd
x <- sumdf2$NA_bhatt_ep_lag1

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(peak_timing_sd ~ NA_bhatt_ep_lag1, analysis(split), family = Gamma(link = "inverse"))
}

boot_models <-
  boots %>%
  mutate(
    model = map(splits, fit_glm_on_bootstrap),
    coef_info = map(model, tidy)
  )

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

m1 <- glm(peak_timing_sd ~ NA_bhatt_ep_lag1, data = sumdf2, family = Gamma(link = "inverse"))
fam <- family(m1)
ilink <- fam$linkinv
coefini <- coef(m1)

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(peak_timing_sd ~ NA_bhatt_ep_lag1, data = .x, start = coefini, family = Gamma(link = "inverse"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(peak_timing_sd, NA_bhatt_ep_lag1, HA_wolf_lag2) %>%
  pivot_longer(cols = c(NA_bhatt_ep_lag1, HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(peak_timing_sd, value)) %>%
  mutate(plot_data = map(model_data, ~ {
    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)

    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model,
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>%
        sample_n(n, TRUE) %>%
        glm(peak_timing_sd ~ value, ., start = coefini, family = Gamma(link = "inverse"), maxit = 150) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      dplyr::summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

na_ep_peak_sd <- ggplot() +
  geom_ribbon(aes(x = peak_timing_sd, ymin = Gamma(link = "inverse")$linkinv(l), ymax = Gamma(link = "inverse")$linkinv(u)),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(peak_timing_sd = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = peak_timing_sd),
    method = "glm", formula = y ~ x,
    method.args = list(family = Gamma(link = "inverse"), start = coefini),
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = peak_timing_sd, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-1)")) +
  ylab("Variation (s.d.) in Peak Timing") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  theme_bw(base_size = 16) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_text(
    y = 45, x = 0.7, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_ep_peak_sd

####################################################
## A/H3 Days from Onset to Peak
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010", "2002-2003")))

scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2, NA_bhatt_ep_lag1), ~ scale_this(.x))

y <- as.numeric(sumdf2$H3.peak.diff.mean)
x <- sumdf2$NA_bhatt_ep_lag1

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(as.numeric(H3.peak.diff.mean) ~ NA_bhatt_ep_lag1, analysis(split), family = Gamma(link = "inverse"))
}

boot_models <-
  boots %>%
  mutate(
    model = map(splits, fit_glm_on_bootstrap),
    coef_info = map(model, tidy)
  )

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

m1 <- glm(as.numeric(H3.peak.diff.mean) ~ NA_bhatt_ep_lag1, data = sumdf2, family = Gamma(link = "inverse"))
fam <- family(m1)
ilink <- fam$linkinv
coefini <- coef(m1)

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.3

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(as.numeric(H3.peak.diff.mean) ~ NA_bhatt_ep_lag1, start = coefini, data = .x, family = Gamma(link = "inverse"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H3.peak.diff.mean, NA_bhatt_ep_lag1, HA_wolf_lag2) %>%
  pivot_longer(cols = c(NA_bhatt_ep_lag1, HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.peak.diff.mean, value)) %>%
  mutate(plot_data = map(model_data, ~ {
    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)

    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model,
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>%
        sample_n(n, TRUE) %>%
        glm(as.numeric(H3.peak.diff.mean) ~ value, ., start = coefini, family = Gamma(link = "inverse")) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      dplyr::summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

sumdf2 %>%
  dplyr::select(season, H3.peak.diff.hiCI) %>%
  arrange(H3.peak.diff.hiCI)

na_ep_peak_diff <- ggplot() +
  geom_ribbon(aes(x = H3.peak.diff.mean, ymin = Gamma(link = "inverse")$linkinv(l), ymax = Gamma(link = "inverse")$linkinv(u)),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(H3.peak.diff.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = H3.peak.diff.mean),
    method = "glm", formula = y ~ x,
    method.args = list(start = coefini, family = Gamma(link = "inverse")),
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2 %>% mutate(H3.peak.diff.lowCI = ifelse(H3.peak.diff.lowCI < 0, 0, H3.peak.diff.lowCI)),
    aes(
      x = NA_bhatt_ep_lag1,
      ymin = H3.peak.diff.lowCI, ymax = H3.peak.diff.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = as.numeric(H3.peak.diff.mean), fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-1)")) +
  ylab("Days between Onset and Peak") +
  theme(legend.position = "bottom") +
  theme_bw(base_size = 16) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_text(
    y = 120, x = 0.75, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_ep_peak_diff

####################################################
## HA epitope distance
####################################################
####################################################
## Onset timing
####################################################

sum_df <- epi_red2 %>%
  dplyr::select(season, dom_type, NA_bhatt_ep_lag1, onset_days_from_Oct1, peak_days_from_Oct1, HA_wolf_lag2) %>%
  distinct()

sum_df <-
  sum_df %>%
  filter(season != "2002-2003") %>%
  group_by(season, dom_type, HA_wolf_lag2, NA_bhatt_ep_lag1) %>%
  summarise(
    H3.onset.mean = ci(as.numeric(onset_days_from_Oct1), na.rm = T)[1],
    H3.onset.lowCI = ci(as.numeric(onset_days_from_Oct1), na.rm = T)[2],
    H3.onset.hiCI = ci(as.numeric(onset_days_from_Oct1), na.rm = T)[3],
    H3.onset.sd = ci(as.numeric(onset_days_from_Oct1), na.rm = T)[4]
  ) %>%
  ungroup()

sum_df <- sum_df %>% tidyr::separate(col = "season", sep = "-", remove = F, into = c("year1", "year2"))
sum_df <- sum_df %>%
  mutate(h1n1_type = if_else(year1 < 2010, "seasonal_h1n1", "pdm_h1n1")) %>%
  mutate(dom_type2 = case_when(
    year1 < 2009 & dom_type == "H1" ~ "H1",
    year1 > 2009 & dom_type == "H1" ~ "H1pdm",
    dom_type == "H3" ~ "H3",
    dom_type == "co-circ" ~ "H3/H1pdm"
  ))

sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2, NA_bhatt_ep_lag1), ~ scale_this(.x))

y <- sumdf2$H3.onset.mean
x <- sumdf2$HA_wolf_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))

model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")
# use LM to be consistent with NA epitope distance

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.onset.mean ~ HA_wolf_lag2, analysis(split), family = gaussian(link = "identity"))
}

boot_models <-
  boots %>%
  mutate(
    model = map(splits, fit_glm_on_bootstrap),
    coef_info = map(model, tidy)
  )

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

m1 <- glm(H3.onset.mean ~ HA_wolf_lag2, data = sumdf2, family = gaussian(link = "identity"))
fam <- family(m1)
ilink <- fam$linkinv
coefini <- coef(m1)

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.onset.mean ~ HA_wolf_lag2, data = .x, start = coefini, family = gaussian(link = "identity"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H3.onset.mean, HA_wolf_lag2, HA_wolf_lag2) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.onset.mean, value)) %>%
  mutate(plot_data = map(model_data, ~ {
    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)

    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model,
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>%
        sample_n(n, TRUE) %>%
        glm(H3.onset.mean ~ value, ., start = coefini, family = gaussian(link = "identity")) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      dplyr::summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

ha_ep_onset <- ggplot() +
  geom_ribbon(aes(x = H3.onset.mean, ymin = gaussian(link = "identity")$linkinv(l), ymax = gaussian(link = "identity")$linkinv(u)),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>% rename(H3.onset.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = H3.onset.mean),
    method = "glm", formula = y ~ x,
    method.args = list(start = coefini, family = gaussian(link = "identity")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_wolf_lag2,
      ymin = H3.onset.lowCI, ymax = H3.onset.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = H3.onset.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab(expression("H3 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("Onset Timing (Days from Oct 1)") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  theme_bw(base_size = 16) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_text(
    y = 95, x = 0.75, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ha_ep_onset

####################################################
## Peak timing
####################################################

sum_df <- epi_red2 %>%
  dplyr::select(season, dom_type, NA_bhatt_ep_lag1, onset_days_from_Oct1, peak_days_from_Oct1, HA_wolf_lag2) %>%
  distinct()

sum_df <-
  sum_df %>%
  # filter(season != "2002-2003") %>%
  group_by(season, dom_type, HA_wolf_lag2, NA_bhatt_ep_lag1) %>%
  summarise(
    H3.peak.mean = ci(as.numeric(peak_days_from_Oct1), na.rm = T)[1],
    H3.peak.lowCI = ci(as.numeric(peak_days_from_Oct1), na.rm = T)[2],
    H3.peak.hiCI = ci(as.numeric(peak_days_from_Oct1), na.rm = T)[3],
    H3.peak.sd = ci(as.numeric(peak_days_from_Oct1), na.rm = T)[4]
  ) %>%
  ungroup()

sum_df <- sum_df %>% tidyr::separate(col = "season", sep = "-", remove = F, into = c("year1", "year2"))

sum_df <- sum_df %>%
  mutate(h1n1_type = if_else(year1 < 2010, "seasonal_h1n1", "pdm_h1n1")) %>%
  mutate(dom_type2 = case_when(
    year1 < 2009 & dom_type == "H1" ~ "H1",
    year1 > 2009 & dom_type == "H1" ~ "H1pdm",
    dom_type == "H3" ~ "H3",
    dom_type == "co-circ" ~ "H3/H1pdm"
  ))

sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2, NA_bhatt_ep_lag1), ~ scale_this(.x))

y <- sumdf2$H3.peak.mean
x <- sumdf2$HA_wolf_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")
## use LM to be consistent with NA epitope distance

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.peak.mean ~ HA_wolf_lag2, analysis(split), family = gaussian(link = "identity"))
}

boot_models <-
  boots %>%
  mutate(
    model = map(splits, fit_glm_on_bootstrap),
    coef_info = map(model, tidy)
  )

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

m1 <- glm(H3.peak.mean ~ HA_wolf_lag2, data = sumdf2, family = gaussian(link = "identity"))
fam <- family(m1)
ilink <- fam$linkinv
coefini <- coef(m1)

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.peak.mean ~ HA_wolf_lag2, data = .x, start = coefini, family = gaussian(link = "identity"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H3.peak.mean, HA_wolf_lag2, HA_wolf_lag2) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.peak.mean, value)) %>%
  mutate(plot_data = map(model_data, ~ {
    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)

    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model,
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>%
        sample_n(n, TRUE) %>%
        glm(H3.peak.mean ~ value, ., start = coefini, family = gaussian(link = "identity")) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      dplyr::summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

ha_ep_peak <- ggplot() +
  geom_ribbon(aes(x = H3.peak.mean, ymin = gaussian(link = "identity")$linkinv(l), ymax = gaussian(link = "identity")$linkinv(u)),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>% rename(H3.peak.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = H3.peak.mean),
    method = "glm", formula = y ~ x,
    method.args = list(start = coefini, family = gaussian(link = "identity")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_wolf_lag2,
      ymin = H3.peak.lowCI, ymax = H3.peak.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = H3.peak.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab(expression("H3 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("Peak Timing (Days from Oct 1)") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  theme_bw(base_size = 16) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_text(
    y = 180, x = 0.7, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ha_ep_peak

# epi_timing_fig <- plot_grid(ha_ep_onset + theme(legend.position = "none"),
#   na_ep_onset + theme(legend.position = "none"),
#   ha_ep_peak + theme(legend.position = "none"),
#   na_ep_peak + theme(legend.position = "none"),
#   rel_widths = c(1, 1),
#   nrow = 2, labels = NULL
# )

epi_timing1 <- plot_grid(ha_ep_onset + theme(legend.position = "none"),
  na_ep_onset + theme(legend.position = "none"),
  rel_widths = c(1, 1),
  nrow = 1, labels = NULL
)

epi_timing2 <- plot_grid(ha_ep_peak + theme(legend.position = "none"),
  na_ep_peak + theme(legend.position = "none"),
  rel_widths = c(1, 1),
  nrow = 1, labels = NULL
)


epi_timing_fig <- plot_grid(epi_timing1, epi_timing2, nrow = 2, labels = "AUTO")
epi_timing_fig

epi_leg <- get_legend(ha_ep_peak + theme(
  legend.direction = "horizontal",
  legend.justification = "center", legend.box.just = "bottom", legend.position = "bottom"
))
epi_timing_fig_leg <- plot_grid(epi_timing_fig, epi_leg, nrow = 2, rel_heights = c(2, 0.1))
epi_timing_fig_leg
save_plot(epi_timing_fig_leg, filename = "figures/peak_and_onset_epi_timing_north_amer_build.png", base_width = 12, base_height = 12)

####################################################
## A/H3 Variation in Peak timing (s.d.)
####################################################
sum_df <- epi_red2 %>%
  dplyr::select(season, dom_type, NA_bhatt_ep_lag1, peak_timing_sd, peak_diff, HA_wolf_lag2) %>%
  distinct()
head(sum_df)

sum_df <-
  sum_df %>%
  # filter(season != "2002-2003") %>%
  group_by(season, dom_type, peak_timing_sd, HA_wolf_lag2, NA_bhatt_ep_lag1) %>%
  summarise(
    H3.peak.diff.mean = ci(as.numeric(peak_diff), na.rm = T)[1],
    H3.peak.diff.lowCI = ci(as.numeric(peak_diff), na.rm = T)[2],
    H3.peak.diff.hiCI = ci(as.numeric(peak_diff), na.rm = T)[3],
    H3.peak.diff.sd = ci(as.numeric(peak_diff), na.rm = T)[4]
  ) %>%
  ungroup()

sum_df <- sum_df %>% tidyr::separate(col = "season", sep = "-", remove = F, into = c("year1", "year2"))

sum_df <- sum_df %>%
  mutate(h1n1_type = if_else(year1 < 2010, "seasonal_h1n1", "pdm_h1n1")) %>%
  mutate(dom_type2 = case_when(
    year1 < 2009 & dom_type == "H1" ~ "H1",
    year1 > 2009 & dom_type == "H1" ~ "H1pdm",
    dom_type == "H3" ~ "H3",
    dom_type == "co-circ" ~ "H3/H1pdm"
  ))

sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2, NA_bhatt_ep_lag1), ~ scale_this(.x))

y <- sumdf2$peak_timing_sd
x <- sumdf2$HA_wolf_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(peak_timing_sd ~ HA_wolf_lag2, analysis(split), family = Gamma(link = "identity"))
}

boot_models <-
  boots %>%
  mutate(
    model = map(splits, fit_glm_on_bootstrap),
    coef_info = map(model, tidy)
  )

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

m1 <- glm(peak_timing_sd ~ HA_wolf_lag2, data = sumdf2, family = Gamma(link = "identity"))
fam <- family(m1)
ilink <- fam$linkinv
coefini <- coef(m1)

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.08

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(peak_timing_sd ~ HA_wolf_lag2, data = .x, start = coefini, family = Gamma(link = "identity"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(peak_timing_sd, HA_wolf_lag2, NA_bhatt_ep_lag1) %>%
  pivot_longer(cols = c(HA_wolf_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(peak_timing_sd, value)) %>%
  mutate(plot_data = map(model_data, ~ {
    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)

    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model,
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>%
        sample_n(n, TRUE) %>%
        glm(peak_timing_sd ~ value, ., start = coefini, family = Gamma(link = "identity")) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      dplyr::summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

ha_ep_peak_sd <- ggplot() +
  geom_ribbon(aes(x = peak_timing_sd, ymin = Gamma(link = "identity")$linkinv(l), ymax = Gamma(link = "identity")$linkinv(u)),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>% rename(peak_timing_sd = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = peak_timing_sd),
    method = "glm", formula = y ~ x,
    method.args = list(start = coefini, family = Gamma(link = "identity")),
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  # geom_line(aes(y = Gamma(link="identity")$linkinv(.fitted), x=HA_wolf_lag2, group = id), alpha = .2, col = "grey",data=boot_aug) +
  geom_point(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = peak_timing_sd, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab(expression("H3 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("Variation (s.d.) in Peak Timing") +
  theme_bw(base_size = 16) +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_text(
    y = 45, x = 0.7, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ha_ep_peak_sd

####################################################
## A/H3 Onset to Peak
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010", "2002-2003")))
head(sumdf2)
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2, NA_bhatt_ep_lag1), ~ scale_this(.x))

y <- as.numeric(sumdf2$H3.peak.diff.mean)
x <- sumdf2$HA_wolf_lag2

linear.model <- glm(y ~ x, family = gaussian())
identity.model <- glm(y ~ x, family = gaussian(link = "identity"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "identity"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(linear.model, identity.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)
names(sumdf2)

fit_glm_on_bootstrap <- function(split) {
  glm(as.numeric(H3.peak.diff.mean) ~ HA_wolf_lag2, analysis(split), family = Gamma(link = "identity"))
}

boot_models <-
  boots %>%
  mutate(
    model = map(splits, fit_glm_on_bootstrap),
    coef_info = map(model, tidy)
  )

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

m1 <- glm(as.numeric(H3.peak.diff.mean) ~ HA_wolf_lag2, data = sumdf2, family = Gamma(link = "identity"))
fam <- family(m1)
ilink <- fam$linkinv
coefini <- coef(m1)

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.14

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(as.numeric(H3.peak.diff.mean) ~ HA_wolf_lag2, start = coefini, data = .x, family = Gamma(link = "identity"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H3.peak.diff.mean, HA_wolf_lag2) %>%
  pivot_longer(cols = c(HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.peak.diff.mean, value)) %>%
  mutate(plot_data = map(model_data, ~ {
    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)

    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model,
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>%
        sample_n(n, TRUE) %>%
        glm(as.numeric(H3.peak.diff.mean) ~ value, ., start = coefini, family = Gamma(link = "identity")) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      dplyr::summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data


ha_ep_peak_diff <- ggplot() +
  geom_ribbon(aes(x = H3.peak.diff.mean, ymin = Gamma(link = "identity")$linkinv(l), ymax = Gamma(link = "identity")$linkinv(u)),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>% rename(H3.peak.diff.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = H3.peak.diff.mean),
    method = "glm", formula = y ~ x,
    method.args = list(start = coefini, family = Gamma(link = "identity")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2 %>% mutate(H3.peak.diff.lowCI = ifelse(H3.peak.diff.lowCI < 0, 0, H3.peak.diff.lowCI)),
    aes(
      x = HA_wolf_lag2,
      ymin = H3.peak.diff.lowCI, ymax = H3.peak.diff.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = as.numeric(H3.peak.diff.mean), fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("H3 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("Days between Onset and Peak") +
  theme(legend.position = "bottom") +
  theme_bw(base_size = 16) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_text(
    y = 120, x = 0.75, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ha_ep_peak_diff

ha_vs_na_epi_fig <- plot_grid(
  na_ep_peak_diff + theme(legend.position = "none"),
  ha_ep_peak_diff + theme(legend.position = "none"),
  rel_widths = c(1, 1),
  nrow = 1, labels = "AUTO"
)
epi_leg <- get_legend(ha_ep_peak_diff + theme(
  legend.direction = "horizontal",
  legend.justification = "center",
  legend.box.just = "bottom", legend.position = "bottom"
))
all_epi_leg_HA_NA <- plot_grid(ha_vs_na_epi_fig, epi_leg, nrow = 2, rel_heights = c(1, 0.1))
all_epi_leg_HA_NA
save_plot(all_epi_leg_HA_NA, filename = "figures/ha_vs_na_onset_to_peak.png", base_width = 12, base_height = 6)

###############################################################################################
## Seasonal Duration
###############################################################################################

###############################################################################################
## H3 LBI vs seasonal duration
###############################################################################################

sum_df <- epi_red2 %>%
  dplyr::select(
    season, dom_type, NA_bhatt_ep_lag2, HA_mean_lbi, NA_mean_lbi, HA_wolf_lag2,
    ha_lbi_shannon, na_lbi_shannon,
    H3_season_duration
  ) %>%
  distinct()
sum_df$H3_season_duration <- as.numeric(sum_df$H3_season_duration)

sum_df <-
  sum_df %>%
  # filter(season != "2002-2003") %>%
  group_by(season, dom_type, NA_bhatt_ep_lag2, HA_wolf_lag2, HA_mean_lbi, NA_mean_lbi, ha_lbi_shannon, na_lbi_shannon) %>%
  summarise(
    H3.duration.mean = ci(as.numeric(H3_season_duration), na.rm = T)[1],
    H3.duration.lowCI = ci(as.numeric(H3_season_duration), na.rm = T)[2],
    H3.duration.hiCI = ci(as.numeric(H3_season_duration), na.rm = T)[3],
    H3.duration.sd = ci(as.numeric(H3_season_duration), na.rm = T)[4]
  ) %>%
  ungroup()

sum_df <- sum_df %>% tidyr::separate(col = "season", sep = "-", remove = F, into = c("year1", "year2"))
sum_df <- sum_df %>%
  mutate(h1n1_type = if_else(year1 < 2010, "seasonal_h1n1", "pdm_h1n1")) %>%
  mutate(dom_type2 = case_when(
    year1 < 2009 & dom_type == "H1" ~ "H1",
    year1 > 2009 & dom_type == "H1" ~ "H1pdm",
    dom_type == "H3" ~ "H3",
    dom_type == "co-circ" ~ "H3/H1pdm"
  ))

sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010")))
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(NA_bhatt_ep_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- as.numeric(sumdf2$H3.duration.mean)
x <- sumdf2$ha_lbi_shannon

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(as.numeric(H3.duration.mean) ~ ha_lbi_shannon, analysis(split), family = gaussian(link = "inverse"))
}

boot_models <-
  boots %>%
  mutate(
    model = map(splits, fit_glm_on_bootstrap),
    coef_info = map(model, tidy)
  )

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

m1 <- glm(as.numeric(H3.duration.mean) ~ ha_lbi_shannon, data = sumdf2, family = gaussian(link = "inverse"))
fam <- family(m1)
ilink <- fam$linkinv
coefini <- coef(m1)

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.54

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(as.numeric(H3.duration.mean) ~ ha_lbi_shannon, start = coefini, data = .x, family = gaussian(link = "inverse"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.3f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H3.duration.mean, ha_lbi_shannon, HA_wolf_lag2) %>%
  pivot_longer(cols = c(ha_lbi_shannon, HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.duration.mean, value)) %>%
  mutate(plot_data = map(model_data, ~ {
    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)

    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model,
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>%
        sample_n(n, TRUE) %>%
        glm(as.numeric(H3.duration.mean) ~ value, ., start = coefini, family = gaussian(link = "inverse")) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      dplyr::summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

ha_lbi_duration <- ggplot() +
  geom_ribbon(aes(x = ha_lbi_shannon, ymin = gaussian(link = "inverse")$linkinv(l), ymax = gaussian(link = "inverse")$linkinv(u)),
    tbl_plot_data %>% filter(name == "ha_lbi_shannon") %>% rename(ha_lbi_shannon = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = ha_lbi_shannon, y = H3.duration.mean),
    method = "glm", formula = y ~ x,
    method.args = list(start = coefini, family = gaussian(link = "inverse")),
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2 %>% mutate(H3.duration.lowCI = ifelse(H3.duration.lowCI < 0, 0, H3.duration.lowCI)),
    aes(
      x = ha_lbi_shannon,
      ymin = H3.duration.lowCI, ymax = H3.duration.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = ha_lbi_shannon, y = as.numeric(H3.duration.mean), fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("H3 LBI Diversity") +
  ylab("Season Duration (Weeks)") +
  theme(legend.position = "bottom") +
  theme_bw(base_size = 16) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_text(
    y = 38, x = -2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ha_lbi_duration

###############################################################################################
## N2 LBI vs seasonal duration
###############################################################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010")))
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(NA_bhatt_ep_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- as.numeric(sumdf2$H3.duration.mean)
x <- sumdf2$na_lbi_shannon

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(as.numeric(H3.duration.mean) ~ na_lbi_shannon, analysis(split), family = gaussian(link = "inverse"))
}

boot_models <-
  boots %>%
  mutate(
    model = map(splits, fit_glm_on_bootstrap),
    coef_info = map(model, tidy)
  )

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

m1 <- glm(as.numeric(H3.duration.mean) ~ na_lbi_shannon, data = sumdf2, family = gaussian(link = "inverse"))
fam <- family(m1)
ilink <- fam$linkinv
coefini <- coef(m1)

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.53

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(as.numeric(H3.duration.mean) ~ na_lbi_shannon, start = coefini, data = .x, family = gaussian(link = "inverse"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.3f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H3.duration.mean, na_lbi_shannon, HA_wolf_lag2) %>%
  pivot_longer(cols = c(na_lbi_shannon, HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.duration.mean, value)) %>%
  mutate(plot_data = map(model_data, ~ {
    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)

    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model,
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>%
        sample_n(n, TRUE) %>%
        glm(as.numeric(H3.duration.mean) ~ value, ., start = coefini, family = gaussian(link = "inverse")) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      dplyr::summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

na_lbi_duration <- ggplot() +
  geom_ribbon(aes(x = na_lbi_shannon, ymin = gaussian(link = "inverse")$linkinv(l), ymax = gaussian(link = "inverse")$linkinv(u)),
    tbl_plot_data %>% filter(name == "na_lbi_shannon") %>% rename(na_lbi_shannon = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = na_lbi_shannon, y = H3.duration.mean),
    method = "glm", formula = y ~ x,
    method.args = list(start = coefini, family = gaussian(link = "inverse")),
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2 %>% mutate(H3.duration.lowCI = ifelse(H3.duration.lowCI < 0, 0, H3.duration.lowCI)),
    aes(
      x = na_lbi_shannon,
      ymin = H3.duration.lowCI, ymax = H3.duration.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = na_lbi_shannon, y = as.numeric(H3.duration.mean), fill = dom_type2), size = 5, pch = 21
  ) +
  # xlab(expression("N2 epitope distance ("~italic(t)~"-2)"))+
  xlab("N2 LBI Diversity") +
  ylab("Season Duration (Weeks)") +
  theme(legend.position = "bottom") +
  theme_bw(base_size = 16) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_text(
    y = 38, x = -2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_lbi_duration


lbi_epi_fig <- plot_grid(ha_lbi_duration + theme(legend.position = "none"),
  na_lbi_duration + theme(legend.position = "none"),
  rel_widths = c(1, 1),
  nrow = 1, labels = c("A", "B")
)
lbi_epi_fig
epi_leg <- get_legend(ha_lbi_duration + theme(
  legend.direction = "horizontal",
  legend.justification = "center", legend.box.just = "bottom", legend.position = "bottom"
))
lbi_leg_NA <- plot_grid(lbi_epi_fig, epi_leg, nrow = 2, rel_heights = c(2, 0.1))
lbi_leg_NA

save_plot(lbi_leg_NA, filename = "figures/ha_vs_na_lbi_season_duration.png", base_width = 12, base_height = 6)
