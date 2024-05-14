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
## Univariate associations between evolutionary indicators and age-specific patterns
########################################################################################################
## load data
load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")
head(epi_red)

load("data/age_distribution_ili_cases.RData")
names(age_dist)
names(age_dist)[grepl("prop",names(age_dist))] <- paste0("ili_", names(age_dist)[grepl("prop",names(age_dist))])
head(age_dist)

epi_red2 <- left_join(
  age_dist %>% dplyr::select(season_description, region, contains("prop")) %>%
    rename(season = season_description) %>%
    filter(!(season %in% c("2019-2020", "2020-2021", "2021-2022", "2009-2010"))),
  epi_red,
  by = c("season", "region")
)
unique(epi_red2$season)

epi_red2 <- epi_red2 %>% tidyr::separate(season, into = c("year1", "year2"), remove = F)

epi_red2 <- epi_red2 %>% filter(!(region == "Region 10" & year1 < 2009))
epi_red2[!complete.cases(epi_red2), ] %>% dplyr::select(season, region, year1)

epi_red2 %>%
  group_by(region) %>%
  tally()

names(epi_red2)

epi_long <- epi_red2 %>%
  filter(!(season %in% c("2009-2010")))%>%
  pivot_longer(
    cols = c(HA_std_lbi:na_lbi_shannon),
    names_to = "evol_metrics",
    values_to = "value"
  ) %>%
  distinct() %>%
  ungroup()

epi_long
sort(unique(epi_long$evol_metrics))

h3_0_4_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ lm(ili_age_0_4_prop ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_0_4_predictors, decreasing = T)[1:5]

h3_5_24_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ lm(ili_age_5_24_prop ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_5_24_predictors, decreasing = T)[1:5]

h3_25_64_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ lm(ili_age_25_64_prop ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_25_64_predictors, decreasing = T)[1:5]

h3_65_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ lm(ili_age_65_prop ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_65_predictors, decreasing = T)[1:5]

epi_red2 <- epi_red2 %>%
  dplyr::select(season, region, dom_type, contains(c(
    "prop",
    "mean", "lag1", "lag2"
  ))) %>%
  dplyr::select(-contains(c("lbi","usa")))%>%
  distinct()
nrow(epi_red2)
head(epi_red2)
epi_red2[!complete.cases(epi_red2), ] %>% dplyr::select(season, region)

###############################################
## Heatmap
###############################################

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

set.seed(27)
boots <- bootstraps(epi_red2, times = 1000, apparent = TRUE)
# boots
datalist <- list()
for (i in boots$id) {
  # i=boots$id[1]
  epi_long <- boots %>%
    filter(id == i) %>%
    pull(splits) %>%
    as.data.frame() %>%
    group_by(season, dom_type) %>%
    summarize_at(vars(
      contains("prop"),
      HA_koel_lag1:NA_krammer_ep_lag2
    ), mean) %>%
    ungroup() %>%
    pivot_longer(
      cols = contains(c("lag")),
      names_to = "evol_metrics",
      values_to = "value"
    )
  
  epi_long = epi_long %>%
    group_by(evol_metrics)%>%
    mutate(value = scale_this(value))

  ili_h3_0_4_predictors <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$ili_age_0_4_prop, .$value, method = "spearman")) %>%
    map_dbl("estimate")

  ili_h3_0_4_pvalues <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$ili_age_0_4_prop, .$value, method = "spearman")) %>%
    map_dbl("p.value")

  ili_h3_5_24_predictors <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$ili_age_5_24_prop, .$value, method = "spearman")) %>%
    map_dbl("estimate")

  ili_h3_5_24_pvalues <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$ili_age_5_24_prop, .$value, method = "spearman")) %>%
    map_dbl("p.value")

  ili_h3_25_64_predictors <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$ili_age_25_64_prop, .$value, method = "spearman")) %>%
    map_dbl("estimate")

  ili_h3_25_64_pvalues <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$ili_age_25_64_prop, .$value, method = "spearman")) %>%
    map_dbl("p.value")

  ili_h3_65_predictors <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$ili_age_65_prop, .$value, method = "spearman")) %>%
    map_dbl("estimate")

  ili_h3_65_pvalues <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$ili_age_65_prop, .$value, method = "spearman")) %>%
    map_dbl("p.value")

  epi_evol_df <- data.frame(
    ili_h3_0_4 = ili_h3_0_4_predictors,
    ili_h3_5_24 = ili_h3_5_24_predictors,
    ili_h3_25_64 = ili_h3_25_64_predictors,
    ili_h3_65 = ili_h3_65_predictors,
    ili_h3_0_4_p = ili_h3_0_4_pvalues,
    ili_h3_5_24_p = ili_h3_5_24_pvalues,
    ili_h3_25_64_p = ili_h3_25_64_pvalues,
    ili_h3_65_p = ili_h3_65_pvalues
  )
  epi_evol_df$evol_metric <- rownames(epi_evol_df)

  epi_evol_df_long <- epi_evol_df %>%
    dplyr::select(-contains("_p")) %>%
    pivot_longer(cols = c(
      ili_h3_0_4,
      ili_h3_5_24,
      ili_h3_25_64,
      ili_h3_65
    ), names_to = "epi_metric", values_to = "cor")

  epi_evol_df_long2 <- epi_evol_df %>%
    dplyr::select(evol_metric, contains("_p")) %>%
    dplyr::rename(c(
      ili_h3_0_4 = ili_h3_0_4_p,
      ili_h3_5_24 = ili_h3_5_24_p,
      ili_h3_25_64 = ili_h3_25_64_p,
      ili_h3_65 = ili_h3_65_p
    )) %>%
    pivot_longer(
      cols = c(
        ili_h3_0_4,
        ili_h3_5_24,
        ili_h3_25_64,
        ili_h3_65
      ),
      names_to = "epi_metric", values_to = "pvalue"
    )

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
epi_evol_df_combined %>% filter(pvalue <= 0.1 & pvalue > 0.05)
epi_evol_df_combined %>% filter(pvalue <= 0.05)

epi_evol_df_combined <-
  epi_evol_df_combined %>%
  mutate(p_symbol = case_when(
    pvalue_round > 0.05 ~ "",
    pvalue_round <= 0.05 & pvalue_round > 0.01 ~ "*",
    pvalue_round <= 0.01 & pvalue_round > 0.001 ~ "**",
    pvalue_round <= 0.001 ~ "***"
  ))

epi_evol_df_combined$evol_metric <- as.factor(epi_evol_df_combined$evol_metric)
levels(epi_evol_df_combined$evol_metric)
levels(epi_evol_df_combined$evol_metric) <- c(
  "H3 RBS distance (t-1)",
  "H3 RBS distance (t-2)",
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
  "N2 epitope distance (N=53) (t-2)"
)
unique(epi_evol_df_combined$evol_metric)
levels(epi_evol_df_combined$evol_metric)
epi_evol_df_combined$evol_metric <- factor(epi_evol_df_combined$evol_metric, levels = c(
  "N2 epitope distance (N=53) (t-2)",
  "N2 epitope distance (N=53) (t-1)",
  "N2 non-epitope distance (t-2)",
  "N2 non-epitope distance (t-1)",
  "N2 epitope distance (N=223) (t-2)",
  "N2 epitope distance (N=223) (t-1)",
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
# epi_evol_df_combined$evol_metric = factor(epi_evol_df_combined$evol_metric,levels=rev(sort(levels(epi_evol_df_combined$evol_metric))))
epi_evol_df_combined$epi_metric <- as.factor(epi_evol_df_combined$epi_metric)
levels(epi_evol_df_combined$epi_metric)
levels(epi_evol_df_combined$epi_metric) <- c(
  "ILI Prop. 0-4 yrs", "ILI Prop. 25-64 yrs",
  "ILI Prop 5-24 yrs", "ILI Prop 65+ yrs"
)
epi_evol_df_combined$epi_metric <- factor(epi_evol_df_combined$epi_metric,
  levels = c(
    "ILI Prop. 0-4 yrs", "ILI Prop 5-24 yrs", "ILI Prop. 25-64 yrs",
    "ILI Prop 65+ yrs"
  )
)

heat_p2 <- ggplot(data = epi_evol_df_combined %>% filter(grepl("ILI", epi_metric)) %>%
  filter(!grepl("non-epitope", evol_metric))) +
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
# save_plot(heat_p2, filename = "figures/Fig6_sup_fig1_age_epi_evol_heatmap_alternate.png", base_width = 12, base_height = 10)
save_plot(heat_p2, filename = "figures/Fig6_sup_fig1_age_epi_evol_heatmap_alternate.pdf",dpi=300, base_width = 12, base_height = 10)

####################################################
## Regression models
####################################################

sum_df <-
  epi_red2 %>%
  dplyr::select(
    contains("ILI"), season, dom_type,
    HA_wolf_lag2, HA_titer_tree_lag2,
    NA_bhatt_ep_lag1,
    NA_bhatt_ep_lag2
  ) %>%
  group_by(
    season, dom_type, HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1,
    NA_bhatt_ep_lag2
  ) %>%
  dplyr::summarise(
    prop.0_4.mean = ci(ili_age_0_4_prop, na.rm = T)[1],
    prop.0_4.lowCI = ci(ili_age_0_4_prop, na.rm = T)[2],
    prop.0_4.hiCI = ci(ili_age_0_4_prop, na.rm = T)[3],
    prop.0_4.sd = ci(ili_age_0_4_prop, na.rm = T)[4],
    prop.5_24.mean = ci(ili_age_5_24_prop, na.rm = T)[1],
    prop.5_24.lowCI = ci(ili_age_5_24_prop, na.rm = T)[2],
    prop.5_24.hiCI = ci(ili_age_5_24_prop, na.rm = T)[3],
    prop.5_24.sd = ci(ili_age_5_24_prop, na.rm = T)[4],
    prop.25_64.mean = ci(ili_age_25_64_prop, na.rm = T)[1],
    prop.25_64.lowCI = ci(ili_age_25_64_prop, na.rm = T)[2],
    prop.25_64.hiCI = ci(ili_age_25_64_prop, na.rm = T)[3],
    prop.25_64.sd = ci(ili_age_25_64_prop, na.rm = T)[4],
    prop.65.mean = ci(ili_age_65_prop, na.rm = T)[1],
    prop.65.lowCI = ci(ili_age_65_prop, na.rm = T)[2],
    prop.65.hiCI = ci(ili_age_65_prop, na.rm = T)[3],
    prop.65.sd = ci(ili_age_65_prop, na.rm = T)[4]
  ) %>%
  ungroup()
head(sum_df)
sum_df <- sum_df %>% tidyr::separate(col = "season", sep = "-", remove = F, into = c("year1", "year2"))
sum_df <- sum_df %>%
  mutate(h1n1_type = if_else(year1 < 2010, "seasonal_h1n1", "pdm_h1n1")) %>%
  mutate(dom_type2 = case_when(
    year1 < 2009 & dom_type == "H1" ~ "H1",
    year1 > 2009 & dom_type == "H1" ~ "H1pdm",
    dom_type == "H3" ~ "H3",
    dom_type == "co-circ" ~ "H3/H1pdm"
  ))

sum_df[!complete.cases(sum_df), ]
unique(sum_df$season)
nrow(sum_df)
####################################################
## HA epitope distance
####################################################

####################################################
## 0-4 yrs
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))
names(sumdf2)
y <- sumdf2$prop.0_4.mean
x <- sumdf2$HA_wolf_lag2
linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))

model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.0_4.mean ~ HA_wolf_lag2, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.0_4.mean ~ HA_wolf_lag2, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.1f", pvalue)
  )
labels

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2,prop.0_4.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.0_4.mean, value)) %>%
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
        betareg(prop.0_4.mean ~ value, .) %>%
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
ep_0_4 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.0_4.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>% rename(prop.0_4.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = prop.0_4.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_wolf_lag2,
      ymin = prop.0_4.lowCI, ymax = prop.0_4.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = prop.0_4.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("H3 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("Prop. 0 - 4 years") +
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
    y = 0.36, x = 0, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_0_4

####################################################
## 5-24 years
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.5_24.mean
x <- sumdf2$HA_wolf_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))

model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.5_24.mean ~ HA_wolf_lag2, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.5_24.mean ~ HA_wolf_lag2, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.1f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.1f", pvalue)
  )
labels

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, prop.5_24.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.5_24.mean, value)) %>%
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
        betareg(prop.5_24.mean ~ value, .) %>%
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


ep_5_24 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.5_24.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>% rename(prop.5_24.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = prop.5_24.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_wolf_lag2,
      ymin = prop.5_24.lowCI, ymax = prop.5_24.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = prop.5_24.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("H3 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("Prop. 5 - 24 years") +
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
    y = 0.45, x = 0, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_5_24

####################################################
## 25-64 years
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.25_64.mean
x <- sumdf2$HA_wolf_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.25_64.mean ~ HA_wolf_lag2, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.25_64.mean ~ HA_wolf_lag2, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.1f", pvalue)
  )
labels

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2,prop.25_64.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.25_64.mean, value)) %>%
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
        betareg(prop.25_64.mean ~ value, .) %>%
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


ep_25_64 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.25_64.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>% rename(prop.25_64.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = prop.25_64.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_wolf_lag2,
      ymin = prop.25_64.lowCI, ymax = prop.25_64.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = prop.25_64.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("H3 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("Prop. 25 - 64 years") +
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
    y = 0.525, x = 0, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_25_64

####################################################
## 65+ years
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.65.mean
x <- sumdf2$HA_wolf_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))

model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.65.mean ~ HA_wolf_lag2, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.65.mean ~ HA_wolf_lag2, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
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
  dplyr::select(HA_wolf_lag2, prop.65.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.65.mean, value)) %>%
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
        betareg(prop.65.mean ~ value, .) %>%
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

ep_65 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.65.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>% rename(prop.65.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = prop.65.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_wolf_lag2,
      ymin = prop.65.lowCI, ymax = prop.65.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = prop.65.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("H3 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("Prop. \u2265 65 years") +
  theme(legend.position = c(0.6, 0.8)) +
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
    y = 0.15, x = 0, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_65


all_epi_wolf_ep <- plot_grid(ep_0_4 + theme(legend.position = "none"),
  ep_5_24 + theme(legend.position = "none"),
  ep_25_64 + theme(legend.position = "none"),
  ep_65 + theme(legend.position = "none"),
  rel_widths = c(1, 1, 1, 1),
  nrow = 1
)
all_epi_wolf_ep

####################################################
## NA bhatt epitope distance (t-2)
####################################################

####################################################
## 0-4 yrs
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.0_4.mean
x <- sumdf2$NA_bhatt_ep_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))

model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)


fit_glm_on_bootstrap <- function(split) {
  betareg(prop.0_4.mean ~ NA_bhatt_ep_lag2, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.0_4.mean ~ NA_bhatt_ep_lag2, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
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
  dplyr::select(NA_bhatt_ep_lag2, prop.0_4.mean) %>%
  pivot_longer(cols = c(NA_bhatt_ep_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.0_4.mean, value)) %>%
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
        betareg(prop.0_4.mean ~ value, .) %>%
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
ep_0_4 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.0_4.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag2") %>% rename(prop.0_4.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag2, y = prop.0_4.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_bhatt_ep_lag2,
      ymin = prop.0_4.lowCI, ymax = prop.0_4.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag2, y = prop.0_4.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("Prop. 0 - 4 years") +
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
    y = 0.36, x = 0, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_0_4

####################################################
## 5-24
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.5_24.mean
x <- sumdf2$NA_bhatt_ep_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))

model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.5_24.mean ~ NA_bhatt_ep_lag2, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.5_24.mean ~ NA_bhatt_ep_lag2, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
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
  dplyr::select(NA_bhatt_ep_lag2, prop.5_24.mean) %>%
  pivot_longer(cols = c(NA_bhatt_ep_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.5_24.mean, value)) %>%
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
        betareg(prop.5_24.mean ~ value, .) %>%
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
ep_5_24 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.5_24.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag2") %>% rename(prop.5_24.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag2, y = prop.5_24.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_bhatt_ep_lag2,
      ymin = prop.5_24.lowCI, ymax = prop.5_24.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag2, y = prop.5_24.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("Prop. 5 - 24 years") +
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
    y = 0.45, x = 0, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_5_24

####################################################
## 25-64
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.25_64.mean
x <- sumdf2$NA_bhatt_ep_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.25_64.mean ~ NA_bhatt_ep_lag2, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.25_64.mean ~ NA_bhatt_ep_lag2, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
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
  dplyr::select(NA_bhatt_ep_lag2, prop.25_64.mean) %>%
  pivot_longer(cols = c(NA_bhatt_ep_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.25_64.mean, value)) %>%
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
        betareg(prop.25_64.mean ~ value, .) %>%
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


ep_25_64 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.25_64.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag2") %>% rename(prop.25_64.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag2, y = prop.25_64.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_bhatt_ep_lag2,
      ymin = prop.25_64.lowCI, ymax = prop.25_64.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag2, y = prop.25_64.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("Prop. 25 - 64 years") +
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
    y = 0.525, x = 0, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_25_64

####################################################
## 65+
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.65.mean
x <- sumdf2$NA_bhatt_ep_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)
names(sumdf2)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.65.mean ~ NA_bhatt_ep_lag2, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.65.mean ~ NA_bhatt_ep_lag2, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
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
  dplyr::select(NA_bhatt_ep_lag2,  prop.65.mean) %>%
  pivot_longer(cols = c(NA_bhatt_ep_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.65.mean, value)) %>%
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
        betareg(prop.65.mean ~ value, .) %>%
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

ep_65 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.65.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag2") %>% rename(prop.65.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag2, y = prop.65.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_bhatt_ep_lag2,
      ymin = prop.65.lowCI, ymax = prop.65.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag2, y = prop.65.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("Prop. \u2265 65 years") +
  theme(legend.position = c(0.6, 0.8)) +
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
    y = 0.15, x = 0, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_65


all_epi_bhatt_ep <- plot_grid(ep_0_4 + theme(legend.position = "none"),
  ep_5_24 + theme(legend.position = "none"),
  ep_25_64 + theme(legend.position = "none"),
  ep_65 + theme(legend.position = "none"),
  rel_widths = c(1, 1, 1, 1),
  nrow = 1
)
all_epi_bhatt_ep

####################################################
## NA bhatt epitope distance (t-1)
####################################################

####################################################
## 0-4 yrs
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.0_4.mean
x <- sumdf2$NA_bhatt_ep_lag1

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.0_4.mean ~ NA_bhatt_ep_lag1, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.0_4.mean ~ NA_bhatt_ep_lag1, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.1f", pvalue)
  )
labels

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(NA_bhatt_ep_lag1, prop.0_4.mean) %>%
  pivot_longer(cols = c(NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.0_4.mean, value)) %>%
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
        betareg(prop.0_4.mean ~ value, .) %>%
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
ep_0_4 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.0_4.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(prop.0_4.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = prop.0_4.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_bhatt_ep_lag1,
      ymin = prop.0_4.lowCI, ymax = prop.0_4.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = prop.0_4.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-1)")) +
  ylab("Prop. 0 - 4 years") +
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
    y = 0.36, x = 0, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_0_4

####################################################
## 5-24
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.5_24.mean
x <- sumdf2$NA_bhatt_ep_lag1

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.5_24.mean ~ NA_bhatt_ep_lag1, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.5_24.mean ~ NA_bhatt_ep_lag1, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
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
  dplyr::select( NA_bhatt_ep_lag1, prop.5_24.mean) %>%
  pivot_longer(cols = c(NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.5_24.mean, value)) %>%
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
        betareg(prop.5_24.mean ~ value, .) %>%
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
ep_5_24 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.5_24.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(prop.5_24.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = prop.5_24.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_bhatt_ep_lag1,
      ymin = prop.5_24.lowCI, ymax = prop.5_24.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = prop.5_24.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-1)")) +
  ylab("Prop. 5 - 24 years") +
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
    y = 0.45, x = 0, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_5_24

####################################################
## 25-64
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.25_64.mean
x <- sumdf2$NA_bhatt_ep_lag1

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.25_64.mean ~ NA_bhatt_ep_lag1, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.25_64.mean ~ NA_bhatt_ep_lag1, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
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
  dplyr::select(NA_bhatt_ep_lag1, HA_titer_tree_lag2, NA_bhatt_ep_lag1, prop.25_64.mean) %>%
  pivot_longer(cols = c(NA_bhatt_ep_lag1, HA_titer_tree_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.25_64.mean, value)) %>%
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
        betareg(prop.25_64.mean ~ value, .) %>%
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


ep_25_64 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.25_64.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(prop.25_64.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = prop.25_64.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_bhatt_ep_lag1,
      ymin = prop.25_64.lowCI, ymax = prop.25_64.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = prop.25_64.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-1)")) +
  ylab("Prop. 25 - 64 years") +
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
    y = 0.525, x = 0, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_25_64

####################################################
## 65+
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.65.mean
x <- sumdf2$NA_bhatt_ep_lag1

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)
names(sumdf2)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.65.mean ~ NA_bhatt_ep_lag1, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.65.mean ~ NA_bhatt_ep_lag1, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
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
  dplyr::select( NA_bhatt_ep_lag1, prop.65.mean) %>%
  pivot_longer(cols = c(NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.65.mean, value)) %>%
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
        betareg(prop.65.mean ~ value, .) %>%
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

ep_65 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.65.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(prop.65.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = prop.65.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_bhatt_ep_lag1,
      ymin = prop.65.lowCI, ymax = prop.65.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = prop.65.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-1)")) +
  ylab("Prop. \u2265 65 years") +
  theme(legend.position = c(0.6, 0.8)) +
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
    y = 0.15, x = 0, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_65

all_epi_bhatt_ep_lag1 <- plot_grid(ep_0_4 + theme(legend.position = "none"),
  ep_5_24 + theme(legend.position = "none"),
  ep_25_64 + theme(legend.position = "none"),
  ep_65 + theme(legend.position = "none"),
  rel_widths = c(1, 1, 1, 1),
  nrow = 1
)
all_epi_bhatt_ep_lag1

####################################################
## HI titer distance
####################################################

####################################################
## 0-4 yrs
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_titer_tree_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))
names(sumdf2)
y <- sumdf2$prop.0_4.mean
x <- sumdf2$HA_titer_tree_lag2
linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))

model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.0_4.mean ~ HA_titer_tree_lag2, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.0_4.mean ~ HA_titer_tree_lag2, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.1f", pvalue)
  )
labels

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_titer_tree_lag2, prop.0_4.mean) %>%
  pivot_longer(cols = c(HA_titer_tree_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.0_4.mean, value)) %>%
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
        betareg(prop.0_4.mean ~ value, .) %>%
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
ep_0_4 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.0_4.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "HA_titer_tree_lag2") %>% rename(prop.0_4.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = prop.0_4.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_titer_tree_lag2,
      ymin = prop.0_4.lowCI, ymax = prop.0_4.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = prop.0_4.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("HI titer distance (" ~ italic(t) ~ "-2)")) +
  ylab("Prop. 0 - 4 years") +
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
    y = 0.36, x = 0, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_0_4

####################################################
## 5-24 years
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_titer_tree_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.5_24.mean
x <- sumdf2$HA_titer_tree_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))

model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.5_24.mean ~ HA_titer_tree_lag2, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.5_24.mean ~ HA_titer_tree_lag2, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.1f", pvalue)
  )
labels

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_titer_tree_lag2,  prop.5_24.mean) %>%
  pivot_longer(cols = c(HA_titer_tree_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.5_24.mean, value)) %>%
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
        betareg(prop.5_24.mean ~ value, .) %>%
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


ep_5_24 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.5_24.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "HA_titer_tree_lag2") %>% rename(prop.5_24.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = prop.5_24.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_titer_tree_lag2,
      ymin = prop.5_24.lowCI, ymax = prop.5_24.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = prop.5_24.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("HI titer distance (" ~ italic(t) ~ "-2)")) +
  ylab("Prop. 5 - 24 years") +
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
    y = 0.45, x = -0.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_5_24

####################################################
## 25-64 years
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_titer_tree_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.25_64.mean
x <- sumdf2$HA_titer_tree_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.25_64.mean ~ HA_titer_tree_lag2, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.25_64.mean ~ HA_titer_tree_lag2, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.1f", pvalue)
  )
labels

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_titer_tree_lag2, prop.25_64.mean) %>%
  pivot_longer(cols = c(HA_titer_tree_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.25_64.mean, value)) %>%
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
        betareg(prop.25_64.mean ~ value, .) %>%
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


ep_25_64 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.25_64.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "HA_titer_tree_lag2") %>% rename(prop.25_64.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = prop.25_64.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_titer_tree_lag2,
      ymin = prop.25_64.lowCI, ymax = prop.25_64.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = prop.25_64.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("HI titer distance (" ~ italic(t) ~ "-2)")) +
  ylab("Prop. 25 - 64 years") +
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
    y = 0.525, x = -0.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_25_64

####################################################
## 65+ years
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_titer_tree_lag2:NA_bhatt_ep_lag2), ~ scale_this(.x))

y <- sumdf2$prop.65.mean
x <- sumdf2$HA_titer_tree_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
beta.model <- betareg(y ~ x, link = c("logit"))

model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, beta.model, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(prop.65.mean ~ HA_titer_tree_lag2, analysis(split))
}

# boot_models <-
#   boots %>%
#   mutate(
#     model = map(splits, fit_glm_on_bootstrap),
#     coef_info = map(model, tidy)
#   )
#
# boot_coefs <-
#   boot_models %>%
#   unnest(coef_info)

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(prop.65.mean ~ HA_titer_tree_lag2, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
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
  dplyr::select(HA_titer_tree_lag2,  prop.65.mean) %>%
  pivot_longer(cols = c(HA_titer_tree_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(prop.65.mean, value)) %>%
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
        betareg(prop.65.mean ~ value, .) %>%
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

ep_65 <- ggplot() +
  geom_ribbon(
    aes(
      x = prop.65.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "HA_titer_tree_lag2") %>% rename(prop.65.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = prop.65.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_titer_tree_lag2,
      ymin = prop.65.lowCI, ymax = prop.65.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = prop.65.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("HI titer distance (" ~ italic(t) ~ "-2)")) +
  ylab("Prop. \u2265 65 years") +
  theme(legend.position = c(0.6, 0.8)) +
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
    y = 0.15, x = -0.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_65


all_epi_HI_titer <- plot_grid(ep_0_4 + theme(legend.position = "none"),
                             ep_5_24 + theme(legend.position = "none"),
                             ep_25_64 + theme(legend.position = "none"),
                             ep_65 + theme(legend.position = "none"),
                             rel_widths = c(1, 1, 1, 1),
                             nrow = 1
)
all_epi_HI_titer


####################################################
## Combine all plots
####################################################

all_measures <- plot_grid(
  # all_epi_HI_titer+ theme(legend.position = "none"),
  all_epi_wolf_ep + theme(legend.position = "none"),
  all_epi_bhatt_ep_lag1 + theme(legend.position = "none"),
  all_epi_bhatt_ep + theme(legend.position = "none"),
  nrow = 3, labels = "AUTO"
)
all_measures
epi_leg <- cowplot::get_plot_component(ep_65+
                                         guides(color = "none") +
                                         theme(
                                           legend.position = "bottom",
                                           legend.direction = "horizontal",
                                           legend.justification = "center",
                                           legend.box.just = "bottom",
                                           legend.text = element_text(size = 14),
                                           legend.title = element_text(size = 16)
                                         ), 
                                       'guide-box-bottom', return_all = TRUE)
cowplot::ggdraw(epi_leg)
all_epi_leg_all_measures <- plot_grid(all_measures, epi_leg, nrow = 2, rel_heights = c(3, 0.2))
all_epi_leg_all_measures
# save_plot(all_epi_leg_all_measures, filename = "figures/Fig6_all_antigenic_measures_vs_age_patterns_north_amer_build_ha_ep.png", base_width = 18, base_height = 14)
save_plot(all_epi_leg_all_measures, filename = "figures/Fig6_all_antigenic_measures_vs_age_patterns_north_amer_build_ha_ep.pdf", dpi = 300, base_width = 18, base_height = 14)

all_measures <- plot_grid(
  all_epi_HI_titer+ theme(legend.position = "none"),
  # all_epi_wolf_ep + theme(legend.position = "none"),
  all_epi_bhatt_ep_lag1 + theme(legend.position = "none"),
  all_epi_bhatt_ep + theme(legend.position = "none"),
  nrow = 3, labels = "AUTO"
)
all_measures
all_epi_leg_all_measures <- plot_grid(all_measures, epi_leg, nrow = 2, rel_heights = c(3, 0.2))
all_epi_leg_all_measures
# save_plot(all_epi_leg_all_measures, filename = "figures/all_antigenic_measures_vs_age_patterns_north_amer_build_hi_titer.png", base_width = 18, base_height = 14)
# save_plot(all_epi_leg_all_measures, filename = "figures/all_antigenic_measures_vs_age_patterns_north_amer_build_hi_titer.pdf", dpi = 300, base_width = 18, base_height = 14)
