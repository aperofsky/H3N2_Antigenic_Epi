## load packages
list.of.packages <- c(
  "dplyr", "ggplot2", "cowplot", "tidyr", "broom", "gmodels", "rsample", "metan",
  "betareg", "purrr", "tidymodels", "MuMIn", "readr", "RColorBrewer"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

theme_set(theme_cowplot(font_size = 16))

########################################################################################################
## Associations between evolutionary indicators and A/H3N2 epidemic metrics
########################################################################################################
## load data
load("data/antigenic_epi_north_amer_build_for_lasso_replicates.Rdata")
names(epi_red)

load("data/subtype_distribution_by_region_season.RData") # subtype_dist
subtype_dist <- subtype_dist %>%
  ungroup() %>%
  mutate(a_total = h3_total + h1_total) %>%
  mutate(
    h3_vs_h1 = h3_total / a_total,
    h3_dom = h3_total / (a_total + b_total)
  ) %>% # prop of h3 out of IAV%>%
  rename(season = season_description)
names(subtype_dist)

epi_red <- left_join(epi_red %>% dplyr::select(-h3_vs_h1, -iva_vs_ivb, -h3_vs_total_flu, -h1_vs_total_flu, -ivb_vs_iva, -total_a),
  subtype_dist %>% dplyr::select(season, region, h3_dom, h3_vs_h1),
  by = c("season", "region")
) %>%
  distinct() %>%
  filter(!(season %in% c("2009-2010", "2019-2020", "2020-2021", "2021-2022")))

epi_long <- epi_red %>%
  pivot_longer(
    cols = c(HA_mean_lbi_lag1:na_lbi_shannon),
    names_to = "evol_metrics",
    values_to = "value"
  ) %>%
  distinct() %>%
  ungroup()

h3_Rt_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ lm(H3_max_Rt ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_Rt_predictors, decreasing = T)[1:5]

h3_peak_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ lm(H3_max_intensity ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_peak_predictors, decreasing = T)[1:5]

h3_size_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ lm(H3_cum_intensity ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_size_predictors, decreasing = T)[1:5]

h3_sh_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ lm(H3_shannon_entropy_res ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_sh_predictors, decreasing = T)[1:5]

####################################################
## Heatmap, Evol indicators vs Epi metrics
####################################################

set.seed(27)
boots <- bootstraps(epi_red, times = 1000, apparent = TRUE)

datalist <- list()
for (i in boots$id) {
  # i=boots$id[1]
  epi_long <- boots %>%
    filter(id == i) %>%
    pull(splits) %>%
    as.data.frame() %>%
    dplyr::select(-contains(c("std", "usa", "shannon_lag2", "lbi_lag2"))) %>%
    group_by(season, dom_type) %>%
    summarize_at(vars(
      h3_vs_h1, h3_dom, H3_max_intensity,
      H3_shannon_entropy_res, H3_cum_intensity, H3_max_Rt,
      HA_mean_lbi_lag1:na_lbi_shannon
    ), mean) %>%
    ungroup() %>%
    pivot_longer(
      cols = contains(c("lag", "lbi")),
      names_to = "evol_metrics",
      values_to = "value"
    ) %>%
    filter(!is.na(value))
  sort(unique(epi_long$evol_metrics))

  h3_dom_predictors <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$h3_dom, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("estimate")

  # h3_dom_predictors = epi_long %>%
  #   # split(f=list(.$replicate,.$evol_metrics))%>%
  #   split(.$evol_metrics) %>% #
  #   map(~ cor.test(.$h3_vs_h1, scale(.$value,center=T)),method="spearman") %>%
  #   # map(summary) %>%
  #   map_dbl("estimate")

  h3_dom_pvalues <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$h3_dom, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("p.value")

  h3_max_predictors <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$H3_max_intensity, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("estimate")

  h3_max_pvalues <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$H3_max_intensity, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("p.value")

  h3_cum_predictors <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$H3_cum_intensity, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("estimate")

  h3_cum_pvalues <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$H3_cum_intensity, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("p.value")

  h3_R0_predictors <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$H3_max_Rt, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("estimate")

  h3_R0_pvalues <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$H3_max_Rt, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("p.value")

  h3_shannon_predictors <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$H3_shannon_entropy_res, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("estimate")

  h3_shannon_pvalues <- epi_long %>%
    split(.$evol_metrics) %>% #
    map(~ cor.test(.$H3_shannon_entropy_res, scale(.$value, center = T)), method = "spearman") %>%
    map_dbl("p.value")

  epi_evol_df <- data.frame(
    h3_dom = h3_dom_predictors,
    h3_max = h3_max_predictors,
    h3_cum = h3_cum_predictors,
    h3_R0 = h3_R0_predictors,
    h3_shannon = h3_shannon_predictors,
    h3_dom_p = h3_dom_pvalues,
    h3_max_p = h3_max_pvalues,
    h3_cum_p = h3_cum_pvalues,
    h3_R0_p = h3_R0_pvalues,
    h3_shannon_p = h3_shannon_pvalues
  )
  epi_evol_df$evol_metric <- rownames(epi_evol_df)

  epi_evol_df_long <- epi_evol_df %>%
    dplyr::select(-c(h3_dom_p, h3_max_p, h3_cum_p, h3_R0_p, h3_shannon_p)) %>%
    pivot_longer(cols = c(h3_dom, h3_max, h3_cum, h3_R0, h3_shannon), names_to = "epi_metric", values_to = "cor")

  epi_evol_df_long2 <- epi_evol_df %>%
    dplyr::select(-c(h3_dom, h3_max, h3_cum, h3_R0, h3_shannon)) %>%
    dplyr::rename(c(h3_dom = h3_dom_p, h3_max = h3_max_p, h3_cum = h3_cum_p, h3_R0 = h3_R0_p, h3_shannon = h3_shannon_p)) %>%
    pivot_longer(cols = c(h3_dom, h3_max, h3_cum, h3_R0, h3_shannon), names_to = "epi_metric", values_to = "pvalue")

  epi_combined <- left_join(epi_evol_df_long, epi_evol_df_long2, by = c("evol_metric", "epi_metric"))
  epi_combined$replicate <- i
  # epi_combined %>% filter(pvalue<0.05)
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
epi_evol_df_combined %>% filter(pvalue <= 0.1 & pvalue > 0.05)

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
unique(epi_evol_df_combined$evol_metric)
levels(epi_evol_df_combined$evol_metric)
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
# epi_evol_df_combined$evol_metric = factor(epi_evol_df_combined$evol_metric,levels=rev(sort(levels(epi_evol_df_combined$evol_metric))))
epi_evol_df_combined$epi_metric <- as.factor(epi_evol_df_combined$epi_metric)
levels(epi_evol_df_combined$epi_metric)
epi_evol_df_combined$epi_metric <- factor(epi_evol_df_combined$epi_metric, levels = c("h3_cum", "h3_max", "h3_R0", "h3_shannon", "h3_dom"))
levels(epi_evol_df_combined$epi_metric) <- c("Epidemic\nSize", "Peak\nIncidence", "Effective Rt", "Epidemic\nIntensity", "A(H3N2) Subtype\nDominance")
levels(epi_evol_df_combined$evol_metric)
head(epi_evol_df_combined)
unique(epi_evol_df_combined$evol_metric)

epi_evol_df_combined$evol_metric <- factor(epi_evol_df_combined$evol_metric, levels = rev(levels(epi_evol_df_combined$evol_metric)))

heat_p2 <- ggplot(data = epi_evol_df_combined) +
  facet_grid(evol_metric ~ epi_metric,
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
save_plot(heat_p2, filename = "figures/epi_evol_heatmap_alternate.png", base_width = 12, base_height = 14)

########################################################################################################
## Evolutionary indicators vs Epi Metrics
########################################################################################################

####################################################
## HI Titer Distance vs Epi Metrics
####################################################

sum_df <-
  epi_red %>%
  tidyr::replace_na(list(
    H1_cum_intensity = 0, H3_cum_intensity = 0,
    H3_cum_intensity = 0, IVB_cum_intensity = 0
  )) %>%
  group_by(
    season, dom_type, HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1,
    NA_mean_lbi, HA_mean_lbi, HA_mean_lbi_lag1,
    ha_lbi_shannon_lag1, ha_lbi_shannon, na_lbi_shannon
  ) %>%
  dplyr::summarise(
    H3.shannon.mean = ci(H3_shannon_entropy_res, na.rm = T)[1],
    H3.shannon.lowCI = ci(H3_shannon_entropy_res, na.rm = T)[2],
    H3.shannon.hiCI = ci(H3_shannon_entropy_res, na.rm = T)[3],
    H3.shannon.sd = ci(H3_shannon_entropy_res, na.rm = T)[4],
    H3.peak.mean = ci(H3_max_intensity, na.rm = T)[1],
    H3.peak.lowCI = ci(H3_max_intensity, na.rm = T)[2],
    H3.peak.hiCI = ci(H3_max_intensity, na.rm = T)[3],
    H3.peak.sd = ci(H3_max_intensity, na.rm = T)[4],
    H3.epi.size.mean = ci(H3_cum_intensity, na.rm = T)[1],
    H3.epi.size.lowCI = ci(H3_cum_intensity, na.rm = T)[2],
    H3.epi.size.hiCI = ci(H3_cum_intensity, na.rm = T)[3],
    H3.epi.size.sd = ci(H3_cum_intensity, na.rm = T)[4],
    H3.R0.mean = ci(H3_max_Rt, na.rm = T)[1],
    H3.R0.lowCI = ci(H3_max_Rt, na.rm = T)[2],
    H3.R0.hiCI = ci(H3_max_Rt, na.rm = T)[3],
    H3.R0.sd = ci(H3_max_Rt, na.rm = T)[4]
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

####################################################
## A/H3 epidemic size
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))
names(sumdf2)
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.epi.size.mean
x <- sumdf2$HA_titer_tree_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, rank = "BIC")
# use Gaussian GLM with identity link across everything so that model results are comparable

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)
names(sumdf2)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.epi.size.mean ~ HA_titer_tree_lag2, analysis(split), family = gaussian(link = "identity"), maxit = 150)
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

m1 <- glm(H3.epi.size.mean ~ HA_titer_tree_lag2, data = sumdf2, family = gaussian(link = "identity"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.17

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.epi.size.mean ~ HA_titer_tree_lag2, data = .x, family = gaussian(link = "identity"), maxit = 150)),
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

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1, H3.epi.size.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.epi.size.mean, value)) %>%
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
        glm(H3.epi.size.mean ~ value, ., family = gaussian(link = "identity"), maxit = 150) %>%
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

hi_epi_size <- ggplot() +
  geom_ribbon(aes(x = H3.epi.size.mean, ymin = gaussian(link = "identity")$linkinv(l), ymax = gaussian(link = "identity")$linkinv(u)),
    tbl_plot_data %>% filter(name == "HA_titer_tree_lag2") %>% rename(H3.epi.size.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = H3.epi.size.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "identity")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_titer_tree_lag2,
      ymin = H3.epi.size.lowCI, ymax = H3.epi.size.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = H3.epi.size.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab(expression("HI titer distance (" ~ italic(t) ~ "-2)")) +
  ylab("A(H3N2) Epidemic Size") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_text(
    y = 77, x = -2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  ) +
  background_grid(major = "xy", minor = "none")

hi_epi_size

####################################################
## A/H3 Scaled Peak Incidence
####################################################

sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.peak.mean
x <- sumdf2$HA_titer_tree_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")
## use Gaussian with identity link so that model results are comparable

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.peak.mean ~ HA_titer_tree_lag2, analysis(split), family = gaussian(link = "identity"), maxit = 150)
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

m1 <- glm(H3.peak.mean ~ HA_titer_tree_lag2, data = sumdf2, family = gaussian(link = "identity"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.17

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.peak.mean ~ HA_titer_tree_lag2, data = .x, family = gaussian(link = "identity"), maxit = 150)),
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

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1, H3.peak.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

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
        glm(H3.peak.mean ~ value, ., family = gaussian(link = "identity"), maxit = 150) %>%
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

hi_peak <- ggplot() +
  geom_ribbon(aes(x = H3.peak.mean, ymin = gaussian(link = "identity")$linkinv(l), ymax = gaussian(link = "identity")$linkinv(u)),
    tbl_plot_data %>% filter(name == "HA_titer_tree_lag2") %>% rename(H3.peak.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = H3.peak.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "identity")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_titer_tree_lag2,
      ymin = H3.peak.lowCI, ymax = H3.peak.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = H3.peak.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("HI titer distance (" ~ italic(t) ~ "-2)")) +
  ylab("A(H3N2) Peak Incidence") +
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
    y = 16.5, x = -2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  ) +
  background_grid(major = "xy", minor = "none")
hi_peak

####################################################
## A/H3 shannon entropy
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2000-2001", "2009-2010")))
## no shannon entropy measurements for 2000-2001 and 2009-2010 (no H3N2 circulation)
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.shannon.mean
x <- sumdf2$HA_titer_tree_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
beta.model <- betareg(y ~ x, link = c("logit"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, beta.model, rank = "BIC")
## Beta fits best

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(H3.shannon.mean ~ HA_titer_tree_lag2, analysis(split))
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

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(H3.shannon.mean ~ HA_titer_tree_lag2, data = .x)),
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

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1, H3.shannon.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.shannon.mean, value)) %>%
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
        betareg(H3.shannon.mean ~ value, .) %>%
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

hi_shannon <- ggplot() +
  geom_ribbon(aes(x = H3.shannon.mean, ymin = l, ymax = u),
    tbl_plot_data %>% filter(name == "HA_titer_tree_lag2") %>% rename(H3.shannon.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = H3.shannon.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_titer_tree_lag2,
      ymin = H3.shannon.lowCI, ymax = H3.shannon.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = H3.shannon.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("HI titer distance (" ~ italic(t) ~ "-2)")) +
  ylab("A(H3N2) Epidemic Intensity") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_text(
    y = 0.8, x = -2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  ) +
  background_grid(major = "xy", minor = "none")

hi_shannon

####################################################
## A/H3N2 Rt
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010")))
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.R0.mean
x <- sumdf2$HA_titer_tree_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")
## gaussian with inverse link

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.R0.mean ~ HA_titer_tree_lag2, analysis(split), family = gaussian(link = "inverse"), maxit = 150)
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

m1 <- glm(H3.R0.mean ~ HA_titer_tree_lag2, data = sumdf2, family = gaussian(link = "inverse"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.21

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.R0.mean ~ HA_titer_tree_lag2, data = .x, family = gaussian(link = "inverse"), maxit = 150)),
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

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1, H3.R0.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.R0.mean, value)) %>%
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
        glm(H3.R0.mean ~ value, ., family = gaussian(link = "inverse"), maxit = 150) %>%
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

hi_R0 <- ggplot() +
  geom_ribbon(aes(x = H3.R0.mean, ymin = gaussian(link = "inverse")$linkinv(l), ymax = gaussian(link = "inverse")$linkinv(u)),
    tbl_plot_data %>% filter(name == "HA_titer_tree_lag2") %>% rename(H3.R0.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = H3.R0.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "inverse")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_titer_tree_lag2,
      ymin = H3.R0.lowCI, ymax = H3.R0.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_titer_tree_lag2, y = H3.R0.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("HI titer distance (" ~ italic(t) ~ "-2)")) +
  ylab("A(H3N2) Effective Rt") +
  theme(legend.position = c(0.6, 0.8)) +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_text(
    y = 2, x = -2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  ) +
  background_grid(major = "xy", minor = "none")
hi_R0

all_hi_titer <- plot_grid(hi_epi_size + theme(legend.position = "none"),
  hi_peak + theme(legend.position = "none"),
  hi_R0 + theme(legend.position = "none"),
  hi_shannon + theme(legend.position = "none"),
  rel_widths = c(1, 1, 1, 1),
  nrow = 1
)
all_hi_titer

####################################################
## HA epitope distance
####################################################

####################################################
## A/H3 Epi Size
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.epi.size.mean
x <- sumdf2$HA_wolf_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, rank = "BIC")
# Gaussian with identity link to make results consistent across evolutionary indicators

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.epi.size.mean ~ HA_wolf_lag2, analysis(split), family = gaussian(link = "identity"), maxit = 150)
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

boot_aug <-
  boot_models %>%
  sample_n(1000) %>%
  mutate(augmented = map(model, augment)) %>%
  unnest(augmented)

m1 <- glm(H3.epi.size.mean ~ HA_wolf_lag2, data = sumdf2, family = gaussian(link = "identity"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.35

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.epi.size.mean ~ HA_wolf_lag2, data = .x, family = gaussian(link = "identity"), maxit = 150)),
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

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1, H3.epi.size.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.epi.size.mean, value)) %>%
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
        glm(H3.epi.size.mean ~ value, ., family = gaussian(link = "identity"), maxit = 150) %>%
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

ep_epi_size <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.epi.size.mean,
      ymin = gaussian(link = "identity")$linkinv(l),
      ymax = gaussian(link = "identity")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>% rename(H3.epi.size.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = H3.epi.size.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "identity")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_wolf_lag2,
      ymin = H3.epi.size.lowCI, ymax = H3.epi.size.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = H3.epi.size.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("H3 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("A(H3N2) Epidemic Size") +
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
    y = 78, x = -1.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_epi_size

####################################################
## A/H3 Scaled Peak Incidence; gaussian, identity link
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010")))
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.peak.mean
x <- sumdf2$HA_wolf_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, rank = "BIC")

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.peak.mean ~ HA_wolf_lag2, analysis(split), family = gaussian(link = "identity"), maxit = 150)
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

m1 <- glm(H3.peak.mean ~ HA_wolf_lag2, data = sumdf2, family = gaussian(link = "identity"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.39

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.peak.mean ~ HA_wolf_lag2, data = .x, family = gaussian(link = "identity"), maxit = 150)),
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

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, NA_bhatt_ep_lag1, H3.peak.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

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
        glm(H3.peak.mean ~ value, ., family = gaussian(link = "identity"), maxit = 150) %>%
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

ep_peak <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.peak.mean,
      ymin = gaussian(link = "identity")$linkinv(l),
      ymax = gaussian(link = "identity")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>% rename(H3.peak.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = H3.peak.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "identity")),
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
  ylab("A(H3N2) Peak Incidence") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
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
    y = 16.5, x = -1.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_peak

####################################################
## A/H3 shannon entropy
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010", "2000-2001")))
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.shannon.mean
x <- sumdf2$HA_wolf_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
beta.model <- betareg(y ~ x)
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(beta.model, linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")
## use Beta to be consistent across evolutionary indicators

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(H3.shannon.mean ~ HA_wolf_lag2, analysis(split))
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

m1 <- betareg(H3.shannon.mean ~ HA_wolf_lag2, data = sumdf2)
fam <- family(m1)
ilink <- fam$linkinv

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(H3.shannon.mean ~ HA_wolf_lag2, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, NA_bhatt_ep_lag1, H3.shannon.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.shannon.mean, value)) %>%
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
        betareg(H3.shannon.mean ~ value, .) %>%
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

ep_shannon <- ggplot() +
  geom_ribbon(aes(x = H3.shannon.mean, ymin = l, ymax = u),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>% rename(H3.shannon.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = H3.shannon.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_wolf_lag2,
      ymin = H3.shannon.lowCI, ymax = H3.shannon.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = H3.shannon.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("H3 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("A(H3N2) Epidemic Intensity") +
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
    y = 0.8, x = -1.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_shannon

####################################################
## A/H3 effective Rt
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010")))
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.R0.mean
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
  glm(H3.R0.mean ~ HA_wolf_lag2, analysis(split), family = gaussian(link = "inverse"), maxit = 150)
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

m1 <- glm(H3.R0.mean ~ HA_wolf_lag2, data = sumdf2, family = gaussian(link = "inverse"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.36

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.R0.mean ~ HA_wolf_lag2, data = .x, family = gaussian(link = "inverse"), maxit = 150)),
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

set.seed(27)
n_boot <- 1000
sumdf2 %>%
  dplyr::select(HA_wolf_lag2, NA_bhatt_ep_lag1, H3.R0.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.R0.mean, value)) %>%
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
        glm(H3.R0.mean ~ value, ., family = gaussian(link = "inverse"), maxit = 150) %>%
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

ep_R0 <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.R0.mean,
      ymin = gaussian(link = "inverse")$linkinv(l),
      ymax = gaussian(link = "inverse")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>% rename(H3.R0.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = H3.R0.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "inverse")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_wolf_lag2,
      ymin = H3.R0.lowCI, ymax = H3.R0.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_wolf_lag2, y = H3.R0.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("H3 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("A(H3N2) Effective Rt") +
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
    y = 2, x = -1.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ep_R0

all_epi_wolf_ep <- plot_grid(ep_epi_size + theme(legend.position = "none"),
  ep_peak + theme(legend.position = "none"),
  ep_R0 + theme(legend.position = "none"),
  ep_shannon + theme(legend.position = "none"),
  rel_widths = c(1, 1, 1, 1),
  nrow = 1
)
all_epi_wolf_ep

####################################################
## NA epitope distance
####################################################

####################################################
## A/H3 Epi Size
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.epi.size.mean
x <- sumdf2$NA_bhatt_ep_lag1

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, rank = "BIC")
## use Gaussian model with identity link so that results are consistent across evolutionary indicators

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)
boots

fit_glm_on_bootstrap <- function(split) {
  glm(H3.epi.size.mean ~ NA_bhatt_ep_lag1, analysis(split), family = gaussian(link = "identity"), maxit = 150)
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

m1 <- glm(H3.epi.size.mean ~ NA_bhatt_ep_lag1, data = sumdf2, family = gaussian(link = "identity"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.24

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.epi.size.mean ~ NA_bhatt_ep_lag1, data = .x, family = gaussian(link = "identity"), maxit = 150)),
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

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1, H3.epi.size.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.epi.size.mean, value)) %>%
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
        glm(H3.epi.size.mean ~ value, ., family = gaussian(link = "identity"), maxit = 150) %>%
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

na_ep_epi_size <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.epi.size.mean,
      ymin = gaussian(link = "identity")$linkinv(l),
      ymax = gaussian(link = "identity")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(H3.epi.size.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = H3.epi.size.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "identity")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_bhatt_ep_lag1,
      ymin = H3.epi.size.lowCI, ymax = H3.epi.size.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = H3.epi.size.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-1)")) +
  ylab("A(H3N2) Epidemic Size") +
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
    y = 78, x = -1.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_ep_epi_size

####################################################
## A/H3 Scaled Peak Incidence
####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.peak.mean
x <- sumdf2$NA_bhatt_ep_lag1

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")
## use Gaussian model with identity link so that results are consistent across evolutionary indicators

set.seed(26)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.peak.mean ~ NA_bhatt_ep_lag1, analysis(split), family = gaussian(link = "identity"), maxit = 150)
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

m1 <- glm(H3.peak.mean ~ NA_bhatt_ep_lag1, data = sumdf2, family = gaussian(link = "identity"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.31

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.peak.mean ~ NA_bhatt_ep_lag1, data = .x, family = gaussian(link = "identity"), maxit = 150)),
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

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1, H3.peak.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

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
        glm(H3.peak.mean ~ value, ., family = gaussian(link = "identity"), maxit = 150) %>%
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
  geom_ribbon(
    aes(
      x = H3.peak.mean,
      ymin = gaussian(link = "identity")$linkinv(l),
      ymax = gaussian(link = "identity")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(H3.peak.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = H3.peak.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "identity")),
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
    aes(x = NA_bhatt_ep_lag1, y = H3.peak.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-1)")) +
  ylab("A(H3N2) Peak Incidence") +
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
    y = 16.5, x = -1.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_ep_peak

####################################################
## A/H3 shannon entropy
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2000-2001", "2009-2010")))
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.shannon.mean
x <- sumdf2$NA_bhatt_ep_lag1

linear.model <- glm(y ~ x, family = gaussian())
beta.model <- betareg(y ~ x)
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, beta.model, rank = "BIC")
# Beta fits best

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(H3.shannon.mean ~ NA_bhatt_ep_lag1, analysis(split))
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

m1 <- betareg(H3.shannon.mean ~ NA_bhatt_ep_lag1, data = sumdf2)
fam <- family(m1)
ilink <- fam$linkinv

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(H3.shannon.mean ~ NA_bhatt_ep_lag1, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1, H3.shannon.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.shannon.mean, value)) %>%
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
        betareg(H3.shannon.mean ~ value, .) %>%
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

na_ep_shannon <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.shannon.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(H3.shannon.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = H3.shannon.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_bhatt_ep_lag1,
      ymin = H3.shannon.lowCI, ymax = H3.shannon.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = H3.shannon.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-1)")) +
  ylab("A(H3N2) Epidemic Intensity") +
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
    y = 0.81, x = -1.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_ep_shannon

####################################################
## A/H3 effective Rt
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010")))
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.R0.mean
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
  glm(H3.R0.mean ~ NA_bhatt_ep_lag1, analysis(split), family = gaussian(link = "inverse"), maxit = 150)
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

m1 <- glm(H3.R0.mean ~ NA_bhatt_ep_lag1, data = sumdf2, family = gaussian(link = "inverse"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.3

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.R0.mean ~ NA_bhatt_ep_lag1, data = .x, family = gaussian(link = "inverse"), maxit = 150)),
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

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1, H3.R0.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_titer_tree_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.R0.mean, value)) %>%
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
        glm(H3.R0.mean ~ value, ., family = gaussian(link = "inverse"), maxit = 150) %>%
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

na_ep_R0 <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.R0.mean,
      ymin = gaussian(link = "inverse")$linkinv(l),
      ymax = gaussian(link = "inverse")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(H3.R0.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = H3.R0.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "inverse"), maxit = 150),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_bhatt_ep_lag1,
      ymin = H3.R0.lowCI, ymax = H3.R0.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_bhatt_ep_lag1, y = H3.R0.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-1)")) +
  ylab("A/H3 Effective Rt") +
  theme(legend.position = "bottom") +
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
    y = 2, x = -1.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_ep_R0


all_epi_bhatt_ep <- plot_grid(na_ep_epi_size + theme(legend.position = "none"),
  na_ep_peak + theme(legend.position = "none"),
  na_ep_R0 + theme(legend.position = "none"),
  na_ep_shannon + theme(legend.position = "none"),
  rel_widths = c(1, 1, 1, 1),
  nrow = 1
)
all_epi_bhatt_ep

########################################################################################################
## combine HA epitope, NA epitope, and HI titer results
########################################################################################################
all_measures <- plot_grid(
  all_epi_wolf_ep + theme(legend.position = "none"),
  all_epi_bhatt_ep + theme(legend.position = "none"),
  all_hi_titer + theme(legend.position = "none"),
  nrow = 3, labels = "AUTO"
)
epi_leg <- get_legend(na_ep_R0 + theme(
  legend.direction = "horizontal",
  legend.justification = "center",
  legend.box.just = "bottom",
  legend.position = "bottom",
  legend.title.align = 0.5
))

all_epi_leg_all_measures <- plot_grid(all_measures, epi_leg, nrow = 2, rel_heights = c(3, 0.1))
all_epi_leg_all_measures

## Figure 3
save_plot(all_epi_leg_all_measures, filename = "figures/all_antigenic_measures_vs_H3_epi_metrics_north_amer_build.png", base_width = 16, base_height = 13)

####################################################
## NA LBI
####################################################

####################################################
## N2 LBI vs A/H3 shannon entropy
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010", "2000-2001")))
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.shannon.mean
x <- sumdf2$na_lbi_shannon

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
beta.model <- betareg(y ~ x)
model.sel(beta.model, linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")
## use Beta for shannon entropy

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  betareg(H3.shannon.mean ~ na_lbi_shannon, analysis(split))
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

m1 <- betareg(H3.shannon.mean ~ na_lbi_shannon, data = sumdf2)
fam <- family(m1)
ilink <- fam$linkinv

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(H3.shannon.mean ~ na_lbi_shannon, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.3f", pvalue)
  )

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, HA_titer_tree_lag2, na_lbi_shannon, H3.shannon.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_titer_tree_lag2, na_lbi_shannon)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.shannon.mean, value)) %>%
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
        betareg(H3.shannon.mean ~ value, .) %>%
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

na_lbi_shannon <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.shannon.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "na_lbi_shannon") %>% rename(H3.shannon.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = na_lbi_shannon, y = H3.shannon.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = na_lbi_shannon,
      ymin = H3.shannon.lowCI, ymax = H3.shannon.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = na_lbi_shannon, y = H3.shannon.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("N2 LBI Diversity") +
  ylab("A(H3N2) Epidemic Intensity") +
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
    y = 0.81, x = -1.7, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_lbi_shannon
####################################################
## N2 LBI vs A/H3 R0
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010")))
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.R0.mean
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
  glm(H3.R0.mean ~ na_lbi_shannon, analysis(split), family = gaussian(link = "identity"))
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

m1 <- glm(H3.R0.mean ~ na_lbi_shannon, data = sumdf2, family = gaussian(link = "identity"))
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.46

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.R0.mean ~ na_lbi_shannon, data = .x, family = gaussian(link = "identity"))),
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

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, HA_titer_tree_lag2, na_lbi_shannon, H3.R0.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_titer_tree_lag2, na_lbi_shannon)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.R0.mean, value)) %>%
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
        glm(H3.R0.mean ~ value, ., family = gaussian(link = "identity")) %>%
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

na_lbi_R0 <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.R0.mean,
      ymin = gaussian(link = "identity")$linkinv(l),
      ymax = gaussian(link = "identity")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "na_lbi_shannon") %>% rename(H3.R0.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = na_lbi_shannon, y = H3.R0.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "identity")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = na_lbi_shannon,
      ymin = H3.R0.lowCI, ymax = H3.R0.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = na_lbi_shannon, y = H3.R0.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab("N2 LBI Diversity") +
  ylab("A/H3 Effective Rt") +
  theme(legend.position = "bottom") +
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
    y = 2, x = -1.7, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_lbi_R0

all_epi_na_lbi <- plot_grid(
  na_lbi_R0 + theme(legend.position = "none"),
  na_lbi_shannon + theme(legend.position = "none"),
  rel_widths = c(1, 1),
  nrow = 1
)
all_epi_na_lbi

####################################################
## HA LBI
####################################################
####################################################
## H3 LBI vs A/H3 shannon entropy
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010", "2000-2001")))
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.shannon.mean
x <- sumdf2$ha_lbi_shannon_lag1

linear.model <- glm(y ~ x, family = gaussian())
beta.model <- betareg(y ~ x)
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
gamma.model3 <- glm(y ~ x, family = Gamma(link = "identity"))
model.sel(beta.model, linear.model, log.model, inv.model, gamma.model, gamma.model2, gamma.model3, rank = "BIC")
## use Beta

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)
boots
names(sumdf2)

fit_glm_on_bootstrap <- function(split) {
  betareg(H3.shannon.mean ~ ha_lbi_shannon_lag1, analysis(split))
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

m1 <- betareg(H3.shannon.mean ~ ha_lbi_shannon_lag1, data = sumdf2)
fam <- family(m1)
ilink <- fam$linkinv
ilink

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(H3.shannon.mean ~ ha_lbi_shannon_lag1, data = .x, link = "logit")),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.3f", pvalue)
  )

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, HA_titer_tree_lag2, ha_lbi_shannon_lag1, H3.shannon.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_titer_tree_lag2, ha_lbi_shannon_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.shannon.mean, value)) %>%
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
        betareg(H3.shannon.mean ~ value, .) %>%
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

ha_lbi_shannon <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.shannon.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "ha_lbi_shannon_lag1") %>% rename(H3.shannon.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = ha_lbi_shannon_lag1, y = H3.shannon.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = ha_lbi_shannon_lag1,
      ymin = H3.shannon.lowCI, ymax = H3.shannon.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = ha_lbi_shannon_lag1, y = H3.shannon.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab(expression("H3 LBI Diversity (" ~ italic(t) ~ "-1)")) +
  ylab("A(H3N2) Epidemic Intensity") +
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
    y = 0.85, x = -2.7, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ha_lbi_shannon

####################################################
## H3 LBI vs A/H3 effective Rt
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010")))
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_wolf_lag2:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.R0.mean
x <- sumdf2$ha_lbi_shannon_lag1

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
  glm(H3.R0.mean ~ ha_lbi_shannon_lag1, analysis(split), family = gaussian(link = "identity"))
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

m1 <- glm(H3.R0.mean ~ ha_lbi_shannon_lag1, data = sumdf2, family = gaussian(link = "identity"))
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.5

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.R0.mean ~ ha_lbi_shannon_lag1, data = .x, family = gaussian(link = "identity"))),
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

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_wolf_lag2, HA_titer_tree_lag2, ha_lbi_shannon_lag1, H3.R0.mean) %>%
  pivot_longer(cols = c(HA_wolf_lag2, HA_titer_tree_lag2, ha_lbi_shannon_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.R0.mean, value)) %>%
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
        glm(H3.R0.mean ~ value, ., family = gaussian(link = "identity")) %>%
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

ha_lbi_R0 <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.R0.mean,
      ymin = gaussian(link = "identity")$linkinv(l),
      ymax = gaussian(link = "identity")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "ha_lbi_shannon_lag1") %>% rename(H3.R0.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = ha_lbi_shannon_lag1, y = H3.R0.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "identity")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = ha_lbi_shannon_lag1,
      ymin = H3.R0.lowCI, ymax = H3.R0.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = ha_lbi_shannon_lag1, y = H3.R0.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab(expression("H3 LBI Diversity (" ~ italic(t) ~ "-1)")) +
  ylab("A/H3 Effective Rt") +
  theme(legend.position = "bottom") +
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
    y = 2, x = -2.7, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ha_lbi_R0

all_epi_ha_lbi <- plot_grid(
  ha_lbi_R0 + theme(legend.position = "none"),
  ha_lbi_shannon + theme(legend.position = "none"),
  rel_widths = c(1, 1),
  nrow = 1
)
all_epi_ha_lbi

####################################################
## combine HA LBI and NA LBI results
####################################################
all_measures <- plot_grid(all_epi_ha_lbi + theme(legend.position = "none"), all_epi_na_lbi + theme(legend.position = "none"), nrow = 2, labels = "AUTO")
all_measures

epi_leg <- get_legend(ha_lbi_R0 + theme(
  legend.direction = "horizontal",
  legend.justification = "center",
  legend.box.just = "bottom",
  legend.position = "bottom",
  legend.title.align = 0.5
))

all_epi_leg_all_measures <- plot_grid(all_measures, epi_leg, nrow = 2, rel_heights = c(3, 0.3))
all_epi_leg_all_measures
## Figure S8
save_plot(all_epi_leg_all_measures, filename = "figures/LBI_vs_H3_epi_metrics_north_amer_build.png", base_width = 10, base_height = 10)
