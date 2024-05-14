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
load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")

epi_red2 <- epi_red %>%
  dplyr::select(
    season, region, dom_type,
    contains(c("duration", "onset", "peak", "mean","lbi","bayes"))
  ) %>%
  dplyr::select(-contains(c("usa")))
names(epi_red2)

epi_long <- epi_red2 %>%
  pivot_longer(
    cols = contains(c("lbi")),
    names_to = "evol_metrics",
    values_to = "value"
  )
head(epi_long)
names(epi_long)
unique(epi_long$season)

## onset days from Oct 1: N/A
## peak days from Oct 1: H3 LBI diversity (t-1)
## peak timing s.d. = N/A
## onset timing s.d. = H3 mean lbi (t-1), N2 mean lbi (t-1), H3 LBI diversity (t-1), N2 LBI diversity (t-1)
## seasonal duration = H3 mean lbi, N2 mean lbi, N2 LBI diversity, H3 LBI diversity, N2 mean lbi (t-1)
## days from onset t0 peak = N/A

h3_onset_predictors <- epi_long %>%
  # filter(!(season %in% c("2009-2010","2000-2001","2002-2003"))) %>% ## no H3 circulation
  split(.$evol_metrics) %>% #
  map(~ lm(as.numeric(onset_days_from_Oct1) ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_onset_predictors, decreasing = T)[1:4]
# ha_lbi_shannon na_lbi_shannon     HA_std_lbi     NA_std_lbi 
# 0.027161158    0.010110784    0.004508245   -0.002604102 

h3_peak_predictors <- epi_long %>%
  # filter(!(season %in% c("2009-2010","2000-2001","2002-2003"))) %>% ## no H3 circulation
  split(.$evol_metrics) %>% #
  map(~ lm(as.numeric(peak_days_from_Oct1) ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_peak_predictors, decreasing = T)[1:4]
# NA_std_lbi na_lbi_shannon ha_lbi_shannon     HA_std_lbi 
# 0.011074819    0.001636952   -0.003875405   -0.004300542 

h3_peak_predictors <- epi_long %>%
  distinct(season,dom_type,peak_timing_sd,evol_metrics,value)%>%
  split(.$evol_metrics) %>% #
  map(~ lm(as.numeric(peak_timing_sd) ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_peak_predictors, decreasing = T)

h3_onset_predictors <- epi_long %>%
  filter(!(season %in% c("2009-2010","2000-2001","2002-2003"))) %>% ## no H3 circulation
  distinct(season,dom_type,onset_timing_sd,evol_metrics,value)%>%
  split(.$evol_metrics) %>% #
  map(~ lm(as.numeric(onset_timing_sd) ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_onset_predictors, decreasing = T)

h3_duration_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ lm(H3_season_duration ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_duration_predictors, decreasing = T)
# na_lbi_shannon     NA_std_lbi ha_lbi_shannon     HA_std_lbi 
# 0.3738007      0.3521435      0.2783960      0.2073845 

h3_peak_diff_predictors <- epi_long %>%
  # filter(!(season %in% c("2009-2010","2000-2001","2002-2003"))) %>% ## no H3 circulation
  split(.$evol_metrics) %>% #
  map(~ lm(as.numeric(peak_diff) ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_peak_diff_predictors, decreasing = T)
# NA_std_lbi na_lbi_shannon     HA_std_lbi ha_lbi_shannon 
# 0.0330539503   0.0301474783   0.0034908543  -0.0006267528

###############################################################################################
## H3 LBI diversity vs seasonal duration
###############################################################################################
sort(names(epi_red2))
sum_df <- epi_red2 %>%
  filter(!(season %in% c("2000-2001","2002-2003")))%>%
  dplyr::select(
    H3_season_duration,
    season, region,dom_type,contains("lbi")
  ) %>%
  distinct()
sum_df$H3_season_duration <- as.numeric(sum_df$H3_season_duration)

sum_df %>%
  group_by(season) %>%
  tally()

names(sum_df)
sum_df <-
  sum_df %>%
  group_by(season, dom_type, 
           HA_std_lbi, NA_std_lbi, 
           ha_lbi_shannon, na_lbi_shannon) %>%
  summarise(
    H3.duration.mean = ci(as.numeric(H3_season_duration), na.rm = T)[1],
    H3.duration.lowCI = ci(as.numeric(H3_season_duration), na.rm = T)[2],
    H3.duration.hiCI = ci(as.numeric(H3_season_duration), na.rm = T)[3],
    H3.duration.sd = ci(as.numeric(H3_season_duration), na.rm = T)[4]
  ) %>%
  ungroup()

sum_df %>% dplyr::select(season,H3.duration.mean) %>% arrange(H3.duration.mean)

ggplot(sum_df)+
  geom_point(aes(x=ha_lbi_shannon,y=H3.duration.mean))

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
names(sum_df)

sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_std_lbi:na_lbi_shannon), ~ scale_this(.x))

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
1 - m1$deviance / m1$null.deviance # 0.28

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
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H3.duration.mean, ha_lbi_shannon, na_lbi_shannon) %>%
  pivot_longer(cols = c(ha_lbi_shannon, na_lbi_shannon)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  filter(!is.na(value))%>%
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

cols <- c("#B24745FF", "#00A1D5FF", "#6A6599FF", "#DF8F44FF")


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
    linewidth = 1, linetype = 2, color = "black"
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
  # xlab(expression("N2 LBI Diversity ("~italic(t)~"-1)"))+
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
    y = 36, x = -1.25, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ha_lbi_duration



###############################################################################################
## N2 LBI vs seasonal duration
###############################################################################################

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
1 - m1$deviance / m1$null.deviance # 0.35

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
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )
labels


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
  # xlab(expression("N2 LBI Diversity ("~italic(t)~"-1)"))+
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
    y = 36, x = -1.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_lbi_duration

lbi_epi_fig <- plot_grid(
  ha_lbi_duration + theme(legend.position = "none"),
  na_lbi_duration + theme(legend.position = "none"),
  rel_widths = c(1, 1),
  nrow = 1, labels = NULL
)
lbi_epi_fig

epi_leg <- get_legend(ha_lbi_duration +
                        guides(color = "none") +
                        theme(
                          legend.direction = "horizontal",
                          legend.justification = "center",
                          legend.box.just = "bottom",
                          legend.text = element_text(size = 12),
                          legend.title = element_text(size = 14)
                        ))

###############################################################################################
## H3 std LBI vs seasonal duration
###############################################################################################
y <- as.numeric(sumdf2$H3.duration.mean)
x <- sumdf2$HA_std_lbi

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
  glm(as.numeric(H3.duration.mean) ~ HA_std_lbi, analysis(split), family = gaussian(link = "identity"))
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

m1 <- glm(as.numeric(H3.duration.mean) ~ HA_std_lbi, data = sumdf2, family = gaussian(link = "identity"))
fam <- family(m1)
ilink <- fam$linkinv
coefini <- coef(m1)

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.27

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(as.numeric(H3.duration.mean) ~ HA_std_lbi, start = coefini, data = .x, family = gaussian(link = "identity"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  dplyr::summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.1f", adj.r.squared),
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H3.duration.mean, HA_std_lbi, NA_std_lbi) %>%
  pivot_longer(cols = c(HA_std_lbi, NA_std_lbi)) -> tbl_mtcars_long

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
        glm(as.numeric(H3.duration.mean) ~ value, ., start = coefini, family = gaussian(link = "identity")) %>%
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

ha_std_lbi_duration <- ggplot() +
  geom_ribbon(aes(x = HA_std_lbi, ymin = gaussian(link = "identity")$linkinv(l), ymax = gaussian(link = "identity")$linkinv(u)),
              tbl_plot_data %>% filter(name == "HA_std_lbi") %>% rename(HA_std_lbi = value),
              alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_std_lbi, y = H3.duration.mean),
    method = "glm", formula = y ~ x,
    method.args = list(start = coefini, family = gaussian(link = "identity")),
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2 %>% mutate(H3.duration.lowCI = ifelse(H3.duration.lowCI < 0, 0, H3.duration.lowCI)),
    aes(
      x = HA_std_lbi,
      ymin = H3.duration.lowCI, ymax = H3.duration.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_std_lbi, y = as.numeric(H3.duration.mean), fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("H3 s.d. LBI") +
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
    y = 34, x = -1.25, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ha_std_lbi_duration


###############################################################################################
## N2 mean LBI vs seasonal duration
###############################################################################################

y <- as.numeric(sumdf2$H3.duration.mean)
x <- sumdf2$NA_std_lbi

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
  glm(as.numeric(H3.duration.mean) ~ NA_std_lbi, analysis(split), family = gaussian(link = "identity"))
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

m1 <- glm(as.numeric(H3.duration.mean) ~ NA_std_lbi, data = sumdf2, family = gaussian(link = "identity"))
fam <- family(m1)
ilink <- fam$linkinv
coefini <- coef(m1)

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.32

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(as.numeric(H3.duration.mean) ~ NA_std_lbi, start = coefini, data = .x, family = gaussian(link = "identity"))),
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

na_std_lbi_duration <- ggplot() +
  geom_ribbon(aes(x = NA_std_lbi, ymin = gaussian(link = "identity")$linkinv(l), ymax = gaussian(link = "identity")$linkinv(u)),
              tbl_plot_data %>% filter(name == "NA_std_lbi") %>% rename(NA_std_lbi = value),
              alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_std_lbi, y = H3.duration.mean),
    method = "glm", formula = y ~ x,
    method.args = list(start = coefini, family = gaussian(link = "identity")),
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2 %>% mutate(H3.duration.lowCI = ifelse(H3.duration.lowCI < 0, 0, H3.duration.lowCI)),
    aes(
      x = NA_std_lbi,
      ymin = H3.duration.lowCI, ymax = H3.duration.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_std_lbi, y = as.numeric(H3.duration.mean), fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("N2 s.d. LBI") +
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
    y = 34, x = -1.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_std_lbi_duration




lbi_epi_fig2 <- plot_grid(
  ha_std_lbi_duration + theme(legend.position = "none"),
  # ha_lbi_duration_lag1 + theme(legend.position = "none"),
  na_std_lbi_duration + theme(legend.position = "none"),
  # na_lbi_duration_lag1 + theme(legend.position = "none"),
  rel_widths = c(1, 1),
  nrow = 1, labels = NULL
)
lbi_epi_fig2

lbi_epi_all = plot_grid(lbi_epi_fig,lbi_epi_fig2,nrow=2,labels="AUTO")
lbi_epi_all


lbi_leg_all <- plot_grid(lbi_epi_all , epi_leg, nrow = 2, rel_heights = c(2, 0.1))
lbi_leg_all

# save_plot(lbi_leg_all, filename = "figures/Fig5_ha_vs_na_lbi_season_duration.png", base_width = 12, base_height = 10)
save_plot(lbi_leg_all, filename = "figures/Fig5_ha_vs_na_lbi_season_duration.pdf", dpi=300, base_width = 12, base_height = 11)

