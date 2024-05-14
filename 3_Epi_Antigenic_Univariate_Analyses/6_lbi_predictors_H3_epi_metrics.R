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
load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")
sort(names(epi_red))

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

names(epi_red)

epi_red <- left_join(epi_red %>% dplyr::select(-contains(c("vs", "total"))),
                     subtype_dist %>% dplyr::select(season, region, h3_dom, h3_vs_h1),
                     by = c("season", "region")) %>%
  distinct() %>%
  filter(!(season %in% c("2009-2010", "2019-2020", "2020-2021", "2021-2022")))

unique(epi_red$season)
sort(names(epi_red))

epi_red = epi_red %>%
  dplyr::select(season, region,
                H3_shannon_entropy_res,H3_max_intensity,
                H3_cum_intensity,H3_max_Rt,
                dom_type,contains(c("lbi")))%>%
  dplyr::select(-contains("lag2"))
unique(epi_red$season)

epi_red %>%
  filter(!is.na(H3_max_Rt))%>%
  group_by(season)%>%
  tally()%>%
  arrange(n)

epi_red %>%
  filter(!is.na(H3_shannon_entropy_res))%>%
  group_by(season)%>%
  tally()%>%
  arrange(n)

names(epi_red)

epi_long <- epi_red %>%
  pivot_longer(
    cols = c(HA_std_lbi:na_lbi_shannon),
    names_to = "evol_metrics",
    values_to = "value"
  ) %>%
  distinct() %>%
  ungroup()

sort(unique(epi_long$evol_metrics))

h3_Rt_predictors <- epi_long %>%
  filter(!(season %in% c("2000-2001","2002-2003")))%>%
  # filter(season!="2000-2001")%>%
  split(.$evol_metrics) %>% #
  map(~ lm(H3_max_Rt ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_Rt_predictors, decreasing = T)[1:4]
# NA_std_lbi na_lbi_shannon     HA_std_lbi ha_lbi_shannon 
# 0.19413410     0.17520536     0.07546195     0.05975110 

h3_sh_predictors <- epi_long %>%
  # filter(!(season %in% c("2000-2001","2002-2003")))%>%
  # filter(season!="2000-2001")%>%
  split(.$evol_metrics) %>% #
  map(~ lm(H3_shannon_entropy_res ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")
sort(h3_sh_predictors, decreasing = T)[1:4]
# NA_std_lbi na_lbi_shannon     HA_std_lbi ha_lbi_shannon 
# 0.2073965      0.1873718      0.1139488      0.1094164

########################################################################################################
## Evolutionary indicators vs Epi Metrics
########################################################################################################

sum_df <-
  epi_red %>%
  dplyr::select(season, region,
                H3_shannon_entropy_res,H3_max_intensity,
                H3_cum_intensity,H3_max_Rt,
                dom_type,contains(c("lbi")))%>%
  dplyr::select(-contains("lag2"))%>%
  group_by(
    season, dom_type,
    # HA_mean_lbi_lag1,NA_mean_lbi_lag1,
    # HA_mean_lbi,NA_mean_lbi,
    HA_std_lbi, NA_std_lbi,
    ha_lbi_shannon, na_lbi_shannon,
    # ha_lbi_shannon_lag1, na_lbi_shannon_lag1
  ) %>%
  dplyr::summarise(
    H3.shannon.mean = ci(H3_shannon_entropy_res, na.rm = T)[1],
    H3.shannon.lowCI = ci(H3_shannon_entropy_res, na.rm = T)[2],
    H3.shannon.hiCI = ci(H3_shannon_entropy_res, na.rm = T)[3],
    H3.shannon.sd = ci(H3_shannon_entropy_res, na.rm = T)[4],
    # H3.peak.mean = ci(H3_max_intensity, na.rm = T)[1],
    # H3.peak.lowCI = ci(H3_max_intensity, na.rm = T)[2],
    # H3.peak.hiCI = ci(H3_max_intensity, na.rm = T)[3],
    # H3.peak.sd = ci(H3_max_intensity, na.rm = T)[4],
    # H3.epi.size.mean = ci(H3_cum_intensity, na.rm = T)[1],
    # H3.epi.size.lowCI = ci(H3_cum_intensity, na.rm = T)[2],
    # H3.epi.size.hiCI = ci(H3_cum_intensity, na.rm = T)[3],
    # H3.epi.size.sd = ci(H3_cum_intensity, na.rm = T)[4],
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
## NA LBI Diversity
####################################################

####################################################
## N2 LBI Diversity (current season) vs A/H3 shannon entropy
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010")))
  # filter(!(season %in% c("2009-2010", "2000-2001","2002-2003")))
scale_this <- function(x) as.vector(scale(x, center = T))
names(sumdf2)
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_std_lbi:na_lbi_shannon), ~ scale_this(.x))

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
    pvalue = sprintf("italic(P) == %.2f", pvalue)
  )
labels

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(ha_lbi_shannon, na_lbi_shannon, H3.shannon.mean) %>%
  pivot_longer(cols = c(ha_lbi_shannon, na_lbi_shannon)) -> tbl_mtcars_long

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

cols <- c("#B24745FF", "#00A1D5FF", "#6A6599FF", "#DF8F44FF")

NA_lbi_shannon <- ggplot() +
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
  # xlab("N2 Mean log(LBI)") +
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
    y = 0.78, x = 0.2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
NA_lbi_shannon

####################################################
## N2 LBI Diversity (current season) vs A/H3 Rt
####################################################
# sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010")))
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010", "2000-2001","2002-2003")))
scale_this <- function(x) as.vector(scale(x, center = T))
names(sumdf2)
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_std_lbi:na_lbi_shannon), ~ scale_this(.x))

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
  glm(H3.R0.mean ~ na_lbi_shannon, analysis(split), family = gaussian(link = "inverse"))
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

m1 <- glm(H3.R0.mean ~ na_lbi_shannon, data = sumdf2, family = gaussian(link = "inverse"))
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.32

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.R0.mean ~ na_lbi_shannon, data = .x, family = gaussian(link = "inverse"))),
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
  dplyr::select(ha_lbi_shannon,na_lbi_shannon, H3.R0.mean) %>%
  pivot_longer(cols = c(ha_lbi_shannon, na_lbi_shannon)) -> tbl_mtcars_long

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
        glm(H3.R0.mean ~ value, ., family = gaussian(link = "inverse")) %>%
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
      ymin = gaussian(link = "inverse")$linkinv(l),
      ymax = gaussian(link = "inverse")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "na_lbi_shannon") %>% rename(H3.R0.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = na_lbi_shannon, y = H3.R0.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "inverse")),
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
  ylab("A(H3N2) Effective Rt") +
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
    y = 1.8, x = 0.2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_lbi_R0

all_epi_na_lbi <- plot_grid(
  na_lbi_R0 + theme(legend.position = "none"),
  NA_lbi_shannon + theme(legend.position = "none"),
  rel_widths = c(1, 1),
  nrow = 1
)
all_epi_na_lbi

####################################################
## HA LBI
####################################################
####################################################
## H3 LBI Diversity (t-1) vs A/H3 shannon entropy
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010")))
# sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010", "2000-2001","2002-2003")))
scale_this <- function(x) as.vector(scale(x, center = T))
names(sumdf2)
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_std_lbi:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.shannon.mean
x <- sumdf2$ha_lbi_shannon

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
  betareg(H3.shannon.mean ~ ha_lbi_shannon, analysis(split))
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

m1 <- betareg(H3.shannon.mean ~ ha_lbi_shannon, data = sumdf2)
fam <- family(m1)
ilink <- fam$linkinv

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(H3.shannon.mean ~ ha_lbi_shannon, data = .x, link = "logit")),
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
labels

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(na_lbi_shannon,ha_lbi_shannon, H3.shannon.mean) %>%
  pivot_longer(cols = c(na_lbi_shannon, ha_lbi_shannon)) -> tbl_mtcars_long


tbl_mtcars_long %>%
  filter(!is.na(value))%>%
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

data.table::data.table(sumdf2[!complete.cases(sumdf2),])

ha_lbi_shannon <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.shannon.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "ha_lbi_shannon") %>% rename(H3.shannon.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = ha_lbi_shannon, y = H3.shannon.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = ha_lbi_shannon,
      ymin = H3.shannon.lowCI, ymax = H3.shannon.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = ha_lbi_shannon, y = H3.shannon.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("H3 LBI Diversity")+
  # xlab(expression("H3 LBI Diversity (" ~ italic(t) ~ "-1)")) +
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
    y = 0.78, x = 0.2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ha_lbi_shannon

####################################################
## H3 LBI Div vs A/H3 effective Rt
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010", "2000-2001","2002-2003")))
scale_this <- function(x) as.vector(scale(x, center = T))
names(sumdf2)
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_std_lbi:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.R0.mean
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
  glm(H3.R0.mean ~ ha_lbi_shannon, analysis(split), family = gaussian(link = "inverse"))
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

m1 <- glm(H3.R0.mean ~ ha_lbi_shannon, data = sumdf2, family = gaussian(link = "inverse"))
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.13

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.R0.mean ~ ha_lbi_shannon, data = .x, family = gaussian(link = "inverse"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
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
  dplyr::select(na_lbi_shannon, ha_lbi_shannon, H3.R0.mean) %>%
  pivot_longer(cols = c(na_lbi_shannon, ha_lbi_shannon)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  filter(!is.na(value))%>%
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
        glm(H3.R0.mean ~ value, ., family = gaussian(link = "inverse")) %>%
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
      ymin = gaussian(link = "inverse")$linkinv(l),
      ymax = gaussian(link = "inverse")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "ha_lbi_shannon") %>% rename(H3.R0.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = ha_lbi_shannon, y = H3.R0.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "inverse")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = ha_lbi_shannon,
      ymin = H3.R0.lowCI, ymax = H3.R0.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = ha_lbi_shannon, y = H3.R0.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab("H3 LBI Diversity")+
  # xlab(expression("H3 LBI Diversity (" ~ italic(t) ~ "-1)")) +
  ylab("A(H3N2) Effective Rt") +
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
    y = 1.81, x = 0.2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ha_lbi_R0

#current season
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
all_measures <- plot_grid(all_epi_ha_lbi + theme(legend.position = "none"), 
                          all_epi_na_lbi + theme(legend.position = "none"), 
                          nrow = 2, labels = "AUTO")
all_measures

epi_leg <- cowplot::get_plot_component(ha_lbi_R0+
                                         guides(color = "none") +
                                         theme(
                                           legend.position = "bottom",
                                           legend.direction = "horizontal",
                                           legend.justification = "center",
                                           legend.box.just = "bottom",
                                           legend.text = element_text(size = 12),
                                           legend.title = element_text(size = 14,hjust = 0.5)
                                         ), 
                                       'guide-box-bottom', return_all = TRUE)
cowplot::ggdraw(epi_leg)

all_epi_leg_all_measures <- plot_grid(all_measures, epi_leg, nrow = 2, rel_heights = c(3, 0.3))
all_epi_leg_all_measures

# save_plot(all_epi_leg_all_measures, filename = "figures/Fig3_sup_fig4_LBI_vs_H3_epi_metrics_north_amer_build.png", base_width = 10, base_height = 10)
save_plot(all_epi_leg_all_measures, filename = "figures/Fig3_sup_fig4_LBI_vs_H3_epi_metrics_north_amer_build.pdf", dpi = 300, base_width = 10, base_height = 10)


####################################################
## NA std LBI
####################################################

####################################################
## N2 Mean LBI vs A/H3 shannon entropy
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010")))
# sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010", "2000-2001","2002-2003")))
scale_this <- function(x) as.vector(scale(x, center = T))
names(sumdf2)
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_std_lbi:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.shannon.mean
x <- sumdf2$NA_std_lbi

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
  betareg(H3.shannon.mean ~ NA_std_lbi, analysis(split))
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

m1 <- betareg(H3.shannon.mean ~ NA_std_lbi, data = sumdf2)
fam <- family(m1)
ilink <- fam$linkinv

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(H3.shannon.mean ~ NA_std_lbi, data = .x)),
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
labels

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(HA_std_lbi,NA_std_lbi, H3.shannon.mean) %>%
  pivot_longer(cols = c(HA_std_lbi,NA_std_lbi)) -> tbl_mtcars_long

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

NA_std_lbi_shannon <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.shannon.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "NA_std_lbi") %>% rename(H3.shannon.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_std_lbi, y = H3.shannon.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_std_lbi,
      ymin = H3.shannon.lowCI, ymax = H3.shannon.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_std_lbi, y = H3.shannon.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("N2 s.d. LBI") +
  # xlab("N2 LBI Diversity") +
  ylab("A(H3N2) Epidemic Intensity") +
  # theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
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
    y = 0.78, x = 0.2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
NA_std_lbi_shannon

####################################################
## N2 LBI vs A/H3 R0
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010", "2000-2001","2002-2003")))
scale_this <- function(x) as.vector(scale(x, center = T))
names(sumdf2)
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_std_lbi:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.R0.mean
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
  glm(H3.R0.mean ~ NA_std_lbi, analysis(split), family = gaussian(link = "inverse"))
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

m1 <- glm(H3.R0.mean ~ NA_std_lbi, data = sumdf2, family = gaussian(link = "inverse"))
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.34

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.R0.mean ~ NA_std_lbi, data = .x, family = gaussian(link = "inverse"))),
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
  dplyr::select(HA_std_lbi, NA_std_lbi, H3.R0.mean) %>%
  pivot_longer(cols = c(HA_std_lbi, NA_std_lbi)) -> tbl_mtcars_long

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
        glm(H3.R0.mean ~ value, ., family = gaussian(link = "inverse")) %>%
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

na_std_lbi_R0 <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.R0.mean,
      ymin = gaussian(link = "inverse")$linkinv(l),
      ymax = gaussian(link = "inverse")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "NA_std_lbi") %>% rename(H3.R0.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = NA_std_lbi, y = H3.R0.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "inverse")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = NA_std_lbi,
      ymin = H3.R0.lowCI, ymax = H3.R0.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = NA_std_lbi, y = H3.R0.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab("N2 s.d. LBI")+
  ylab("A(H3N2) Effective Rt") +
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
    y = 1.81, x = 0.2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
na_std_lbi_R0

all_epi_na_lbi <- plot_grid(
  na_std_lbi_R0 + theme(legend.position = "none"),
  NA_std_lbi_shannon + theme(legend.position = "none"),
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
unique(sum_df$season)
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010")))
scale_this <- function(x) as.vector(scale(x, center = T))
names(sumdf2)
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_std_lbi:na_lbi_shannon), ~ scale_this(.x))

y <- sumdf2$H3.shannon.mean
x <- sumdf2$HA_std_lbi

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
  betareg(H3.shannon.mean ~ HA_std_lbi, analysis(split))
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

m1 <- betareg(H3.shannon.mean ~ HA_std_lbi, data = sumdf2)
fam <- family(m1)
ilink <- fam$linkinv
ilink

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ betareg(H3.shannon.mean ~ HA_std_lbi, data = .x, link = "logit")),
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
labels

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(NA_std_lbi, HA_std_lbi, H3.shannon.mean) %>%
  pivot_longer(cols = c(NA_std_lbi,  HA_std_lbi)) -> tbl_mtcars_long

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

HA_std_lbi_shannon <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.shannon.mean,
      ymin = l,
      ymax = u
    ),
    tbl_plot_data %>% filter(name == "HA_std_lbi") %>% rename(H3.shannon.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_std_lbi, y = H3.shannon.mean),
    method = "betareg", formula = y ~ x,
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_std_lbi,
      ymin = H3.shannon.lowCI, ymax = H3.shannon.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_std_lbi, y = H3.shannon.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("H3 s.d. LBI")+
  # xlab(expression("H3 mean LBI (" ~ italic(t) ~ "-1)")) +
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
    y = 0.78, x = 0.2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
HA_std_lbi_shannon

####################################################
## H3 LBI vs A/H3 effective Rt
####################################################
sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010", "2000-2001","2002-2003")))
scale_this <- function(x) as.vector(scale(x, center = T))
sumdf2 <- sumdf2 %>%
  mutate_at(vars(HA_std_lbi:na_lbi_shannon), ~ scale_this(.x))


y <- sumdf2$H3.R0.mean
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
  glm(H3.R0.mean ~ HA_std_lbi, analysis(split), family = gaussian(link = "inverse"))
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

m1 <- glm(H3.R0.mean ~ HA_std_lbi, data = sumdf2, family = gaussian(link = "inverse"))
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.41

labels <- boots %>%
  mutate(
    model = purrr::map(splits, ~ glm(H3.R0.mean ~ HA_std_lbi, data = .x, family = gaussian(link = "inverse"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
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
  dplyr::select(NA_std_lbi, HA_std_lbi, H3.R0.mean) %>%
  pivot_longer(cols = c(NA_std_lbi, HA_std_lbi)) -> tbl_mtcars_long

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
        glm(H3.R0.mean ~ value, ., family = gaussian(link = "inverse")) %>%
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
      ymin = gaussian(link = "inverse")$linkinv(l),
      ymax = gaussian(link = "inverse")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "HA_std_lbi") %>% rename(H3.R0.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = HA_std_lbi, y = H3.R0.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "inverse")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = HA_std_lbi,
      ymin = H3.R0.lowCI, ymax = H3.R0.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = HA_std_lbi, y = H3.R0.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab("H3 s.d. LBI")+
  # xlab(expression("H3 mean LBI (" ~ italic(t) ~ "-1)")) +
  ylab("A(H3N2) Effective Rt") +
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
    y = 1.81, x = 0.2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
ha_lbi_R0

all_epi_ha_lbi <- plot_grid(
  ha_lbi_R0 + theme(legend.position = "none"),
  HA_std_lbi_shannon + theme(legend.position = "none"),
  rel_widths = c(1, 1),
  nrow = 1
)
all_epi_ha_lbi

####################################################
## combine HA LBI and NA LBI results
####################################################
all_measures <- plot_grid(all_epi_ha_lbi + theme(legend.position = "none"), 
                          all_epi_na_lbi + theme(legend.position = "none"), 
                          nrow = 2, labels = "AUTO")
all_measures

epi_leg <- cowplot::get_plot_component(ha_lbi_R0+
                                         guides(color = "none") +
                                         theme(
                                           legend.position = "bottom",
                                           legend.direction = "horizontal",
                                           legend.justification = "center",
                                           legend.box.just = "bottom",
                                           legend.text = element_text(size = 12),
                                           legend.title = element_text(size = 14,hjust = 0.5)
                                         ), 
                                       'guide-box-bottom', return_all = TRUE)
cowplot::ggdraw(epi_leg)

all_epi_leg_all_measures <- plot_grid(all_measures, epi_leg, nrow = 2, rel_heights = c(3, 0.3))
all_epi_leg_all_measures

# save_plot(all_epi_leg_all_measures, filename = "figures/Fig3_sup_fig3_std_lbi_vs_H3_epi_metrics_north_amer_build.png", base_width = 10, base_height = 10)
save_plot(all_epi_leg_all_measures, filename = "figures/Fig3_sup_fig3_std_lbi_vs_H3_epi_metrics_north_amer_build.pdf", dpi=300,base_width = 10, base_height = 10)
