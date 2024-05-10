# # ## load packages
# # # you will need your own API key to use ggmap
list.of.packages <- c(
  "dplyr", "ggplot2", "cowplot", "tidyr", "vegan", "broom", "gmodels", "gghighlight",
  "betareg", "purrr", "tidymodels", "MuMIn", "ggmap", "scatterpie","ggrepel"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

theme_set(theme_cowplot(font_size = 16))
##################################################################################
## Associations between A/H3N2 subtype dominance and evolutionary indicators
##################################################################################

## load data
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

load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")
head(epi_red)
names(epi_red)

combined <- left_join(subtype_dist,
  epi_red %>% dplyr::select(
    -h3_vs_h1, -iva_vs_ivb, -h3_vs_total_flu,
    -h1_vs_total_flu, -ivb_vs_iva, -total_a
  ),
  by = c("season", "region")
) %>%
  distinct() %>%
  filter(!(season %in% c("2009-2010", "2019-2020", "2020-2021", "2021-2022"))) %>%
  tidyr::separate(season, into = c("year1", "year2"), sep = "-", remove = F)
names(combined)

combined %>%
  mutate(pre_pandemic = if_else(year1 < 2010, "pre", "post")) %>%
  group_by(dom_type, pre_pandemic) %>%
  summarize(
    mean_dom = mean(h3_vs_flu_samples),
    mean_iva_dom = mean(h3_vs_h1)
  )

names(combined)
combined %>%
  # dplyr::select(season, region, contains(c("total", "vs", "samples", "dom", "lbi", "mean", "lag1", "lag2"))) %>%
  # dplyr::select(-contains(c("lag0", "sd", "usa", "inv", "sub", "std", "HA_stem_lag"))) %>%
  dplyr::select(season, region, contains(c("vs", "dom", "lbi", "mean", "lag1", "lag2"))) %>%
  names()
names(combined)

## which evolutionary indicators have strongest correlations with subtype dominance?
epi_long <- combined %>%
  dplyr::select(season, region, contains(c("vs", "dom", "lbi", "mean", "lag1", "lag2"))) %>%
  dplyr::select(-contains(c("lag0", "usa", "inv", "sub", "HA_stem_lag"))) %>%
  group_by(season, dom_type) %>%
  summarize_at(vars(h3_vs_h1, h3_dom, HA_std_lbi:NA_krammer_ep_lag2), mean) %>%
  ungroup() %>%
  pivot_longer(
    cols = contains(c("lag", "lbi")),
    names_to = "evol_metrics",
    values_to = "value"
  )

unique(epi_long$evol_metrics)

h3_dom_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ betareg(h3_vs_h1 ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("pseudo.r.squared")
sort(h3_dom_predictors, decreasing = T)[1:5]

h3_dom_predictors <- epi_long %>%
  split(.$evol_metrics) %>% #
  map(~ betareg(h3_dom ~ scale(value), data = .)) %>%
  map(summary) %>%
  map_dbl("pseudo.r.squared")
sort(h3_dom_predictors, decreasing = T)[1:5]

## dataset for regression models
epi_red2 <- combined %>%
  group_by(season, year1, year2, dom_type, HA_wolf_lag2, NA_bhatt_ep_lag1, HA_koel_lag2) %>%
  dplyr::summarize(
    h3_dom_mean = gmodels::ci(h3_dom, na.rm = T)[1],
    h3_dom_lwr = gmodels::ci(h3_dom, na.rm = T)[2],
    h3_dom_upper = gmodels::ci(h3_dom, na.rm = T)[3],
    h3_dom_se = gmodels::ci(h3_dom, na.rm = T)[4]
  ) %>%
  ungroup() %>%
  filter(!(season %in% c("1995-1996", "1996-1997", "2009-2010")))

### scale predictor variables
scale_this <- function(x) as.vector(scale(x, center = T))
epi_red2 <- epi_red2 %>%
  mutate_at(vars(HA_wolf_lag2:HA_koel_lag2), ~ scale_this(.x))

ggplot(epi_red2) +
  geom_point(aes(x = as.numeric(year2), y = h3_dom_mean, color = dom_type))

epi_red2 %>%
  mutate(pre_pandemic = if_else(as.numeric(year1) < 2010, "pre", "post")) %>%
  group_by(dom_type, pre_pandemic) %>%
  summarize(mean_dom = mean(h3_dom_mean))

##################################################################################
## NA epitope distance (t-1) vs H3 dominance
##################################################################################

## compare fit of different models
y <- epi_red2$h3_dom_mean
x <- epi_red2$NA_bhatt_ep_lag1

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
quasibinomial.model <- glm(y ~ x, family = quasibinomial(link = "logit"))
beta.model <- betareg(y ~ x, link = c("logit"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, quasibinomial.model, beta.model, rank = "BIC")
# beta model fits best

#### create bootstrap samples of dataset
set.seed(27)
boots <- bootstraps(epi_red2, times = 1000, apparent = TRUE)

fit_lm_on_bootstrap <- function(split) {
  betareg(h3_dom_mean ~ NA_bhatt_ep_lag1, analysis(split))
}

boot_models <- boots %>%
  mutate(
    model = purrr::map(splits, fit_lm_on_bootstrap),
    coef_inf = map(model, tidy)
  )

# boot_aug <-
#   boot_models %>%
#   sample_n(200) %>%
#   mutate(augmented = map(model, augment)) %>%
#   unnest(augmented)

# boot_models <- boots %>%
#   mutate(
#     model = purrr::map(splits, ~ betareg(h3_dom_mean ~ NA_bhatt_ep_lag1, data = .)),
#     coef_inf = map(model, tidy)
#   )

boot_coefs <- boot_models %>%
  unnest(coef_inf)

percentile_intervals <- int_pctl(boot_models, coef_inf)
percentile_intervals

## R^2 and P-values for scatterplots
labels <- boots %>%
  mutate(
    model = map(splits, fit_lm_on_bootstrap),
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

### fit model to each bootstrap sample
set.seed(27)
n_boot <- 1000

# epi_red2 %>%
#   dplyr::select(-c(season, dom_type, h3_dom_lwr:h3_dom_se)) %>%
#   pivot_longer(cols = c(HA_wolf_lag2, NA_bhatt_ep_lag1)) -> tbl_mtcars_long

epi_red2 %>%
  dplyr::select(-c(season, dom_type, h3_dom_lwr:h3_dom_se, year1, year2)) %>%
  pivot_longer(cols = c(NA_bhatt_ep_lag1, HA_wolf_lag2, HA_koel_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  dplyr::filter(!is.na(value)) %>%
  nest(model_data = c(h3_dom_mean, value)) %>%
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
        sample_n(n, replace = TRUE) %>%
        betareg(h3_dom_mean ~ value, .) %>%
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
        median = quantile(.fitted, 0.5),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data


epi_red2 <- epi_red2 %>% tidyr::separate(col = "season", sep = "-", remove = F, into = c("year1", "year2"))

## differentiate pre and post-2009 H1
epi_red2 <- epi_red2 %>%
  mutate(dom_type2 = case_when(
    year1 < 2009 & dom_type == "H1" ~ "H1",
    year1 > 2009 & dom_type == "H1" ~ "H1pdm",
    dom_type == "H3" ~ "H3",
    dom_type == "co-circ" ~ "H3/H1pdm"
  ))


## plot
cols <- c("#B24745FF", "#00A1D5FF", "#6A6599FF", "#DF8F44FF")

NA_ep_dom <- ggplot(
  epi_red2,
  aes(x = NA_bhatt_ep_lag1, y = h3_dom_mean)
) +
  geom_ribbon(aes(x = h3_dom_mean, ymin = l, ymax = u),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(h3_dom_mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_line(aes(x = h3_dom_mean, y = median),
    tbl_plot_data %>% filter(name == "NA_bhatt_ep_lag1") %>% rename(h3_dom_mean = value),
    linewidth = 1, linetype = 2, color = "black"
  ) +
  # geom_smooth(
  #   data = epi_red2,
  #   aes(y = h3_dom_mean, x = NA_bhatt_ep_lag1),
  #   method = "betareg", formula = y ~ x,
  #   se = F,
  #   linewidth = 1, linetype = 2, color = "red"
  # ) +
  geom_point(
    data = subset(epi_red2, season %in% c("2003-2004", "2007-2008")),
    aes(x = NA_bhatt_ep_lag1, y = h3_dom_mean), fill = "green", pch = 21, size = 8
  ) +
  geom_errorbar(
    data = epi_red2,
    aes(
      x = NA_bhatt_ep_lag1,
      ymin = h3_dom_lwr, ymax = h3_dom_upper,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = epi_red2,
    aes(x = NA_bhatt_ep_lag1, y = h3_dom_mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab(expression("N2 epitope distance (" ~ italic(t) ~ "-1)")) +
  ylab("A(H3N2) Dominance") +
  # theme(legend.title.align = 0.5) +
  theme(legend.title = element_text(0.5)) +
  theme(legend.position = "right") +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant\nIAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant\nIAV"
  ) +
  geom_text(
    y = 0, x = 0.2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  ) +
  geom_text(
    data = subset(epi_red2, season %in% c("2003-2004", "2007-2008")), nudge_y = -0.08, nudge_x = 0.2,
    aes(x = NA_bhatt_ep_lag1, y = h3_dom_mean, label = season), cex = 3.5
  )

NA_ep_dom

##################################################################################
## HA epitope distance (t-2) vs H3 dominance
##################################################################################

y <- epi_red2$h3_dom_mean
x <- epi_red2$HA_wolf_lag2

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
quasibinomial.model <- glm(y ~ x, family = quasibinomial(link = "logit"))
beta.model <- betareg(y ~ x, link = c("logit"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, quasibinomial.model, beta.model, rank = "BIC")
## beta model fits best

## bootstrap samples of the dataset
# set.seed(27)
# boots <- bootstraps(epi_red2, times = 1000, apparent = TRUE)
#
# boot_models <- boots %>%
#   mutate(
#     model = map(splits, ~ betareg(h3_dom_mean ~ HA_wolf_lag2, data = .)),
#     coef_inf = map(model, tidy)
#   )

fit_lm_on_bootstrap <- function(split) {
  betareg(h3_dom_mean ~ HA_wolf_lag2, analysis(split))
}

boot_models <- boots %>%
  mutate(
    model = purrr::map(splits, fit_lm_on_bootstrap),
    coef_inf = map(model, tidy)
  )

boot_coefs <- boot_models %>% unnest(coef_inf)

percentile_intervals <- int_pctl(boot_models, coef_inf)
percentile_intervals

## R^2 and P-values for scatterplot
labels <- boots %>%
  mutate(
    model = map(splits, fit_lm_on_bootstrap),
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

## plot
HA_ep_dom <- ggplot(
  epi_red2,
  aes(x = HA_wolf_lag2, y = h3_dom_mean)
) +
  geom_ribbon(aes(x = h3_dom_mean, ymin = l, ymax = u),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>%
      rename(h3_dom_mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_line(aes(x = h3_dom_mean, y = median),
    tbl_plot_data %>% filter(name == "HA_wolf_lag2") %>% rename(h3_dom_mean = value),
    linewidth = 1, linetype = 2, color = "black"
  ) +
  # geom_smooth(
  #   data = epi_red2,
  #   aes(y = h3_dom_mean, x = HA_wolf_lag2),
  #   method = "betareg", formula = y ~ x,
  #   se = F,
  #   linewidth = 1, linetype = 2, color = "black"
  # ) +
  geom_point(
    data = subset(epi_red2, season %in% c("2003-2004", "2007-2008")),
    aes(x = HA_wolf_lag2, y = h3_dom_mean), fill = "green", pch = 21, size = 8
  ) +
  geom_errorbar(
    data = epi_red2,
    aes(
      x = HA_wolf_lag2,
      ymin = h3_dom_lwr, ymax = h3_dom_upper,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = epi_red2,
    aes(x = HA_wolf_lag2, y = h3_dom_mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab(expression("H3 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab("A(H3N2) Dominance") +
  theme(legend.title.align = 0.5) +
  theme(legend.position = "right") +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant\nIAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant\nIAV"
  ) +
  geom_text(
    y = 0, x = 0.2, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  ) +
  geom_text(
    data = subset(epi_red2, season %in% c("2003-2004", "2007-2008")), nudge_y = -0.08, nudge_x = 0.2,
    aes(x = HA_wolf_lag2, y = h3_dom_mean, label = season), cex = 3.5
  )
HA_ep_dom

both <- plot_grid(HA_ep_dom + theme(legend.position = "none"),
  NA_ep_dom + theme(legend.position = "none"),
  labels = "AUTO", nrow = 2
)

epi_leg <- get_legend(HA_ep_dom +
  guides(color = "none") +
  theme(
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box.just = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  ))
all_epi_leg <- plot_grid(both, epi_leg, nrow = 2, rel_heights = c(2, 0.1))
all_epi_leg

##################################################################################
## Pie charts of subtype composition for each season/region
##################################################################################
load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")

# weekly incidence data
load("data/hhs_division_level_ILI_and_virology_interp_smoothed_sampling_effort.RData") # regionflu_ili_vir_adj

# season-level H3 data
load("data/region_level_flu_metrics_H3.RData") # region_flu_metrics_H3

seasons <- unique(regionflu_ili_vir_adj$season_description)
seasons # check

states_to_keep <- regionflu_ili_vir_adj %>%
  mutate(prop_positive = percent_positive / 100) %>%
  dplyr::select(season_description, month.day, region, prop_h1, prop_h3, prop_b, prop_a, prop_positive, total_specimens) %>%
  group_by(season_description, region) %>%
  filter(month.day >= "10-01" | month.day <= "06-01") %>%
  summarise(
    non_na_count = sum(!is.na(prop_positive)),
    sample_count = sum(total_specimens)
  ) %>%
  filter(non_na_count > 21) %>%
  ungroup()

sort(names(epi_red))
subtype_dist <- left_join(states_to_keep,
  subtype_dist %>% dplyr::select(region, season, h3_total, h1_total, b_total, flu_samples, resp_samples),
  by = c("region", "season_description" = "season")
)
names(subtype_dist)
subtype_dist <- subtype_dist %>%
  mutate(pos_prop = flu_samples / resp_samples) %>%
  rename(season = season_description)

### you need your own key to access Google API and download maps
# ggmap::register_google(key = "yourkey",write=T)
cities <- c("Boston", "Syracuse", "Philadelphia", "Atlanta", "Chicago", "Dallas", "Kansas City", "Denver", "San Francisco", "Seattle")
cities_df <- as.data.frame(cities, place = cities)
cities_df$cities <- as.character(cities_df$cities)
locations_df <- mutate_geocode(cities_df, cities)
locations_df$region <- c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5", "Region 6", "Region 7", "Region 8", "Region 9", "Region 10")

pie_df <- left_join(subtype_dist, locations_df, by = "region") %>%
  mutate(untyped = flu_samples - h3_total - h1_total - b_total)

label_df = pie_df %>% filter(season == "2010-2011") %>%
  dplyr::select(lon, lat, region)%>%
  mutate(region2 = substr(region,start=8,stop=9))

region1 <- c("Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont")
region2 <- c("New Jersey", "New York", "Puerto Rico", "Virgin Islands")
region3 <- c("Delaware", "District of Columbia", "Maryland", "Pennsylvania", "Virginia", "West Virginia")
region4 <- c("Alabama", "Florida", "Georgia", "Kentucky", "Mississippi", "North Carolina", "South Carolina", "Tennessee")
region5 <- c("Illinois", "Indiana", "Michigan", "Minnesota", "Ohio", "Wisconsin")
region6 <- c("Arkansas", "Louisiana", "New Mexico", "Oklahoma", "Texas")
region7 <- c("Iowa", "Kansas", "Missouri", "Nebraska")
region8 <- c("Colorado", "Montana", "North Dakota", "South Dakota", "Utah", "Wyoming")
region9 <- c("Arizona", "California", "Hawaii", "Nevada")
region10 <- c("Alaska", "Idaho", "Oregon", "Washington")

hhs.list <- list(
  Region_1 = region1, Region_2 = region2, Region_3 = region3,
  Region_4 = region4, Region_5 = region5, Region_6 = region6,
  Region_7 = region7, Region_8 = region8, Region_9 = region9,
  Region_10 = region10
)

hhsdf <- plyr::ldply(hhs.list, data.frame)
colnames(hhsdf) <- c("hhs_region", "region")

hhsdf$region <- tolower(hhsdf$region)

us <- map_data("state")
us

us <- left_join(us, hhsdf, by = "region")
us

cols <- c(
  "#40abf0",
  "#dd3b3f",
  "#00a100"
)

fujian_pie <- ggplot(us, aes(long, lat)) +
  geom_map(map = us, aes(map_id = region, fill = hhs_region), color = "grey", alpha = 0.3) +
  geom_scatterpie(
    data = pie_df %>% filter(season == "2003-2004"),
    aes(lon, lat, r = sqrt(pos_prop * 100) / 2.5),
    cols = c("h3_total", "h1_total", "b_total")
  ) +
  coord_fixed() +
  scale_fill_manual(
    name = "Flu type",
    breaks = c("h3_total", "h1_total", "b_total"),
    labels = c("A(H3N2)", "A(H1N1)", "B"),
    values = c(
      "h3_total" = "#B24745FF",
      "h1_total" = "#00A1D5FF",
      "b_total" = "#00a100",
      # "untyped" = "white",
      "Region_1" = "#8dd3c7",
      "Region_2" = "#ffffb3",
      "Region_3" = "#bebada",
      "Region_4" = "#fb8072",
      "Region_5" = "#80b1d3",
      "Region_6" = "#fdb462",
      "Region_7" = "#b3de69",
      "Region_8" = "#fccde5",
      "Region_9" = "#d9d9d9",
      "Region_10" = "#bc80bd"
    )
  ) +
  geom_label_repel(data=label_df,
                   aes(x = lon, y = lat, label = region2),
                   fontface="bold",
                   size = 4,
                   nudge_x = 0.2,
                   label.r = 0.1,
                   # label.padding = 0.1,
                   point.padding = 3, # additional padding around each point
                   min.segment.length = Inf) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing = unit(c(0, 0, 0, 0), "cm"),
    plot.background = element_rect(fill = "transparent", colour = NA),
    plot.margin = unit(c(-1, 0, -1, 0), "cm"),
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_blank()
  ) +
  ggtitle("2003-2004")
fujian_pie

cols <- c("#DF8F44FF", "#00A1D5FF", "#B24745FF")
brisbane_pie <- ggplot(us, aes(long, lat)) +
  geom_map(map = us, aes(map_id = region, fill = hhs_region), color = "grey", alpha = 0.3) +
  geom_scatterpie(
    data = pie_df %>% filter(season == "2007-2008"),
    aes(lon, lat, r = sqrt(pos_prop * 100) / 2.5),
    cols = c("h3_total", "h1_total", "b_total")
  ) +
  coord_fixed() +
  scale_fill_manual(
    name = "Flu type",
    breaks = c("h3_total", "h1_total", "b_total"),
    labels = c("A(H3N2)", "A(H1N1)", "B"),
    values = c(
      "h3_total" = "#B24745FF",
      "h1_total" = "#00A1D5FF",
      "b_total" = "#00a100",
      # "h3_total" = "#80cdc1",
      # "h1_total" = "#dfc27d",
      # "b_total" ="#D67236",
      # "untyped" = "white",
      "Region_1" = "#8dd3c7",
      "Region_2" = "#ffffb3",
      "Region_3" = "#bebada",
      "Region_4" = "#fb8072",
      "Region_5" = "#80b1d3",
      "Region_6" = "#fdb462",
      "Region_7" = "#b3de69",
      "Region_8" = "#fccde5",
      "Region_9" = "#d9d9d9",
      "Region_10" = "#bc80bd"
    )
  ) +
  geom_label_repel(data=label_df,
                   aes(x = lon, y = lat, label = region2),
                   fontface="bold",
                   size = 4,
                   nudge_x = 0.2,
                   label.r = 0.1,
                   # label.padding = 0.1,
                   point.padding = 3, # additional padding around each point
                   min.segment.length = Inf) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing = unit(c(0, 0, 0, 0), "cm"),
    plot.background = element_rect(fill = "transparent", colour = NA),
    plot.margin = unit(c(-1, 0, -1, 0), "cm"),
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_blank()
  ) +
  ggtitle("2007-2008")
brisbane_pie
fujian_pie

legend_b <- get_legend(
  brisbane_pie +
    theme(
      legend.direction = "horizontal",
      legend.justification = "center", legend.box.just = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12)
    )
)

scatter_grid <- plot_grid(fujian_pie + theme(legend.position = "none"),
  brisbane_pie + theme(legend.position = "none"),
  align = "vh",
  labels = c("C", "D"),
  hjust = -1,
  nrow = 2
)
maps_and_leg <- plot_grid(scatter_grid, legend_b, nrow = 2, rel_heights = c(2, 0.2))

maps_and_leg

scatter_and_maps <- plot_grid(all_epi_leg, maps_and_leg, nrow = 1, rel_widths = c(1, 1))
scatter_and_maps
save_plot(scatter_and_maps, filename = "figures/fig_4_h3_dominance_scatter_and_pie_maps.png", base_width = 10, base_height = 9)
save_plot(scatter_and_maps, filename = "figures/fig_4_h3_dominance_scatter_and_pie_maps.pdf", dpi=300, base_width = 10, base_height = 9)

##################################################################################
## Pie charts for all seasons for supplement
##################################################################################
unique(pie_df$season)
pie_df <- pie_df %>% filter(!(season %in% c("2019-2020", "2020-2021")))

all_seasons <- ggplot(us, aes(long, lat)) +
  geom_map(map = us, aes(map_id = region, fill = hhs_region), color = "grey", alpha = 0.3) +
  geom_scatterpie(
    data = pie_df[complete.cases(pie_df), ],
    aes(lon, lat, r = 3),
    cols = c("h3_total", "h1_total", "b_total")
  ) +
  coord_fixed() +
  facet_wrap(~season) +
  scale_fill_manual(
    name = "Flu type",
    breaks = c("h3_total", "h1_total", "b_total"),
    labels = c("A(H3N2)", "A(H1N1) or A(H1N1)pdm09", "B"),
    values = c(
      "h3_total" = "#B24745FF",
      "h1_total" = "#00A1D5FF",
      "b_total" = "#00a100",
      # "h3_total" = "#80cdc1",
      # "h1_total" = "#dfc27d",
      # "b_total" ="#D67236",
      # "untyped" = "white",
      "Region_1" = "#8dd3c7",
      "Region_2" = "#ffffb3",
      "Region_3" = "#bebada",
      "Region_4" = "#fb8072",
      "Region_5" = "#80b1d3",
      "Region_6" = "#fdb462",
      "Region_7" = "#b3de69",
      "Region_8" = "#fccde5",
      "Region_9" = "#d9d9d9",
      "Region_10" = "#bc80bd"
    )
  ) +
  theme_bw(base_size = 16) +
  geom_label_repel(data=label_df,
    aes(x = lon, y = lat, label = region2),
    fontface="bold",
    size = 3,
    nudge_x = 0.1,
    label.r = 0.1,
    label.padding = 0.1,
    # point.padding = 3, # additional padding around each point
    min.segment.length = Inf) +
  theme(
    legend.position = "bottom",
    # legend.justification = c(1, 0),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing = unit(c(0, 0, 0, 0), "cm"),
    plot.background = element_rect(fill = "transparent", colour = NA),
    # plot.margin = unit(c(-1, 0, -1, 0), "cm"),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    legend.title.align = 0.5,
    strip.background = element_blank()
  )
all_seasons
save_plot(all_seasons, filename = "figures/fig_s10_all_seasons_flu_type_hhs_region_pies.png", base_width = 12, base_height = 8)
save_plot(all_seasons, filename = "figures/fig_s10_all_seasons_flu_type_hhs_region_pies.pdf", dpi=300,base_width = 12, base_height = 8)

