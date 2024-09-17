## load packages
list.of.packages <- c(
  "dplyr", "ggplot2", "cowplot", "tidyr", "broom", "gmodels", "rsample",
  "betareg", "purrr", "tidymodels", "MuMIn", "cdcfluview","rstatix"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

theme_set(theme_cowplot(font_size = 16))
########################################################################################################
## Associations between H1 or B epidemic size and H3 epidemic metrics
########################################################################################################

## load data
load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")
head(epi_red)

timing_df = epi_red %>% 
  filter(season!="2009-2010" & dom_type %in% c("H1","H3"))%>%
  dplyr::select(season,region,dom_type,onset_days_from_Oct1,peak_days_from_Oct1,
                          onset_timing_sd,peak_timing_sd,H3_season_duration)

timing_df %>%
  group_by(dom_type) %>%
  get_summary_stats(onset_days_from_Oct1,peak_days_from_Oct1, type = "common")
# dom_type variable                 n   min   max median   iqr  mean    sd    se    ci
# <chr>    <fct>                <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 H1       onset_days_from_Oct1    38    35   147    77   26.2  79.0  25.7  4.16  8.44
# 2 H1       peak_days_from_Oct1     56    28   203   144.  35   145.   33.9  4.52  9.07
# 3 H3       onset_days_from_Oct1   122    14    98    56   35    53.5  21.8  1.97  3.90
# 4 H3       peak_days_from_Oct1    122    56   196   119   49   114.   28.2  2.55  5.05

## difference in onsets 
77-56 #21 days (median)
78-53.5 #24.5 days (mean)

## difference in peaks
144-119 #25 days (median)
145-114 #31 days (mean)

onset.stat.test <- timing_df %>% 
  wilcox_test(onset_days_from_Oct1 ~ dom_type, paired = FALSE) 
onset.stat.test

wilcox.test(onset_days_from_Oct1~dom_type,data=timing_df,paired=F)

peak.stat.test <- timing_df %>% 
  wilcox_test(peak_days_from_Oct1 ~ dom_type, paired = FALSE) 
peak.stat.test

wilcox.test(peak_days_from_Oct1~dom_type,data=timing_df,paired=F)

timing_df %>%
  group_by(dom_type) %>%
  get_summary_stats(onset_timing_sd,peak_timing_sd, type = "common")%>%
  arrange(variable)
# dom_type variable            n   min   max median   iqr  mean    sd    se    ci
# <chr>    <fct>           <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 H1       onset_timing_sd    47  7.79  38.7  16.3   5.17  20.5 10.5  1.54   3.10
# 2 H3       onset_timing_sd   122  4.21  25.0   9.59 10.0   12.3  6.50 0.588  1.16
# 3 H1       peak_timing_sd     56 18.7   47.5  22.6  18.7   27.8 11.3  1.51   3.02
# 4 H3       peak_timing_sd    122  4.95  31.3  12.0   7.63  14.1  6.84 0.62   1.23

onset.stat.test <- timing_df %>% 
  wilcox_test(onset_timing_sd ~ dom_type, paired = FALSE,detailed = TRUE) 
onset.stat.test
wilcox.test(onset_timing_sd~dom_type,data=timing_df,paired=F)

peak.stat.test <- timing_df %>% 
  wilcox_test(peak_timing_sd ~ dom_type, paired = FALSE,detailed = TRUE) 
peak.stat.test
wilcox.test(peak_timing_sd~dom_type,data=timing_df,paired=F)
t = wilcox.test(peak_timing_sd~dom_type,data=timing_df,paired=F)
t$p.value #6.425816e-18

timing_df %>%
  group_by(dom_type) %>%
  get_summary_stats(H3_season_duration, type = "common")
# dom_type variable               n   min   max median   iqr  mean    sd    se    ci
# <chr>    <fct>              <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 H1       H3_season_duration    56     1    33   21.5  17.2  20.2 10.0  1.34  2.68 
# 2 H3       H3_season_duration   122    14    34   28     9    27.2  5.29 0.479 0.948

wilcox.test(H3_season_duration~dom_type,data=timing_df,paired=F)

sum_df <-
  epi_red %>%
  tidyr::replace_na(list(H1_cum_intensity = 0, H3_cum_intensity = 0, IVB_cum_intensity = 0)) %>%
  group_by(season, dom_type) %>%
  summarise(
    H1.epi.size.mean = ci(H1_cum_intensity, na.rm = T)[1],
    H1.epi.size.lowCI = ci(H1_cum_intensity, na.rm = T)[2],
    H1.epi.size.hiCI = ci(H1_cum_intensity, na.rm = T)[3],
    H1.epi.size.sd = ci(H1_cum_intensity, na.rm = T)[4],
    IVB.epi.size.mean = ci(IVB_cum_intensity, na.rm = T)[1],
    IVB.epi.size.lowCI = ci(IVB_cum_intensity, na.rm = T)[2],
    IVB.epi.size.hiCI = ci(IVB_cum_intensity, na.rm = T)[3],
    IVB.epi.size.sd = ci(IVB_cum_intensity, na.rm = T)[4],
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
    H3.R0.sd = ci(H3_max_Rt, na.rm = T)[4],
    H3.shannon.mean = ci(H3_shannon_entropy_res, na.rm = T)[1],
    H3.shannon.lowCI = ci(H3_shannon_entropy_res, na.rm = T)[2],
    H3.shannon.hiCI = ci(H3_shannon_entropy_res, na.rm = T)[3],
    H3.shannon.sd = ci(H3_shannon_entropy_res, na.rm = T)[4]
  ) %>%
  ungroup()
sum_df <- sum_df %>% tidyr::separate(col = "season", sep = "-", remove = F, into = c("year1", "year2"))

sum_df <- sum_df %>%
  mutate(dom_type2 = case_when(
    year1 < 2009 & dom_type == "H1" ~ "H1",
    year1 > 2009 & dom_type == "H1" ~ "H1pdm",
    dom_type == "H3" ~ "H3",
    dom_type == "co-circ" ~ "H3/H1pdm"
  ))

sumdf2 <- sum_df %>% filter(!(season %in% c("2009-2010", "1995-1996", "1996-1997")))

#####################################################
## H1 epi size vs H3 peak incidence
#####################################################

y <- sumdf2$H3.peak.mean
x <- sumdf2$H1.epi.size.mean

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, rank = "BIC")
## Gamma with inverse link produces wonky results; use Gaussian with inverse link instead

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.peak.mean ~ H1.epi.size.mean, analysis(split), family = gaussian(link = "inverse"), maxit = 150)
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

m1 <- glm(H3.peak.mean ~ H1.epi.size.mean, data = sumdf2, family = gaussian(link = "inverse"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.64

labels <- boots %>%
  mutate(
    model = map(splits, ~ glm(H3.peak.mean ~ H1.epi.size.mean, data = .x, family = gaussian(link = "inverse"), maxit = 150)),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(p) == %.2f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H1.epi.size.mean, IVB.epi.size.mean, H3.peak.mean) %>%
  pivot_longer(cols = c(H1.epi.size.mean, IVB.epi.size.mean)) -> tbl_mtcars_long

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
        glm(H3.peak.mean ~ value, ., family = gaussian(link = "inverse"), maxit = 200) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data


cols <- c("#B24745FF", "#00A1D5FF", "#6A6599FF", "#DF8F44FF")

h3_vs_h1_max_peak <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.peak.mean, ymin = gaussian(link = "inverse")$linkinv(l),
      ymax = gaussian(link = "inverse")$linkinv(u)
    ),
    data = tbl_plot_data %>% filter(name == "H1.epi.size.mean") %>% rename(H3.peak.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = H1.epi.size.mean,
      ymin = H3.peak.lowCI, ymax = H3.peak.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.peak.mean,
      xmin = H1.epi.size.lowCI, xmax = H1.epi.size.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.peak.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("A(H1N1) Epidemic Size") +
  ylab("A(H3N2) Peak Incidence") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.peak.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "inverse")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 17, x = 1, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_h1_max_peak

#####################################################
## B epi size vs H3 peak incidence
#####################################################

y <- sumdf2$H3.peak.mean
x <- sumdf2$IVB.epi.size.mean

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, rank = "BIC")
## Gamma with inverse link produces wonky results; use Gaussian with inverse link instead

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

boot_models <-
  boots %>%
  mutate(
    model = map(splits, ~ glm(H3.peak.mean ~ IVB.epi.size.mean, data = ., family = gaussian(link = "inverse"), maxit = 150)),
    coef_info = map(model, tidy)
  )

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

labels <- boots %>%
  mutate(
    model = map(splits, ~ glm(H3.peak.mean ~ IVB.epi.size.mean, data = .x, family = gaussian(link = "inverse"), maxit = 150)),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.1f", adj.r.squared),
    pvalue = sprintf("italic(p) == %.2f", pvalue)
  )
labels

m1 <- glm(H3.peak.mean ~ IVB.epi.size.mean, data = sumdf2, family = gaussian(link = "inverse"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv
1 - m1$deviance / m1$null.deviance # 0.29

h3_vs_ivb_max_peak <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.peak.mean, ymin = gaussian(link = "inverse")$linkinv(l),
      ymax = gaussian(link = "inverse")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "IVB.epi.size.mean") %>% rename(H3.peak.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = IVB.epi.size.mean,
      ymin = H3.peak.lowCI, ymax = H3.peak.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.peak.mean,
      xmin = IVB.epi.size.lowCI, xmax = IVB.epi.size.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = IVB.epi.size.mean, y = H3.peak.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("B Epidemic Size") +
  ylab("A(H3N2) Peak Incidence") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = IVB.epi.size.mean, y = H3.peak.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "inverse")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 17, x = 1, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_ivb_max_peak

#####################################################
## H1 epi size vs H3 epi size
#####################################################
y <- sumdf2$H3.epi.size.mean
x <- sumdf2$H1.epi.size.mean

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))

model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, rank = "BIC")

fit_glm_on_bootstrap <- function(split) {
  glm(H3.epi.size.mean ~ H1.epi.size.mean, analysis(split), family = gaussian(link = "inverse"), maxit = 150)
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

m1 <- glm(H3.epi.size.mean ~ H1.epi.size.mean, data = sumdf2, family = gaussian(link = "inverse"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv
1 - m1$deviance / m1$null.deviance # 0.65

labels <- boots %>%
  mutate(
    model = map(splits, ~ glm(H3.epi.size.mean ~ H1.epi.size.mean, data = .x, family = gaussian(link = "inverse"), maxit = 150)),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(p) == %.2f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H1.epi.size.mean, IVB.epi.size.mean, H3.epi.size.mean) %>%
  pivot_longer(cols = c(H1.epi.size.mean, IVB.epi.size.mean)) -> tbl_mtcars_long

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
        glm(H3.epi.size.mean ~ value, ., family = gaussian(link = "inverse"), maxit = 150) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

h3_vs_h1_epi_size <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.epi.size.mean, ymin = gaussian(link = "inverse")$linkinv(l),
      ymax = gaussian(link = "inverse")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "H1.epi.size.mean") %>% rename(H3.epi.size.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = H1.epi.size.mean,
      ymin = H3.epi.size.lowCI, ymax = H3.epi.size.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.epi.size.mean,
      xmin = H1.epi.size.lowCI, xmax = H1.epi.size.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.epi.size.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("A(H1N1) Epidemic Size") +
  ylab("A(H3N2) Epidemic Size") +
  theme(
    legend.title = element_blank(),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.position = "bottom"
  ) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.epi.size.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "inverse")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 80, x = 1, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_h1_epi_size

#####################################################
## B epi size vs H3 epi size
#####################################################

y <- sumdf2$H3.epi.size.mean
x <- sumdf2$IVB.epi.size.mean

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, rank = "BIC")

fit_glm_on_bootstrap <- function(split) {
  glm(H3.epi.size.mean ~ IVB.epi.size.mean, analysis(split), family = gaussian(link = "inverse"), maxit = 150)
}

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

boot_models <-
  boots %>%
  mutate(
    model = map(splits, fit_glm_on_bootstrap),
    coef_info = map(model, tidy)
  )

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

m1 <- glm(H3.epi.size.mean ~ IVB.epi.size.mean, data = sumdf2, family = gaussian(link = "inverse"))
fam <- family(m1)
ilink <- fam$linkinv
1 - m1$deviance / m1$null.deviance # 0.05

labels <- boots %>%
  mutate(
    model = map(splits, ~ glm(H3.epi.size.mean ~ IVB.epi.size.mean, data = .x, family = gaussian(link = "inverse"), maxit = 150)),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(p) == %.1f", pvalue)
  )
labels

h3_vs_ivb_epi_size <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.epi.size.mean,
      ymin = gaussian(link = "inverse")$linkinv(l), ymax = gaussian(link = "inverse")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "IVB.epi.size.mean") %>% rename(H3.epi.size.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = IVB.epi.size.mean,
      ymin = H3.epi.size.lowCI, ymax = H3.epi.size.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.epi.size.mean,
      xmin = IVB.epi.size.lowCI, xmax = IVB.epi.size.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = IVB.epi.size.mean, y = H3.epi.size.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("B Epidemic Size") +
  ylab("A(H3N2) Epidemic Size") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = IVB.epi.size.mean, y = H3.epi.size.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "inverse")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 80, x = 1, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_ivb_epi_size

#####################################################
## H1 epi size vs H3 Rt
#####################################################
df <- sumdf2 %>%
  dplyr::select(season, H3.R0.mean, H1.epi.size.mean, IVB.epi.size.mean) %>%
  filter(season != "2000-2001") ### no H3 Rt for 2000-2001

y <- df$H3.R0.mean
x <- df$H1.epi.size.mean

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))

model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, rank = "BIC")
# : Gamma, log link (very little difference across models, "inverse" link produces wonky results)

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.R0.mean ~ H1.epi.size.mean, analysis(split), family = Gamma(link = "log"), maxit = 150)
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

m1 <- glm(H3.R0.mean ~ H1.epi.size.mean, data = sumdf2, family = Gamma(link = "log"))
fam <- family(m1)
ilink <- fam$linkinv
1 - m1$deviance / m1$null.deviance # 0.45

labels <- boots %>%
  mutate(
    model = map(splits, ~ glm(H3.R0.mean ~ H1.epi.size.mean, data = .x, family = Gamma(link = "log"), maxit = 150)),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(p) == %.2f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H1.epi.size.mean, IVB.epi.size.mean, H3.R0.mean) %>%
  pivot_longer(cols = c(H1.epi.size.mean, IVB.epi.size.mean)) -> tbl_mtcars_long

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
    # the _sumdf2r means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>%
        sample_n(n, TRUE) %>%
        glm(H3.R0.mean ~ value, ., family = Gamma(link = "log"), maxit = 150) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

## no H3 R0 for 2000-2001
h3_vs_h1_r0 <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.R0.mean, ymin = Gamma(link = "log")$linkinv(l),
      ymax = Gamma(link = "log")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "H1.epi.size.mean") %>% rename(H3.R0.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = H1.epi.size.mean,
      ymin = H3.R0.lowCI, ymax = H3.R0.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.R0.mean,
      xmin = H1.epi.size.lowCI, xmax = H1.epi.size.hiCI,
      color = dom_type2
    ), width = .0001
  ) +
  geom_point(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.R0.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("A(H1N1) Epidemic Size") +
  ylab("A(H3N2) Effective Rt") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.R0.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = Gamma(link = "log")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 2, x = 0.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_h1_r0

#####################################################
## B epi size vs H3 Rt
#####################################################
sumdf2 %>%
  filter(season == "2000-2001") %>%
  dplyr::select(H3.R0.mean, H3.shannon.mean)

df <- sumdf2 %>% filter(season != "2000-2001")

y <- df$H3.R0.mean
x <- df$IVB.epi.size.mean

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, rank = "BIC")
# no difference between gaussian log, inv, and identity link

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.R0.mean ~ IVB.epi.size.mean, analysis(split), family = Gamma(link = "log"), maxit = 150)
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

m1 <- glm(H3.R0.mean ~ IVB.epi.size.mean, data = sumdf2, family = Gamma(link = "log"))
fam <- family(m1)
ilink <- fam$linkinv
1 - m1$deviance / m1$null.deviance # 0.004

labels <- boots %>%
  mutate(
    model = map(splits, ~ glm(H3.R0.mean ~ IVB.epi.size.mean, data = .x, family = Gamma(link = "log"), maxit = 150)),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(p) == %.1f", pvalue)
  )
labels

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
        glm(H3.R0.mean ~ value, ., family = Gamma(link = "log"), maxit = 150) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted values
      summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

h3_vs_ivb_r0 <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.R0.mean,
      ymin = Gamma(link = "log")$linkinv(l), ymax = Gamma(link = "log")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "IVB.epi.size.mean") %>% rename(H3.R0.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = IVB.epi.size.mean,
      ymin = H3.R0.lowCI, ymax = H3.R0.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.R0.mean,
      xmin = IVB.epi.size.lowCI, xmax = IVB.epi.size.hiCI,
      color = dom_type2
    ), width = .001
  ) +
  geom_point(
    data = sumdf2,
    aes(x = IVB.epi.size.mean, y = H3.R0.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("B Epidemic Size") +
  ylab("A(H3N2) Effective Rt") +
  theme(legend.position = "bottom", legend.justification = "center") +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = IVB.epi.size.mean, y = H3.R0.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = Gamma(link = "log")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 2, x = 0.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )

h3_vs_ivb_r0

#####################################################
## H1 epi size vs H3 epidemic intensity (inverse shannon entropy)
#####################################################

y <- sumdf2$H3.shannon.mean
range(y, na.rm = T)
x <- sumdf2$H1.epi.size.mean

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
quasibinomial.model <- glm(y ~ x, family = quasibinomial(link = "logit"))
beta.model <- betareg(y ~ x, link = c("logit"))

model.sel(
  linear.model, log.model, inv.model, gamma.model, gamma.model2,
  quasibinomial.model, beta.model,
  rank = "BIC"
)
## beta fits best

fit_glm_on_bootstrap <- function(split) {
  betareg(H3.shannon.mean ~ H1.epi.size.mean, analysis(split))
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

m1 <- betareg(H3.shannon.mean ~ H1.epi.size.mean, data = sumdf2)
fam <- family(m1)
ilink <- fam$linkinv

labels <- boots %>%
  mutate(
    model = map(splits, ~ betareg(H3.shannon.mean ~ H1.epi.size.mean, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(p) == %.2f", pvalue)
  )
labels

set.seed(27)
n_boot <- 1000

sumdf2 %>%
  dplyr::select(H1.epi.size.mean, IVB.epi.size.mean, H3.shannon.mean) %>%
  pivot_longer(cols = c(H1.epi.size.mean, IVB.epi.size.mean)) -> tbl_mtcars_long

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
      summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

h3_vs_h1_shannon <- ggplot() +
  geom_ribbon(aes(x = H3.shannon.mean, ymin = l, ymax = u),
    tbl_plot_data %>% filter(name == "H1.epi.size.mean") %>% rename(H3.shannon.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = H1.epi.size.mean,
      ymin = H3.shannon.lowCI, ymax = H3.shannon.hiCI,
      color = dom_type2
    ), width = .003
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.shannon.mean,
      xmin = H1.epi.size.lowCI, xmax = H1.epi.size.hiCI,
      color = dom_type2
    ), width = .003
  ) +
  geom_point(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.shannon.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  xlab("A(H1N1) Epidemic Size") +
  ylab("A(H3N2) Epidemic Intensity") +
  theme(
    legend.title = element_blank(),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.position = "bottom"
  ) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.shannon.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = Gamma(link = "log")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 0.8, x = 1, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_h1_shannon

#####################################################
## B epi size vs H3 epidemic intensity
#####################################################

y <- sumdf2$H3.shannon.mean
x <- sumdf2$IVB.epi.size.mean

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
quasibinomial.model <- glm(y ~ x, family = quasibinomial(link = "logit"))
beta.model <- betareg(y ~ x, link = c("logit"))

model.sel(
  linear.model, log.model, inv.model, gamma.model, gamma.model2,
  quasibinomial.model, beta.model,
  rank = "BIC"
)

fit_glm_on_bootstrap <- function(split) {
  betareg(H3.shannon.mean ~ IVB.epi.size.mean, analysis(split))
}

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

boot_models <-
  boots %>%
  mutate(
    model = map(splits, fit_glm_on_bootstrap),
    coef_info = map(model, tidy)
  )

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

m1 <- betareg(H3.shannon.mean ~ IVB.epi.size.mean, data = sumdf2)
fam <- family(m1)
ilink <- fam$linkinv

labels <- boots %>%
  mutate(
    model = map(splits, ~ betareg(H3.shannon.mean ~ IVB.epi.size.mean, data = .x)),
    adj.r.squared = map_dbl(model, ~ signif(summary(.x)$pseudo.r.squared, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients$mean[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(p) == %.1f", pvalue)
  )
labels

h3_vs_ivb_shannon <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.shannon.mean,
      ymin = l, ymax = u
    ),
    tbl_plot_data %>% filter(name == "IVB.epi.size.mean") %>%
      rename(H3.shannon.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = IVB.epi.size.mean,
      ymin = H3.shannon.lowCI, ymax = H3.shannon.hiCI,
      color = dom_type2
    ), width = .06
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.shannon.mean,
      xmin = IVB.epi.size.lowCI, xmax = IVB.epi.size.hiCI,
      color = dom_type2
    ), width = .005
  ) +
  geom_point(
    data = sumdf2,
    aes(x = IVB.epi.size.mean, y = H3.shannon.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab("B Epidemic Size") +
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
  geom_smooth(
    data = sumdf2,
    aes(x = IVB.epi.size.mean, y = H3.shannon.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = Gamma(link = "inverse")),
    se = F,
    linewidth = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 0.8, x = 1, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_ivb_shannon

combined1 <- plot_grid(h3_vs_h1_epi_size + theme(legend.position = "none"),
  h3_vs_h1_max_peak + theme(legend.position = "none"),
  h3_vs_h1_r0 + theme(legend.position = "none"),
  h3_vs_h1_shannon + theme(legend.position = "none"),
  nrow = 1, labels = NULL
)

combined2 <- plot_grid(h3_vs_ivb_epi_size + theme(legend.position = "none"),
  h3_vs_ivb_max_peak + theme(legend.position = "none"),
  h3_vs_ivb_r0 + theme(legend.position = "none"),
  h3_vs_ivb_shannon + theme(legend.position = "none"),
  labels = NULL, nrow = 1
)
combined3 <- plot_grid(combined1, combined2, nrow = 2, labels = "AUTO")
leg <- cowplot::get_plot_component(h3_vs_ivb_shannon +
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
cowplot::ggdraw(leg)

combined3 <- plot_grid(combined3, leg, nrow = 2, rel_heights = c(2, 0.1))
combined3

save_plot(combined3,
  filename = "figures/Fig7_h1_or_ivb_epi_size_vs_h3_parameters.png",
  dpi = 300, bg = "white",
  base_width = 20, base_height = 10,
)

# save_plot(combined3,
#           filename = "figures/Fig7_h1_or_ivb_epi_size_vs_h3_paramters.pdf",
#           dpi = 300,
#           base_width = 20, base_height = 10
# )

#####################################################
## Wavelet Results
#####################################################
##### phase difference: lead and lag
### load wavelet results
path_file <- "data"
filenames_list <- list.files(path = path_file, full.names = TRUE, pattern = "h3_vs_h1.csv")
head(filenames_list)
merging_manifests <- lapply(filenames_list, function(filename) {
  print(paste("Merging", filename, sep = " "))
  read.csv(filename)
})
master_file <- do.call(plyr::rbind.fill, merging_manifests)
head(master_file)
master_file$comparison <- "h3_vs_h1"

path_file <- "data"
filenames_list <- list.files(path = path_file, full.names = TRUE, pattern = "h3_vs_ivb.csv")
head(filenames_list)
merging_manifests <- lapply(filenames_list, function(filename) {
  print(paste("Merging", filename, sep = " "))
  read.csv(filename)
})
master_file2 <- do.call(plyr::rbind.fill, merging_manifests)
head(master_file2)
master_file2$comparison <- "h3_vs_ivb"

all_lag_df <- bind_rows(master_file, master_file2)
mwk <- mmwr_week(as.Date(all_lag_df$date))[2]
all_lag_df$epi_week <- mwk$mmwr_week

lag_summary <- all_lag_df %>%
  mutate(flu_season = if_else(epi_week >= 40 | epi_week <= 20, "flu", "not_flu")) %>%
  filter(flu_season == "flu") %>%
  group_by(season, comparison, region) %>%
  summarize(
    mean_lag_radians = mean(diff),
    mean_lag_weeks = mean(diff_weeks)
  ) %>%
  ungroup()

sum_df2 <- left_join(lag_summary, epi_red %>% dplyr::select(
  region, season, dom_type, H1_cum_intensity, H3_cum_intensity,
  IVB_cum_intensity
), by = c("region", "season"))

sum_df3 <-
  sum_df2 %>%
  filter(season != "2009-2010") %>%
  tidyr::replace_na(list(
    H1_cum_intensity = 0,
    H3_cum_intensity = 0, IVB_cum_intensity = 0
  )) %>%
  group_by(season, dom_type, comparison) %>%
  summarise(
    H1.epi.size.mean = ci(H1_cum_intensity, na.rm = T)[1],
    H1.epi.size.lowCI = ci(H1_cum_intensity, na.rm = T)[2],
    H1.epi.size.hiCI = ci(H1_cum_intensity, na.rm = T)[3],
    H1.epi.size.sd = ci(H1_cum_intensity, na.rm = T)[4],
    IVB.epi.size.mean = ci(IVB_cum_intensity, na.rm = T)[1],
    IVB.epi.size.lowCI = ci(IVB_cum_intensity, na.rm = T)[2],
    IVB.epi.size.hiCI = ci(IVB_cum_intensity, na.rm = T)[3],
    IVB.epi.size.sd = ci(IVB_cum_intensity, na.rm = T)[4],
    H3.epi.size.mean = ci(H3_cum_intensity, na.rm = T)[1],
    H3.epi.size.lowCI = ci(H3_cum_intensity, na.rm = T)[2],
    H3.epi.size.hiCI = ci(H3_cum_intensity, na.rm = T)[3],
    H3.epi.size.sd = ci(H3_cum_intensity, na.rm = T)[4],
    lag.weeks.mean = ci(mean_lag_weeks, na.rm = T)[1],
    lag.weeks.lowCI = ci(mean_lag_weeks, na.rm = T)[2],
    lag.weeks.hiCI = ci(mean_lag_weeks, na.rm = T)[3],
    lag.weeks.sd = ci(mean_lag_weeks, na.rm = T)[4],
    lag.radians.mean = ci(mean_lag_radians, na.rm = T)[1],
    lag.radians.lowCI = ci(mean_lag_radians, na.rm = T)[2],
    lag.radians.hiCI = ci(mean_lag_radians, na.rm = T)[3],
    lag.radians.sd = ci(mean_lag_radians, na.rm = T)[4]
  ) %>%
  ungroup()
sum_df3 <- sum_df3 %>% tidyr::separate(col = "season", sep = "-", remove = F, into = c("year1", "year2"))

sum_df3 <- sum_df3 %>%
  mutate(dom_type2 = case_when(
    year1 < 2009 & dom_type == "H1" ~ "H1",
    year1 > 2009 & dom_type == "H1" ~ "H1pdm",
    dom_type == "H3" ~ "H3",
    dom_type == "co-circ" ~ "H3/H1pdm"
  ))

#####################################################
# B epi size vs H3/B timing
#####################################################
h3_lag_vs_ivb <- sum_df3 %>% filter(comparison == "h3_vs_ivb")

y <- h3_lag_vs_ivb$lag.weeks.mean
x <- h3_lag_vs_ivb$IVB.epi.size.mean

linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2, rank = "BIC")

fit_glm_on_bootstrap <- function(split) {
  glm(lag.weeks.mean ~ IVB.epi.size.mean, analysis(split), family = gaussian(link = "identity"))
}

set.seed(27)
boots <- bootstraps(h3_lag_vs_ivb, times = 1000, apparent = TRUE)
boots

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
    model = map(splits, ~ glm(lag.weeks.mean ~ IVB.epi.size.mean, data = .x, family = gaussian(link = "identity"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(p) == %.1f", pvalue)
  )
labels

# a data frame with all the annotation info
annotation <- data.frame(
  x = c(15, 15),
  y = c(8, -8),
  label = c("H3 leads", "B leads")
)

h3_lag_vs_ivb_epi_size_fig <- ggplot() +
  geom_point(
    data = h3_lag_vs_ivb %>% filter(season != "2009-2010"),
    aes(x = IVB.epi.size.mean, y = lag.weeks.mean, fill = dom_type2), size = 5, pch = 21
  ) +
  geom_errorbar(
    data = h3_lag_vs_ivb %>% filter(season != "2009-2010"),
    aes(
      x = IVB.epi.size.mean,
      ymin = lag.weeks.lowCI, ymax = lag.weeks.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_errorbar(
    data = h3_lag_vs_ivb %>% filter(season != "2009-2010"),
    aes(
      y = lag.weeks.mean,
      xmin = IVB.epi.size.lowCI, xmax = IVB.epi.size.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  xlab("B Epidemic Size") +
  ylab("Phase Lag (weeks)") +
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
  geom_hline(yintercept = 0, lty = "dashed", lwd = 1) +
  ylim(c(-8, 8)) +
  geom_text(data = annotation, aes(x = x, y = y, label = label), size = 6) +
  geom_text(
    y = 6, x = 20, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_lag_vs_ivb_epi_size_fig

#####################################################
# H1 epi size vs H3/H1 timing
#####################################################

h3_lag_vs_h1 <- sum_df3 %>% filter(comparison == "h3_vs_h1")

y <- h3_lag_vs_h1$lag.weeks.mean
x <- h3_lag_vs_h1$H1.epi.size.mean
linear.model <- glm(y ~ x, family = gaussian())
# log.model = glm(y~x,family=gaussian(link="log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
# gamma.model = glm(y~x,family=Gamma(link="log"))
# gamma.model2 = glm(y~x,family=Gamma(link="inverse"))
# gamma.model3 = glm(y~x,family=Gamma(link="identity"))
model.sel(linear.model, inv.model, rank = "BIC")

fit_glm_on_bootstrap <- function(split) {
  glm(lag.weeks.mean ~ H1.epi.size.mean, analysis(split), family = gaussian(link = "identity"))
}

set.seed(27)
boots <- bootstraps(h3_lag_vs_h1, times = 1000, apparent = TRUE)

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
    model = map(splits, ~ glm(lag.weeks.mean ~ H1.epi.size.mean, data = .x, family = gaussian(link = "identity"))),
    adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance / .x$null.deviance, 5)),
    pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))
  ) %>%
  summarize(
    adj.r.squared = mean(adj.r.squared),
    pvalue = mean(pvalue)
  ) %>%
  mutate(
    adj.r.squared = sprintf("italic(R^2) == %.2f", adj.r.squared),
    pvalue = sprintf("italic(p) == %.1f", pvalue)
  )
labels

# a data frame with all the annotation info
annotation2 <- data.frame(
  x = c(20, 20),
  y = c(8, -8),
  label = c("H3 leads", "H1 leads")
)

h3_lag_vs_h1_epi_size_fig <- ggplot() +
  geom_errorbar(
    data = h3_lag_vs_h1,
    aes(
      x = H1.epi.size.mean,
      ymin = lag.weeks.lowCI, ymax = lag.weeks.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_errorbar(
    data = h3_lag_vs_h1,
    aes(
      y = lag.weeks.mean,
      xmin = H1.epi.size.lowCI, xmax = H1.epi.size.hiCI,
      color = dom_type2
    ), width = .05
  ) +
  geom_point(
    data = h3_lag_vs_h1,
    aes(x = H1.epi.size.mean, y = lag.weeks.mean, fill = dom_type2), pch = 21, size = 5
  ) +
  xlab("A(H1N1) Epidemic Size") +
  ylab("Phase Lag (weeks)") +
  theme(legend.position = "bottom", legend.justification = "center") +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "H1pdm", "H3/H1pdm"),
    labels = c("H3", "H1", "H1pdm", "H3/H1pdm"), values = cols, name = "Dominant IAV"
  ) +
  geom_hline(yintercept = 0, lty = "dashed", lwd = 1) +
  ylim(c(-8, 8)) +
  geom_text(data = annotation2, aes(x = x, y = y, label = label), size = 6) +
  geom_text(
    y = 6, x = 25, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_lag_vs_h1_epi_size_fig

combined <- plot_grid(h3_lag_vs_h1_epi_size_fig + theme(legend.position = "none"),
  h3_lag_vs_ivb_epi_size_fig + theme(legend.position = "none"),
  labels = "AUTO", nrow = 1
)
combined

leg <- cowplot::get_plot_component(h3_lag_vs_h1_epi_size_fig+
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
cowplot::ggdraw(leg)

# leg <- get_legend(h3_lag_vs_h1_epi_size_fig)

combined2 <- plot_grid(combined, leg, rel_heights = c(1, 0.1), nrow = 2)
combined2
save_plot(combined2, filename = "figures/Fig7_sup_fig3_h3_vs_h1_or_ivb_epi_size_wavelet.png", 
          base_width = 12, base_height = 6, dpi = 300, bg="white")
# save_plot(combined2, filename = "figures/Fig7_sup_fig3_h3_vs_h1_or_ivb_epi_size_wavelet.pdf", dpi = 300, base_width = 12, base_height = 6)
