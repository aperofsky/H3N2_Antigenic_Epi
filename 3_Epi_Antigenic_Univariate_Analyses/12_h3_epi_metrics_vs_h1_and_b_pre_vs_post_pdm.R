## load packages
list.of.packages <- c(
  "dplyr", "ggplot2", "cowplot", "tidyr", "broom", "gmodels", "rsample",
  "betareg", "purrr", "tidymodels", "MuMIn", "cdcfluview"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

theme_set(theme_cowplot(font_size = 16))

#######################################################################################################
## Associations between H1 or B epidemic size and H3 epidemic metrics, pre- and post-2009 pandemic
########################################################################################################

## load data
load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")

sum_df <-
  epi_red %>%
  filter(season != "2009-2010") %>%
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

#####################################################
### ALL SEASONS
#####################################################
#####################################################
## H1 epi size vs H3 peak incidence
#####################################################
sumdf2 <- sum_df %>% filter(season != "2009-2010")

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
    # legend.title = element_blank(),
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
    y = 0.81, x = 1, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_h1_shannon


combined_all_seasons <- plot_grid(h3_vs_h1_epi_size + theme(legend.position = "none"),
  h3_vs_h1_max_peak + theme(legend.position = "none"),
  h3_vs_h1_r0 + theme(legend.position = "none"),
  h3_vs_h1_shannon + theme(legend.position = "none"),
  labels = NULL, nrow = 1
)
combined_all_seasons

leg <- cowplot::get_plot_component(
  h3_vs_h1_shannon +
    guides(color = "none") +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.box.just = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    ),
  "guide-box-bottom",
  return_all = TRUE
)


#####################################################
## Peak incidence: Gamma distribution, inverse link
#####################################################
sumdf2 <- sum_df %>% filter(h1n1_type == "seasonal_h1n1")
length(unique(sumdf2$season)) # 12 seasons

y <- sumdf2$H3.peak.mean
x <- sumdf2$H1.epi.size.mean
linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
# gamma.model3 = glm(y~x,family=Gamma(link="identity"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2)

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
1 - m1$deviance / m1$null.deviance # 0.7

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
  # for mpg and value observations within each level of name (e.g., disp, hp, ...)
  mutate(plot_data = map(model_data, ~ {
    # calculate information about the observed mpg and value observations
    # within each level of name to be used in each bootstrap sample
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
        glm(H3.peak.mean ~ value, ., family = gaussian(link = "inverse")) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted mpgs
      summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data


cols <- c("#DF8F44FF", "#00A1D5FF", "#B24745FF")
unique(sumdf2$dom_type)

h3_vs_h1_max_peak_pre_pdm <- ggplot() +
  geom_ribbon(aes(x = H3.peak.mean, ymin = gaussian(link = "inverse")$linkinv(l), ymax = gaussian(link = "inverse")$linkinv(u)),
    tbl_plot_data %>% filter(name == "H1.epi.size.mean") %>% rename(H3.peak.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = H1.epi.size.mean,
      ymin = H3.peak.lowCI, ymax = H3.peak.hiCI,
      color = dom_type
    ), width = .05
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.peak.mean,
      xmin = H1.epi.size.lowCI, xmax = H1.epi.size.hiCI,
      color = dom_type
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.peak.mean, fill = dom_type), size = 5, pch = 21
  ) +
  xlab("A(H1N1) Epidemic Size") +
  ylab("A(H3N2) Peak Incidence") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(labels = c("H3/H1", "H1", "H3"), breaks = c("H3/H1", "H1", "H3"), values = cols) +
  scale_fill_manual(labels = c("H3/H1", "H1", "H3"), breaks = c("H3/H1", "H1", "H3"), values = cols) +
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
h3_vs_h1_max_peak_pre_pdm

#####################################################
## Epi size, Gaussian with inverse link
#####################################################
sumdf2 <- sum_df %>% filter(h1n1_type == "seasonal_h1n1")
y <- sumdf2$H3.epi.size.mean
x <- sumdf2$H1.epi.size.mean
linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2)

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
  # for mpg and value observations within each level of name (e.g., disp, hp, ...)
  mutate(plot_data = map(model_data, ~ {
    # calculate information about the observed mpg and value observations
    # within each level of name to be used in each bootstrap sample
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
      # to estimate the lower and upper 95% quantiles of predicted mpgs
      summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

h3_vs_h1_epi_size_pre_pdm <- ggplot() +
  geom_ribbon(aes(x = H3.epi.size.mean, ymin = gaussian(link = "inverse")$linkinv(l), ymax = gaussian(link = "inverse")$linkinv(u)),
    tbl_plot_data %>% filter(name == "H1.epi.size.mean") %>% rename(H3.epi.size.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  # geom_line(aes(y = gaussian(link="inverse")$linkinv(.fitted), x=H1.epi.size.mean, group = id), alpha = .2, col = "grey",data=boot_aug) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = H1.epi.size.mean,
      ymin = H3.epi.size.lowCI, ymax = H3.epi.size.hiCI,
      color = dom_type
    ), width = .05
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.epi.size.mean,
      xmin = H1.epi.size.lowCI, xmax = H1.epi.size.hiCI,
      color = dom_type
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.epi.size.mean, fill = dom_type), size = 5, pch = 21
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
  scale_color_manual(labels = c("H3/H1", "H1", "H3"), breaks = c("H3/H1", "H1", "H3"), values = cols) +
  scale_fill_manual(labels = c("H3/H1", "H1", "H3"), breaks = c("H3/H1", "H1", "H3"), values = cols) +
  geom_smooth(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.epi.size.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "inverse")),
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 80, x = 1, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_h1_epi_size_pre_pdm

#####################################################
## R0: Gamma, log link (very little difference across models, "inverse" link produces wonky results)
#####################################################
df <- sumdf2 %>%
  dplyr::select(season, H3.R0.mean, H1.epi.size.mean, IVB.epi.size.mean) %>%
  filter(season != "2000-2001") ### no H3 Rt for 2000-2001

unique(df$season)

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
1 - m1$deviance / m1$null.deviance # 0.61

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
h3_vs_h1_r0_pre_pdm <- ggplot() +
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
  scale_color_manual(labels = c("H3/H1", "H1", "H3"), breaks = c("H3/H1", "H1", "H3"), values = cols) +
  scale_fill_manual(labels = c("H3/H1", "H1", "H3"), breaks = c("H3/H1", "H1", "H3"), values = cols) +
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
h3_vs_h1_r0_pre_pdm


#####################################################
## Epidemic Intensity, Beta distribution
#####################################################

sumdf2 <- sumdf2 %>% filter(!(season %in% c("2000-2001", "2009-2010")))
range(sumdf2$H3.shannon.mean)

y <- sumdf2$H3.shannon.mean
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
  quasibinomial.model, beta.model
)

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
ilink

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
    pvalue = sprintf("italic(p) == %.1f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H1.epi.size.mean, IVB.epi.size.mean, H3.shannon.mean) %>%
  pivot_longer(cols = c(H1.epi.size.mean, IVB.epi.size.mean)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.shannon.mean, value)) %>%
  # for mpg and value observations within each level of name (e.g., disp, hp, ...)
  mutate(plot_data = map(model_data, ~ {
    # calculate information about the observed mpg and value observations
    # within each level of name to be used in each bootstrap sample
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
      # to estimate the lower and upper 95% quantiles of predicted mpgs
      summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

h3_vs_h1_shannon_pre_pdm <- ggplot() +
  geom_ribbon(aes(x = H3.shannon.mean, ymin = l, ymax = u),
    tbl_plot_data %>% filter(name == "H1.epi.size.mean") %>% rename(H3.shannon.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  # geom_line(aes(y = gaussian(link="inverse")$linkinv(.fitted), x=H1.epi.size.mean, group = id), alpha = .2, col = "grey",data=boot_aug) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = H1.epi.size.mean,
      ymin = H3.shannon.lowCI, ymax = H3.shannon.hiCI,
      color = dom_type
    ), width = .003
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.shannon.mean,
      xmin = H1.epi.size.lowCI, xmax = H1.epi.size.hiCI,
      color = dom_type
    ), width = .003
  ) +
  geom_point(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.shannon.mean, fill = dom_type), size = 5, pch = 21
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
  scale_color_manual(labels = c("H3/H1", "H1", "H3"), breaks = c("H3/H1", "H1", "H3"), values = cols) +
  scale_fill_manual(labels = c("H3/H1", "H1", "H3"), breaks = c("H3/H1", "H1", "H3"), values = cols) +
  geom_smooth(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.shannon.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = Gamma(link = "log")),
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 0.82, x = 1, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_h1_shannon_pre_pdm

combined_pre_pdm <- plot_grid(h3_vs_h1_epi_size_pre_pdm + theme(legend.position = "none"),
  h3_vs_h1_max_peak_pre_pdm + theme(legend.position = "none"),
  h3_vs_h1_r0_pre_pdm + theme(legend.position = "none"),
  h3_vs_h1_shannon_pre_pdm + theme(legend.position = "none"),
  labels = NULL, nrow = 1
)
combined_pre_pdm

###############################################
## post pandemic
sumdf2 <- sum_df %>% filter(h1n1_type == "pdm_h1n1")
length(unique(sumdf2$season)) # 9 seasons

y <- sumdf2$H3.peak.mean
x <- sumdf2$H1.epi.size.mean
linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
# gamma.model3 = glm(y~x,family=Gamma(link="identity"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2)

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.peak.mean ~ H1.epi.size.mean, analysis(split), family = gaussian(link = "log"), maxit = 200)
}

boot_models <-
  boots %>%
  mutate(
    model = map(splits, fit_glm_on_bootstrap),
    coef_info = map(model, tidy)
  )

# boot_models <-
#   boots %>%
#   mutate(model = map(splits, ~glm(H3.peak.mean ~ H1.epi.size.mean,data=.,family=gaussian(link="log"))),
#          coef_info = map(model, tidy))

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

# boot_aug <-
#   boot_models %>%
#   sample_n(200) %>%
#   mutate(augmented = map(model, augment)) %>%
#   unnest(augmented)

# library(RColorBrewer)
# col_pal = brewer.pal(11,"RdBu")
# col_pal
# plot(1:length(col_pal), 1:length(col_pal), col = col_pal, cex = 5, pch=19)
# cols = col_pal[c(4,2,9)]

m1 <- glm(H3.peak.mean ~ H1.epi.size.mean, data = sumdf2, family = gaussian(link = "log"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance / m1$null.deviance # 0.47

labels <- boots %>%
  mutate(
    model = map(splits, ~ glm(H3.peak.mean ~ H1.epi.size.mean, data = .x, family = gaussian(link = "log"), maxit = 150)),
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
  # for mpg and value observations within each level of name (e.g., disp, hp, ...)
  mutate(plot_data = map(model_data, ~ {
    # calculate information about the observed mpg and value observations
    # within each level of name to be used in each bootstrap sample
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
        glm(H3.peak.mean ~ value, ., family = gaussian(link = "log"), maxit = 150) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted mpgs
      summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

cols <- c("#B24745FF", "#6A6599FF", "#DF8F44FF")

h3_vs_h1_max_peak_post_pdm <- ggplot() +
  geom_ribbon(aes(x = H3.peak.mean, ymin = gaussian(link = "log")$linkinv(l), ymax = gaussian(link = "log")$linkinv(u)),
    tbl_plot_data %>% filter(name == "H1.epi.size.mean") %>% rename(H3.peak.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = H1.epi.size.mean,
      ymin = H3.peak.lowCI, ymax = H3.peak.hiCI,
      color = dom_type
    ), width = .05
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.peak.mean,
      xmin = H1.epi.size.lowCI, xmax = H1.epi.size.hiCI,
      color = dom_type
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.peak.mean, fill = dom_type), size = 5, pch = 21
  ) +
  xlab("A(H1N1) Epidemic Size") +
  ylab("A(H3N2) Peak Incidence") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "co-circ"),
    labels = c("H3", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "co-circ"),
    labels = c("H3", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  # scale_color_manual(labels=c("H3/H1","H1","H3"),breaks=c("co-circ","H1","H3"),values=cols)+
  # scale_fill_manual(labels=c("H3/H1","H1","H3"),breaks=c("co-circ","H1","H3"),values=cols)+
  geom_smooth(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.peak.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "log")),
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 10, x = 1, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_h1_max_peak_post_pdm

#####################################################
## Epi size, Gaussian with inverse link
#####################################################
y <- sumdf2$H3.epi.size.mean
x <- sumdf2$H1.epi.size.mean
linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2)

fit_glm_on_bootstrap <- function(split) {
  glm(H3.epi.size.mean ~ H1.epi.size.mean, analysis(split), family = gaussian(link = "log"), maxit = 150)
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

m1 <- glm(H3.epi.size.mean ~ H1.epi.size.mean, data = sumdf2, family = gaussian(link = "log"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv

labels <- boots %>%
  mutate(
    model = map(splits, ~ glm(H3.epi.size.mean ~ H1.epi.size.mean, data = .x, family = gaussian(link = "log"), maxit = 150)),
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
  # for mpg and value observations within each level of name (e.g., disp, hp, ...)
  mutate(plot_data = map(model_data, ~ {
    # calculate information about the observed mpg and value observations
    # within each level of name to be used in each bootstrap sample
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
        glm(H3.epi.size.mean ~ value, ., family = gaussian(link = "log"), maxit = 150) %>%
        # suppress augment() warnings about dropping columns
        {
          suppressWarnings(augment(., newdata = tibble(value = pred_x)))
        }
    }) %>%
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>%
      # to estimate the lower and upper 95% quantiles of predicted mpgs
      summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

h3_vs_h1_epi_size_post_pdm <- ggplot() +
  geom_ribbon(aes(x = H3.epi.size.mean, ymin = gaussian(link = "log")$linkinv(l), ymax = gaussian(link = "log")$linkinv(u)),
    tbl_plot_data %>% filter(name == "H1.epi.size.mean") %>% rename(H3.epi.size.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = H1.epi.size.mean,
      ymin = H3.epi.size.lowCI, ymax = H3.epi.size.hiCI,
      color = dom_type
    ), width = .05
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.epi.size.mean,
      xmin = H1.epi.size.lowCI, xmax = H1.epi.size.hiCI,
      color = dom_type
    ), width = .05
  ) +
  geom_point(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.epi.size.mean, fill = dom_type), size = 5, pch = 21
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
    breaks = c("H3", "H1", "co-circ"),
    labels = c("H3", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "co-circ"),
    labels = c("H3", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  # scale_color_manual(labels=c("H3/H1","H1","H3"),breaks=c("co-circ","H1","H3"),values=cols)+
  # scale_fill_manual(labels=c("H3/H1","H1","H3"),breaks=c("co-circ","H1","H3"),values=cols)+
  geom_smooth(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.epi.size.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = gaussian(link = "log")),
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 57, x = 1, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_h1_epi_size_post_pdm

#####################################################
## R0: Gamma, log link (very little difference across models, "inverse" link produces wonky results)
#####################################################
df <- sumdf2 %>%
  dplyr::select(season, H3.R0.mean, H1.epi.size.mean) %>%
  filter(season != "2000-2001")

y <- df$H3.R0.mean
x <- df$H1.epi.size.mean
linear.model <- glm(y ~ x, family = gaussian())
log.model <- glm(y ~ x, family = gaussian(link = "log"))
inv.model <- glm(y ~ x, family = gaussian(link = "inverse"))
gamma.model <- glm(y ~ x, family = Gamma(link = "log"))
gamma.model2 <- glm(y ~ x, family = Gamma(link = "inverse"))
model.sel(linear.model, log.model, inv.model, gamma.model, gamma.model2)

set.seed(27)
boots <- bootstraps(sumdf2, times = 1000, apparent = TRUE)
boots

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

m1 <- glm(H3.R0.mean ~ H1.epi.size.mean, data = sumdf2, family = Gamma(link = "log"), maxit = 150)
fam <- family(m1)
ilink <- fam$linkinv

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
    pvalue = sprintf("italic(p) == %.1f", pvalue)
  )
labels

set.seed(27)

n_boot <- 1000

sumdf2 %>%
  dplyr::select(H1.epi.size.mean, IVB.epi.size.mean, H3.R0.mean) %>%
  pivot_longer(cols = c(H1.epi.size.mean, IVB.epi.size.mean)) -> tbl_mtcars_long

tbl_mtcars_long %>%
  nest(model_data = c(H3.R0.mean, value)) %>%
  # for mpg and value observations within each level of name (e.g., disp, hp, ...)
  mutate(plot_data = map(model_data, ~ {
    # calculate information about the observed mpg and value observations
    # within each level of name to be used in each bootstrap sample
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
      # to estimate the lower and upper 95% quantiles of predicted mpgs
      summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

h3_vs_h1_r0_post_pdm <- ggplot() +
  geom_ribbon(
    aes(
      x = H3.R0.mean,
      ymin = Gamma(link = "log")$linkinv(l), ymax = Gamma(link = "log")$linkinv(u)
    ),
    tbl_plot_data %>% filter(name == "H1.epi.size.mean") %>% rename(H3.R0.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = H1.epi.size.mean,
      ymin = H3.R0.lowCI, ymax = H3.R0.hiCI,
      color = dom_type
    ), width = .05
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.R0.mean,
      xmin = H1.epi.size.lowCI, xmax = H1.epi.size.hiCI,
      color = dom_type
    ), width = .0001
  ) +
  geom_point(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.R0.mean, fill = dom_type), size = 5, pch = 21
  ) +
  xlab("A(H1N1) Epidemic Size") +
  ylab("A(H3N2) Effective Rt") +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(
    breaks = c("H3", "H1", "co-circ"),
    labels = c("H3", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "co-circ"),
    labels = c("H3", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  # scale_color_manual(labels=c("H3/H1","H1","H3"),breaks=c("co-circ","H1","H3"),values=cols)+
  # scale_fill_manual(labels=c("H3/H1","H1","H3"),breaks=c("co-circ","H1","H3"),values=cols)+
  geom_smooth(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.R0.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = Gamma(link = "log")),
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 1.5, x = 0.5, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_h1_r0_post_pdm

#####################################################
## Epidemic Intensity, Beta distribution
#####################################################

sumdf2 <- sumdf2 %>% filter(!(season %in% c("2000-2001", "2009-2010")))
y <- sumdf2$H3.shannon.mean
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
  quasibinomial.model, beta.model
)

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
  # for mpg and value observations within each level of name (e.g., disp, hp, ...)
  mutate(plot_data = map(model_data, ~ {
    # calculate information about the observed mpg and value observations
    # within each level of name to be used in each bootstrap sample
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
      # to estimate the lower and upper 95% quantiles of predicted mpgs
      summarize(
        l = quantile(.fitted, .025),
        u = quantile(.fitted, .975),
        .groups = "drop"
      ) %>%
      arrange(value)
  })) %>%
  dplyr::select(-model_data) %>%
  unnest(plot_data) -> tbl_plot_data

h3_vs_h1_shannon_post_pdm <- ggplot() +
  geom_ribbon(aes(x = H3.shannon.mean, ymin = l, ymax = u),
    tbl_plot_data %>% filter(name == "H1.epi.size.mean") %>% rename(H3.shannon.mean = value),
    alpha = 0.3, fill = "grey"
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      x = H1.epi.size.mean,
      ymin = H3.shannon.lowCI, ymax = H3.shannon.hiCI,
      color = dom_type
    ), width = .003
  ) +
  geom_errorbar(
    data = sumdf2,
    aes(
      y = H3.shannon.mean,
      xmin = H1.epi.size.lowCI, xmax = H1.epi.size.hiCI,
      color = dom_type
    ), width = .003
  ) +
  geom_point(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.shannon.mean, fill = dom_type), size = 5, pch = 21
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
    breaks = c("H3", "H1", "co-circ"),
    labels = c("H3", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  scale_fill_manual(
    breaks = c("H3", "H1", "co-circ"),
    labels = c("H3", "H1pdm", "H3/H1pdm"), values = cols
  ) +
  geom_smooth(
    data = sumdf2,
    aes(x = H1.epi.size.mean, y = H3.shannon.mean),
    method = "glm", formula = y ~ x,
    method.args = list(family = Gamma(link = "log")),
    se = F,
    size = 1, linetype = 2, color = "black"
  ) +
  geom_text(
    y = 0.43, x = 1, aes(label = paste(adj.r.squared, pvalue, sep = "~`,`~")), size = 5,
    data = labels, parse = TRUE, hjust = 0
  )
h3_vs_h1_shannon_post_pdm

combined_post_pdm <- plot_grid(h3_vs_h1_epi_size_post_pdm + theme(legend.position = "none"),
  h3_vs_h1_max_peak_post_pdm + theme(legend.position = "none"),
  h3_vs_h1_r0_post_pdm + theme(legend.position = "none"),
  h3_vs_h1_shannon_post_pdm + theme(legend.position = "none"),
  labels = NULL, nrow = 1
)
combined_post_pdm

combined_all_seasons

title <- ggdraw() +
  draw_label(
    "All seasons",
    fontface = "bold",
    x = 0,
    hjust = 0,
    size = 18
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 6)
  )

combined_all_seasons_title <-
  plot_grid(
    title, combined_all_seasons,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )


combined_pre_pdm

title <- ggdraw() +
  draw_label(
    "Pre-2009 seasons",
    fontface = "bold",
    x = 0,
    hjust = 0,
    size = 18
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 6)
  )

combined_pre_pdm_title <-
  plot_grid(
    title, combined_pre_pdm,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )

title <- ggdraw() +
  draw_label(
    "Post-2009 seasons",
    fontface = "bold",
    x = 0,
    hjust = 0,
    size = 18
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 6)
  )

combined_post_pdm_title <-
  plot_grid(
    title, combined_post_pdm,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )

all_plots <- plot_grid(combined_all_seasons_title, combined_pre_pdm_title, combined_post_pdm_title, nrow = 3, labels = NULL)

all_plots2 <- plot_grid(all_plots, leg, nrow = 2, rel_heights = c(3, 0.1))
all_plots2
save_plot(all_plots2,
  filename = "figures/Fig7_sup_fig2_h1_vs_h3_parameters_pre_post_pdm.png",
  base_width = 18, base_height = 15, dpi = 300, bg = 'white'
)
save_plot(all_plots2,
  filename = "figures/Fig7_sup_fig2_h1_vs_h3_parameters_pre_post_pdm.pdf",
  base_width = 18, base_height = 15, dpi = 300, bg = 'white'
)
