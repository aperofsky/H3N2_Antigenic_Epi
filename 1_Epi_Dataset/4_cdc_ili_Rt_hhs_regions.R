# install RStan before installing epidemia https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
# After installing RStan, run:
# devtools::install_github("ImperialCollegeLondon/epidemia")

## load packages
list.of.packages <- c("cdcfluview", "dplyr", "ggplot2", "epidemia", "tempdisagg", "R0", "lognorm", "EpiNow2")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

## load data
load("data/hhs_division_level_ILI_and_virology_interp_smoothed_sampling_effort.RData")
load("data/region_level_flu_metrics_H3.RData")

range(regionflu_ili_vir_adj$wk_date) # ""1997-09-28" "2021-09-26"

regions <- unique(regionflu_ili_vir_adj$region)
seasons <- unique(regionflu_ili_vir_adj$season_description)

regions
seasons
## warning: this takes a long time to run on a laptop. this is better suited to run on a cluster. ~20 seasons x 9-10 regions
r0_eg_list <- list()
for (k in seasons) {
  for (j in regions) {
    if (k == "2009-2010" | k == "2000-2001") next
    data1 <- regionflu_ili_vir_adj[regionflu_ili_vir_adj$season_description == k & regionflu_ili_vir_adj$region == j, ]
    colnames(data1)[names(data1) == "season_description"] <- "season"
    data1 <- data1[order(data1$wk_date), ]
    data1$ili_pos2 <- ifelse(is.na(data1$ili_h3_st), 0, data1$ili_h3_st)
    epi_dates <- region_flu_metrics_H3[region_flu_metrics_H3$season == k & region_flu_metrics_H3$region == j, ]
    if (dim(epi_dates)[1] == 0) next
    region <- data1[, c("wk_date", "ili_pos2")] # subset fluhosp to one region
    start <- epi_dates$onset_week ## epidemic onset
    peak <- epi_dates$peak_week ## peak day
    region2 <- region[region$wk_date >= start - 21 & region$wk_date <= peak + 21, ]

    if (all(is.na(region2) == TRUE)) next
    if (dim(region2)[1] < 3) next
    if (dim(region2[region2$ili_pos2 > 0, ])[1] < 3) next
    if (max(with(rle(region2$ili_pos2 > 0), lengths[values == TRUE])) < 3) next
    if (sum(region2$ili_pos2) < 0.09) next
    ind2 <- region2[, 2]
    ind2 <- ceiling(ind2 * 1000)

    data <- data.frame(date = region2[, 1]$wk_date, cases = ind2$ili_pos2)

    m.d.noind <- td(data ~ 1, to = "daily", method = "fast", conversion = "sum")
    gdp.d.noind <- predict(m.d.noind)

    # ts_plot(
    #   ts_scale(
    #     ts_c(gdp.d.noind, data)
    #   ),
    #   title = "Daily disaggregated ILI",
    #   subtitle = "no indicator"
    # )

    gdp.d.noind$value <- gdp.d.noind$value + abs(min(gdp.d.noind$value))

    # vector of dates for incidence data
    date <- (first(gdp.d.noind$time) - 1) + seq(0, along.with = c(NA, gdp.d.noind$value))

    data <- data.frame(region = j, cases = c(NA, round(gdp.d.noind$value, digits = 0)), date = date)

    ## influenza incubation period
    mu1 <- EpiNow2::convert_to_logmean(1.4, 1.51)
    sigma1 <- EpiNow2::convert_to_logsd(1.4, 1.51)

    ## influenza symptom onset to report
    mu2 <- EpiNow2::convert_to_logmean(2, 1.4858)
    sigma2 <- EpiNow2::convert_to_logsd(2, 1.4858)

    ## infection to symptom onset + symptom onset to care seeking
    combined_dist <- estimateSumLognormal(c(mu1, mu2), c(sigma1, sigma2))

    ## code borrowed from R0 package
    tmax <- 13
    t.scale <- c(0, 0.5 + c(0:tmax))
    meanlog <- combined_dist[[1]]
    sdlog <- combined_dist[[2]]
    inf_to_case <- diff(plnorm(t.scale, meanlog = meanlog, sdlog = sdlog))
    inf_to_case_adj <- inf_to_case / sum(inf_to_case) # make sure sums to 1

    ## influenza generation time
    GT.flu <- R0::generation.time("weibull", c(3.6, 1.6))

    # transmission
    ## Rt in a given group at the given date
    ## simple random walk model
    rt <- epirt(
      formula = R(region, date) ~ 1 + rw(time = date, prior_scale = 0.01),
      prior_intercept = normal(log(1.3), 0.2), link = "log"
    )

    ## percent who seek care for ILI = 45% source: doi: 10.1093/infdis/jiu224
    # observations are a function of latent infections in the population
    data$off <- 0.45
    obs <- epiobs(
      formula = cases ~ 0 + offset(off),
      family = "neg_binom",
      link = "logit",
      prior_aux = normal(location = 10, scale = 5),
      i2o = inf_to_case_adj
    )

    inf <- epiinf(gen = GT.flu$GT)

    # treat infections as latent parameters
    inf_ext <- epiinf(gen = GT.flu$GT, latent = TRUE)

    args <- list(rt = rt, obs = obs, inf = inf, data = data, iter = 1e4, warmup = 2e3, seed = 12345)
    args_ext <- args
    args_ext$inf <- inf_ext

    options(mc.cores = parallel::detectCores())
    fm2 <- do.call(epim, args_ext)

    post_rt <- posterior_rt(fm2)

    post_rt_median <- data.frame(
      date = post_rt$time,
      median = apply(post_rt$draws, 2, function(x) quantile(x, 0.5)),
      group = post_rt$group
    )

    post_rt_median$year <- mmwr_week(post_rt_median$date)[, 1]$mmwr_year
    post_rt_median$week <- mmwr_week(post_rt_median$date)[, 2]$mmwr_week
    post_rt_median$cdc_date <- mmwr_week_to_date(post_rt_median$year, post_rt_median$week)

    initial_rt <- post_rt_median %>%
      group_by(group, cdc_date) %>%
      summarize(rt_weekly = mean(median)) %>%
      ungroup() %>%
      filter(cdc_date >= start & cdc_date <= peak) %>%
      summarize(rt_weekly = mean(rt_weekly)) %>%
      pull(rt_weekly)

    max_rt <- post_rt_median %>%
      group_by(group, cdc_date) %>%
      summarize(rt_weekly = mean(median)) %>%
      ungroup() %>%
      slice_max(n = 2, order_by = rt_weekly) %>%
      summarize(rt_weekly = mean(rt_weekly)) %>%
      pull(rt_weekly)


    R0_values <- data.frame(region = j, season = k, initial_Rt = initial_rt, max_Rt = max_rt)
    print(R0_values)
    r0_eg_list[[length(r0_eg_list) + 1]] <- R0_values
  }
}
R0_eg_values <- do.call(rbind.data.frame, r0_eg_list)
head(R0_eg_values)
table(R0_eg_values$region)
table(R0_eg_values$season)
R0_eg_values %>% arrange(-max_Rt)
save(R0_eg_values, file = "data/epidemic_Rt_estimates.RData")
