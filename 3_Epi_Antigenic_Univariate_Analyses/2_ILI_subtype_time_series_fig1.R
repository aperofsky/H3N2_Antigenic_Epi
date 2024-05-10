## load packages
# # Hmisc package is needed for stat_summary function to work inside ggplot but shouldn't be loaded
list.of.packages <- c("dplyr", "ggplot2", "cowplot", "ggExtra", "tidyr", "ggsci", "RColorBrewer", "scales")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

theme_set(theme_cowplot(font_size = 16))

####################################################################################
## Plot ILI x percent positive for types/subtypes
####################################################################################
## load data
load("data/hhs_division_level_ILI_and_virology_interp_smoothed_sampling_effort.RData")
head(regionflu_ili_vir_adj)
names(regionflu_ili_vir_adj)

df <- regionflu_ili_vir_adj %>%
  filter(!(season_description %in% c("1996-1997", "2019-2020", "2020-2021", "2021-2022"))) %>%
  mutate(
    ili_h3_st = ifelse(region == "Region 10" & wk_date < as.Date("2009-10-01"), NA, ili_h3_st),
    ili_ivb_st = ifelse(region == "Region 10" & wk_date < as.Date("2009-10-01"), NA, ili_ivb_st),
    ili_ivb_st = ifelse(region == "Region 10" & wk_date < as.Date("2009-10-01"), NA, ili_ivb_st)
  )
df$region <- as.factor(df$region)
df$region <- factor(df$region, levels = c(
  "Region 1", "Region 2", "Region 3", "Region 4", "Region 5",
  "Region 6", "Region 7", "Region 8", "Region 9", "Region 10"
))
range(df$wk_date) # "1997-09-28" "2019-06-30"

df <- df %>%
  group_by(region) %>%
  complete(wk_date = seq.Date(range(df$wk_date)[1], range(df$wk_date)[2], by = "week")) %>%
  ungroup()

df %>%
  dplyr::select(wk_date, region, ili_h3_st) %>%
  filter(is.na(ili_h3_st)) %>%
  group_by(region) %>%
  tally()

df <- df %>% mutate_at(c("ili_h3_st"), ~ (scale(.) %>% as.vector()))

myPalette <- colorRampPalette(brewer.pal(9, "Reds"), space = "Lab")

ylabels <- rev(c("1-Boston", "2-NYC", "3-DC", "4-Atlanta", "5-Chicago", "6-Dallas", "7-Kansas City", "8-Denver", "9-San Francisco", "10-Seattle"))
labels.wrap <- lapply(strwrap(ylabels, 30, simplify = F), paste, collapse = "") # word wrap

p <- ggplot(df, aes(x = wk_date, y = rev(region), fill = ili_h3_st)) +
  # geom_tile(color = "grey", linewidth = 0.00001) +
  geom_tile()+
  scale_fill_gradientn(
    colors = myPalette(100), name = "ILI x % A(H3N2)",
    na.value = "white"
  ) +
  scale_x_date(
    expand = c(0, 0), date_breaks = "3 years", date_labels = "%Y",
    limits = c(as.Date("1997-06-01"), as.Date("2019-06-01"))
  ) +
  theme_minimal(base_size = 14) +
  labs(x = "Date", y = "HHS Region") +
  scale_y_discrete(labels = rev(seq(1:10)), breaks = c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5", "Region 6", "Region 7", "Region 8", "Region 9", "Region 10"))
p

p2 <- p + theme(legend.position = "bottom") +
  theme(strip.background = element_rect(colour = "white")) +
  theme(axis.ticks = element_blank()) +
  theme(legend.title = element_text(size = 10)) +
  theme(legend.text = element_text(size = 8), legend.title.align = 0.5) +
  theme(axis.text.y = element_text(hjust = 0)) +
  removeGrid() # ggExtra
p2

long_regionflu <- regionflu_ili_vir_adj %>%
  filter(!(season_description %in%
    c("1996-1997", "2019-2020", "2020-2021", "2021-2022"))) %>%
  dplyr::select(wk_date, season_description, region, ili_h3_st:ili_ivb_st) %>%
  gather(type, incidence, ili_h3_st:ili_ivb_st)

years <- data.frame(y = seq(as.Date("1998-01-01"), as.Date("2020-01-01"), "year"))

long_regionflu[!complete.cases(long_regionflu), ] %>% filter(region == "Region 5")

## note: need to install Hmisc for stat_summary to work
### ILI x % virus type average across regions
incidence <- ggplot(
  long_regionflu,
  aes(x = wk_date, y = incidence, group = factor(type), color = factor(type), fill = factor(type))
) +
  stat_summary(geom = "ribbon", fun.data = "mean_cl_boot", alpha = 0.6) +
  geom_vline(xintercept = as.numeric(years$y), linetype = "dashed", col = "grey", alpha = 0.9) +
  scale_x_date(
    date_breaks = "3 years", date_labels = "%Y",
    limits = c(as.Date("1997-06-01"), as.Date("2019-06-01")), expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Scaled ILI x % positive") +
  theme_minimal(base_size = 14) +
  removeGrid() # ggExtra
incidence

cols <- c(
  "#40abf0",
  "#dd3b3f",
  "#00a100"
)

incidence2 <- incidence +
  scale_fill_manual(guide = guide_legend(), values = cols, labels = c("ILI x % A(H1N1)", "ILI x % A(H3N2)", "ILI x % B")) +
  scale_color_manual(guide = guide_legend(), values = cols, labels = c("ILI x % A(H1N1)", "ILI x % A(H3N2)", "ILI x % B")) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(),
    legend.justification = "center"
  )
incidence2

#####################################################
## Plot seasonal antigenic and genetic distances
#####################################################

all_measures_HA <- readr::read_tsv("2_Phylo_Dataset/distance_tables/north-america_summarized_mean_seasonal_distances_h3n2_ha_21y.tsv")

all_measures_HA2 <- all_measures_HA %>%
  mutate(
    time_diff = as.numeric(format(as.Date(other_season_start), "%Y")) - as.numeric(format(as.Date(current_season_start), "%Y")),
    season1 = paste(as.numeric(format(as.Date(current_season_start), "%Y")), as.numeric(format(as.Date(current_season_end), "%Y")), sep = "-"),
    season2 = paste(as.numeric(format(as.Date(other_season_start), "%Y")), as.numeric(format(as.Date(other_season_end), "%Y")), sep = "-")
  )

HA_one_season_lag <- all_measures_HA2 %>%
  filter(time_diff == 1) %>%
  dplyr::rename(season = season2) %>%
  dplyr::select(season, count, mean, distance_map, std)

all_measures_NA <- readr::read_tsv("2_Phylo_Dataset/distance_tables/north-america_summarized_mean_seasonal_distances_h3n2_na_21y.tsv")
all_measures_NA2 <- all_measures_NA %>%
  mutate(
    time_diff = as.numeric(format(as.Date(other_season_start), "%Y")) - as.numeric(format(as.Date(current_season_start), "%Y")),
    season1 = paste(as.numeric(format(as.Date(current_season_start), "%Y")), as.numeric(format(as.Date(current_season_end), "%Y")), sep = "-"),
    season2 = paste(as.numeric(format(as.Date(other_season_start), "%Y")), as.numeric(format(as.Date(other_season_end), "%Y")), sep = "-")
  )

NA_one_season_lag <- all_measures_NA2 %>%
  filter(time_diff == 1) %>%
  dplyr::rename(season = season2) %>%
  dplyr::select(season, count, mean, distance_map, std)


load("data/north_amer_build_season_h3n2_replicates_HA_direct_ag_distances.RData")
names(seas.ag.HA)

# HA_lbi = seas.ag.HA %>%
#   dplyr::select(season,replicate,HA_mean_lbi,HA_std_lbi)%>%
#   group_by(season)%>% 
#   summarize(mean = mean(HA_mean_lbi),
#             std = mean(HA_mean_lbi))%>%
#   mutate(distance_map = "HA mean LBI",
#          count = 5)

HA_lbi_std = seas.ag.HA %>%
  dplyr::select(season,replicate,HA_std_lbi)%>%
  group_by(season)%>% 
  summarize(mean = mean(HA_std_lbi),
            std = mean(HA_std_lbi))%>%
  mutate(distance_map = "HA s.d. LBI",
         count = 5)

load("data/north_amer_build_season_h3n2_replicates_NA_direct_ag_distances.RData")
names(seas.ag.NA)

# NA_lbi = seas.ag.NA %>%
#   dplyr::select(season,replicate,NA_mean_lbi,NA_std_lbi)%>%
#   group_by(season)%>% 
#   summarize(mean = mean(NA_mean_lbi),
#             std = mean(NA_mean_lbi))%>%
#   mutate(distance_map = "NA mean LBI",
#          count = 5)

NA_lbi_std = seas.ag.NA %>%
  dplyr::select(season,replicate,NA_std_lbi)%>%
  group_by(season)%>% 
  summarize(mean = mean(NA_std_lbi),
            std = mean(NA_std_lbi))%>%
  mutate(distance_map = "NA s.d. LBI",
         count = 5)

# ha_lbi <- readr::read_tsv("2_Phylo_Dataset/distance_tables/north-america_mean_seasonal_lbi_h3n2_ha_21y.tsv")
# head(ha_lbi)
# ha_lbi <- ha_lbi %>%
#   filter(replicate == 0) %>%
#   mutate(
#     season1 = as.numeric(format(as.Date(season_start), "%Y")),
#     season2 = as.numeric(format(as.Date(season_end), "%Y"))
#   ) %>%
#   mutate(season = paste(season1, season2, sep = "-")) %>%
#   mutate(distance_map = "HA LBI") %>%
#   dplyr::select(season, count, mean, distance_map, std)
# 
# 
# na_lbi <- readr::read_tsv("2_Phylo_Dataset/distance_tables/north-america_mean_seasonal_lbi_h3n2_na_21y.txt")
# na_lbi <- na_lbi %>%
#   filter(replicate == 0) %>%
#   mutate(
#     season1 = as.numeric(format(as.Date(season_start), "%Y")),
#     season2 = as.numeric(format(as.Date(season_end), "%Y"))
#   ) %>%
#   mutate(season = paste(season1, season2, sep = "-")) %>%
#   mutate(distance_map = "NA LBI") %>%
#   dplyr::select(season, count, mean, distance_map, std)

lbi_div <- readr::read_rds("data/LBI_diversity_by_season.rds") %>%
  dplyr::select(season, replicate, contains("lbi_shannon")) %>%
  group_by(season) %>%
  summarize_at(c("ha_lbi_shannon", "na_lbi_shannon"), ~ mean(.x))

lbi_div_long <- lbi_div %>% pivot_longer(cols = c(ha_lbi_shannon, na_lbi_shannon), names_to = "distance_map", values_to = "mean")
lbi_div_long$count <- 5

lbi_div_std <- readr::read_rds("data/LBI_diversity_by_season.rds") %>%
  dplyr::select(season, replicate, contains("lbi_shannon")) %>%
  group_by(season) %>%
  summarize_at(c("ha_lbi_shannon", "na_lbi_shannon"), ~ sd(.x))

lbi_div_std_long <- lbi_div_std %>%
  pivot_longer(cols = c(ha_lbi_shannon, na_lbi_shannon), names_to = "distance_map", values_to = "std")

lbi_div_all <- full_join(lbi_div_long, lbi_div_std_long, by = c("season", "distance_map"))

all_measures <-
  bind_rows(
    HA_one_season_lag, NA_one_season_lag,
    # ha_lbi, na_lbi, 
    # HA_lbi,NA_lbi,
    HA_lbi_std,NA_lbi_std,
    lbi_div_all
  ) %>%
  filter(!(season %in% c("1995-1996", "1996-1997", "2019-2020"))) %>%
  mutate(margin_of_error = qt(0.975, df = count - 1) * std / sqrt(count)) %>%
  mutate(
    lwr = mean - margin_of_error,
    upr = mean + margin_of_error
  )

all_measures <- all_measures %>% tidyr::separate(season, into = c("year1", "year2"), remove = F)
all_measures$year.new <- as.Date(paste(all_measures$year2, "01", "01", sep = "-"))
unique(all_measures$distance_map)

all_measures$distance_map <- factor(all_measures$distance_map, levels = c(
  "wolf", 
  "wolf_nonepitope",
  "koel", 
  "stem", "stem_ep",
  "titer_substitution_model", "titer_tree_model",
  "bhatt", "bhatt_nonepitope",
  "krammer", 
  # "HA mean LBI","NA mean LBI",
  "HA s.d. LBI","NA s.d. LBI",
  "ha_lbi_shannon", "na_lbi_shannon"
))
unique(all_measures$distance_map)

ep_plot <- ggplot(data = all_measures %>%
  filter(!(distance_map %in% c(
    "titer_substitution_model", "titer_tree_model",
    "bhatt_nonepitope", "wolf_nonepitope", "stem",
    "ha_lbi_shannon", "na_lbi_shannon",
    # "HA mean LBI","NA mean LBI",
    "HA s.d. LBI","NA s.d. LBI",
    "ha_lbi_shannon", "na_lbi_shannon"
  )))) +
  geom_point(
    pch = 21,
    aes(x = year.new, y = mean, color = distance_map),
    alpha = 0.8, size = 4
  ) +
  geom_line(aes(x = year.new, y = mean, color = distance_map),
    alpha = 0.8, linewidth = 2
  ) +
  geom_errorbar(aes(
    x = year.new,
    ymin = lwr, ymax = upr,
    color = distance_map
  ), width = 0.2) +
  xlab("Season") +
  ylab("Mean epitope distance\nfrom the prior season") +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.66, 0.9),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.direction = "horizontal",
    legend.key = element_rect(colour = "transparent", fill = "transparent")
  ) +
  scale_fill_locuszoom(
    name = "Indicator",
    labels = c(
      "HA epitope (N= 129)",
      "HA RBS (N = 7)",
      "HA stalk footprint (N = 34)",
      "NA epitope (N = 223)",
      "NA epitope (N = 53)"
    )
  ) +
  scale_color_locuszoom(
    name = "Indicator",
    labels = c(
      "HA epitope (N= 129)",
      "HA RBS (N = 7)",
      "HA stalk footprint (N = 34)",
      "NA epitope (N = 223)",
      "NA epitope (N = 53)"
    )
  ) +
  scale_x_date(expand = c(0.01, 0.01), date_breaks = "3 years", date_labels = "%Y") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(fill = NA)))
ep_plot

ag_plot <- ggplot(data = all_measures %>%
  filter(distance_map %in% c("titer_tree_model"))) +
  geom_point(
    pch = 21,
    aes(x = year.new, y = mean, color = distance_map),
    alpha = 0.8, size = 4
  ) +
  geom_line(aes(x = year.new, y = mean, color = distance_map),
    alpha = 0.8, linewidth = 2
  ) +
  geom_errorbar(aes(
    x = year.new,
    ymin = lwr, ymax = upr,
    color = distance_map
  ), width = 0.2) +
  xlab("Season") +
  ylab(expression(atop(paste("Mean HI log"[2], "titer distance"), paste("from the prior season")))) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    legend.title = element_blank(),
    legend.direction = "horizontal"
  ) +
  scale_color_locuszoom(
    name = "Indicator",
    labels = c("HI titer (substitution model)", "HI titer (tree model)")
  ) +
  scale_x_date(expand = c(0.01, 0.01), date_breaks = "3 years", date_labels = "%Y")
ag_plot

lbi_div_plot <- ggplot(data = all_measures %>%
  filter(distance_map %in% c("ha_lbi_shannon", "na_lbi_shannon"))) +
  geom_point(
    pch = 21,
    aes(x = year.new, y = mean, color = distance_map),
    alpha = 0.8, size = 4
  ) +
  geom_line(aes(x = year.new, y = mean, color = distance_map),
    alpha = 0.8, linewidth = 2
  ) +
  geom_errorbar(aes(
    x = year.new,
    ymin = lwr, ymax = upr,
    color = distance_map
  ), width = 0.2) +
  xlab("Season") +
  ylab("LBI Diversity") +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.8, 0.1), legend.title = element_blank(),
    legend.direction = "horizontal"
  ) +
  scale_fill_npg(
    name = "Indicator",
    labels = c("HA LBI", "NA LBI")
  ) +
  scale_color_npg(
    name = "Indicator",
    labels = c("HA LBI", "NA LBI")
  ) +
  scale_x_date(expand = c(0.01, 0.01), date_breaks = "3 years", date_labels = "%Y") +
  scale_y_continuous(expand=c(0.01,0.01),limits = c(1,9),n.breaks = 10)+
  guides(color = guide_legend(override.aes = list(fill = NA)))
lbi_div_plot

lbi_mean_plot <- ggplot(data = all_measures %>%
                         filter(distance_map %in% c("HA mean LBI","NA mean LBI"))) +
  geom_point(
    pch = 21,
    aes(x = year.new, y = mean, color = distance_map),
    alpha = 0.8, size = 4
  ) +
  geom_line(aes(x = year.new, y = mean, color = distance_map),
            alpha = 0.8, linewidth = 2
  ) +
  # geom_errorbar(aes(
  #   x = year.new,
  #   ymin = lwr, ymax = upr,
  #   color = distance_map
  # ), width = 0.2) +
  xlab("Season") +
  ylab("Mean LBI") +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.2, 0.9), legend.title = element_blank(),
    legend.direction = "horizontal"
  ) +
  scale_fill_npg(
    name = "Indicator",
    labels = c("HA LBI", "NA LBI")
  ) +
  scale_color_npg(
    name = "Indicator",
    labels = c("HA LBI", "NA LBI")
  ) +
  scale_x_date(expand = c(0.01, 0.01), date_breaks = "3 years", date_labels = "%Y") +
  guides(color = guide_legend(override.aes = list(fill = NA)))
lbi_mean_plot

# log_lbi_std_plot <- ggplot(data = all_measures %>%
#                               filter(distance_map %in% c("HA log(LBI) diversity","NA log(LBI) diversity"))) +
#   geom_point(
#     pch = 21,
#     aes(x = year.new, y = mean, color = distance_map),
#     alpha = 0.8, size = 4
#   ) +
#   geom_line(aes(x = year.new, y = mean, color = distance_map),
#             alpha = 0.8, linewidth = 2
#   ) +
#   geom_errorbar(aes(
#     x = year.new,
#     ymin = lwr, ymax = upr,
#     color = distance_map
#   ), width = 0.2) +
#   xlab("Season") +
#   ylab("S.D. log LBI") +
#   theme_bw(base_size = 16) +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     legend.position = c(0.2, 0.9), legend.title = element_blank(),
#     legend.direction = "horizontal"
#   ) +
#   scale_fill_npg(
#     name = "Indicator",
#     labels = c("HA LBI", "NA LBI")
#   ) +
#   scale_color_npg(
#     name = "Indicator",
#     labels = c("HA LBI", "NA LBI")
#   ) +
#   scale_x_date(expand = c(0.01, 0.01), date_breaks = "3 years", date_labels = "%Y") +
#   guides(color = guide_legend(override.aes = list(fill = NA)))
# log_lbi_std_plot

inc_combined <- plot_grid(incidence2, p2, nrow = 2, rel_heights = c(1, 1), align = "hv", labels = "AUTO")
inc_combined
save_plot(inc_combined, filename = "figures/ILI_time_series.png", base_width = 10, base_height = 10)
save_plot(inc_combined, filename = "figures/ILI_time_series.pdf", dpi = 300, base_width = 10, base_height = 10)

ag_combined <- plot_grid(ep_plot, ag_plot, lbi_div_plot,
  nrow = 3,
  rel_heights = c(1, 1, 1), align = "hv", labels = c("C", "D", "E")
)
ag_combined
save_plot(ag_combined, filename = "figures/all_evol_time_series.png", base_width = 12, base_height = 16)
save_plot(ag_combined, filename = "figures/all_evol_time_series.pdf", dpi = 300, base_width = 12, base_height = 16)

#####################################################
## H1N1 Incidence
#####################################################

ylabels <- rev(c("1-Boston", "2-NYC", "3-DC", "4-Atlanta", "5-Chicago", "6-Dallas", "7-Kansas City", "8-Denver", "9-San Francisco", "10-Seattle"))
labels.wrap <- lapply(strwrap(ylabels, 30, simplify = F), paste, collapse = "") # word wrap

h1_p <- ggplot(df, aes(x = wk_date, y = rev(region), fill = ili_h1_st)) +
  geom_tile()+
  # geom_tile(color = "white", linewidth = 0.00001) +
  scale_fill_gradientn(
    colors = myPalette(100), name = "ILI x % A(H1N1)",
    na.value = "white",
    rescaler = function(x, to = c(0, 1), from = NULL) {
      ifelse(x < 12,
        scales::rescale(x,
          to = to,
          from = c(min(x, na.rm = TRUE), 12)
        ),
        1
      )
    }
  ) +
  scale_x_date(
    expand = c(0, 0), date_breaks = "3 years", date_labels = "%Y",
    limits = c(as.Date("1997-06-01"), as.Date("2019-06-01"))
  ) +
  theme_minimal(base_size = 14) +
  labs(x = "Date", y = "HHS Region") +
  scale_y_discrete(labels = labels.wrap, breaks = c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5", "Region 6", "Region 7", "Region 8", "Region 9", "Region 10"))
h1_p

h1_p2 <- h1_p + theme(legend.position = "bottom") +
  theme(strip.background = element_rect(colour = "white")) +
  theme(axis.ticks = element_blank()) +
  theme(legend.title = element_text(size = 10)) +
  theme(legend.text = element_text(size = 8), legend.title.align = 0.5) +
  theme(axis.text.y = element_text(hjust = 0)) +
  removeGrid() # ggExtra
h1_p2

ivb_p <- ggplot(df, aes(x = wk_date, y = rev(region), fill = ili_ivb_st)) +
  geom_tile()+
  # geom_tile(color = "white", linewidth = 0.00001) +
  scale_fill_gradientn(
    colors = myPalette(100), name = "ILI x % B",
    na.value = "white"
  #   rescaler = function(x, to = c(0, 1), from = NULL) {
  #     ifelse(x < 12,
  #       scales::rescale(x,
  #         to = to,
  #         from = c(min(x, na.rm = TRUE), 12)
  #       ),
  #       1
  #     )
  #   }
  ) +
  scale_x_date(
    expand = c(0, 0), date_breaks = "3 years", date_labels = "%Y",
    limits = c(as.Date("1997-06-01"), as.Date("2019-06-01"))
  ) +
  theme_minimal(base_size = 14) +
  labs(x = "Date", y = "HHS Region") +
  scale_y_discrete(labels = labels.wrap, breaks = c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5", "Region 6", "Region 7", "Region 8", "Region 9", "Region 10"))
ivb_p

ivb_p2 <- ivb_p + theme(legend.position = "bottom") +
  theme(strip.background = element_rect(colour = "white")) +
  theme(axis.ticks = element_blank()) +
  theme(legend.title = element_text(size = 10)) +
  theme(legend.text = element_text(size = 8), legend.title.align = 0.5) +
  theme(axis.text.y = element_text(hjust = 0)) +
  removeGrid() # ggExtra
ivb_p2

combined <- plot_grid(h1_p2, ivb_p2, labels = "AUTO", nrow = 2)
combined
save_plot(combined, filename = "figures/ILI_H1_and_IVB_timing.png", base_width = 12, base_height = 12)
save_plot(combined, filename = "figures/ILI_H1_and_IVB_timing.pdf", dpi=300,base_width = 12, base_height = 12)
