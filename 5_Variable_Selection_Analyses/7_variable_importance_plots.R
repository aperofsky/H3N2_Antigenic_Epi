list.of.packages <- c("dplyr", "ggplot2", "cowplot", "viridis", "grid", "gridExtra", "ggpubr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

######################################################################################
## Plot random forest and LASSO variable rankings
######################################################################################

######################################################################################
## Random forest variable importance
######################################################################################
## load variable importance results
h3_peak <- read.csv("data/H3_peak_cforest_variable_importance.csv")
h3_peak$epi_measure <- "peak"

h3_size <- read.csv("data/H3_epi_size_cforest_variable_importance.csv")
h3_size$epi_measure <- "cum_size"

h3_r0 <- read.csv("data/H3_R0_cforest_variable_importance.csv")
h3_r0$epi_measure <- "R0"

h3_shannon <- read.csv("data/H3_shannon_entropy_cforest_variable_importance.csv")
h3_shannon$epi_measure <- "shannon_entropy"

h3_dom <- read.csv("data/h3_dom_cforest_variable_importance.csv")
h3_dom$epi_measure <- "dominance"

all_cforest <- bind_rows(h3_peak, h3_size, h3_dom, h3_r0, h3_shannon)
range(all_cforest$CI_upper)

all_cforest$epi_measure <- as.factor(all_cforest$epi_measure)
levels(all_cforest$epi_measure)
levels(all_cforest$epi_measure) <- c(
  "Epidemic Size", "Subtype Dominance",
  "Peak Incidence", "Effective Rt", "Epidemic Intensity"
)
head(all_cforest)

sort(unique(all_cforest$var))

label_vec <- c(
  "A(H1N1) epidemic size",
  "A(H1N1) epidemic size (t-1)",
  "A(H3N2) epidemic size (t-1)",
  "H3 LBI Diversity",
  "H3 LBI Diversity (t-1)",
  "HI titer (t-2)",
  "H3 epitope (t-2)",
  "B epidemic size",
  "B epidemic size (t-1)",
  "N2 epitope (t-1)",
  "N2 LBI Diversity",
  "N2 LBI Diversity (t-1)",
  "Dominant IAV (t-1)",
  "N2 distance to vaccine",
  "H3 distance to vaccine",
  "Vaccine coverage x VE",
  "Vaccine coverage x VE (t-1)",
  "Vaccine coverage",
  "Vaccine coverage (t-1)",
  "A(H3N2) VE",
  "A(H3N2) VE (t-1)"
)

length(label_vec)

label_df <- data.frame(var = sort(unique(all_cforest$var)), label = label_vec)
label_df
all_cforest <- left_join(all_cforest, label_df, by = "var")

head(all_cforest)
length(unique(all_cforest$label))
cforest_limited <- all_cforest %>%
  group_by(epi_measure) %>%
  arrange(epi_measure, desc(Estimate)) %>%
  slice_max(order_by = Estimate, n = 21) %>%
  ungroup()

######################################################################################
# epidemic size
######################################################################################
cforest_plot_epi_size <- ggplot(
  cforest_limited %>% filter(epi_measure == "Epidemic Size"),
  aes(
    x = reorder(label, Estimate), y = Estimate,
    fill = Estimate
  )
) +
  geom_bar(stat = "identity", position = "dodge", color = "black", lwd = 0.1) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, x = reorder(label, Estimate)), color = "black", width = 0.5, lwd = 0.3) +
  ylab("") +
  xlab("") +
  coord_flip() +
  theme_bw(base_size = 16) +
  scale_fill_viridis(option = "rocket", direction = -1, begin = 0.5) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  ggtitle("Epidemic Size") +
  ylim(0, 0.33)
cforest_plot_epi_size

######################################################################################
## Peak incidence
######################################################################################
cforest_plot_peak <- ggplot(
  cforest_limited %>% filter(epi_measure == "Peak Incidence"),
  aes(x = reorder(label, Estimate), y = Estimate, fill = Estimate)
) +
  geom_bar(stat = "identity", position = "dodge", color = "black", lwd = 0.1) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, x = reorder(label, Estimate)), color = "black", width = 0.5, lwd = 0.3) +
  coord_flip() +
  ylab("") +
  xlab("") +
  theme_bw(base_size = 16) +
  scale_fill_viridis(option = "rocket", direction = -1, begin = 0.5) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  ggtitle("Peak Incidence") +
  ylim(0, 0.33)
cforest_plot_peak

######################################################################################
## subtype dominance
######################################################################################
cforest_plot_dom <- ggplot(
  cforest_limited %>% filter(epi_measure == "Subtype Dominance"),
  aes(x = reorder(label, Estimate), y = Estimate, fill = Estimate)
) +
  geom_bar(stat = "identity", position = "dodge", color = "black", lwd = 0.1) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, x = reorder(label, Estimate)), color = "black", width = 0.5, lwd = 0.3) +
  ylab("") +
  xlab("") +
  coord_flip() +
  theme_bw(base_size = 16) +
  scale_fill_viridis(option = "rocket", direction = -1, begin = 0.5) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  ggtitle("Subtype Dominance") +
  ylim(0, 0.33)
cforest_plot_dom

######################################################################################
## effective Rt
######################################################################################
cforest_plot_r0 <- ggplot(
  cforest_limited %>% filter(epi_measure == "Effective Rt"),
  aes(x = reorder(label, Estimate), y = Estimate, fill = Estimate)
) +
  geom_bar(stat = "identity", position = "dodge", color = "black", lwd = 0.1) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, x = reorder(label, Estimate)), color = "black", width = 0.5, lwd = 0.3) +
  coord_flip() +
  ylab("") +
  xlab("") +
  theme_bw(base_size = 16) +
  scale_fill_viridis(option = "rocket", direction = -1, begin = 0.5) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  ggtitle("Effective Rt") +
  ylim(0, 0.33)
cforest_plot_r0

######################################################################################
## inverse shannon entropy
######################################################################################
cforest_plot_SE <- ggplot(
  cforest_limited %>% filter(epi_measure == "Epidemic Intensity"),
  aes(x = reorder(label, Estimate), y = Estimate, fill = Estimate)
) +
  geom_bar(stat = "identity", position = "dodge", color = "black", lwd = 0.1) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, x = reorder(label, Estimate)), color = "black", width = 0.5, lwd = 0.3) +
  coord_flip() +
  ylab("") +
  xlab("") +
  theme_bw(base_size = 16) +
  scale_fill_viridis(option = "rocket", direction = -1, begin = 0.5) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  ggtitle("Epidemic Intensity") +
  ylim(0, 0.33)
cforest_plot_SE

combined_cforest <- plot_grid(cforest_plot_epi_size,
  cforest_plot_peak,
  cforest_plot_r0,
  cforest_plot_SE,
  cforest_plot_dom,
  labels = "AUTO", nrow = 2
)
combined_cforest

x.grob <- textGrob("Conditional Permutation Importance",
  gp = gpar(col = "black", fontsize = 18), hjust = 0.2
)

combined_cforest <- plot_grid(combined_cforest, x.grob, rel_heights = c(2, 0.1), nrow = 2)
combined_cforest
save_plot(combined_cforest, filename = "figures/cforest_variable_importance_H3_epi_measures.png", base_width = 18, base_height = 10)

######################################################################################
## LASSO variable importance
######################################################################################
h3_peak <- read.csv("data/H3_peak_ML_variable_importance.csv") %>% filter(model == "Lasso")
h3_peak$epi_measure <- "peak"

h3_size <- read.csv("data/H3_epi_size_ML_variable_importance.csv") %>% filter(model == "Lasso")
h3_size$epi_measure <- "cum_size"

h3_r0 <- read.csv("data/H3_R0_ML_variable_importance.csv") %>% filter(model == "Lasso")
h3_r0$epi_measure <- "R0"

h3_shannon <- read.csv("data/H3_shannon_entropy_ML_variable_importance.csv") %>% filter(model == "Lasso")
h3_shannon$epi_measure <- "shannon_entropy"

h3_dom <- read.csv("data/h3_dom_ML_variable_importance.csv") %>% filter(model == "Lasso")
h3_dom$epi_measure <- "dominance"

all_lasso <- bind_rows(h3_peak, h3_size, h3_dom, h3_r0, h3_shannon)
names(all_lasso)[1] <- "var"
all_lasso$epi_measure <- as.factor(all_lasso$epi_measure)
levels(all_lasso$epi_measure)
levels(all_lasso$epi_measure) <- c(
  "Epidemic Size", "Subtype Dominance",
  "Peak Incidence", "Effective Rt", "Epidemic Intensity"
)

label_vec <- c(
  "A(H1N1) epidemic size",
  "A(H1N1) epidemic size (t-1)",
  "A(H3N2) epidemic size (t-1)",
  "H3 LBI Diversity",
  "H3 LBI Diversity (t-1)",
  "HI titer (t-2)",
  "H3 epitope (t-2)",
  "B epidemic size",
  "B epidemic size (t-1)",
  "N2 epitope (t-1)",
  "N2 LBI Diversity",
  "N2 LBI Diversity (t-1)",
  "Dominant IAV (t-1)",
  "N2 distance to vaccine",
  "H3 distance to vaccine",
  "Vaccine coverage x VE",
  "Vaccine coverage x VE (t-1)",
  "Vaccine coverage",
  "Vaccine coverage (t-1)",
  "A(H3N2) VE",
  "A(H3N2) VE (t-1)"
)

label_df <- data.frame(var = sort(unique(all_lasso$var)), label = label_vec)
all_lasso <- left_join(all_lasso, label_df, by = "var")

lasso_limited <- all_lasso %>%
  group_by(epi_measure) %>%
  arrange(epi_measure, desc(importance)) %>%
  slice_max(order_by = importance, n = 21) %>%
  ungroup()

######################################################################################
## Epidemic size
######################################################################################

lasso_plot_epi_size <- ggplot(
  lasso_limited %>% filter(epi_measure == "Epidemic Size"),
  aes(
    x = reorder(label, importance), y = importance,
    fill = importance
  )
) +
  geom_bar(stat = "identity", position = "dodge", color = "black", lwd = 0.1) +
  ylab("") +
  xlab("") +
  coord_flip() +
  theme_bw(base_size = 16) +
  scale_fill_viridis(option = "rocket", direction = -1, begin = 0.5) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  ggtitle("Epidemic Size") +
  ylim(0, 100)
lasso_plot_epi_size

######################################################################################
## peak incidence
######################################################################################

lasso_plot_peak <- ggplot(
  lasso_limited %>% filter(epi_measure == "Peak Incidence"),
  aes(x = reorder(label, importance), y = importance, fill = importance)
) +
  geom_bar(stat = "identity", position = "dodge", color = "black", lwd = 0.1) +
  coord_flip() +
  ylab("") +
  xlab("") +
  theme_bw(base_size = 16) +
  scale_fill_viridis(option = "rocket", direction = -1, begin = 0.5) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  ggtitle("Peak Incidence") +
  ylim(0, 100)
lasso_plot_peak

######################################################################################
## subtype dominance
######################################################################################

lasso_plot_dom <- ggplot(
  lasso_limited %>% filter(epi_measure == "Subtype Dominance"),
  aes(x = reorder(label, importance), y = importance, fill = importance)
) +
  geom_bar(stat = "identity", position = "dodge", color = "black", lwd = 0.1) +
  ylab("") +
  xlab("") +
  coord_flip() +
  theme_bw(base_size = 16) +
  scale_fill_viridis(option = "rocket", direction = -1, begin = 0.5) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  ggtitle("Subtype Dominance") +
  ylim(0, 100)
lasso_plot_dom

######################################################################################
## effective Rt
######################################################################################
lasso_plot_r0 <- ggplot(
  lasso_limited %>% filter(epi_measure == "Effective Rt"),
  aes(x = reorder(label, importance), y = importance, fill = importance)
) +
  geom_bar(stat = "identity", position = "dodge", color = "black", lwd = 0.1) +
  coord_flip() +
  ylab("") +
  xlab("") +
  theme_bw(base_size = 16) +
  scale_fill_viridis(option = "rocket", direction = -1, begin = 0.5) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  ggtitle("Effective Rt") +
  ylim(0, 100)
lasso_plot_r0

######################################################################################
## inverse shannon entropy
######################################################################################
lasso_plot_SE <- ggplot(
  lasso_limited %>% filter(epi_measure == "Epidemic Intensity"),
  aes(x = reorder(label, importance), y = importance, fill = importance)
) +
  geom_bar(stat = "identity", position = "dodge", color = "black", lwd = 0.1) +
  coord_flip() +
  ylab("") +
  xlab("") +
  theme_bw(base_size = 16) +
  scale_fill_viridis(option = "rocket", direction = -1, begin = 0.5) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  ggtitle("Epidemic Intensity") +
  ylim(0, 100)
lasso_plot_SE

combined_lasso <- plot_grid(lasso_plot_epi_size,
  lasso_plot_peak,
  lasso_plot_r0,
  lasso_plot_SE,
  lasso_plot_dom,
  labels = "AUTO", nrow = 2
)
combined_lasso

x.grob <- textGrob("Variable Importance",
  gp = gpar(col = "black", fontsize = 18), hjust = 0.2
)

combined_lasso <- plot_grid(combined_lasso, x.grob, rel_heights = c(2, 0.1), nrow = 2)
combined_lasso
save_plot(combined_lasso, filename = "figures/lasso_variable_importance_H3_epi_measures.png", base_width = 18, base_height = 10)

#############################################
# Residuals
#############################################

epi_size_resid <- read.csv("data/epi_size_incidence_obs_vs_predicted.csv") %>% filter(object == "party RF")
epi_size_resid$metric <- "epi size"
epi_size_resid$resid <- epi_size_resid$H3_cum_intensity - epi_size_resid$pred
epi_size_resid$obs <- epi_size_resid$H3_cum_intensity
epi_size_resid$rmse <- sqrt(mean((epi_size_resid$resid)^2))

peak_resid <- read.csv("data/peak_incidence_obs_vs_predicted.csv") %>% filter(object == "party RF")
peak_resid$metric <- "peak inc"
peak_resid$resid <- peak_resid$H3_max_intensity - peak_resid$pred
peak_resid$obs <- peak_resid$H3_max_intensity
peak_resid$rmse <- sqrt(mean((peak_resid$resid)^2))

shannon_resid <- read.csv("data/shannon_entropy_obs_vs_predicted.csv") %>% filter(object == "party RF")
shannon_resid$metric <- "shannon"
shannon_resid$resid <- shannon_resid$H3_shannon_entropy_res - shannon_resid$pred
shannon_resid$obs <- shannon_resid$H3_shannon_entropy_res
shannon_resid$rmse <- sqrt(mean((shannon_resid$resid)^2))

r0_resid <- read.csv("data/r0_obs_vs_predicted.csv") %>% filter(object == "party RF")
r0_resid$metric <- "effective R"
r0_resid$resid <- r0_resid$H3_max_Rt - r0_resid$pred
r0_resid$obs <- r0_resid$H3_max_Rt
r0_resid$rmse <- sqrt(mean((r0_resid$resid)^2))

dom_resid <- read.csv("data/dom_obs_vs_predicted.csv") %>% filter(object == "party RF")
dom_resid$metric <- "dom"
dom_resid$resid <- dom_resid$h3_dom - dom_resid$pred
dom_resid$obs <- dom_resid$h3_dom
dom_resid$rmse <- sqrt(mean((dom_resid$resid)^2))

combined_resid <- bind_rows(epi_size_resid, peak_resid, dom_resid, shannon_resid, r0_resid) %>%
  dplyr::select(
    metric, season, region, pred, obs, resid, rmse, model, object, HA_wolf_lag2, NA_bhatt_ep_lag1,
    HA_titer_tree_lag2
  )
combined_resid$region <- factor(combined_resid$region,
  levels = c(
    "Region 1", "Region 2", "Region 3", "Region 4", "Region 5", "Region 6",
    "Region 7", "Region 8", "Region 9", "Region 10"
  )
)
levels(combined_resid$region) <- c(
  "1-Boston", "2-NYC", "3-DC", "4-Atlanta",
  "5-Chicago", "6-Dallas", "7-Kansas City",
  "8-Denver", "9-San Francisco", "10-Seattle"
)
unique(combined_resid$region)

######################################################################################
### epidemic size
######################################################################################
combined_resid %>%
  filter(metric == "epi size") %>%
  pull(obs) %>%
  max()
combined_resid %>%
  filter(metric == "epi size") %>%
  pull(pred) %>%
  max()
combined_resid %>%
  filter(metric == "epi size") %>%
  pull(obs) %>%
  min()
combined_resid %>%
  filter(metric == "epi size") %>%
  pull(pred) %>%
  min()

ggplot(combined_resid %>% filter(metric == "epi size"), aes(x = pred, y = obs)) +
  geom_point(aes(fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  scale_fill_viridis_c(name = "H3 epitope distance (t-2)") +
  scale_size_continuous(name = "H3 epitope distance (t-2)", range = c(1, 5)) +
  ylab("Observed Epidemic Size") +
  xlab("Predicted Epidemic Size") +
  coord_equal() +
  theme_bw() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cex = 3, method = "spearman") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(size = "none") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 99)
  ) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 99)
  )

combined_resid$HA_wolf_lag2 <- as.vector(scale(combined_resid$HA_wolf_lag2))
sort(unique(combined_resid$HA_wolf_lag2))

epi_size_plot <- ggplot(data = combined_resid %>% filter(metric == "epi size"), aes(x = pred, y = obs)) +
  geom_point(
    data = combined_resid %>% filter(metric == "epi size") %>% filter(HA_wolf_lag2 < 1.5),
    aes(x = pred, y = obs, fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21
  ) +
  geom_point(
    data = combined_resid %>% filter(metric == "epi size") %>% filter(HA_wolf_lag2 > 1.5),
    aes(x = pred, y = obs, fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21
  ) +
  facet_wrap(
    nrow = 1,
    ~ factor(region, levels = c(
      "1-Boston", "2-NYC", "3-DC", "4-Atlanta",
      "5-Chicago", "6-Dallas", "7-Kansas City",
      "8-Denver", "9-San Francisco", "10-Seattle"
    ))
  ) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  scale_fill_viridis_c(name = "H3 epitope distance (t-2)") +
  scale_size_continuous(name = "H3 epitope distance (t-2)", range = c(1, 5)) +
  ylab("Observed Epidemic Size") +
  xlab("Predicted Epidemic Size") +
  coord_equal() +
  theme_bw() +
  ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cex = 3, method = "spearman") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(size = "none") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 99)
  ) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 99)
  )
epi_size_plot

######################################################################################
### peak incidence
######################################################################################
combined_resid %>%
  filter(metric == "peak inc") %>%
  pull(obs) %>%
  max()
combined_resid %>%
  filter(metric == "peak inc") %>%
  pull(pred) %>%
  max()
combined_resid %>%
  filter(metric == "peak inc") %>%
  pull(obs) %>%
  min()
combined_resid %>%
  filter(metric == "peak inc") %>%
  pull(pred) %>%
  min()

ggplot(combined_resid %>% filter(metric == "peak inc"), aes(x = pred, y = obs)) +
  geom_point(aes(fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  scale_fill_viridis_c(name = "H3 epitope distance (t-2)") +
  scale_size_continuous(name = "H3 epitope distance (t-2)", range = c(1, 5)) +
  ylab("Observed Peak Incidence") +
  xlab("Predicted Peak Incidence") +
  coord_equal() +
  theme_bw() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cex = 3, method = "spearman") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(size = "none") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 22)
  ) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 22)
  )

peak_plot <- ggplot(data = combined_resid %>% filter(metric == "peak inc"), aes(x = pred, y = obs)) +
  geom_point(
    data = combined_resid %>% filter(metric == "peak inc") %>% filter(HA_wolf_lag2 < 1.5),
    aes(x = pred, y = obs, fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21
  ) +
  geom_point(
    data = combined_resid %>% filter(metric == "peak inc") %>% filter(HA_wolf_lag2 > 1.5),
    aes(x = pred, y = obs, fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21
  ) +
  facet_wrap(
    nrow = 1,
    ~ factor(region, levels = c(
      "1-Boston", "2-NYC", "3-DC", "4-Atlanta",
      "5-Chicago", "6-Dallas", "7-Kansas City",
      "8-Denver", "9-San Francisco", "10-Seattle"
    ))
  ) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  scale_fill_viridis_c(name = "H3 epitope distance (t-2)") +
  scale_size_continuous(name = "H3 epitope distance (t-2)", range = c(1, 5)) +
  ylab("Observed Peak Incidence") +
  xlab("Predicted Peak Incidence") +
  coord_equal() +
  theme_bw() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cex = 3, method = "spearman") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(size = "none") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 22)
  ) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 22)
  )
peak_plot

######################################################################################
## subtype dominance
######################################################################################
combined_resid %>%
  filter(metric == "dom") %>%
  pull(obs) %>%
  max()
combined_resid %>%
  filter(metric == "dom") %>%
  pull(pred) %>%
  max()
combined_resid %>%
  filter(metric == "dom") %>%
  pull(obs) %>%
  min()
combined_resid %>%
  filter(metric == "dom") %>%
  pull(pred) %>%
  min()

ggplot(combined_resid %>% filter(metric == "dom"), aes(x = pred, y = obs)) +
  geom_point(aes(fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  scale_fill_viridis_c(name = "H3 epitope distance (t-2)") +
  scale_size_continuous(name = "H3 epitope distance (t-2)", range = c(1, 5)) +
  # scale_color_manual(values = mycolors)+
  ylab("Observed A/H3 Dominance") +
  xlab("Predicted A/H3 Dominance") +
  coord_equal() +
  theme_bw() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cex = 3, method = "spearman") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(size = "none") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 1)
  ) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 1)
  )

dom_plot <- ggplot(data = combined_resid %>% filter(metric == "dom"), aes(x = pred, y = obs)) +
  geom_point(
    data = combined_resid %>% filter(metric == "dom") %>% filter(HA_wolf_lag2 < 1.5),
    aes(x = pred, y = obs, fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21
  ) +
  geom_point(
    data = combined_resid %>% filter(metric == "dom") %>% filter(HA_wolf_lag2 > 1.5),
    aes(x = pred, y = obs, fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21
  ) +
  facet_wrap(
    nrow = 1,
    ~ factor(region, levels = c(
      "1-Boston", "2-NYC", "3-DC", "4-Atlanta",
      "5-Chicago", "6-Dallas", "7-Kansas City",
      "8-Denver", "9-San Francisco", "10-Seattle"
    ))
  ) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  scale_fill_viridis_c(name = "H3 epitope distance (t-2)") +
  scale_size_continuous(name = "H3 epitope distance (t-2)", range = c(1, 5)) +
  ylab("Observed Dominance") +
  xlab("Predicted Dominance") +
  coord_equal() +
  theme_bw() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cex = 3, method = "spearman") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(size = "none") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 1)
  ) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 1)
  )
dom_plot

######################################################################################
#### inverse shannon entropy
######################################################################################
combined_resid %>%
  filter(metric == "shannon") %>%
  pull(obs) %>%
  max()
combined_resid %>%
  filter(metric == "shannon") %>%
  pull(pred) %>%
  max()
combined_resid %>%
  filter(metric == "shannon") %>%
  pull(obs) %>%
  min()
combined_resid %>%
  filter(metric == "shannon") %>%
  pull(pred) %>%
  min()


ggplot(
  data = combined_resid %>% filter(metric == "shannon"),
  aes(x = pred, y = obs)
) +
  geom_point(
    data = combined_resid %>% filter(metric == "shannon") %>% filter(HA_wolf_lag2 < 9),
    aes(x = pred, y = obs, fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21
  ) +
  geom_point(
    data = combined_resid %>% filter(metric == "shannon") %>% filter(HA_wolf_lag2 > 9),
    aes(x = pred, y = obs, fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21
  ) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  scale_fill_viridis_c(name = "H3 epitope distance (t-2)") +
  scale_size_continuous(name = "H3 epitope distance (t-2)", range = c(1, 5)) +
  ylab("Observed Epidemic Intensity") +
  xlab("Predicted Epidemic Intensity") +
  coord_equal() +
  theme_bw() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cex = 3, method = "spearman") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(size = "none") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 1)
  ) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 1)
  )

shannon_plot <- ggplot(combined_resid %>% filter(metric == "shannon"), aes(x = pred, y = obs)) +
  geom_point(
    data = combined_resid %>% filter(metric == "shannon") %>% filter(HA_wolf_lag2 < 1.5),
    aes(x = pred, y = obs, fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21
  ) +
  geom_point(
    data = combined_resid %>% filter(metric == "shannon") %>% filter(HA_wolf_lag2 > 1.5),
    aes(x = pred, y = obs, fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21
  ) +
  facet_wrap(
    nrow = 1,
    ~ factor(region, levels = c(
      "1-Boston", "2-NYC", "3-DC", "4-Atlanta",
      "5-Chicago", "6-Dallas", "7-Kansas City",
      "8-Denver", "9-San Francisco", "10-Seattle"
    ))
  ) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  scale_fill_viridis_c(name = "H3 epitope distance (t-2)") +
  scale_size_continuous(name = "H3 epitope distance (t-2)", range = c(1, 5)) +
  ylab("Observed Epidemic Intensity") +
  xlab("Predicted Epidemic Intensity") +
  coord_equal() +
  theme_bw() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cex = 3, method = "spearman") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(size = "none") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 1)
  ) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0, 1)
  )
shannon_plot

######################################################################################
##### Effective Rt
######################################################################################
combined_resid %>%
  filter(metric == "effective R") %>%
  pull(obs) %>%
  max()
combined_resid %>%
  filter(metric == "effective R") %>%
  pull(pred) %>%
  max()
combined_resid %>%
  filter(metric == "effective R") %>%
  pull(obs) %>%
  min()
combined_resid %>%
  filter(metric == "effective R") %>%
  pull(pred) %>%
  min()

ggplot(combined_resid %>% filter(metric == "effective R"), aes(x = pred, y = obs)) +
  geom_point(aes(fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  scale_fill_viridis_c(name = "H3 epitope distance (t-2)") +
  scale_size_continuous(name = "H3 epitope distance (t-2)", range = c(1, 5)) +
  ylab("Observed Rt") +
  xlab("Predicted Rt") +
  coord_equal() +
  theme_bw() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cex = 3, method = "spearman") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(size = "none") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0.86, 2.2)
  ) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0.86, 2.2)
  )

r0_plot <- ggplot(data = combined_resid %>% filter(metric == "effective R"), aes(x = pred, y = obs)) +
  geom_point(
    data = combined_resid %>% filter(metric == "effective R") %>% filter(HA_wolf_lag2 < 9),
    aes(x = pred, y = obs, fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21
  ) +
  geom_point(
    data = combined_resid %>% filter(metric == "effective R") %>% filter(HA_wolf_lag2 > 9),
    aes(x = pred, y = obs, fill = HA_wolf_lag2, size = HA_wolf_lag2), alpha = 0.8, pch = 21
  ) +
  facet_wrap(
    nrow = 1,
    ~ factor(region, levels = c(
      "1-Boston", "2-NYC", "3-DC", "4-Atlanta",
      "5-Chicago", "6-Dallas", "7-Kansas City",
      "8-Denver", "9-San Francisco", "10-Seattle"
    ))
  ) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  scale_fill_viridis_c(name = "H3 epitope distance (t-2)") +
  scale_size_continuous(name = "H3 epitope distance (t-2)", range = c(1, 5)) +
  ylab("Observed Rt") +
  xlab("Predicted Rt") +
  coord_equal() +
  theme_bw() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cex = 3, method = "spearman") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(size = "none") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0.86, 2.2)
  ) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.1), limits = c(0.86, 2.2)
  )
r0_plot

all_resid_plot <- plot_grid(epi_size_plot + theme(legend.position = "none"),
  peak_plot + theme(legend.position = "none"),
  r0_plot + theme(legend.position = "none"),
  shannon_plot + theme(legend.position = "none"),
  dom_plot + theme(legend.position = "none"),
  nrow = 5, labels = "AUTO"
)
leg <- get_legend(epi_size_plot + theme(
  legend.direction = "horizontal",
  legend.justification = "center",
  legend.box.just = "bottom",
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 12)
))
all_resid_plot2 <- plot_grid(all_resid_plot, leg, rel_heights = c(5, 0.4), nrow = 2)
all_resid_plot2
save_plot(all_resid_plot2, filename = "figures/cond_inf_forest_residual_plots_by_region_and_epi_metric.png", base_width = 15, base_height = 14)

combined_resid %>%
  group_by(season, metric, HA_wolf_lag2, NA_bhatt_ep_lag1) %>%
  summarize(residual = sum(resid)) %>%
  ungroup() %>%
  arrange(metric, -residual) %>%
  group_by(metric) %>%
  slice_max(residual, n = 2)

combined_resid %>%
  filter(metric %in% c("peak inc")) %>%
  dplyr::select(HA_wolf_lag2, season, region, metric, obs) %>%
  arrange(metric, region, -obs) %>%
  distinct()

######################################################################################
#### Seasonal RMSE
######################################################################################
seasonal_plots <- combined_resid %>%
  group_by(season, metric, HA_wolf_lag2, NA_bhatt_ep_lag1) %>%
  summarize(rmse = sqrt(mean((resid)^2))) %>%
  ungroup()
seasonal_plots$metric <- as.factor(seasonal_plots$metric)
levels(seasonal_plots$metric)
levels(seasonal_plots$metric) <- c("Subtype Dominance", "Effective Rt", "Epidemic Size", "Peak Incidence", "Epidemic Intensity")
seasonal_plots$metric <- factor(seasonal_plots$metric, levels = c("Epidemic Size", "Peak Incidence", "Effective Rt", "Epidemic Intensity", "Subtype Dominance"))

seasonal_plots %>%
  group_by(metric) %>%
  slice_max(rmse, n = 2)

p <- ggplot(seasonal_plots, aes(x = scale(HA_wolf_lag2), y = rmse)) +
  geom_point(aes(fill = HA_wolf_lag2), alpha = 0.8, pch = 21, size = 5) +
  scale_fill_viridis_c(name = "H3 epitope distance (t-2)") +
  geom_smooth(se = F, lty = "dashed", method = "lm", color = "black") +
  # scale_size_continuous(name="H3 epitope distance (t-2)",range=c(1,5))+
  facet_wrap(~metric, scales = "free") +
  stat_cor(p.accuracy = 0.01, r.accuracy = 0.01, cex = 4, method = "spearman") +
  xlab("H3 epitope distance (t-2)") +
  ylab("Seasonal RMSE") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none", legend.direction = "horizontal")
p
save_plot(p, filename = "figures/model_total_rmse_by_metric.png", base_width = 10, base_height = 8)

p <- ggplot(seasonal_plots, aes(x = scale(NA_bhatt_ep_lag1), y = rmse)) +
  geom_point(aes(fill = NA_bhatt_ep_lag1), alpha = 0.8, pch = 21, size = 5) +
  scale_fill_viridis_c(name = "N2 epitope distance (t-1)") +
  geom_smooth(se = F, lty = "dashed", method = "lm", color = "black") +
  facet_wrap(~metric, scales = "free") +
  stat_cor(p.accuracy = 0.01, r.accuracy = 0.01, cex = 4, method = "spearman") +
  xlab("N2 epitope distance (t-1)") +
  ylab("Seasonal RMSE") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none", legend.direction = "horizontal")
p
save_plot(p, filename = "figures/model_total_rmse_by_metric_NA.png", base_width = 10, base_height = 8)
