## load packages
list.of.packages <- c(
  "dplyr", "ggplot2", "cowplot", "tidyr", "viridis", "RColorBrewer", "ggrepel", "corrplot", "rstatix", "ggpubr"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

theme_set(theme_cowplot(font_size = 16))

############################################
## HI titer vs H3 epitope
############################################
## load
load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")

ag_df <- epi_red %>%
  dplyr::select(season, dom_type, HA_wolf_lag1, HA_wolf_lag2, contains("titer")) %>%
  dplyr::select(-contains("lag0")) %>%
  distinct()
unique(ag_df$season)

range(scale(ag_df$HA_titer_tree_lag1))
range(scale(ag_df$HA_wolf_lag1))

p <- ggplot(data = ag_df, aes(x = scale(HA_wolf_lag1), y = scale(HA_titer_tree_lag1))) +
  geom_point(aes(x = scale(HA_wolf_lag1), y = scale(HA_titer_tree_lag1), fill = season), pch = 21, size = 8) +
  xlab(expression("H3 epitope distance ("~italic(t) ~ "-1)")) +
  ylab(expression("HI titer distance ("~italic(t) ~ "-1)")) +
  theme(legend.position = "right") +
  scale_fill_viridis(option = "magma", discrete = T, name = "Season") +
  theme_cowplot(font_size = 16) +
  background_grid(major = "xy") +
  stat_cor(
    method = "spearman",
    p.digits = 1,
    label.x = 1.7,
    label.y = 2.7,
    size = 4,
    cor.coef.name = "rho"
  ) +
  coord_equal() +
  xlim(c(-2.01, 2.8)) +
  ylim(c(-2.01, 2.8))
p
p <- p + ggrepel::geom_label_repel(
  data = ag_df,
  mapping = aes(x = scale(HA_wolf_lag1), y = scale(HA_titer_tree_lag1), label = season),
  point.padding = 0.25, size = 3, min.segment.length = 0
) + theme(legend.position = "none")
p

range(scale(ag_df$HA_titer_tree_lag2))
range(scale(ag_df$HA_wolf_lag2))
q <- ggplot(data = ag_df, aes(x = scale(HA_wolf_lag2), y = scale(HA_titer_tree_lag2))) +
  geom_point(aes(x = scale(HA_wolf_lag2), y = scale(HA_titer_tree_lag2), fill = season), pch = 21, size = 8) +
  xlab(expression("H3 epitope distance ("~italic(t) ~ "-2)")) +
  ylab(expression("HI titer distance ("~italic(t) ~ "-2)")) +
  theme(legend.position = "right") +
  scale_fill_viridis(option = "magma", discrete = T, name = "Season") +
  theme_cowplot(font_size = 16) +
  background_grid(major = "xy") +
  stat_cor(
    method = "spearman",
    p.digits = 1,
    label.x = 1.7,
    label.y = 2.7,
    size = 4,
    cor.coef.name = "rho"
  ) +
  coord_equal() +
  xlim(c(-2.01, 2.8)) +
  ylim(c(-2.01, 2.8))
q
q <- q + ggrepel::geom_label_repel(
  data = ag_df,
  mapping = aes(x = scale(HA_wolf_lag2), y = scale(HA_titer_tree_lag2), label = season),
  point.padding = 0.25, size = 3, min.segment.length = 0
) + theme(legend.position = "none")
q

com <- plot_grid(p, q, labels = "AUTO")
com

save_plot(com, filename = "figures/H3_ep_vs_HI_tree_distance_scatterplot_both_lags.png", base_width = 12, base_height = 6)
save_plot(com, filename = "figures/H3_ep_vs_HI_tree_distance_scatterplot_both_lags.pdf", dpi = 300, base_width = 12, base_height = 6)

########################################################################################
## H3 vs N2
########################################################################################

ag_df <- epi_red %>%
  dplyr::select(season, dom_type, HA_wolf_lag1, HA_wolf_lag2, NA_bhatt_ep_lag1, NA_bhatt_ep_lag2) %>%
  distinct()
head(ag_df)

range(scale(ag_df$HA_wolf_lag1))
range(scale(ag_df$HA_wolf_lag2))
range(scale(ag_df$NA_bhatt_ep_lag1))
range(scale(ag_df$NA_bhatt_ep_lag2))

p <- ggplot(data = ag_df, aes(x = scale(HA_wolf_lag1), y = scale(NA_bhatt_ep_lag1))) +
  geom_point(aes(x = scale(HA_wolf_lag1), y = scale(NA_bhatt_ep_lag1), fill = season), pch = 21, size = 8) +
  xlab(expression("H3 epitope distance ("~italic(t) ~ "-1)")) +
  ylab(expression("N2 epitope distance ("~italic(t) ~ "-1)")) +
  theme(legend.position = "right") +
  scale_fill_viridis(option = "magma", discrete = T, name = "Season") +
  theme_cowplot(font_size = 16) +
  background_grid(major = "xy") +
  stat_cor(
    method = "spearman",
    label.x = 1.7,
    label.y = 2.7,
    size = 4,
    p.digits = 1,
    cor.coef.name = "rho"
  ) +
  coord_equal() +
  xlim(c(-1.5, 2.8)) +
  ylim(c(-1.5, 2.8))

p <- p + ggrepel::geom_label_repel(
  data = ag_df,
  mapping = aes(x = scale(HA_wolf_lag1), y = scale(NA_bhatt_ep_lag1), label = season),
  point.padding = 0.25, size = 3, min.segment.length = 0
) + theme(legend.position = "none")
p

q <- ggplot(data = ag_df, aes(x = scale(HA_wolf_lag2), y = scale(NA_bhatt_ep_lag2))) +
  geom_point(aes(x = scale(HA_wolf_lag2), y = scale(NA_bhatt_ep_lag2), fill = season), pch = 21, size = 8) +
  xlab(expression("H3 epitope distance (" ~ italic(t) ~ "-2)")) +
  ylab(expression("N2 epitope distance (" ~ italic(t) ~ "-2)")) +
  theme(legend.position = "right") +
  scale_fill_viridis(option = "magma", discrete = T, name = "Season") +
  theme_cowplot(font_size = 16) +
  background_grid(major = "xy") +
  stat_cor(
    method = "spearman",
    label.x = 1.7,
    label.y = 2.7,
    p.digits = 1,
    size = 4,
    cor.coef.name = "rho"
  ) +
  coord_equal() +
  xlim(c(-1.5, 2.8)) +
  ylim(c(-1.5, 2.8))

q <- q + ggrepel::geom_label_repel(
  data = ag_df,
  mapping = aes(x = scale(HA_wolf_lag2), y = scale(NA_bhatt_ep_lag2), label = season),
  point.padding = 0.25, size = 3, min.segment.length = 0
) + theme(legend.position = "none")
q

comb <- plot_grid(p, q, nrow = 1, labels = "AUTO")
comb

save_plot(comb, filename = "figures/H3_vs_N2_epitope_distance_scatterplot_both_lags.png", base_width = 12, base_height = 6)
save_plot(comb, filename = "figures/H3_vs_N2_epitope_distance_scatterplot_both_lags.pdf",dpi=300, base_width = 12, base_height = 6)

########################################################################################
## Evolutionary Indicators Scatterplot
## One season lag (t-1)
########################################################################################
ag_df <- epi_red %>%
  dplyr::select(season, dom_type, contains(c("lag1","lbi"))) %>%
  distinct()
head(ag_df)

ag_df <- ag_df %>%
  dplyr::select(season, contains(c("koel", "wolf", "wolf_nonepitope", "stem_ep", "titer", "bhatt", "krammer", "lbi"))) %>%
  dplyr::select(-contains("sub"))

names(ag_df)
corr <- round(cor(
  ag_df %>% dplyr::select(-season) %>% mutate_at(vars(HA_koel_lag1:na_lbi_shannon), as.numeric),
  use = "pairwise.complete.obs", method = "spearman"
), 2)

p.mat <- cor_pmat(
  ag_df %>% dplyr::select(-season) %>% mutate_at(vars(HA_koel_lag1:na_lbi_shannon), as.numeric),
  method = "spearman"
) %>%
  dplyr::select(-rowname)

new_mat <- matrix(p.adjust(as.vector(as.matrix(p.mat)), method = "BH"), ncol = ncol(p.mat))

colnames(corr)

colnames(corr) <- c(
  "H3 RBS Epitope (t-1)",
  "H3 Epitope (t-1)",
  "H3 Non-Epitope (t-1)",
  "H3 Stalk Footprint (t-1)",
  "HI Titer (t-1)",
  "N2 Epitope (N=223) (t-1)",
  "N2 Non-Epitope (t-1)",
  "N2 Epitope (N=53) (t-1)",
  "H3 s.d. LBI t",
  "N2 s.d. LBI t",
  "H3 LBI Diversity t",
  "N2 LBI Diversity t"
)

rownames(corr) <- c(
  "H3 RBS Epitope (t-1)",
  "H3 Epitope (t-1)",
  "H3 Non-Epitope (t-1)",
  "H3 Stalk Footprint (t-1)",
  "HI Titer (t-1)",
  "N2 Epitope (N=223) (t-1)",
  "N2 Non-Epitope (t-1)",
  "N2 Epitope (N=53) (t-1)",
  "H3 s.d. LBI t",
  "N2 s.d. LBI t",
  "H3 LBI Diversity t",
  "N2 LBI Diversity t"
)

colnames(new_mat) <- colnames(corr)
rownames(new_mat) <- rownames(corr)

corrplot(corr,
  type = "lower", diag = F, outline = T, col = rev(COL2("RdBu", 200)),
  order = "original", p.mat = new_mat, method = "circle",
  insig = "label_sig", tl.col = "black"
)

pdf("figures/evol_metrics_pairwise_correlations_lag1.pdf", height = 8, width = 10, paper = "letter")
## Create a graphical object g here
corrplot(corr,
  type = "lower", diag = F, outline = T, col = rev(COL2("RdBu", 200)),
  order = "original", p.mat = new_mat, method = "circle",
  insig = "label_sig", tl.col = "black"
)
## Stop writing to the PDF file
dev.off()

########################################################################################
## Evolutionary Indicators Scatterplot
## Two season lag (t-2)
########################################################################################
ag_df <- epi_red %>%
  dplyr::select(season, dom_type, contains(c("lag2","lbi"))) %>%
  distinct()

ag_df <- ag_df %>%
  dplyr::select(season, contains(c("koel", "wolf", "wolf_nonepitope", "stem_ep", "titer", "bhatt", "krammer", "lbi"))) %>%
  dplyr::select(-contains("sub"))
names(ag_df)

corr <- round(cor(
  ag_df %>% dplyr::select(-season) %>%mutate_at(vars(HA_koel_lag2:na_lbi_shannon), as.numeric),
  method = "spearman", 
  use = "pairwise.complete.obs"
), 2)

p.mat <- cor_pmat(
  ag_df %>%
    dplyr::select(-season) %>%
    mutate_at(vars(HA_koel_lag2:na_lbi_shannon), as.numeric),
  method = "spearman"
) %>%
  dplyr::select(-rowname)

new_mat <- matrix(p.adjust(as.vector(as.matrix(p.mat)), method = "BH"), ncol = ncol(p.mat))

colnames(corr) 
colnames(corr) <- c(
  "H3 RBS Epitope (t-2)",
  "H3 Epitope (t-2)",
  "H3 Non-Epitope (t-2)",
  "H3 Stalk Footprint (t-2)",
  "HI Titer (t-2)",
  "N2 Epitope (N=223) (t-2)",
  "N2 Non-Epitope (t-2)",
  "N2 Epitope (N=53) (t-2)",
  "H3 s.d. LBI t",
  "N2 s.d. LBI t",
  "H3 LBI Diversity t",
  "N2 LBI Diversity t"
)

rownames(corr) <- c(
  "H3 RBS Epitope (t-2)",
  "H3 Epitope (t-2)",
  "H3 Non-Epitope (t-2)",
  "H3 Stalk Footprint (t-2)",
  "HI Titer (t-2)",
  "N2 Epitope (N=223) (t-2)",
  "N2 Non-Epitope (t-2)",
  "N2 Epitope (N=53) (t-2)",
  "H3 s.d. LBI t",
  "N2 s.d. LBI t",
  "H3 LBI Diversity t",
  "N2 LBI Diversity t"
)

colnames(new_mat) <- colnames(corr)
rownames(new_mat) <- rownames(corr)

corrplot(corr,
  type = "lower", diag = F, outline = T, col = rev(COL2("RdBu", 200)),
  order = "original", p.mat = new_mat, method = "circle",
  insig = "label_sig", tl.col = "black"
)

pdf("figures/evol_metrics_pairwise_correlations_lag2.pdf", height = 8, width = 10, paper = "letter")
## Create a graphical object g here
corrplot(corr,
  type = "lower", diag = F, outline = T, col = rev(COL2("RdBu", 200)),
  order = "original", p.mat = new_mat, method = "circle",
  insig = "label_sig", tl.col = "black"
)
## Stop writing to the PDF file
dev.off()
