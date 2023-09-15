## load packages
list.of.packages <- c(
  "dplyr", "ggplot2", "cowplot", "tidyr", "viridis", "RColorBrewer", "corrplot", "rstatix"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

theme_set(theme_cowplot(font_size = 16))
########################################################################################
## Epi Metrics Scatterplot
########################################################################################
## load
load("data/antigenic_epi_north_amer_build_for_lasso_replicates.Rdata")
load("data/subtype_distribution_by_region_season.RData")

subtype_dist <- subtype_dist %>%
  ungroup() %>%
  mutate(a_total = h3_total + h1_total) %>%
  mutate(
    h3_vs_h1 = h3_total / a_total, # prop of h3 out of IAV
    h3_dom = h3_total / (a_total + b_total) ## prop of h3 out of all flu positives
  ) %>%
  rename(season = season_description)
names(subtype_dist)

epi_red2 <- left_join(subtype_dist,
  epi_red %>% dplyr::select(
    -h3_vs_h1, -iva_vs_ivb, -h3_vs_total_flu,
    -h1_vs_total_flu, -ivb_vs_iva, -total_a
  ),
  by = c("season", "region")
) %>%
  distinct() %>%
  filter(!(season %in% c("2009-2010", "2019-2020", "2020-2021", "2021-2022")))

ag_df <- epi_red2 %>%
  dplyr::select(
    season, region,
    H3_season_duration,
    onset_days_from_Oct1,
    onset_timing_sd, peak_timing_sd,
    peak_days_from_Oct1, peak_diff,
    H3_shannon_entropy_res,
    H3_max_Rt, h3_dom,
    H3_cum_intensity,
    H3_max_intensity,
    H1_cum_intensity, IVB_cum_intensity
  ) %>%
  distinct()

corr <- round(cor(
  ag_df %>%
    dplyr::select(-season, -region) %>%
    mutate_at(vars(H3_season_duration:IVB_cum_intensity), as.numeric),
  method = "spearman", use = "pairwise.complete.obs"
), 2)

p.mat <- cor_pmat(
  ag_df %>%
    dplyr::select(-season, -region) %>%
    mutate_at(vars(H3_season_duration:IVB_cum_intensity), as.numeric),
  method = "spearman"
) %>%
  dplyr::select(-rowname)

new_mat <- matrix(p.adjust(as.vector(as.matrix(p.mat)), method = "BH"), ncol = ncol(p.mat))

colnames(corr) <- c(
  "H3N2 Season Duration",
  "H3N2 Onset Week",
  "H3N2 Onset Timing s.d.",
  "H3N2 Peak Timing s.d.",
  "H3N2 Peak Week",
  "Days from Onset to Peak",
  "H3N2 Epidemic Intensity",
  "H3N2 Effective Rt",
  "H3N2 Subtype Dominance",
  "H3N2 Epidemic Size",
  "H3N2 Peak Incidence",
  "H1N1 Epidemic Size",
  "B Epidemic Size"
)

rownames(corr) <- c(
  "H3N2 Season Duration",
  "H3N2 Onset Week",
  "H3N2 Onset Timing s.d.",
  "H3N2 Peak Timing s.d.",
  "H3N2 Peak Week",
  "Days from Onset to Peak",
  "H3N2 Epidemic Intensity",
  "H3N2 Effective Rt",
  "H3N2 Subtype Dominance",
  "H3N2 Epidemic Size",
  "H3N2 Peak Incidence",
  "H1N1 Epidemic Size",
  "B Epidemic Size"
)

colnames(new_mat) <- colnames(corr)
rownames(new_mat) <- rownames(corr)

corrplot(corr,
  type = "lower", diag = F, outline = T, col = rev(COL2("RdBu", 200)),
  order = "original", p.mat = new_mat, method = "circle",
  insig = "label_sig", tl.col = "black"
)

pdf("figures/epi_metrics_pairwise_correlations.pdf", height = 8, width = 8, paper = "letter")
## Create a graphical object g here
corrplot(corr,
  type = "lower", diag = F, outline = T, col = rev(COL2("RdBu", 200)),
  order = "original", p.mat = new_mat, method = "circle",
  insig = "label_sig", tl.col = "black"
)
## Stop writing to the PDF file
dev.off()
