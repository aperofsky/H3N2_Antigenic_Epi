## load packages
list.of.packages <- c(
  "dplyr", "ggplot2", "cowplot", "tidyr", "viridis", "RColorBrewer", "ggrepel", "corrplot", "rstatix", "ggpubr"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

theme_set(theme_cowplot(font_size = 16))

# load data
load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")
unique(epi_red$season)

#### correlations among all predictors
ag_df <- epi_red %>%
  dplyr::select(season, dom_type, contains(c("lag","lbi"))) %>%
  distinct()
head(ag_df)
names(ag_df)

corr <- round(cor(
  ag_df %>% dplyr::select(-season,-dom_type) %>% mutate_at(vars(HA_koel_lag1:na_lbi_shannon), as.numeric),
  use = "pairwise.complete.obs", method = "spearman"
), 2)

p.mat <- rstatix::cor_pmat(
  ag_df %>% dplyr::select(-season,-dom_type) %>% mutate_at(vars(HA_koel_lag1:na_lbi_shannon), as.numeric),
  method = "spearman"
) %>%
  dplyr::select(-rowname)

new_mat <- matrix(p.adjust(as.vector(as.matrix(p.mat)), method = "BH"), ncol = ncol(p.mat))
dim(corr)
dim(new_mat)

colnames(corr) <- c(
  "H3 RBS Epitope (t-1)",
  "H3 Stalk Footprint (t-1)",
  "HI Titer (t-1)",
  "H3 Epitope (t-1)",
  "H3 Non-Epitope (t-1)",
  "N2 Epitope (N=223) (t-1)",
  "N2 Non-Epitope (t-1)",
  "N2 Epitope (N=53) (t-1)",
  "H3 RBS Epitope (t-2)",
  "H3 Stalk Footprint (t-2)",
  "HI Titer (t-2)",
  "H3 Epitope (t-2)",
  "H3 Non-Epitope (t-2)",
  "N2 Epitope (N=223) (t-2)",
  "N2 Non-Epitope (t-2)",
  "N2 Epitope (N=53) (t-2)",
  "H3 s.d. LBI t",
  "N2 s.d. LBI t",
  "H3 LBI Diversity t",
  "N2 LBI Diversity t"
)

rownames(corr) <- c(
  "H3 RBS Epitope (t-1)",
  "H3 Stalk Footprint (t-1)",
  "HI Titer (t-1)",
  "H3 Epitope (t-1)",
  "H3 Non-Epitope (t-1)",
  "N2 Epitope (N=223) (t-1)",
  "N2 Non-Epitope (t-1)",
  "N2 Epitope (N=53) (t-1)",
  "H3 RBS Epitope (t-2)",
  "H3 Stalk Footprint (t-2)",
  "HI Titer (t-2)",
  "H3 Epitope (t-2)",
  "H3 Non-Epitope (t-2)",
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

# pdf("figures/Fig2_sup_fig6_evol_metrics_pairwise_correlations_all_lags.pdf", height = 8, width = 10, paper = "letter")
png("figures/Fig2_sup_fig6_evol_metrics_pairwise_correlations_all_lags.png", height = 8, width = 10, units="in",res=300,bg = "white")
## Create a graphical object g here
corrplot(corr,
         type = "lower", diag = F, outline = T, col = rev(COL2("RdBu", 200)),
         order = "original", p.mat = new_mat, method = "circle",
         insig = "label_sig", tl.col = "black"
)
dev.off()
