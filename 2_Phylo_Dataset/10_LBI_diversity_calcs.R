## load packages
list.of.packages <- c("dplyr", "ggplot2", "readr", "tidyr", "vegan", "lubridate")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

############################################################################
## Calculate seasonal Shannon entropy of LBI values
############################################################################

lbi_tb <- read_tsv("2_Phylo_Dataset/distance_tables/north-america_strain_seasonal_lbi_h3n2_ha_21y.tsv.gz")

lbi_tb$year <- lubridate::year(lbi_tb$season_end)
range(lbi_tb$lbi)

breaks <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)

tags <- c("(0-2]", "(2-4]", "(4-6]", "(6-8]", "(8-10]", "(10-12]", "(12-14]", "(14-16]", "(16-18]", "(18-20]")

lbi_tb$lbi_bin <- cut(lbi_tb$lbi,
  breaks = breaks,
  include.lowest = TRUE,
  right = T,
  labels = tags
)

lbi_tb_sum <- lbi_tb %>%
  group_by(year, lbi_bin, replicate) %>%
  tally() %>%
  ungroup()

ggplot(data = lbi_tb_sum, aes(x = year, fill = lbi_bin)) +
  geom_bar(position = "fill") +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  scale_x_continuous(breaks = seq(from = 1997, to = 2020, 2), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~replicate)

lbi_tb_sum_wide <- lbi_tb_sum %>% pivot_wider(values_from = n, names_from = lbi_bin)
lbi_tb_sum_wide[is.na(lbi_tb_sum_wide)] <- 0

lbi_tb_sum_wide$ha_lbi_shannon <- diversity(lbi_tb_sum_wide, index = "shannon", MARGIN = 1)

# lbi_tb_sum_wide$ha_lbi_inv_simpson <- diversity(lbi_tb_sum_wide, index = "invsimpson", MARGIN = 1)

lbi_tb_na <- read_tsv("2_Phylo_Dataset/distance_tables/north-america_strain_seasonal_lbi_h3n2_na_21y.tsv.gz")
lbi_tb_na$year <- lubridate::year(lbi_tb_na$season_end)
range(lbi_tb_na$lbi)

breaks <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)

tags <- c("(0-2]", "(2-4]", "(4-6]", "(6-8]", "(8-10]", "(10-12]", "(12-14]", "(14-16]", "(16-18]", "(18-20]")

lbi_tb_na$lbi_bin <- cut(lbi_tb_na$lbi,
  breaks = breaks,
  include.lowest = TRUE,
  right = T,
  labels = tags
)

lbi_tb_na_sum <- lbi_tb_na %>%
  group_by(year, replicate, lbi_bin) %>%
  tally() %>%
  ungroup()

ggplot(data = lbi_tb_na_sum, aes(x = year, fill = lbi_bin)) +
  geom_bar(position = "fill") +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  scale_x_continuous(breaks = seq(from = 1997, to = 2020, 2), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~replicate)

lbi_tb_na_sum_wide <- lbi_tb_na_sum %>% pivot_wider(values_from = n, names_from = lbi_bin)
lbi_tb_na_sum_wide[is.na(lbi_tb_na_sum_wide)] <- 0

lbi_tb_na_sum_wide$na_lbi_shannon <- diversity(lbi_tb_na_sum_wide, index = "shannon", MARGIN = 1)

# lbi_tb_na_sum_wide$na_lbi_inv_simpson <- diversity(lbi_tb_na_sum_wide, index = "invsimpson", MARGIN = 1)

lbi_eco_div <- full_join(
  lbi_tb_sum_wide %>% dplyr::select(year, replicate, contains("ha")),
  lbi_tb_na_sum_wide %>% dplyr::select(year, replicate, contains("na"))
)

lbi_eco_div <- lbi_eco_div %>%
  mutate(year1 = year - 1) %>%
  mutate(season = paste(year1, year, sep = "-"))

write_rds(lbi_eco_div, file = "data/LBI_diversity_by_season.rds")
