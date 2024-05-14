##################################################
## Sequence Volume by Month and Season
##################################################

## load packages
list.of.packages <- c("dplyr", "tidyr", "readr", "ggplot2", "zoo", "jcolors", "cowplot")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

##################################################
## H3 plots and summary table
##################################################

HA_north_amer_build <- readr::read_tsv("2_Phylo_Dataset/auspice_tables/flu_seasonal_h3n2_ha_21y_north-america.tsv")
head(HA_north_amer_build)
HA_north_amer_build <- HA_north_amer_build %>%
  mutate(
    month.day = format(date, "%m-%d"),
    month = lubridate::month(date),
    year = lubridate::year(date),
    month.yr = as.Date(as.yearmon(date, "%b-%Y"))
  ) %>%
  mutate(season = ifelse(month.day < "07-01", sprintf("%d-%d", year - 1, year), # season = july 1 to june 30
    sprintf("%d-%d", year, year + 1)
  ))

unique(HA_north_amer_build$season)

HA_seq_table <- HA_north_amer_build %>%
  filter(!(season %in% c("1995-1996", "1996-1997", "2019-2020"))) %>%
  distinct() %>%
  group_by(replicate, season, month.yr, region) %>%
  tally() %>%
  pivot_wider(names_from = region, values_from = n) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(replicate, season, month.yr, north_america, south_america, europe, africa, oceania, china, southeast_asia, south_asia, west_asia, japan_korea) %>%
  arrange(season, month.yr, replicate) %>%
  rowwise() %>%
  mutate(other_regions = sum(c_across(south_america:japan_korea))) %>%
  mutate(total_sequences = north_america + other_regions) %>%
  dplyr::select(replicate, season, month.yr, total_sequences, north_america, other_regions, south_america:japan_korea)

HA_seq_table$replicate <- as.factor(HA_seq_table$replicate)
unique(HA_seq_table$replicate)
levels(HA_seq_table$replicate) <- c("0", "1", "2", "3", "4")
head(HA_seq_table)
HA_seq_table <- HA_seq_table %>% tidyr::separate(season, into = c("year1", "year2"), remove = F)

range(HA_seq_table$north_america)
range(HA_seq_table$other_regions)
p <- ggplot(HA_seq_table) +
  geom_vline(aes(xintercept = as.Date(paste(year2, "01", "01", sep = "-"))), lty = "dashed", alpha = 0.2, linewidth = 0.5) +
  geom_hline(yintercept = 25, lty = "dashed", alpha = 0.5, linewidth = 0.5) +
  geom_line(aes(x = month.yr, y = north_america, color = "North America", group = replicate), lwd = 0.7, alpha = 0.4) +
  geom_line(aes(x = month.yr, y = other_regions, color = "Other world regions", group = replicate), lwd = 0.7, alpha = 0.4) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 10, limits = c(0, 50)) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 years", date_labels = "%Y") +
  theme_bw(base_size = 16) +
  ylab("Number of sequences") +
  xlab("Month (Collection Date)") +
  ggtitle("A/H3 sequences by region") +
  scale_color_jcolors(palette = "pal8", name = NULL) +
  theme(legend.position = c(0.8, 0.85), legend.box.background = element_blank(), legend.text = element_text(size = 12))
p

range(HA_seq_table$total_sequences)
q <- ggplot(HA_seq_table) +
  geom_vline(aes(xintercept = as.Date(paste(year2, "01", "01", sep = "-"))), lty = "dashed", alpha = 0.2, linewidth = 0.5) +
  geom_hline(yintercept = 50, lty = "dashed", alpha = 0.5, linewidth = 0.5) +
  geom_line(aes(x = month.yr, y = total_sequences, group = replicate), lwd = 0.7, alpha = 0.4) +
  # geom_line(aes(x=month.yr,y=north_america,color="North America"),lwd=1)+
  # geom_line(aes(x=month.yr,y=other_regions,color="Other regions"),lwd=1)+
  scale_y_continuous(expand = c(0, 0), n.breaks = 10, limits = c(0, 65)) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 years", date_labels = "%Y") +
  theme_bw(base_size = 16) +
  ylab("Number of sequences") +
  xlab("Month (Collection Date)") +
  ggtitle("Total A/H3 sequences")
q

h3_plot_by_month <- plot_grid(p, q, nrow = 2, labels = "AUTO")
h3_plot_by_month

HA_seq_table_season <- HA_north_amer_build %>%
  filter(!(season %in% c("1995-1996", "1996-1997", "2019-2020"))) %>%
  distinct() %>%
  group_by(replicate, season, region) %>%
  tally() %>%
  pivot_wider(names_from = region, values_from = n) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(replicate, season, north_america, africa:europe, oceania:west_asia) %>%
  arrange(season, replicate) %>%
  rowwise() %>%
  mutate(other_regions = sum(c_across(africa:west_asia))) %>%
  mutate(total_sequences = sum(c_across(north_america:west_asia))) %>%
  dplyr::select(replicate, season, total_sequences, north_america, other_regions, africa:west_asia)

HA_seq_table_season$replicate <- as.factor(HA_seq_table_season$replicate)
unique(HA_seq_table_season$replicate)
levels(HA_seq_table_season$replicate) <- c("0", "1", "2", "3", "4")
head(HA_seq_table_season)
HA_seq_table_season <- HA_seq_table_season %>% tidyr::separate(season, into = c("year1", "year2"), remove = F)
HA_seq_table_season$year2 <- as.numeric(HA_seq_table_season$year2)
HA_seq_table_season$year2 <- as.Date(paste(HA_seq_table_season$year2, "01", "01", sep = "-"))

p <- ggplot(HA_seq_table_season) +
  geom_vline(aes(xintercept = year2), lty = "dashed", alpha = 0.2, linewidth = 0.5) +
  # geom_hline(yintercept=25,lty="dashed",alpha=0.5,linewidth=0.5)+
  geom_line(aes(x = year2, y = north_america, color = "North America", group = replicate), lwd = 0.7, alpha = 0.4) +
  geom_line(aes(x = year2, y = other_regions, color = "Other world regions", group = replicate), lwd = 0.7, alpha = 0.4) +
  geom_point(aes(x = year2, y = north_america, fill = "North America", group = replicate), pch = 21, size = 3, alpha = 0.4) +
  geom_point(aes(x = year2, y = other_regions, fill = "Other world regions", group = replicate), pch = 21, size = 3, alpha = 0.4) +
  scale_y_continuous(expand = c(0.04, 0.04), n.breaks = 10) +
  scale_x_date(expand = c(0.04, 0.04), date_breaks = "2 years", date_labels = "%Y") +
  theme_bw(base_size = 16) +
  ylab("Number of sequences") +
  xlab("Season") +
  ggtitle("A/H3 sequences by region") +
  scale_color_jcolors(palette = "pal8", name = NULL) +
  scale_fill_jcolors(palette = "pal8", name = NULL) +
  theme(legend.position = c(0.8, 0.2), legend.box.background = element_blank(), legend.text = element_text(size = 12))
p

range(HA_seq_table_season$total_sequences)
q <- ggplot(HA_seq_table_season) +
  geom_vline(aes(xintercept = year2), lty = "dashed", alpha = 0.2, linewidth = 0.5) +
  # geom_hline(yintercept=25,lty="dashed",alpha=0.5,linewidth=0.5)+
  geom_line(aes(x = year2, y = total_sequences, group = replicate), lwd = 0.7, alpha = 0.4) +
  geom_point(aes(x = year2, y = total_sequences, group = replicate), pch = 21, fill = "black", size = 3, alpha = 0.4) +
  scale_y_continuous(expand = c(0.04, 0.04), n.breaks = 10, limits = c(0, 650)) +
  scale_x_date(expand = c(0.04, 0.04), date_breaks = "2 years", date_labels = "%Y") +
  theme_bw(base_size = 16) +
  ylab("Number of sequences") +
  xlab("Season") +
  ggtitle("Total A/H3 sequences") +
  scale_color_jcolors(palette = "pal8", name = NULL) +
  scale_fill_jcolors(palette = "pal8", name = NULL)
# theme(legend.position = c(0.9, 0.1), legend.box.background = element_blank())
q
h3_plot_by_season <- plot_grid(p, q, nrow = 2, labels = c("C", "D"))
h3_plot_by_season

h3_month_and_season_plots <- plot_grid(h3_plot_by_month, h3_plot_by_season, ncol = 2)
h3_month_and_season_plots
save_plot(h3_month_and_season_plots, file = "figures/h3_sequence_volume.png", base_width = 20, base_height = 8)
save_plot(h3_month_and_season_plots, file = "figures/h3_sequence_volume.pdf", dpi = 300, base_width = 20, base_height = 8)

summary_ha <- HA_seq_table_season %>%
  dplyr::select(replicate, season, total_sequences:west_asia) %>%
  dplyr::select(replicate, season, total_sequences:other_regions, china, southeast_asia, west_asia, japan_korea, south_asia, oceania, europe, south_america, africa)
summary_ha
summary_ha <- summary_ha %>%
  mutate(percent_north_america = round(100 * north_america / total_sequences, 2))
write_csv(summary_ha, file = "data/h3_sequence_summary_table_by_region_and_season.csv")

all_seasons_ha <- summary_ha %>%
  group_by(replicate) %>%
  summarize_at(vars(total_sequences:africa), ~ sum(.x)) %>%
  mutate(percent_north_america = round(100 * north_america / total_sequences, 2))
all_seasons_ha
write_csv(all_seasons_ha, file = "data/h3_sequence_summary_table_by_replicate_all_seasons_combined.csv")

##################################################
## N2 plots and summary table
##################################################
NA_north_amer_build <- readr::read_tsv("2_Phylo_Dataset/auspice_tables/flu_seasonal_h3n2_na_21y_north-america.tsv")
head(NA_north_amer_build)

NA_north_amer_build <- NA_north_amer_build %>%
  mutate(
    month.day = format(date, "%m-%d"),
    month = lubridate::month(date),
    year = lubridate::year(date),
    month.yr = as.Date(as.yearmon(date, "%b-%Y"))
  ) %>%
  mutate(season = ifelse(month.day < "07-01", sprintf("%d-%d", year - 1, year), # season = july 1 to june 30
    sprintf("%d-%d", year, year + 1)
  ))

unique(NA_north_amer_build$season)

NA_seq_table <- NA_north_amer_build %>%
  filter(!(season %in% c("1995-1996", "1996-1997", "2019-2020", "1994-1995"))) %>%
  distinct() %>%
  dplyr::group_by(replicate, season, month.yr, region) %>%
  tally() %>%
  pivot_wider(names_from = region, values_from = n) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(replicate, season, month.yr, north_america, south_america, europe, africa, oceania, china, southeast_asia, south_asia, west_asia, japan_korea) %>%
  arrange(season, month.yr, replicate) %>%
  rowwise() %>%
  mutate(other_regions = sum(c_across(south_america:japan_korea))) %>%
  mutate(total_sequences = north_america + other_regions) %>%
  dplyr::select(replicate, season, month.yr, total_sequences, north_america, other_regions, south_america:japan_korea)

NA_seq_table$replicate <- as.factor(NA_seq_table$replicate)
unique(NA_seq_table$replicate)
levels(NA_seq_table$replicate) <- c("0", "1", "2", "3", "4")
head(NA_seq_table)
NA_seq_table <- NA_seq_table %>% tidyr::separate(season, into = c("year1", "year2"), remove = F)


p <- ggplot(NA_seq_table) +
  geom_vline(aes(xintercept = as.Date(paste(year2, "01", "01", sep = "-"))), lty = "dashed", alpha = 0.2, linewidth = 0.5) +
  geom_hline(yintercept = 25, lty = "dashed", alpha = 0.5, linewidth = 0.5) +
  geom_line(aes(x = month.yr, y = north_america, color = "North America", group = replicate), lwd = 0.7, alpha = 0.4) +
  geom_line(aes(x = month.yr, y = other_regions, color = "Other world regions", group = replicate), lwd = 0.7, alpha = 0.4) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 10, limits = c(0, 50)) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 years", date_labels = "%Y") +
  theme_bw(base_size = 16) +
  ylab("Number of sequences") +
  xlab("Month (Collection Date)") +
  ggtitle("A/N2 sequences by region") +
  scale_color_jcolors(palette = "pal8", name = NULL) +
  theme(legend.position = c(0.8, 0.85), legend.box.background = element_blank(), legend.text = element_text(size = 12))
p

range(NA_seq_table$total_sequences)
q <- ggplot(NA_seq_table) +
  geom_vline(aes(xintercept = as.Date(paste(year2, "01", "01", sep = "-"))), lty = "dashed", alpha = 0.2, linewidth = 0.5) +
  geom_hline(yintercept = 50, lty = "dashed", alpha = 0.5, linewidth = 0.5) +
  geom_line(aes(x = month.yr, y = total_sequences, group = replicate), lwd = 0.7, alpha = 0.4) +
  # geom_line(aes(x=month.yr,y=north_america,color="North America"),lwd=1)+
  # geom_line(aes(x=month.yr,y=other_regions,color="Other regions"),lwd=1)+
  scale_y_continuous(expand = c(0, 0), n.breaks = 10, limits = c(0, 65)) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 years", date_labels = "%Y") +
  theme_bw(base_size = 16) +
  ylab("Number of sequences") +
  xlab("Month (Collection Date)") +
  ggtitle("Total A/N2 sequences")
q

n2_plot_by_month <- plot_grid(p, q, nrow = 2, labels = "AUTO")
n2_plot_by_month

NA_seq_table_season <- NA_north_amer_build %>%
  filter(!(season %in% c("1995-1996", "1996-1997", "2019-2020", "1994-1995"))) %>%
  distinct() %>%
  group_by(replicate, season, region) %>%
  tally() %>%
  pivot_wider(names_from = region, values_from = n) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(replicate, season, north_america, china:europe, oceania:south_asia) %>%
  arrange(season, replicate) %>%
  rowwise() %>%
  mutate(other_regions = sum(c_across(china:south_asia))) %>%
  mutate(total_sequences = sum(c_across(north_america:south_asia))) %>%
  dplyr::select(replicate, season, total_sequences, north_america, other_regions, china:south_asia)

NA_seq_table_season$replicate <- as.factor(NA_seq_table_season$replicate)
unique(NA_seq_table_season$replicate)
levels(NA_seq_table_season$replicate) <- c("0", "1", "2", "3", "4")
head(NA_seq_table_season)
NA_seq_table_season <- NA_seq_table_season %>% tidyr::separate(season, into = c("year1", "year2"), remove = F)
NA_seq_table_season$year2 <- as.numeric(NA_seq_table_season$year2)
NA_seq_table_season$year2 <- as.Date(paste(NA_seq_table_season$year2, "01", "01", sep = "-"))

p <- ggplot(NA_seq_table_season) +
  geom_vline(aes(xintercept = year2), lty = "dashed", alpha = 0.2, linewidth = 0.5) +
  # geom_hline(yintercept=25,lty="dashed",alpha=0.5,linewidth=0.5)+
  geom_line(aes(x = year2, y = north_america, color = "North America", group = replicate), lwd = 0.7, alpha = 0.4) +
  geom_line(aes(x = year2, y = other_regions, color = "Other world regions", group = replicate), lwd = 0.7, alpha = 0.4) +
  geom_point(aes(x = year2, y = north_america, fill = "North America", group = replicate), pch = 21, size = 3, alpha = 0.4) +
  geom_point(aes(x = year2, y = other_regions, fill = "Other world regions", group = replicate), pch = 21, size = 3, alpha = 0.4) +
  scale_y_continuous(expand = c(0.04, 0.04), n.breaks = 10) +
  scale_x_date(expand = c(0.04, 0.04), date_breaks = "2 years", date_labels = "%Y") +
  theme_bw(base_size = 16) +
  ylab("Number of sequences") +
  xlab("Season") +
  ggtitle("A/N2 sequences by region") +
  scale_color_jcolors(palette = "pal8", name = NULL) +
  scale_fill_jcolors(palette = "pal8", name = NULL) +
  theme(legend.position = c(0.8, 0.2), legend.box.background = element_blank(), legend.text = element_text(size = 12))
p

range(NA_seq_table_season$total_sequences)
q <- ggplot(NA_seq_table_season) +
  geom_vline(aes(xintercept = year2), lty = "dashed", alpha = 0.2, linewidth = 0.5) +
  # geom_hline(yintercept=25,lty="dashed",alpha=0.5,linewidth=0.5)+
  geom_line(aes(x = year2, y = total_sequences, group = replicate), lwd = 0.7, alpha = 0.4) +
  geom_point(aes(x = year2, y = total_sequences, group = replicate), pch = 21, fill = "black", size = 3, alpha = 0.4) +
  scale_y_continuous(expand = c(0.04, 0.04), n.breaks = 10, limits = c(0, 650)) +
  scale_x_date(expand = c(0.04, 0.04), date_breaks = "2 years", date_labels = "%Y") +
  theme_bw(base_size = 16) +
  ylab("Number of sequences") +
  xlab("Season") +
  ggtitle("Total A/N2 sequences") +
  scale_color_jcolors(palette = "pal8", name = NULL) +
  scale_fill_jcolors(palette = "pal8", name = NULL) +
  theme(legend.position = c(0.9, 0.1), legend.box.background = element_blank())
q
n2_plot_by_season <- plot_grid(p, q, nrow = 2, labels = c("C", "D"))
n2_plot_by_season

n2_month_and_season_plots <- plot_grid(n2_plot_by_month, n2_plot_by_season, ncol = 2)
save_plot(n2_month_and_season_plots, file = "figures/n2_sequence_volume.png", base_width = 20, base_height = 8)
save_plot(n2_month_and_season_plots, file = "figures/n2_sequence_volume.pdf", dpi = 300, base_width = 20, base_height = 8)

summary_na <- NA_seq_table_season %>%
  dplyr::select(replicate, season, total_sequences:south_asia) %>%
  dplyr::select(replicate, season, total_sequences:other_regions, china, southeast_asia, west_asia, japan_korea, south_asia, oceania, europe, south_america, africa)
summary_na <- summary_na %>%
  mutate(percent_north_america = round(100 * north_america / total_sequences, 2))
write_csv(summary_na, file = "data/n2_sequence_summary_table_by_region_and_season.csv")

all_seasons_na <- summary_na %>%
  group_by(replicate) %>%
  summarize_at(vars(total_sequences:africa), ~ sum(.x)) %>%
  mutate(percent_north_america = round(100 * north_america / total_sequences, 2))
write_csv(all_seasons_na, file = "data/n2_sequence_summary_table_by_replicate_all_seasons_combined.csv")
