## load packages
list.of.packages <- c("dplyr", "ggplot2", "readr", "tidyr", "vegan", "lubridate", "lsr", "iNEXT","RColorBrewer")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

############################################################################
## Calculate seasonal Shannon entropy of LBI values
############################################################################

############################################################################
## HA
############################################################################

lbi_tb_ha <- read_tsv("2_Phylo_Dataset/distance_tables/north-america_strain_seasonal_lbi_h3n2_ha_21y.tsv.gz")
lbi_tb_ha$year <- lubridate::year(lbi_tb_ha$season_end)
lbi_tb_ha <- lbi_tb_ha %>% filter(year < 2020 & year > 1996)
range(lbi_tb_ha$lbi)
# breaks <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
# length(breaks)
# tags <- c("(0-2]", "(2-4]", "(4-6]", "(6-8]", "(8-10]", "(10-12]", "(12-14]", "(14-16]", "(16-18]", "(18-20]")
# length(tags)

breaks <- seq(from = 0, to = 19, by = 1)
breaks
length(breaks)
tags <- c(
  "(0-1]", "(1-2]", "(2-3]", "(3-4]", "(4-5]", "(5-6]", "(6-7]", "(7-8]","(8-9]", "(9-10]", "(10-11]",
  "(11-12]", "(12-13]", "(13-14]", "(14-15]", "(15-16]", "(16-17]", "(17-18]", "(18-19]"
)
length(tags)

lbi_tb_ha$lbi_bin <- cut(lbi_tb_ha$lbi,
  breaks = breaks,
  include.lowest = TRUE,
  right = T,
  labels = tags
)
table(lbi_tb_ha$lbi_bin)

lbi_tb_ha_sum <- lbi_tb_ha %>%
  group_by(year, replicate, lbi_bin) %>%
  tally() %>%
  ungroup()

# Define the number of colors you want
nb.cols <- 19
cols <- colorRampPalette(brewer.pal(9, "Spectral"))(nb.cols)

lbi_tb_ha_sum$replicate = as.factor(lbi_tb_ha_sum$replicate)
levels(lbi_tb_ha_sum$replicate) <- paste("Replicate",unique(lbi_tb_ha_sum$replicate),sep = " ")


ha_bins = ggplot(data = lbi_tb_ha_sum, aes(x = year, fill = lbi_bin)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols,name="LBI Bin") +
  scale_x_continuous(breaks = seq(from = 1997, to = 2020, 4), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~replicate)+
  theme(strip.background = element_blank())+
  xlab("season")+
  ylab("Proportion")+
  ggtitle("A/H3 LBI Bins by Season")
ha_bins

sp_tb <- lbi_tb_ha_sum %>%
  mutate(rep_yr = paste(year, replicate, sep = "-")) %>%
  dplyr::select(-year, -replicate) %>%
  pivot_wider(names_from = rep_yr, values_from = n) %>%
  as.data.frame()

rownames(sp_tb) <- sp_tb$lbi_bin
sp_tb <- sp_tb %>% dplyr::select(-lbi_bin)
sp_tb[is.na(sp_tb)] <- 0
sp_tb %>% tail()

#iNEXT package: iNterpolation and EXTrapolation of Hill number with order q
# q = 1 is Shannon diversity index (geometric mean of proportional abundance)
out <- iNEXT(sp_tb, q = 1, datatype = "abundance")
plot(out)

div_df <- out$iNextEst$size_based
div_df <- div_df %>% tidyr::separate(Assemblage, into = c("year", "replicate"), sep = "-", remove = F)
div_df

# max # sequences for extrapolated samples in each season-replicate pair
div_df %>%
  filter(SC>0.99 & m >=344)%>%
  group_by(year) %>%
  slice_min(m)%>%
  arrange(year)%>%
  dplyr::select(year,m)%>%
  distinct()%>%
  print(n=30)

div_df %>%
  group_by(replicate, year) %>%
  slice_max(m)%>%
  arrange(year)

div_df %>%
  filter(Method!="Extrapolation")%>%
  group_by(replicate, year) %>%
  slice_max(m)%>%
  arrange(year)

max_m = div_df %>%
  filter(year!="1997")%>%
  group_by(replicate, year) %>%
  slice_max(m) %>%
  pull(m)
min(max_m) #344

ggplot(div_df) +
  geom_line(aes(x = m, y = qD, group = year, color = year)) +
  geom_vline(xintercept = 88,lty="dashed")+
  geom_vline(xintercept = 344,lty="dashed")+
  facet_wrap(~replicate)

div_df %>%
  filter(SC>0.995 & m >=344)%>%
  group_by(Assemblage,replicate, year) %>%
  slice_min(m)%>%
  print(n=40)

ha_lbi_shannon_rarefied <- div_df %>%
  filter(m <= 360) %>%
  # filter(m <= min(max_m)) %>%
  group_by(replicate, year) %>%
  slice_max(m) %>%
  dplyr::select(year, replicate,m,qD,SC) %>%
  rename(ha_lbi_shannon = qD) %>%
  ungroup()

ha_lbi_shannon_rarefied %>% dplyr::select(replicate,year,m,ha_lbi_shannon) %>% distinct() %>% print(n=55)

ggplot(ha_lbi_shannon_rarefied) +
  geom_line(aes(x = year, y = ha_lbi_shannon, group = replicate, color = replicate))

############################################################################
## NA
############################################################################

lbi_tb_na <- read_tsv("2_Phylo_Dataset/distance_tables/north-america_strain_seasonal_lbi_h3n2_na_21y.tsv.gz")
lbi_tb_na$year <- lubridate::year(lbi_tb_na$season_end)
lbi_tb_na <- lbi_tb_na %>% filter(year < 2020 & year > 1996)
range(lbi_tb_na$lbi)

# breaks <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18)
# tags <- c("(0-2]", "(2-4]", "(4-6]", "(6-8]", "(8-10]", "(10-12]", "(12-14]", "(14-16]", "(16-18]")
# length(tags)

breaks <- seq(from = 0, to = 17, by = 1)
breaks
length(breaks)
tags <- c(
  "(0-1]", "(1-2]", "(2-3]", "(3-4]", "(4-5]", "(5-6]", "(6-7]", "(7-8]","(8-9]", "(9-10]", "(10-11]",
  "(11-12]", "(12-13]", "(13-14]", "(14-15]", "(15-16]", "(16-17]"
)
length(tags)

lbi_tb_na$lbi_bin <- cut(lbi_tb_na$lbi,
  breaks = breaks,
  include.lowest = TRUE,
  right = T,
  labels = tags
)
table(lbi_tb_na$lbi_bin)

lbi_tb_na_sum <- lbi_tb_na %>%
  group_by(year, replicate, lbi_bin) %>%
  tally() %>%
  ungroup()

lbi_tb_na_sum$replicate = as.factor(lbi_tb_na_sum$replicate)
levels(lbi_tb_na_sum$replicate) <- paste("Replicate",unique(lbi_tb_na_sum$replicate),sep = " ")

## higher LBI values in later seasons = beginning in 2012 more large, dense clades compared to early seasons
na_bins = ggplot(data = lbi_tb_na_sum, aes(x = year, fill = lbi_bin)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols,name="LBI Bin") +
  scale_x_continuous(breaks = seq(from = 1997, to = 2020, 4), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~replicate)+
  theme(strip.background = element_blank())+
  xlab("season")+
  ylab("Proportion")+
  ggtitle("A/N2 LBI Bins by Season")
na_bins

sp_tb <- lbi_tb_na_sum %>%
  mutate(rep_yr = paste(year, replicate, sep = "-")) %>%
  dplyr::select(-year, -replicate) %>%
  pivot_wider(names_from = rep_yr, values_from = n) %>%
  as.data.frame()
rownames(sp_tb) <- sp_tb$lbi_bin
sp_tb <- sp_tb %>% dplyr::select(-lbi_bin)
sp_tb[is.na(sp_tb)] <- 0

#iNEXT package: iNterpolation and EXTrapolation of Hill number with order q
# q = 1 is Shannon diversity index (geometric mean of proportional abundance)
out <- iNEXT(sp_tb, q = 1, datatype = "abundance")
plot(out)

div_df <- out$iNextEst$size_based
div_df <- div_df %>% tidyr::separate(Assemblage, into = c("year", "replicate"), sep = "-", remove = F)

ggplot(div_df) +
  geom_line(aes(x = m, y = qD, group = year,color=year)) +
  facet_wrap(~replicate)

div_df %>%
  group_by(replicate, year) %>%
  slice_max(m)

max_m = div_df %>%
  filter(year>1997)%>%
  group_by(replicate, year) %>%
  slice_max(m) %>%
  pull(m)
max_m
sort(max_m)
min(max_m) #216

na_lbi_shannon_rarefied <- div_df %>%
  filter(m <= 230) %>%
  # filter(m <= min(max_m)) %>%
  group_by(replicate, year) %>%
  slice_max(m) %>%
  dplyr::select(year, replicate,m, qD) %>%
  rename(na_lbi_shannon = qD) %>%
  ungroup()

na_lbi_shannon_rarefied %>%
  distinct(year,m)%>%
  print(n=55)

ggplot(na_lbi_shannon_rarefied) +
  geom_line(aes(x = year, y = na_lbi_shannon, group = replicate, color = replicate))

lbi_eco_div <- full_join(
  na_lbi_shannon_rarefied %>% dplyr::select(year,replicate,na_lbi_shannon),
  ha_lbi_shannon_rarefied %>% dplyr::select(year,replicate,ha_lbi_shannon)
)

lbi_eco_div <- lbi_eco_div %>%
  mutate(year1 = as.numeric(year) - 1) %>%
  mutate(season = paste(year1, year, sep = "-"))

ggplot(lbi_eco_div) +
  geom_boxplot(aes(x = as.numeric(year), y = ha_lbi_shannon, group = season))
ggplot(lbi_eco_div) +
  geom_boxplot(aes(x = as.numeric(year), y = na_lbi_shannon, group = season))

write_rds(lbi_eco_div, file = "data/LBI_diversity_by_season.rds")

comb_plot = plot_grid(ha_bins,na_bins,nrow=2,labels="AUTO")
comb_plot
save_plot(comb_plot,filename="figures/H3_and_N2_LBI_Bins_by_Season.pdf",base_width = 16, base_height = 16, dpi = 300)
