##################################################
## manually calculate HA Koel epitope distances
##################################################

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("Biostrings")

## load packages
list.of.packages <- c("dplyr", "stringi", "Biostrings", "data.table", "tidytable", "reshape2", "tidyr", "readr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

##################################################
## HA epitope sites
##################################################
sites <- c(145, 155, 156, 158, 159, 189, 193)

sites <- sort(sites)
length(sites) # 7
vec <- rep(0, 329)
vec <- paste(vec, collapse = "")

stri_sub_all(vec, from = sites, length = 1) <- "1"
vec
HA_mask <- vec
##################################################
## HA epitope sites - North American build
##################################################
## calculate direct epitope distance just for 1997-1998 (missing from John's dataset)

HA_north_amer_build <- readr::read_tsv("2_Phylo_Dataset/auspice_tables/flu_seasonal_h3n2_ha_21y_north-america.tsv")
head(HA_north_amer_build)
HA_north_amer_build <- HA_north_amer_build %>%
  mutate(
    month.day = format(date, "%m-%d"),
    year = as.numeric(format(date, "%Y"))
  ) %>%
  mutate(season = ifelse(month.day < "07-01", sprintf("%d-%d", year - 1, year), # season = july 1 to june 30
    sprintf("%d-%d", year, year + 1)
  ))
sort(unique(HA_north_amer_build$season))
head(HA_north_amer_build)
nchar(HA_north_amer_build$HA1)[1]
HA_north_amer_build$seq <- HA_north_amer_build$HA1
HA_north_amer_build$seq[grepl("-", HA_north_amer_build$seq)]

HA_north_amer_build %>%
  group_by(season) %>%
  summarize(mean(ep))

replicates <- unique(HA_north_amer_build$replicate)
datalist <- list()
for (i in replicates) {
  HA_north_amer_build_red <- HA_north_amer_build %>%
    filter(season %in% c("1995-1996", "1996-1997", "1997-1998") & replicate == i)

  HA_substrings <- stri_sub_all(HA_north_amer_build_red$seq, from = sites, length = 1)

  collapse_fun <- function(HA_list) paste0(HA_list, collapse = "")

  HA_substrings <- lapply(HA_substrings, FUN = collapse_fun)

  HA_substrings[grepl("-", HA_substrings)]
  HA_substrings_df <- data.frame(matrix(unlist(HA_substrings), nrow = length(HA_substrings), byrow = T), stringsAsFactors = FALSE)
  names(HA_substrings_df) <- "HA_epitope_sites"
  HA_north_amer_build2 <- bind_cols(HA_north_amer_build_red, HA_substrings_df)


  ### putative HA epitope sites
  HA_north_amer_AA_epitope <- AAStringSet(x = HA_north_amer_build2$HA_epitope_sites, start = NA, end = NA, width = NA, use.names = TRUE)
  names(HA_north_amer_AA_epitope) <- paste(HA_north_amer_build2$name, HA_north_amer_build2$date, sep = "|")

  HA_north_amer_dist_epitope <- stringDist(HA_north_amer_AA_epitope, method = "hamming")

  HA_north_amer_dist_epitope_EL <- reshape2::melt(as.matrix(HA_north_amer_dist_epitope)) %>% filter(Var1 != Var2)
  HA_north_amer_dist_epitope_EL$Var1 <- as.character(HA_north_amer_dist_epitope_EL$Var1)
  HA_north_amer_dist_epitope_EL$Var2 <- as.character(HA_north_amer_dist_epitope_EL$Var2)

  HA_north_amer_dist_epitope_EL.raw <- as.data.table(HA_north_amer_dist_epitope_EL)
  HA_north_amer_dist_epitope_EL.raw[, c("Var1_1", "Var1_2") := tstrsplit(Var1, "[|]", fixed = F)]
  HA_north_amer_dist_epitope_EL.raw[, c("Var2_1", "Var2_2") := tstrsplit(Var2, "[|]", fixed = F)]

  split_df2 <- HA_north_amer_dist_epitope_EL.raw %>%
    mutate.(
      Var1_date = as.Date(Var1_2),
      Var2_date = as.Date(Var2_2)
    ) %>%
    mutate.(
      Var1.month.day = format(Var1_date, "%m-%d"),
      Var1.year = as.numeric(format(Var1_date, "%Y")),
      Var2.month.day = format(Var2_date, "%m-%d"),
      Var2.year = as.numeric(format(Var2_date, "%Y"))
    ) %>%
    mutate.(
      Var1.season = ifelse(Var1.month.day < "07-01", sprintf("%d-%d", Var1.year - 1, Var1.year),
        # season = july 1 to june 30
        sprintf("%d-%d", Var1.year, Var1.year + 1)
      ),
      Var2.season = ifelse(Var2.month.day < "07-01", sprintf("%d-%d", Var2.year - 1, Var2.year),
        # season = july 1 to june 30
        sprintf("%d-%d", Var2.year, Var2.year + 1)
      )
    ) %>%
    mutate.(
      Var1.year.new = as.numeric(substr(Var1.season, 6, 9)),
      Var2.year.new = as.numeric(substr(Var2.season, 6, 9))
    ) %>%
    mutate.(season_diff = abs(Var2.year.new - Var1.year.new)) %>%
    tidytable::filter(season_diff %in% c(0, 1, 2))

  split_df2 <- as.data.table(split_df2)
  split_df2[, year_diff := paste(sort(c(Var1.year.new, Var2.year.new)), collapse = "-"), by = 1:nrow(split_df2)]

  HA_north_amer_dist_epitope_summary <-
    split_df2 %>%
    group_by(year_diff, season_diff) %>%
    summarize(
      HA.ep.global.mean = mean(value),
      HA.ep.global.sd = sd(value),
      comparisons = n()
    ) %>%
    ungroup()
  HA_north_amer_dist_epitope_summary$replicate <- i
  print(as.data.frame(HA_north_amer_dist_epitope_summary))
  datalist[[i]] <- HA_north_amer_dist_epitope_summary
}
df <- do.call("rbind", datalist)
HA_north_amer_dist_epitope_summary <- df %>% distinct()
HA_north_amer_dist_epitope_summary <- HA_north_amer_dist_epitope_summary %>% tidyr::separate(year_diff, into = c("year1", "year2"), remove = F)
HA_north_amer_dist_epitope_summary$season <- paste(as.numeric(HA_north_amer_dist_epitope_summary$year2) - 1, HA_north_amer_dist_epitope_summary$year2, sep = "-")
HA_north_amer_dist_epitope_summary %>%
  arrange(season, season_diff) %>%
  filter(season != "1995-1996")
save(HA_north_amer_dist_epitope_summary, file = "data/HA_mean_koel_change_north_amer_build_summary.Rdata")
