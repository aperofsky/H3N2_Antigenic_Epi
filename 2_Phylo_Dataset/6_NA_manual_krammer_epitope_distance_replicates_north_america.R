##################################################
## manually calculate NA epitope distances (Krammer mask)
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
## NA epitope sites
##################################################

## krammer mask
sites <- c(
  147, 150, 152, 153, 154, 155, 197, 198, 199, 218, 219, 220, 221, 222, 224, 245, 246, 247, 248, 249,
  251, 253, 328, 329, 330, 331, 332, 333, 334, 335, 336, 339, 340, 341, 342, 343, 344, 345, 346, 347,
  367, 368, 369, 370, 400, 401, 402, 403, 429, 431, 432, 433, 434
)

sites <- sort(sites)
length(sites) # 53
vec <- rep(0, 469)
vec <- paste(vec, collapse = "")

stri_sub_all(vec, from = sites, length = 1) <- "1"
vec

NA_mask <- vec
##################################################
## NA epitope sites - North American build
##################################################
## calculate direct epitope distance just for 1997-1998 (missing from John's dataset)

NA_north_amer_build <- readr::read_tsv("2_Phylo_Dataset//auspice_tables/flu_seasonal_h3n2_na_21y_north-america.tsv")
head(NA_north_amer_build)
NA_north_amer_build <- NA_north_amer_build %>%
  mutate(
    month.day = format(date, "%m-%d"),
    year = as.numeric(format(date, "%Y"))
  ) %>%
  mutate(season = ifelse(month.day < "07-01", sprintf("%d-%d", year - 1, year), # season = july 1 to june 30
    sprintf("%d-%d", year, year + 1)
  ))
sort(unique(NA_north_amer_build$season))
head(NA_north_amer_build)
nchar(NA_north_amer_build$`NA`)[1]
NA_north_amer_build$seq <- NA_north_amer_build$`NA`
NA_north_amer_build$seq[grepl("-", NA_north_amer_build$seq)]

NA_north_amer_build %>%
  group_by(season) %>%
  summarize(mean(ep))

replicates <- unique(NA_north_amer_build$replicate)
datalist <- list()
for (i in replicates) {
  NA_north_amer_build_red <- NA_north_amer_build %>% filter(season %in% c("1995-1996", "1996-1997", "1997-1998") & replicate == i)

  NA_substrings <- stri_sub_all(NA_north_amer_build_red$seq, from = sites, length = 1)

  collapse_fun <- function(NA_list) paste0(NA_list, collapse = "")

  NA_substrings <- lapply(NA_substrings, FUN = collapse_fun)

  NA_substrings[grepl("-", NA_substrings)]
  NA_substrings_df <- data.frame(matrix(unlist(NA_substrings), nrow = length(NA_substrings), byrow = T), stringsAsFactors = FALSE)
  names(NA_substrings_df) <- "NA_epitope_sites"
  NA_north_amer_build2 <- bind_cols(NA_north_amer_build_red, NA_substrings_df)

  ### putative NA epitope sites
  NA_north_amer_AA_epitope <- AAStringSet(x = NA_north_amer_build2$NA_epitope_sites, start = NA, end = NA, width = NA, use.names = TRUE)
  names(NA_north_amer_AA_epitope) <- paste(NA_north_amer_build2$name, NA_north_amer_build2$date, sep = "|")

  NA_north_amer_dist_epitope <- stringDist(NA_north_amer_AA_epitope, method = "hamming")

  NA_north_amer_dist_epitope_EL <- reshape2::melt(as.matrix(NA_north_amer_dist_epitope)) %>% filter(Var1 != Var2)
  NA_north_amer_dist_epitope_EL$Var1 <- as.character(NA_north_amer_dist_epitope_EL$Var1)
  NA_north_amer_dist_epitope_EL$Var2 <- as.character(NA_north_amer_dist_epitope_EL$Var2)

  NA_north_amer_dist_epitope_EL.raw <- as.data.frame(NA_north_amer_dist_epitope_EL)
  NA_north_amer_dist_epitope_EL.raw <- as.data.table(NA_north_amer_dist_epitope_EL)
  NA_north_amer_dist_epitope_EL.raw[, c("Var1_1", "Var1_2") := tstrsplit(Var1, "[|]", fixed = F)]
  NA_north_amer_dist_epitope_EL.raw[, c("Var2_1", "Var2_2") := tstrsplit(Var2, "[|]", fixed = F)]

  split_df2 <- NA_north_amer_dist_epitope_EL.raw %>%
    tidytable::mutate(
      Var1_date = as.Date(Var1_2),
      Var2_date = as.Date(Var2_2)
    ) %>%
    tidytable::mutate(
      Var1.month.day = format(Var1_date, "%m-%d"),
      Var1.year = as.numeric(format(Var1_date, "%Y")),
      Var2.month.day = format(Var2_date, "%m-%d"),
      Var2.year = as.numeric(format(Var2_date, "%Y"))
    ) %>%
    tidytable::mutate(
      Var1.season = ifelse(Var1.month.day < "07-01", sprintf("%d-%d", Var1.year - 1, Var1.year),
        # season = july 1 to june 30
        sprintf("%d-%d", Var1.year, Var1.year + 1)
      ),
      Var2.season = ifelse(Var2.month.day < "07-01", sprintf("%d-%d", Var2.year - 1, Var2.year),
        # season = july 1 to june 30
        sprintf("%d-%d", Var2.year, Var2.year + 1)
      )
    ) %>%
    tidytable::mutate(
      Var1.year.new = as.numeric(substr(Var1.season, 6, 9)),
      Var2.year.new = as.numeric(substr(Var2.season, 6, 9))
    ) %>%
    tidytable::mutate(season_diff = abs(Var2.year.new - Var1.year.new)) %>%
    tidytable::filter(season_diff %in% c(0, 1, 2))

  split_df2 <- as.data.table(split_df2)
  split_df2[, year_diff := paste(sort(c(Var1.year.new, Var2.year.new)), collapse = "-"), by = 1:nrow(split_df2)]


  NA_north_amer_dist_epitope_summary <-
    split_df2 %>%
    group_by(year_diff, season_diff) %>%
    summarize(
      NA.ep.global.mean = mean(value),
      NA.ep.global.sd = sd(value),
      comparisons = n()
    ) %>%
    ungroup()
  NA_north_amer_dist_epitope_summary$replicate <- i
  datalist[[i]] <- NA_north_amer_dist_epitope_summary
}
df <- do.call("rbind", datalist)
head(df)
NA_north_amer_dist_epitope_summary <- df %>% distinct()

NA_north_amer_dist_epitope_summary <- NA_north_amer_dist_epitope_summary %>% tidyr::separate(year_diff, into = c("year1", "year2"), remove = F)
NA_north_amer_dist_epitope_summary$season <- paste(as.numeric(NA_north_amer_dist_epitope_summary$year2) - 1, NA_north_amer_dist_epitope_summary$year2, sep = "-")

NA_north_amer_dist_epitope_summary %>%
  arrange(season, season_diff) %>%
  filter(season != "1995-1996")
head(NA_north_amer_dist_epitope_summary)

save(NA_north_amer_dist_epitope_summary, file = "data/NA_mean_krammer_epitope_change_north_amer_build_summary.Rdata")
