##################################################
## manually calculate HA non-epitope distances
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
## HA1
sites1 <- c(
  1,
  10,
  100,
  101,
  104,
  105,
  106,
  107,
  108,
  11,
  110,
  111,
  112,
  113,
  114,
  115,
  116,
  118,
  119,
  12,
  120,
  123,
  125,
  127,
  13,
  134,
  136,
  139,
  14,
  141,
  147,
  148,
  149,
  15,
  151,
  153,
  154,
  16,
  160,
  161,
  162,
  164,
  166,
  169,
  17,
  178,
  18,
  180,
  181,
  183,
  184,
  185,
  19,
  191,
  195,
  199,
  2,
  20,
  200,
  202,
  204,
  205,
  206,
  21,
  210,
  211,
  22,
  220,
  221,
  222,
  223,
  224,
  225,
  23,
  231,
  232,
  233,
  234,
  235,
  236,
  237,
  239,
  24,
  241,
  243,
  245,
  249,
  25,
  250,
  251,
  252,
  253,
  254,
  255,
  256,
  257,
  258,
  259,
  26,
  263,
  264,
  266,
  267,
  268,
  269,
  27,
  270,
  271,
  272,
  274,
  277,
  28,
  281,
  282,
  283,
  284,
  285,
  286,
  287,
  288,
  289,
  29,
  290,
  291,
  292,
  293,
  295,
  296,
  298,
  3,
  30,
  301,
  302,
  303,
  306,
  31,
  313,
  314,
  315,
  316,
  317,
  318,
  319,
  32,
  320,
  321,
  322,
  323,
  324,
  325,
  326,
  327,
  328,
  329,
  33,
  34,
  35,
  36,
  37,
  38,
  39,
  4,
  40,
  41,
  42,
  43,
  49,
  5,
  52,
  55,
  56,
  58,
  6,
  60,
  61,
  64,
  65,
  66,
  68,
  69,
  7,
  70,
  71,
  72,
  73,
  74,
  76,
  77,
  79,
  8,
  84,
  85,
  89,
  9,
  90,
  93,
  95,
  97,
  98,
  99
)

sites1 <- sort(sites1)
length(sites1) # 200
vec1 <- rep(0, 329)
vec1 <- paste(vec1, collapse = "")

stri_sub_all(vec1, from = sites1, length = 1) <- "1"
vec1

##################################################
## HA2
sites2 <- c(
  1,
  10,
  100,
  101,
  102,
  103,
  104,
  105,
  106,
  107,
  108,
  109,
  11,
  110,
  111,
  112,
  113,
  114,
  115,
  116,
  117,
  118,
  119,
  12,
  120,
  121,
  122,
  123,
  124,
  125,
  126,
  127,
  128,
  129,
  13,
  130,
  131,
  132,
  133,
  134,
  135,
  136,
  137,
  138,
  139,
  14,
  140,
  141,
  142,
  143,
  144,
  145,
  146,
  147,
  148,
  149,
  15,
  150,
  151,
  152,
  153,
  154,
  155,
  156,
  157,
  158,
  159,
  16,
  160,
  161,
  162,
  163,
  164,
  165,
  166,
  167,
  168,
  169,
  17,
  170,
  171,
  172,
  173,
  174,
  175,
  176,
  177,
  178,
  179,
  18,
  180,
  181,
  182,
  183,
  184,
  185,
  186,
  187,
  188,
  189,
  19,
  190,
  191,
  192,
  193,
  194,
  195,
  196,
  197,
  198,
  199,
  2,
  20,
  200,
  201,
  202,
  203,
  204,
  205,
  206,
  207,
  208,
  209,
  21,
  210,
  211,
  212,
  213,
  214,
  215,
  216,
  217,
  218,
  219,
  22,
  220,
  221,
  23,
  24,
  25,
  26,
  27,
  28,
  29,
  3,
  30,
  31,
  32,
  33,
  34,
  35,
  36,
  37,
  38,
  39,
  4,
  40,
  41,
  42,
  43,
  44,
  45,
  46,
  47,
  48,
  49,
  5,
  50,
  51,
  52,
  53,
  54,
  55,
  56,
  57,
  58,
  59,
  6,
  60,
  61,
  62,
  63,
  64,
  65,
  66,
  67,
  68,
  69,
  7,
  70,
  71,
  72,
  73,
  74,
  75,
  76,
  77,
  78,
  79,
  8,
  80,
  81,
  82,
  83,
  84,
  85,
  86,
  87,
  88,
  89,
  9,
  90,
  91,
  92,
  93,
  94,
  95,
  96,
  97,
  98,
  99
)

sites2 <- sort(sites2)
length(sites2)
vec2 <- rep(0, 221)
vec2 <- paste(vec2, collapse = "")

stri_sub_all(vec2, from = sites2, length = 1) <- "1"
vec2

##################################################
## SigPep
sites3 <- c(
  1,
  10,
  11,
  12,
  13,
  14,
  15,
  16,
  2,
  3,
  4,
  5,
  6,
  7,
  8,
  9
)

sites3 <- sort(sites3)
length(sites3)
vec3 <- rep(0, 16)
vec3 <- paste(vec3, collapse = "")

stri_sub_all(vec3, from = sites3, length = 1) <- "1"
vec3

##################################################
## combined HA1 and HA2
vec <- paste0(vec3, vec1, vec2)
HA_mask <- vec
sites <- unlist(gregexpr(pattern = "1", HA_mask))
##################################################
## HA epitope sites - north_amer build
##################################################
## calculate direct epitope distances for all seasons

HA_north_amer_build <- readr::read_tsv("2_Phylo_Dataset/auspice_tables/flu_seasonal_h3n2_ha_21y_north-america.tsv")
unique(nchar(HA_north_amer_build$HA1)) # 329
unique(nchar(HA_north_amer_build$HA2)) # 221
unique(nchar(HA_north_amer_build$SigPep)) # 16

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
HA_north_amer_build$seq <- paste0(HA_north_amer_build$SigPep, HA_north_amer_build$HA1, HA_north_amer_build$HA2)
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
    mutate.(season_diff = abs(Var2.year.new - Var1.year.new)) %>%
    tidytable::filter(season_diff %in% c(0, 1, 2))

  split_df2 <- as.data.table(split_df2)
  system.time(split_df2[, year_diff := paste(sort(c(Var1.year.new, Var2.year.new)), collapse = "-"), by = 1:nrow(split_df2)])

  HA_north_amer_dist_epitope_summary <-
    split_df2 %>%
    group_by(year_diff, season_diff) %>%
    summarize(
      HA.ep.north_amer.mean = mean(value),
      HA.ep.north_amer.sd = sd(value),
      comparisons = n()
    ) %>%
    ungroup()
  HA_north_amer_dist_epitope_summary$replicate <- i
  datalist[[i]] <- HA_north_amer_dist_epitope_summary
}
df <- do.call("rbind", datalist)
HA_north_amer_dist_epitope_summary <- df %>% distinct()
HA_north_amer_dist_epitope_summary <- HA_north_amer_dist_epitope_summary %>% tidyr::separate(year_diff, into = c("year1", "year2"), remove = F)
HA_north_amer_dist_epitope_summary$season <- paste(as.numeric(HA_north_amer_dist_epitope_summary$year2) - 1, HA_north_amer_dist_epitope_summary$year2, sep = "-")

HA_north_amer_dist_epitope_summary %>%
  arrange(season, season_diff) %>%
  filter(season != "1995-1996")

save(HA_north_amer_dist_epitope_summary, file = "data/HA_mean_non_epitope_change_north_amer_build_summary.Rdata")
