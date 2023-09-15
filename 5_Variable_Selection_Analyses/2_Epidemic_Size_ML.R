## load packages
list.of.packages <- c(
  "dplyr", "ggplot2", "cowplot", "tidyr", "caret", "party", "glmnet", "foreach",
  "gmodels", "moreparty", "permimp", "partykit", "rsample", "forcats", "tibble",
  "doParallel", "parallel"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list = ls(name = env), pos = env)
}

###########################################################################
## compile dataset
###########################################################################

## load data
load("data/antigenic_epi_north_amer_build_for_lasso_replicates.Rdata")

epi_red <- epi_red %>% replace_na(list(
  H3_max_Rt = 0, H3_max_intensity = 0, H3_cum_intensity = 0,
  H3_season_duration = 0, H3_epi_size_prior = 0
))

epi_red <- epi_red %>%
  mutate(
    vac_combined = (adult_18_49_vac_cov / 100) * (adult_65_vac_cov / 100) * weighted_VE,
    vac_cov_combined = (adult_18_49_vac_cov / 100) * (adult_65_vac_cov / 100)
  )

vac_prior <- epi_red %>%
  distinct(season, adult_18_49_vac_cov, adult_65_vac_cov, vac_combined, vac_cov_combined, weighted_VE) %>%
  arrange(season) %>%
  mutate(
    adult_18_49_vac_cov_prior = lag(adult_18_49_vac_cov, n = 1),
    adult_65_vac_cov_prior = lag(adult_65_vac_cov, n = 1),
    vac_combined_prior = lag(vac_combined, n = 1),
    vac_cov_combined_prior = lag(vac_cov_combined, n = 1),
    weighted_VE_prior_season = lag(weighted_VE, n = 1)
  ) %>%
  ungroup()
vac_prior

epi_red <- left_join(epi_red %>% dplyr::select(-contains(c("vac", "VE"))), vac_prior, by = "season")

unique(epi_red$region)

vac_df <- epi_red %>%
  dplyr::select(H3_max_intensity, season, region, contains(c("vac_cov", "VE", "vac_combined"))) %>%
  distinct()
head(vac_df)

vac_df_mean <- vac_df %>%
  dplyr::select(-region) %>%
  group_by(season) %>%
  summarize_at(vars(H3_max_intensity, adult_18_49_vac_cov:vac_combined_prior), mean, na.rm = T) %>%
  distinct()

vac_df <- vac_df %>%
  mutate_at(vars(adult_18_49_vac_cov:vac_combined_prior), scale)

epi_red2 <- epi_red %>% dplyr::select(
  region, H3_cum_intensity, season, year.new,
  HA_titer_tree_lag2, HA_wolf_lag2, NA_bhatt_ep_lag1, # antigenic drift
  ha_lbi_shannon, na_lbi_shannon, # LBI diversity
  ha_lbi_shannon_lag1, na_lbi_shannon_lag1,
  prior_dom_type_national, # prior dominant IAV
  H1_cum_intensity, IVB_cum_intensity, # H1N1 and IBV epidemic size
  H3_epi_size_prior, H1_epi_size_prior, # prior epidemic sizes
  IVB_epi_size_prior,
  vac_combined, vac_cov_combined, weighted_VE, # vaccine parameters
  vac_combined_prior, vac_cov_combined_prior,
  weighted_VE_prior_season,
  usa_bhatt_ep_mean, usa_wolf_ep_mean # distance between circulating strains and US vaccine strain
) 

cum_df <- epi_red2 %>%
  tidyr::replace_na(list(
    IVB_max_intensity = 0, H1_max_intensity = 0,
    IVB_cum_intensity = 0, H1_cum_intensity = 0
  )) %>%
  mutate(
    IVB_epi_size_prior = ifelse(year.new > 1997 & is.na(IVB_epi_size_prior), 0, IVB_epi_size_prior),
    H1_epi_size_prior = ifelse(year.new > 1997 & is.na(H1_epi_size_prior), 0, H1_epi_size_prior)
  ) %>%
  distinct()

cum_df_cleaned <- cum_df %>%
  filter(!(season %in% c("2009-2010", "1995-1996", "1996-1997"))) %>%
  dplyr::select(-year.new) %>%
  drop_na()
cum_df_cleaned$region <- gsub(" ", "_", cum_df_cleaned$region)

cum_df_cleaned <- cum_df_cleaned %>% mutate(prior_dom_type_national = case_when(
  prior_dom_type_national == "H3" ~ 1,
  prior_dom_type_national == "H1" ~ 0,
  prior_dom_type_national == "co-circ" ~ 0.5
))
cum_df_cleaned$prior_dom_type_national <- as.numeric(cum_df_cleaned$prior_dom_type_national)

###########################################################################
## prepare group CV folds
###########################################################################

cum_df_cleaned %>%
  group_by(season) %>%
  group_indices() -> indeces

set.seed(825)
# note: data set was ultimately too small to do a training and test set
split1 <- createDataPartition(as.factor(indeces), p = 1, list = FALSE)
length(split1)

table(cum_df_cleaned[split1, 1])
table(cum_df_cleaned[-split1, 1])

train_df <- cum_df_cleaned[split1, ]
nrow(train_df) # 189
test_df <- cum_df_cleaned[-split1, ]
nrow(test_df)

split_by_season <- rsample::group_vfold_cv(cum_df_cleaned, group = season)
folds <- rsample2caret(split_by_season, data = c("analysis", "assessment"))
folds$index
folds2 <- lapply(folds$index, function(x) lapply(1:10, function(i) sample(x, size = length(x), replace = TRUE)))
folds2 <- unlist(folds2, recursive = FALSE, use.names = TRUE)

folds_out <- folds$indexOut
folds_out <- rep(folds_out, each = 10)
names(folds_out) <- names(folds2)
setdiff(names(folds_out), names(folds2)) # sanity check
setdiff(names(folds2), names(folds_out))

# control parameters for model training
group_fit_control <- trainControl( ## use grouped CV folds
  index = folds2,
  indexOut = folds_out,
  method = "boot",
  allowParallel = TRUE
)

###########################################################################
## cforest
###########################################################################
set.seed(825)
cforest_grid <- expand.grid(mtry = seq(from = 6, to = 14, by = 2)) # normally start with 2 to 20

# Calculate the number of cores
no_cores <- detectCores() - 1
no_cores

# create the cluster for caret to use
cl <- makePSOCKcluster(no_cores)
registerDoParallel(cl)

rf_fit <- train(H3_cum_intensity ~ .,
  data = train_df %>% dplyr::select(-region, -season),
  method = "cforest",
  preProcess = c("center", "scale"),
  tuneGrid = cforest_grid,
  trControl = group_fit_control,
  controls = cforest_unbiased(ntree = 3000)
)
stopCluster(cl)
unregister_dopar()

rf_fit$bestTune # mtry 10

## conditional permutation importance
vi <- replicate(50, permimp(rf_fit$finalModel, conditional = TRUE, progressBar = T, scaled = T, nperm = 1)$values)

sum_df <- apply(as.matrix(vi), 1, function(x) ci(x))
str(sum_df)
sum_df <- sum_df %>%
  t() %>%
  as.data.frame()
vec <- rownames(sum_df)
str(sum_df)
sum_df$var <- rownames(sum_df)
colnames(sum_df)[2:4] <- c("CI_lower", "CI_upper", "std_error")

sum_df2 <- sum_df %>% mutate(var = fct_reorder(var, Estimate))

ggplot(sum_df2, aes(x = Estimate, y = var)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5, fill = "blue") +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper, y = var), color = "red") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme_bw()
write.csv(sum_df2, file = "data/H3_epi_size_cforest_variable_importance.csv", row.names = F)

imp_vec <- sum_df2$Estimate
names(imp_vec) <- sum_df2$var

rf_robustvarimp <-
  imp_vec %>%
  as.data.frame() %>%
  rownames_to_column()
names(rf_robustvarimp) <- c("variable", "importance")
rf_robustvarimp %>% arrange(-importance)
rf_robustvarimp$model <- "cforest"

###############################################################
## LASSO
###############################################################

myGrid <- expand.grid(
  alpha = 1,
  lambda = seq(0.001, 10, length = 101)
)

glm_fit <- train(H3_cum_intensity ~ .,
  data = dplyr::select(train_df, -region, -season),
  method = "glmnet",
  trControl = group_fit_control,
  tuneGrid = myGrid,
  preProcess = c("center", "scale")
)

head(glm_fit)
glm_fit$bestTune$lambda
ggplot(glm_fit) +
  labs(title = "Lasso Regression Parameter Tuning", x = "lambda") +
  theme_bw()
lasso_varimp <- caret::varImp(glm_fit)
plot(lasso_varimp)
lasso_varimp <- lasso_varimp$importance %>%
  as.data.frame() %>%
  rownames_to_column()
names(lasso_varimp) <- c("variable", "importance")
lasso_varimp$model <- "lasso"

########################################################################################
## Extract predictions
########################################################################################

model_list <- list("party RF" = rf_fit, "lasso" = glm_fit)

predVals <- caret::extractPrediction(models = model_list)

ggplot(predVals %>% filter(object == "lasso")) +
  geom_point(aes(x = obs, y = pred))

range(predVals$pred)
predVals %>% filter(pred < 0)
predVals$pred <- if_else(predVals$pred < 0, 0, predVals$pred)
pred_df <- full_join(epi_red, predVals, by = c("H3_cum_intensity" = "obs"))

pred_df %>% filter(!is.na(model) & pred < 0)

plot_df <- pred_df %>% filter(!is.na(model))
plot_df$model <- as.factor(plot_df$model)
levels(plot_df$model)
levels(plot_df$model) <- c("Conditional\nRandom Forest", "Lasso")

plot_df$resid <- plot_df$H3_cum_intensity - plot_df$pred
write.csv(plot_df, "data/epi_size_incidence_obs_vs_predicted.csv", row.names = F)

# save variable importance output
all_imp_measures <- bind_rows(rf_robustvarimp, lasso_varimp)
all_imp_measures$model <- as.factor(all_imp_measures$model)
levels(all_imp_measures$model) <- c("Cond Random Forest", "Lasso")
all_imp_measures$model <- factor(all_imp_measures$model, levels = c("Cond Random Forest", "Lasso"))
unique(all_imp_measures$variable)
write.csv(all_imp_measures, file = "data/H3_epi_size_ML_variable_importance.csv", row.names = F)
