## load packages
list.of.packages <- c(
  "dplyr", "ggplot2", "cowplot", "tidyr", "readr", "caret", "party", "glmnet", "foreach","readr",
  "gmodels", "moreparty", "permimp", "partykit", "rsample", "forcats", "tibble",
  "doParallel", "parallel"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

# stop the cluster
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list = ls(name = env), pos = env)
}
###########################################################################
## compile dataset
###########################################################################
# load data
load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")

pred_df <- read_rds("data/predictor_df_for_var_sel.rds")
sort(names(pred_df))

max_df <- left_join(pred_df, 
                    epi_red %>% dplyr::select(region, season, year.new, H3_max_intensity), 
                    by = c("region", "season", "year.new"))

max_df_cleaned <- max_df %>%
  filter(!(season %in% c("2009-2010", "1995-1996", "1996-1997"))) %>%
  dplyr::select(-year.new) %>%
  drop_na()
max_df_cleaned$region <- gsub(" ", "_", max_df_cleaned$region)

unique(max_df_cleaned$prior_dom_type_national)

max_df %>% filter(region=="Region 6")

###########################################################################
## prepare group CV folds
###########################################################################

max_df_cleaned %>%
  group_by(season) %>%
  group_indices() -> indeces
length(indeces)

# note: data set was ultimately too small to do a training and test set
set.seed(825)
split1 <- createDataPartition(as.factor(indeces), p = 1, list = FALSE)
length(split1)

unique(max_df_cleaned$season)
table(max_df_cleaned[split1, 1])
table(max_df_cleaned[-split1, 1])

train_df <- max_df_cleaned[split1, ]
nrow(train_df) # 189
test_df <- max_df_cleaned[-split1, ]
nrow(test_df) # 0

unique(max_df_cleaned$region)
unique(max_df_cleaned$season)

max_df_cleaned %>%
  group_by(season) %>%
  tally()

length(unique(max_df_cleaned$season)) # 20 seasons
split_by_season <- rsample::group_vfold_cv(max_df_cleaned, group = season) # 20 partitions
folds <- rsample2caret(split_by_season, data = c("analysis", "assessment"))
folds$index
# 10 repeats
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
  # returnResamp = "all"
)

###########################################################################
## cforest
###########################################################################
train_df %>% dplyr::select(-region, -season,-H3_max_intensity) %>% ncol()
set.seed(825)
# hyperparameter tuning
# mtry = number of randomly pre-selected variables
cforest_grid <- expand.grid(mtry = seq(from = 2, to = 21, by = 1)) # normally start with 2 to 20 (top limit of total features)

# Calculate the number of cores
no_cores <- detectCores() - 1
no_cores

# create the cluster for caret
cl <- makePSOCKcluster(no_cores)
registerDoParallel(cl)

rf_fit <- train(H3_max_intensity ~ .,
  data = train_df %>% dplyr::select(-region, -season),
  method = "cforest",
  preProcess = c("center", "scale"),
  tuneGrid = cforest_grid,
  trControl = group_fit_control,
  controls = cforest_unbiased(ntree = 3000)
)
stopCluster(cl)
unregister_dopar()

rf_fit$bestTune # mtry 12
head(rf_fit$resample)

## conditional permutation importance
set.seed(290875)
vi <- replicate(50, permimp(rf_fit$finalModel, conditional = TRUE, progressBar = T, scaled = T, nperm = 1)$values)

sum_df <- apply(as.matrix(vi), 1, function(x) ci(x))
sum_df <- sum_df %>%
  t() %>%
  as.data.frame()
vec <- rownames(sum_df)
str(sum_df)
sum_df$var <- rownames(sum_df)
colnames(sum_df)
colnames(sum_df)[2:4] <- c("CI_lower", "CI_upper", "std_error")

sum_df2 <- sum_df %>%
  mutate(var = fct_reorder(var, Estimate)) %>%
  arrange(-Estimate)
sum_df2
nrow(sum_df2)
write.csv(sum_df2, file = "data/H3_peak_cforest_variable_importance.csv", row.names = F)

ggplot(sum_df2, aes(x = Estimate, y = var)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5, fill = "blue") +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper, y = var), color = "red") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme_bw()

imp_vec <- sum_df2$Estimate
names(imp_vec) <- sum_df2$var

rf_robustvarimp <-
  imp_vec %>%
  as.data.frame() %>%
  rownames_to_column()
names(rf_robustvarimp) <- c("variable", "importance")
rf_robustvarimp %>% arrange(-importance)
rf_robustvarimp$model <- "cforest"

########################################################################################
## LASSO
########################################################################################
myGrid <- expand.grid(
  alpha = 1,
  lambda = seq(0.1, 2.5, length = 101)
)

set.seed(825)
glm_fit <- train(H3_max_intensity ~ .,
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
range(predVals$pred)
predVals %>% filter(pred < 0) # one observation from Lasso model
predVals$pred <- if_else(predVals$pred < 0, 0, predVals$pred)
ggplot(predVals %>% filter(object == "lasso")) +
  geom_point(aes(x = obs, y = pred))

pred_df <- full_join(epi_red, predVals, by = c("H3_max_intensity" = "obs"))

pred_df %>% filter(!is.na(model) & pred < 0)

plot_df <- pred_df %>% filter(!is.na(model))
plot_df$model <- as.factor(plot_df$model)
levels(plot_df$model) <- c("Conditional\nRandom Forest", "Lasso")

plot_df$resid <- plot_df$H3_max_intensity - plot_df$pred

write.csv(plot_df, "data/peak_incidence_obs_vs_predicted.csv", row.names = F)

# save variable importance output
all_imp_measures <- bind_rows(rf_robustvarimp, lasso_varimp)
all_imp_measures$model <- as.factor(all_imp_measures$model)
levels(all_imp_measures$model) <- c("Cond Random Forest", "Lasso")
write.csv(all_imp_measures, file = "data/H3_peak_ML_variable_importance.csv", row.names = F)
