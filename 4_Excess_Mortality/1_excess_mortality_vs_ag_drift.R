## load packages
list.of.packages <- c(
  "dplyr", "ggplot2", "cowplot", "tidyr", "lubridate","zoo","magicfor","purrr",
  "ggpubr","tidymodels", "MuMIn","performance","rsample"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

theme_set(theme_cowplot(font_size = 16))

######################################################################################
## Excess mortality attributable to H3N2 vs antigenic drift
######################################################################################
#load data
load("data/hhs_division_level_ILI_and_virology_interp_smoothed_sampling_effort.RData")
head(regionflu_ili_vir_adj)
names(regionflu_ili_vir_adj)

viral_surv = regionflu_ili_vir_adj %>% 
  dplyr::select(region,season_description,year,week,wk_date,total_specimens,percent_positive,adjusted_h3,prop_h3)%>%
  mutate(adjusted_h3 = na.approx(adjusted_h3,na.rm=4,maxgap=4),
         total_flu = total_specimens * (percent_positive/100))%>%
  group_by(season_description)%>%
  summarize(total_flu = sum(total_flu,na.rm = T),
            total_h3 = sum(adjusted_h3,na.rm = T))%>%
  mutate(prop_h3 = total_h3/total_flu)


load("data/antigenic_epi_north_amer_build_for_ML_replicates.Rdata")
head(epi_red)

ag_df <- 
  epi_red %>%
  dplyr::select(season,dom_type,contains(c("lag1","lag2","lbi"))) %>% 
  dplyr::select(-contains("sd"))%>%
  distinct()
head(ag_df)

# Age groups are (1=<1; 2=1-4; 3=5-49; 4=50-64; 5=65+)
path_file =  "4_Excess_Mortality/excess_mortality_estimates/"

filenames_list <- list.files(path = path_file, full.names=TRUE,pattern="csv")
head(filenames_list)
merging_manifests <- lapply(filenames_list, function(filename){
  print(paste("Merging",filename,sep = " "))
  read.csv(filename)
})
master_file <- do.call(plyr::rbind.fill, merging_manifests)
head(master_file)
unique(master_file$var)
master_file$season = gsub("/","-",master_file$season)
head(master_file)
range(master_file$predicted)
master_file %>%filter(age==5)

flu_deaths = master_file[,c("date","season","flu_ah3","flu_upper_ah3","flu_lower_ah3", "age")]
range(flu_deaths$flu_ah3)
unique(flu_deaths$age)
head(master_file)
# Age groups are (1=<1; 2=1-4; 3=5-49; 4=50-64; 5=65+)

flu_deaths = flu_deaths %>%
  mutate(age_grp = case_when(age==1~"<1",
                             age==2~"1-4",
                             age==3~"5-49",
                             age==4~"50-64",
                             age==5~"65+"))%>%
  mutate(date=as.Date(date))%>%
  mutate(year = lubridate::year(date))%>%
  mutate(month = lubridate::month(date))%>%
  mutate(month.day = format(date,"%m-%d")) %>%
  mutate(season=ifelse(month.day < "07-01",sprintf("%d-%d", year-1, year), #season = july 1 to june 30
                       sprintf("%d-%d", year, year+1))) %>%
  mutate(year.new=as.numeric(substr(season,1,4)))
unique(flu_deaths$season)

flu_deaths %>% filter(age==5)

flu_deaths_agg = 
  flu_deaths %>%
  group_by(season,age,age_grp)%>%
  summarize(h3_deaths=sum(flu_ah3),
            h3_deaths_upr=sum(flu_upper_ah3),
            h3_deaths_lwr=sum(flu_lower_ah3))%>%
  ungroup()

mortality_ag_df <- left_join(flu_deaths_agg, ag_df,by="season") %>%
  filter(!(season %in% c("2009-2010","2019-2020")))

mortality_long = mortality_ag_df %>% pivot_longer(cols=contains(c("lag1","lag2","lbi","titer")),
                                                  names_to = "evol_metrics",
                                                  values_to="evol_value")
head(mortality_long)
names(mortality_long)

age_grp = unique(mortality_long$age_grp)
age_grp = age_grp[!is.na(age_grp)]
age_grp

magic_for(func=put,temporary = T)
for (i in age_grp){
  mortality_predictors = mortality_long %>%
    filter(age_grp == i) %>%
    dplyr::select(season,evol_metrics,evol_value,h3_deaths)%>%
    split(.$evol_metrics) %>% #
    map(~ lm(h3_deaths ~ evol_value, data = .)) %>%
    map(summary) %>%
    map_dbl("adj.r.squared")
  top = sort(mortality_predictors,decreasing = T)[1:2]
  df = data.frame(outcome = i ,r2 = top, metric = names(top))
  put(df)
}
result = magic_result()
result_df = bind_rows(result)%>% arrange(-r2) %>% as.data.frame()
rownames(result_df) <- NULL
head(result_df)
result_df

unique(mortality_long$dom_type)

magic_for(func=put,temporary = T)
for (i in age_grp){
  mortality_predictors = mortality_long %>%
    filter(age_grp == i & dom_type %in% c("H3","co-circ")) %>%
    dplyr::select(season,evol_metrics,evol_value,h3_deaths)%>%
    split(.$evol_metrics) %>% #
    # split(f=list(.$agegrp,.$evol_metrics))%>%
    map(~ lm(h3_deaths ~ evol_value, data = .)) %>%
    map(summary) %>%
    map_dbl("adj.r.squared")
  top = sort(mortality_predictors,decreasing = T)[1:2]
  df = data.frame(outcome = i ,r2 = top, metric = names(top))
  put(df)
}
result = magic_result()
result_df = bind_rows(result)%>% arrange(-r2) %>% as.data.frame()
rownames(result_df) <- NULL
head(result_df)
result_df

names(mortality_long)
unique(mortality_long$evol_metrics)

red = mortality_long %>%
  dplyr::select(season,dom_type,h3_deaths,h3_deaths_lwr,h3_deaths_upr, 
                evol_metrics,evol_value,age,age_grp)%>%
  filter(evol_metrics %in% c("HA_wolf_lag2","NA_bhatt_ep_lag1","HA_titer_tree_lag2"))%>%
  filter(!is.na(age))
head(red)

red = red %>% tidyr::separate(col="season",sep="-",remove=F,into=c("year1","year2"))

red = red %>%
  mutate(dom_type2 = case_when(year1<2009&dom_type=="H1"~"H1",
                               year1>2009&dom_type=="H1"~"H1pdm",
                               dom_type=="H3"~"H3",
                               dom_type=="co-circ"~"H3/H1pdm"))

ggplot(red %>% filter(dom_type %in% c("H3","co-circ")),aes(x=evol_value,y=h3_deaths))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(age_grp~evol_metrics,scales = "free")

unique(mortality_long$age_grp)

######################################################################################
## all ages combined
######################################################################################

red2 = mortality_long %>%
  dplyr::select(season,dom_type,h3_deaths,evol_metrics,h3_deaths_lwr,h3_deaths_upr,evol_value,age,age_grp)%>%
  filter(!is.na(age))%>%
  group_by(season,evol_metrics,evol_value,dom_type)%>%
  summarise(h3_deaths=sum(h3_deaths),
            h3_deaths_lwr=sum(h3_deaths_lwr),
            h3_deaths_upr=sum(h3_deaths_upr) )%>%
  ungroup()

red2 = red2 %>% tidyr::separate(col="season",sep="-",remove=F,into=c("year1","year2"))

red2 = red2 %>%
  mutate(dom_type2 = case_when(year1<2009&dom_type=="H1"~"H1",
                               year1>2009&dom_type=="H1"~"H1pdm",
                               dom_type=="H3"~"H3",
                               dom_type=="co-circ"~"H3/H1pdm"))

ggplot(red2,aes(x=evol_value,y=h3_deaths))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~evol_metrics,scales = "free")

red2%>%
  split(.$evol_metrics) %>% #
  map(~ lm(h3_deaths ~ evol_value, data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")%>%
  sort(decreasing = T)

######################################################################################
## age 5 is 65>
######################################################################################

red3 = mortality_long %>%
  dplyr::select(season,dom_type,h3_deaths,h3_deaths_lwr,h3_deaths_upr,evol_metrics,evol_value,age,age_grp)%>%
  filter(age==5)%>%
  group_by(season,evol_metrics,evol_value,dom_type)%>%
  summarise(h3_deaths=sum(h3_deaths),
            h3_deaths_lwr=sum(h3_deaths_lwr),
            h3_deaths_upr=sum(h3_deaths_upr))%>%
  ungroup()
head(red3)

red3 = red3 %>% tidyr::separate(col="season",sep="-",remove=F,into=c("year1","year2"))

red3 = red3 %>%
  mutate(dom_type2 = case_when(year1<2009&dom_type=="H1"~"H1",
                               year1>2009&dom_type=="H1"~"H1pdm",
                               dom_type=="H3"~"H3",
                               dom_type=="co-circ"~"H3/H1pdm"))

red3%>%
  split(.$evol_metrics) %>% #
  map(~ lm(h3_deaths ~ scale(evol_value), data = .)) %>%
  map(summary) %>%
  map_dbl("adj.r.squared")%>%
  sort(decreasing = T)

cols= c("#B24745FF","#00A1D5FF","#6A6599FF","#DF8F44FF")
# red$evol_metrics = as.factor(red$evol_metrics)
# levels(red$evol_metrics)
# levels(red$evol_metrics)<- c("HI titer distance (t-2)", "H3 epitope distance (t-2)","N2 epitope distance (t-1)")
# 
# red$age = as.factor(red$age)
# red$age_grp = as.factor(red$age_grp)
# unique(red$dom_type)

######################################################################################
## HA epitope distance
######################################################################################

######################################################################################
## all ages combined
######################################################################################

red_wide = red2 %>%
  pivot_wider(names_from = evol_metrics,values_from = evol_value)

names(red_wide)
scale_this <- function(x) as.vector(scale(x,center=T))
red_wide = red_wide %>%
  mutate_at(vars(HA_koel_lag1:na_lbi_shannon),~scale_this(.x))

set.seed(27)
boots <- bootstraps(red_wide, times = 1000, apparent = TRUE)

y =red_wide$h3_deaths
x= red_wide$HA_wolf_lag2

linear.model = glm(y~x,family=gaussian())
log.model = glm(y~x,family=gaussian(link="log"))
inv.model = glm(y~x,family=gaussian(link="inverse"))
gamma.model = glm(y~x,family=Gamma(link="log"))
gamma.model2 = glm(y~x,family=Gamma(link="inverse"))
# gamma.model3 = glm(y~x,family=Gamma(link="identity"))
model.sel(linear.model,log.model,inv.model,gamma.model,gamma.model2)

compare_performance(linear.model,log.model,inv.model,gamma.model,gamma.model2,rank = T)

set.seed(27)
boots <- bootstraps(red_wide, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(h3_deaths ~ HA_wolf_lag2, analysis(split), family=gaussian(link="log"),maxit=150)
}

boot_models <-
  boots %>% 
  mutate(model = map(splits, fit_glm_on_bootstrap),
         coef_info = map(model, tidy))

boot_coefs <- 
  boot_models %>% 
  unnest(coef_info)

m1 = glm(h3_deaths ~ HA_wolf_lag2,data=red_wide,family=gaussian(link="log"),maxit=150)
summary(m1)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance/m1$null.deviance #0.23


labels = boots %>%
  mutate(model = purrr::map(splits, ~ glm(h3_deaths ~ HA_wolf_lag2, data = .x,family=gaussian(link="log"),maxit=150)), 
         adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance/.x$null.deviance, 5)),
         pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))) %>%
  summarize(adj.r.squared = mean(adj.r.squared),
            pvalue = mean(pvalue))%>%
  mutate(adj.r.squared =  sprintf("italic(R^2) == %.2f", adj.r.squared),
         pvalue =  sprintf("italic(p) == %.2f", pvalue))
labels

set.seed(27)
n_boot <- 1000
red_wide %>% 
  dplyr::select(h3_deaths,NA_bhatt_ep_lag1,HA_wolf_lag2) %>% 
  pivot_longer(cols = c(NA_bhatt_ep_lag1,HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>% 
  nest(model_data = c(h3_deaths, value)) %>% 
  mutate(plot_data = map(model_data, ~ {

    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)
    
    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model, 
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>% 
        sample_n(n, TRUE) %>% 
        glm(h3_deaths ~ value, .,family=gaussian(link="log"),maxit=150) %>% 
        # suppress augment() warnings about dropping columns
        { suppressWarnings(augment(., newdata = tibble(value = pred_x))) }
    }) %>% 
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>% 
      # to estimate the lower and upper 95% quantiles
      summarize(l = quantile(.fitted, .025),
                u = quantile(.fitted, .975),
                .groups = "drop"
      ) %>% 
      arrange(value)
  })) %>% 
  dplyr::select(-model_data) %>% 
  unnest(plot_data) -> tbl_plot_data

cols= c("#B24745FF","#00A1D5FF","#6A6599FF","#DF8F44FF")

ep_deaths = ggplot()+
  geom_ribbon(aes(x=h3_deaths, ymin = gaussian(link="log")$linkinv(l), ymax = gaussian(link="log")$linkinv(u)), 
              tbl_plot_data %>% filter(name=="HA_wolf_lag2") %>% rename(h3_deaths=value), 
              alpha = 0.3, fill = "grey") +
  geom_smooth(data=red_wide%>% filter(season!="2009-2010"),
              aes(x=HA_wolf_lag2,y=h3_deaths),
              method = "glm", formula = y~x,
              method.args = list(family = gaussian(link = 'log')),
              se=F,
              linewidth = 1, linetype = 2,color="black")+
  geom_point(data=red_wide%>% filter(season!="2009-2010"),
             aes(x=HA_wolf_lag2,y=h3_deaths,fill=dom_type2),size=4,pch=21)+
  geom_errorbar(data=red_wide%>% filter(season!="2009-2010"),
                aes(x=HA_wolf_lag2,
                    ymin=h3_deaths_lwr, ymax=h3_deaths_upr,
                    color=dom_type2), width=.05)+
  theme_bw(base_size = 16)+
  xlab(expression("H3 epitope distance ("~italic(t)~"-2)"))+
  ylab("Excess A(H3N2) Mortality")+
  theme(legend.position = c(0.6,0.8),legend.title = element_blank())+
  background_grid(major = "xy", minor = "none")+
  scale_color_manual(breaks=c("H3","H1","H1pdm","H3/H1pdm"),
                     labels=c("H3","H1","H1pdm","H3/H1pdm"),values=cols)+
  scale_fill_manual(breaks=c("H3","H1","H1pdm","H3/H1pdm"),
                    labels=c("H3","H1","H1pdm","H3/H1pdm"),values=cols)+
  geom_text(y = 26, x =-1.5, aes(label = paste(adj.r.squared,pvalue,sep = "~`,`~")), size=5,
            data = labels, parse = TRUE, hjust = 0)+
  ggtitle("All Ages")+
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5,face = "plain"))+
  ylim(NA,29)
ep_deaths

######################################################################################
##### over 65
######################################################################################

red_wide = red3 %>%
  pivot_wider(names_from = evol_metrics,values_from = evol_value)

scale_this <- function(x) as.vector(scale(x,center=T))
red_wide = red_wide %>%
  mutate_at(vars(HA_koel_lag1:na_lbi_shannon),~scale_this(.x))

red3 %>% filter(dom_type=="H1") %>% distinct(season,h3_deaths,h3_deaths_lwr,h3_deaths_upr)

set.seed(27)
boots <- bootstraps(red_wide, times = 1000, apparent = TRUE)

y =red_wide$h3_deaths
x= red_wide$HA_wolf_lag2

linear.model = glm(y~x,family=gaussian())
log.model = glm(y~x,family=gaussian(link="log"))
inv.model = glm(y~x,family=gaussian(link="inverse"))
gamma.model = glm(y~x,family=Gamma(link="log"))
gamma.model2 = glm(y~x,family=Gamma(link="inverse"))
# gamma.model3 = glm(y~x,family=Gamma(link="identity"))

model.sel(linear.model,log.model,inv.model,gamma.model,gamma.model2)

set.seed(27)
boots <- bootstraps(red_wide, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(h3_deaths ~ HA_wolf_lag2, analysis(split), family=gaussian(link="log"),maxit=150)
}

boot_models <-
  boots %>% 
  mutate(model = map(splits, fit_glm_on_bootstrap),
         coef_info = map(model, tidy))

boot_coefs <- 
  boot_models %>% 
  unnest(coef_info)

m1 = glm(h3_deaths ~ HA_wolf_lag2,data=red_wide,family=gaussian(link="log"),maxit=150)
summary(m1)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance/m1$null.deviance

labels = boots %>%
  mutate(model = purrr::map(splits, ~ glm(h3_deaths ~ HA_wolf_lag2, data = .x,family=gaussian(link="log"),maxit=150)), 
         adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance/.x$null.deviance, 5)),
         pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))) %>%
  summarize(adj.r.squared = mean(adj.r.squared),
            pvalue = mean(pvalue))%>%
  mutate(adj.r.squared =  sprintf("italic(R^2) == %.2f", adj.r.squared),
         pvalue =  sprintf("italic(p) == %.2f", pvalue))
labels

set.seed(27)
n_boot <- 1000
red_wide %>% 
  dplyr::select(h3_deaths,NA_bhatt_ep_lag1,HA_wolf_lag2) %>% 
  pivot_longer(cols = c(NA_bhatt_ep_lag1,HA_wolf_lag2)) -> tbl_mtcars_long

tbl_mtcars_long %>% 
  nest(model_data = c(h3_deaths, value)) %>% 
  mutate(plot_data = map(model_data, ~ {

    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)
    
    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model, 
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>% 
        sample_n(n, TRUE) %>% 
        glm(h3_deaths ~ value, .,family=gaussian(link="log"),maxit=150) %>% 
        # suppress augment() warnings about dropping columns
        { suppressWarnings(augment(., newdata = tibble(value = pred_x))) }
    }) %>% 
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>% 
      # to estimate the lower and upper 95% quantiles of predicted mpgs
      summarize(l = quantile(.fitted, .025),
                u = quantile(.fitted, .975),
                .groups = "drop"
      ) %>% 
      arrange(value)
  })) %>% 
  dplyr::select(-model_data) %>% 
  unnest(plot_data) -> tbl_plot_data

ep_deaths_over65 = ggplot()+
  geom_ribbon(aes(x=h3_deaths, ymin = gaussian(link="log")$linkinv(l), ymax = gaussian(link="log")$linkinv(u)), 
              tbl_plot_data %>% filter(name=="HA_wolf_lag2") %>% rename(h3_deaths=value), 
              alpha = 0.3, fill = "grey") +
  geom_point(data=red_wide%>% filter(season!="2009-2010"),
             aes(x=HA_wolf_lag2,y=h3_deaths,fill=dom_type2),size=4,pch=21)+
  geom_errorbar(data=red_wide%>% filter(season!="2009-2010"),
                aes(x=HA_wolf_lag2,
                    ymin=h3_deaths_lwr, ymax=h3_deaths_upr,
                    color=dom_type2), width=.05)+
  theme_bw(base_size = 16)+
  xlab(expression("H3 epitope distance ("~italic(t)~"-2)"))+
  ylab("Excess A(H3N2) Mortality")+
  theme(legend.position = "bottom")+
  background_grid(major = "xy", minor = "none")+
  scale_color_manual(breaks=c("H3","H1","H1pdm","H3/H1pdm"),
                     labels=c("H3","H1","H1pdm","H3/H1pdm"),values=cols,name="Dominant IAV")+
  scale_fill_manual(breaks=c("H3","H1","H1pdm","H3/H1pdm"),
                    labels=c("H3","H1","H1pdm","H3/H1pdm"),values=cols,name="Dominant IAV")+
  geom_smooth(data=red_wide%>% filter(season!="2009-2010"),
              aes(x=HA_wolf_lag2,y=h3_deaths),
              method = "glm", formula = y~x,
              method.args = list(family = gaussian(link = 'log')),
              se=F,
              size = 1, linetype = 2,color="black")+
  geom_text(y = 26, x =-1.5, aes(label = paste(adj.r.squared,pvalue,sep = "~`,`~")), size=5,
            data = labels, parse = TRUE, hjust = 0)+
  ggtitle(expression(paste('\u2265',"65 years of age")))+
  theme(plot.title = element_text(hjust = 0.5,face = "plain"))+
  ylim(NA,29)
ep_deaths_over65

HA_excess_plot = plot_grid(ep_deaths+theme(legend.position = "none"),
                           ep_deaths_over65+theme(legend.position = "none"),
                           nrow=1,labels = "AUTO")
leg = get_legend(ep_deaths_over65+
                   theme(legend.direction = "horizontal",legend.justification="center" ,
                         legend.position = "bottom"))
both = plot_grid(HA_excess_plot,leg,nrow=2,rel_heights = c(1,0.1))
both

######################################################################################
## NA epitope distance
######################################################################################

######################################################################################
## all ages combined
######################################################################################

red_wide = red2 %>%
  pivot_wider(names_from = evol_metrics,values_from = evol_value)

scale_this <- function(x) as.vector(scale(x,center=T))
red_wide = red_wide %>%
  mutate_at(vars(HA_koel_lag1:na_lbi_shannon),~scale_this(.x))

set.seed(27)
boots <- bootstraps(red_wide, times = 1000, apparent = TRUE)

y =red_wide$h3_deaths
x= red_wide$NA_bhatt_ep_lag1

linear.model = glm(y~x,family=gaussian())
log.model = glm(y~x,family=gaussian(link="log"))
inv.model = glm(y~x,family=gaussian(link="inverse"))
gamma.model = glm(y~x,family=Gamma(link="log"))
gamma.model2 = glm(y~x,family=Gamma(link="inverse"))
# gamma.model3 = glm(y~x,family=Gamma(link="identity"))

model.sel(linear.model,log.model,inv.model,gamma.model,gamma.model2)

compare_performance(linear.model,log.model,inv.model,gamma.model,gamma.model2,rank = T)

set.seed(27)
boots <- bootstraps(red_wide, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(h3_deaths ~ NA_bhatt_ep_lag1, analysis(split), family=gaussian(link="log"),maxit=150)
}

boot_models <-
  boots %>% 
  mutate(model = map(splits, fit_glm_on_bootstrap),
         coef_info = map(model, tidy))

boot_coefs <- 
  boot_models %>% 
  unnest(coef_info)

m1 = glm(h3_deaths ~ NA_bhatt_ep_lag1,data=red_wide,family=gaussian(link="log"),maxit=150)
summary(m1)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance/m1$null.deviance #0.07

labels = boots %>%
  mutate(model = purrr::map(splits, ~ glm(h3_deaths ~ NA_bhatt_ep_lag1, data = .x,family=gaussian(link="log"),maxit=150)), 
         adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance/.x$null.deviance, 5)),
         pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))) %>%
  summarize(adj.r.squared = mean(adj.r.squared),
            pvalue = mean(pvalue))%>%
  mutate(adj.r.squared =  sprintf("italic(R^2) == %.1f", adj.r.squared),
         pvalue =  sprintf("italic(p) == %.1f", pvalue))
labels

set.seed(27)
n_boot <- 1000
red_wide %>% 
  dplyr::select(h3_deaths,NA_bhatt_ep_lag1,NA_bhatt_ep_lag1) %>% 
  pivot_longer(cols = c(NA_bhatt_ep_lag1,NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>% 
  nest(model_data = c(h3_deaths, value)) %>% 
  mutate(plot_data = map(model_data, ~ {
    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)
    
    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model, 
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>% 
        sample_n(n, TRUE) %>% 
        glm(h3_deaths ~ value, .,family=gaussian(link="log"),maxit=150) %>% 
        # suppress augment() warnings about dropping columns
        { suppressWarnings(augment(., newdata = tibble(value = pred_x))) }
    }) %>% 
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>% 
      # to estimate the lower and upper 95% quantiles of predicted mpgs
      summarize(l = quantile(.fitted, .025),
                u = quantile(.fitted, .975),
                .groups = "drop"
      ) %>% 
      arrange(value)
  })) %>% 
  dplyr::select(-model_data) %>% 
  unnest(plot_data) -> tbl_plot_data

ep_deaths = ggplot()+
  geom_ribbon(aes(x=h3_deaths, ymin = gaussian(link="log")$linkinv(l), ymax = gaussian(link="log")$linkinv(u)), 
              tbl_plot_data %>% filter(name=="NA_bhatt_ep_lag1") %>% rename(h3_deaths=value), 
              alpha = 0.3, fill = "grey") +
  geom_smooth(data=red_wide%>% filter(season!="2009-2010"),
              aes(x=NA_bhatt_ep_lag1,y=h3_deaths),
              method = "glm", formula = y~x,
              method.args = list(family = gaussian(link = 'log')),
              se=F,
              size = 1, linetype = 2,color="black")+
  geom_point(data=red_wide%>% filter(season!="2009-2010"),
             aes(x=NA_bhatt_ep_lag1,y=h3_deaths,fill=dom_type2),size=4,pch=21)+
  geom_errorbar(data=red_wide%>% filter(season!="2009-2010"),
                aes(x=NA_bhatt_ep_lag1,
                    ymin=h3_deaths_lwr, ymax=h3_deaths_upr,
                    color=dom_type2), width=.05)+
  theme_bw(base_size = 16)+
  xlab(expression("N2 epitope distance ("~italic(t)~"-1)"))+
  ylab("Excess A(H3N2) Mortality")+
  theme(legend.position = c(0.6,0.8),legend.title = element_blank())+
  background_grid(major = "xy", minor = "none")+
  scale_color_manual(breaks=c("H3","H1","H1pdm","H3/H1pdm"),
                     labels=c("H3","H1","H1pdm","H3/H1pdm"),values=cols,name="Dominant IAV")+
  scale_fill_manual(breaks=c("H3","H1","H1pdm","H3/H1pdm"),
                    labels=c("H3","H1","H1pdm","H3/H1pdm"),values=cols,name="Dominant IAV")+  
  geom_text(y = 26, x =-1.2, aes(label = paste(adj.r.squared,pvalue,sep = "~`,`~")), size=5,
            data = labels, parse = TRUE, hjust = 0)+
  ggtitle("All Ages")+
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5,face = "plain"))+
  ylim(NA,29)
ep_deaths

######################################################################################
## over 65
######################################################################################

red_wide = red3 %>%
  pivot_wider(names_from = evol_metrics,values_from = evol_value)

scale_this <- function(x) as.vector(scale(x,center=T))
red_wide = red_wide %>%
  mutate_at(vars(HA_koel_lag1:na_lbi_shannon),~scale_this(.x))

red3 %>% filter(dom_type=="H1") %>% distinct(season,h3_deaths,h3_deaths_lwr,h3_deaths_upr)

set.seed(27)
boots <- bootstraps(red_wide, times = 1000, apparent = TRUE)
boots

y =red_wide$h3_deaths
x= red_wide$NA_bhatt_ep_lag1

linear.model = glm(y~x,family=gaussian())
log.model = glm(y~x,family=gaussian(link="log"))
inv.model = glm(y~x,family=gaussian(link="inverse"))
gamma.model = glm(y~x,family=Gamma(link="log"))
gamma.model2 = glm(y~x,family=Gamma(link="inverse"))
# gamma.model3 = glm(y~x,family=Gamma(link="identity"))

model.sel(linear.model,log.model,inv.model,gamma.model,gamma.model2)

set.seed(27)
boots <- bootstraps(red_wide, times = 1000, apparent = TRUE)

fit_glm_on_bootstrap <- function(split) {
  glm(h3_deaths ~ NA_bhatt_ep_lag1, analysis(split), family=gaussian(link="log"),maxit=150)
}

boot_models <-
  boots %>% 
  mutate(model = map(splits, fit_glm_on_bootstrap),
         coef_info = map(model, tidy))

boot_coefs <- 
  boot_models %>% 
  unnest(coef_info)

m1 = glm(h3_deaths ~ NA_bhatt_ep_lag1,data=red_wide,family=gaussian(link="log"),maxit=150)
summary(m1)
fam <- family(m1)
ilink <- fam$linkinv

summary(m1)$coefficients[8]
1 - m1$deviance/m1$null.deviance #0.07

labels = boots %>%
  mutate(model = purrr::map(splits, ~ glm(h3_deaths ~ NA_bhatt_ep_lag1, data = .x,family=gaussian(link="log"),maxit=150)), 
         adj.r.squared = map_dbl(model, ~ signif(1 - .x$deviance/.x$null.deviance, 5)),
         pvalue = map_dbl(model, ~ signif(summary(.x)$coefficients[8], 5))) %>%
  summarize(adj.r.squared = mean(adj.r.squared),
            pvalue = mean(pvalue))%>%
  mutate(adj.r.squared =  sprintf("italic(R^2) == %.1f", adj.r.squared),
         pvalue =  sprintf("italic(p) == %.1f", pvalue))
labels

set.seed(27)
n_boot <- 1000
red_wide %>% 
  dplyr::select(h3_deaths,NA_bhatt_ep_lag1,NA_bhatt_ep_lag1) %>% 
  pivot_longer(cols = c(NA_bhatt_ep_lag1,NA_bhatt_ep_lag1)) -> tbl_mtcars_long

tbl_mtcars_long %>% 
  nest(model_data = c(h3_deaths, value)) %>% 
  mutate(plot_data = map(model_data, ~ {

    submodel_data <- .x
    n <- nrow(submodel_data)
    min_x <- min(submodel_data$value)
    max_x <- max(submodel_data$value)
    pred_x <- seq(min_x, max_x, length.out = 100)
    
    # do the bootstrapping by
    # 1) repeatedly sampling samples of size n with replacement n_boot times,
    # 2) for each bootstrap sample, fit a model, 
    # 3) and make a tibble of predictions
    # the _dfr means to stack each tibble of predictions on top of one another
    map_dfr(1:n_boot, ~ {
      submodel_data %>% 
        sample_n(n, TRUE) %>% 
        glm(h3_deaths ~ value, .,family=gaussian(link="log"),maxit=150) %>% 
        # suppress augment() warnings about dropping columns
        { suppressWarnings(augment(., newdata = tibble(value = pred_x))) }
    }) %>% 
      # the bootstrapping is finished at this point
      # now work across bootstrap samples at each value
      group_by(value) %>% 
      # to estimate the lower and upper 95% quantiles of predicted mpgs
      summarize(l = quantile(.fitted, .025),
                u = quantile(.fitted, .975),
                .groups = "drop"
      ) %>% 
      arrange(value)
  })) %>% 
  dplyr::select(-model_data) %>% 
  unnest(plot_data) -> tbl_plot_data

ep_deaths_over65 = ggplot()+
  geom_ribbon(aes(x=h3_deaths, ymin = gaussian(link="log")$linkinv(l), ymax = gaussian(link="log")$linkinv(u)), 
              tbl_plot_data %>% filter(name=="NA_bhatt_ep_lag1") %>% rename(h3_deaths=value), 
              alpha = 0.3, fill = "grey") +
  geom_point(data=red_wide%>% filter(season!="2009-2010"),
             aes(x=NA_bhatt_ep_lag1,y=h3_deaths,fill=dom_type2),size=4,pch=21)+
  geom_errorbar(data=red_wide%>% filter(season!="2009-2010"),
                aes(x=NA_bhatt_ep_lag1,
                    ymin=h3_deaths_lwr, ymax=h3_deaths_upr,
                    color=dom_type2), width=.05)+
  theme_bw(base_size = 16)+
  xlab(expression("N2 epitope distance ("~italic(t)~"-1)"))+
  ylab("Excess A(H3N2) Mortality")+
  theme(legend.position = "bottom")+
  background_grid(major = "xy", minor = "none")+
  scale_color_manual(breaks=c("H3","H1","H1pdm","H3/H1pdm"),
                     labels=c("H3","H1","H1pdm","H3/H1pdm"),values=cols,name="Dominant IAV")+
  scale_fill_manual(breaks=c("H3","H1","H1pdm","H3/H1pdm"),
                    labels=c("H3","H1","H1pdm","H3/H1pdm"),values=cols,name="Dominant IAV")+
  geom_smooth(data=red_wide%>% filter(season!="2009-2010"),
              aes(x=NA_bhatt_ep_lag1,y=h3_deaths),
              method = "glm", formula = y~x,
              method.args = list(family = gaussian(link = 'log')),
              se=F,
              size = 1, linetype = 2,color="black")+
  geom_text(y = 26, x =-1.2, aes(label = paste(adj.r.squared,pvalue,sep = "~`,`~")), size=5,
            data = labels, parse = TRUE, hjust = 0)+
  ggtitle(expression(paste('\u2265',"65 years of age")))+
  theme(plot.title = element_text(hjust = 0.5,face = "plain"))+
  ylim(0,29)
ep_deaths_over65


NA_excess_plot = plot_grid(ep_deaths+theme(legend.position = "none"),
                           ep_deaths_over65+theme(legend.position = "none"),
                           nrow=1,labels = c("C","D"))
leg = get_legend(ep_deaths_over65+
                   theme(legend.direction = "horizontal",legend.justification="center" ,
                         legend.position = "bottom"))
both = plot_grid(NA_excess_plot,leg,nrow=2,rel_heights = c(1,0.1))
both

######################################################################################
## combine HA and NA figures
######################################################################################

all_plots = plot_grid(HA_excess_plot,NA_excess_plot,leg,rel_heights = c(1,1,0.1),nrow=3)
all_plots
save_plot(all_plots,filename = "figures/Fig3_sup_fig2_excess_mortality_all_ages_and_elderly_HA_and_NA_ep_distance.png",
          dpi=300,base_width = 10,base_height = 10,bg="white")

