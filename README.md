# H3N2_Antigenic_Epi
Code and data to reproduce the results and figures in Perofsky _et al._ 2023. "Antigenic drift and subtype interference shape A(H3N2) epidemic dynamics in the United States" [_medRxiv_ pre-print](https://doi.org/10.1101/2023.10.02.23296453)

## Abstract

Influenza viruses continually evolve new antigenic variants, through mutations in epitopes of their major surface proteins, hemagglutinin (HA) and neuraminidase (NA). Antigenic drift potentiates the reinfection of previously infected individuals, but the contribution of this process to variability in annual epidemics is not well understood. Here we link influenza A(H3N2) virus evolution to regional epidemic dynamics in the United States during 1997—2019. We integrate phenotypic measures of HA antigenic drift and sequence-based measures of HA and NA fitness to infer antigenic and genetic distances between viruses circulating in successive seasons. We estimate the magnitude, severity, timing, transmission rate, age-specific patterns, and subtype dominance of each regional outbreak and find that genetic distance based on broad sets of epitope sites is the strongest evolutionary predictor of A(H3N2) virus epidemiology. Increased HA and NA epitope distance between seasons correlates with larger, more intense epidemics, higher transmission, greater A(H3N2) subtype dominance, and a greater proportion of cases in adults relative to children, consistent with increased population susceptibility. Based on random forest models, A(H1N1) incidence impacts A(H3N2) epidemics to a greater extent than viral evolution, suggesting that subtype interference is a major driver of influenza A virus infection dynamics, presumably via heterosubtypic cross-immunity.

Data processing and statistical analyses are performed with the statistical computing software [R](https://www.r-project.org/) (version 4.3.0). The phylogenetics workflow can be found at the GitHub repository [blab/perofsky-ili-antigenicity](https://github.com/blab/perofsky-ili-antigenicity). Key outputs from [blab/perofsky-ili-antigenicity](https://github.com/blab/perofsky-ili-antigenicity) are in the `2_Phylo_Dataset` folder, so it is not necessary to run the phylogenetic analysis before running the code in this repository.


## R scripts are split into 5 chunks:

* `1_Epi_Dataset` folder
  * `1_cdc_hhs_level_ili_viral_surv_df.R`:
    * Pull and compile syndromic and virologic data from [CDC FluView](https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html) to estimate HHS region-level type and subtype-specific incidences for seasons 1997-1998 to 2018-2019.
    * Incidences are calculated my multiplying the proportion of outpatient encounters for influenza-like illness (weighted by regional population size) by the proportion of respiratory specimens testing positive for influenza A(H3N2), A(H1N1), or B.
  * `2_cdc_virology_surv_interp_smooth_and_onset_estimates.R`:
    * Interpolate missing values and smooth incidence time series.
    * Estimate the timing of epidemic onsets by fitting piecewise linear models to smoothed incidence curves.
  * `3_cdc_ili_burden_metrics_hhs_regions.R`:
    * Correct for differences in sampling effort between pre- and post-2009 seasons and across HHS regions.
    * For each region and season, estimate type/subtype distribution and age-specific ILI case patterns.
    * For each region and season, estimate type/subtype-specific epidemic size, peak incidence, epidemic intensity (inverse Shannon entropy of the weekly incidence distribution), peak week, and seasonal duration.
  * `4_cdc_ili_Rt_hhs_regions.R`:
    * Use the [Epidemia R package](https://imperialcollegelondon.github.io/epidemia/index.html) to estimate regional A(H3N2) time-varying effective reproduction numbers (effective Rt) during each season.
    * Calculates the maximum Rt and initial Rt (mean Rt across the weeks spanning epidemic onset to epidemic peak) for each season.

* `2_Phylo_Dataset` folder
  * `HA_manual` and `NA_manual` scripts (1-9) estimate antigenic and genetic distances between viruses circulating during 1997-1998 and 1996-1997 (one-season lag) and during 1997-1998 and 1995-1996 (two-season lag), which are not included in the output of the phylogenetics workflow [blab/perofsky-ili-antigenicity](https://github.com/blab/perofsky-ili-antigenicity).
  * `10_LBI_diversity_calcs.R`: Calculate the Shannon entropy of HA and NA local branching index (LBI) values during each season.
  * Folder `auspice_tables`: Sequence-level evolutionary fitness measurements from [blab/perofsky-ili-antigenicity](https://github.com/blab/perofsky-ili-antigenicity).
  * Folder `distance_tables`: Mean antigenic and genetic distances between viruses circulating in successive seasons from [blab/perofsky-ili-antigenicity](https://github.com/blab/perofsky-ili-antigenicity).

* `3_Epi_Antigenic_Univariate_Analyses` folder
  * `1_make_ili_and_antigenic_dataset.R`: Combine seasonal epidemic metrics and evolutionary indicators into one dataset.
  * `2_ILI_subtype_time_series_fig1.R`: Make Figure 1 showing regional influenza type and subtype specific incidences from 1997 to 2019.
  * `3_predictors_H3_subtype_dom.R`: Correlations between A(H3N2) viral evolution and A(H3N2) subtype dominance.
  * `4_predictors_H3_epi_metrics.R`: Correlations between A(H3N2) viral evolution and A(H3N2) epidemic size, peak incidence, epidemic intensity, and effective Rt.
  * `5_predictors_H3_age_prop.R`: Correlations between A(H3N2) viral evolution and age-specific ILI case patterns.
  * `6_predictors_H3_epi_timing.R`: Correlations between A(H3N2) viral evolution and A(H3N2) epidemic onset and peak timing, epidemic speed (e.g., seasonal duration, days from onset to peak), and spatiotemporal synchrony across regions.
  * `7_h3_epi_metrics_vs_h1_and_b.R`: Correlations between A(H3N2) epidemic metrics and A(H1N1) and B epidemic size.
  * `8_evol_indicators_correlations.R`: Pairwise correlations between A(H3N2) evolutionary indicators.
  * `9_epi_metric_correlations.R`: Pairwise correlations between A(H3N2), A(H1N1), and B epidemic metrics.

* `4_Wavelet_Analysis` folder
  * `h3_vs_h1_wavelet_coherence.R`: Supplementary wavelet analysis that compares the relative timing of influenza A(H3N2), A(H1N1), and B epidemics during each season. This script sources functions in `WaveletPkg.R` located in the `Wavelets` subfolder.

* `5_Variable_Selection_Analyses` folder
  * Scripts 1-5 run conditional inference random forest models and LASSO regression models predicting A(H3N2) epidemic metrics in each region and season, including peak incidence, epidemic size, epidemic intensity, effective Rt, and subtype dominance. Model features include viral evolutionary indicators, co-circulation of other influenza types/subtypes, proxies for prior immunity, and vaccine-related parameters.
  * `6_top_features_regression_model.R`: Extract the top 10 ranked predictors from random forest models and use model selection (BIC) to determine the best fit "minimal" linear regression model, allowing each candidate model to include up to 3 covariates.
  * `7_variable_importance_plots.R`: Plot variable importance rankings from random forest and LASSO models. Plot observed versus predicted values for each region/season and estimate seasonal RMSE of random forest predictions.

## Other folders

* The `figures` folder contains figures created with the analysis scripts.

* The `data` folder contains all data inputs and outputs, with the exception of evolutionary fitness measurements, which are in the `2_Phylo_Dataset` folder.
