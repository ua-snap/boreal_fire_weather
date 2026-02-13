# retrieve estimates of beta coefficient for each ecoregion and each model
# Anna Spiers 
# Oct 2025

library(dplyr)
library(MCMCvis)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)
library(purrr) #map_dfr

### Start with 1 ecoregion and feedback type, then exapnd
# Load data
ecos_size_path <- "data/ancillary/ecoregion_size.csv"
ecos_size <- read.csv(ecos_size_path)
unique_ecos <- ecos_size$ecos

feedback_types <- c("feedback", "no-feedback")

# Best models 
regression_forms <- read.csv("data/ancillary/regression_predictors.csv")

param_est_df <- data.frame()
for (eco in unique_ecos) {
    for (fdbk in feedback_types) {
        # load jags results
        jags_out_path <- file.path("data/model_results/jags-results",
                              paste0("jags-results_Gamma_",fdbk,
                                     "-model_",eco,".RData"))
        load(jags_out_path)
        
        # Visualize traceplots as sanity check
        MCMCtrace(jags_output, 
                  params = c('beta0', 'beta1', 'beta2',
                             "aab_init","shape"), 
                  ISB=FALSE, exact=TRUE, pdf=F)
        
        # look at parameter estaimtes
        jags_summ <- MCMCsummary(jags_output, round = 2)
        # param estimate for aab_init, shape, betas, offsets][i]
        # aab_init: log-normal estimate of AAB in 1980
        # beta0: intercept
        # beta1: BUI coeff
        # beta2: ISI coeff
        # shape: shape parameter for AAB gamma distribution
        # (for feedback only) offset: 
        
        # Summarize estimates
        offsets = MCMCsummary(jags_output, params = 'offset')
        offset_mean = ifelse(fdbk=="feedback",round(mean(offsets$'50%'),digits=2), 0)
        offset_sd = ifelse(fdbk=="feedback",round(sd(offsets$'50%'),digits=2), 0)
        param_est_df = rbind(param_est_df,
                             cbind("eco"=eco,
                             "mdl" = regression_forms %>% 
                                 filter(ecos==eco) %>% pull(mdl),
                             "fdbk"=fdbk,
                             "b0"=jags_summ[row.names(jags_summ)=="beta0",'50%'],
                             "b1_BUI"=jags_summ[row.names(jags_summ)=="beta1",'50%'],
                             "b2_ISI"=jags_summ[row.names(jags_summ)=="beta2",'50%'],
                             "offset_mean"=offset_mean,
                             "offset_sd"=offset_sd,
                             "dev"=jags_summ[row.names(jags_summ)=="deviance",'50%'],
                             "aab_init"=jags_summ[row.names(jags_summ)=="aab_init",'50%'], 
                             "shape"=jags_summ[row.names(jags_summ)=="shape",'50%']))
    }
}

# save as csv
write.csv(param_est_df,"data/dataframes/beta_estimates.csv",
          row.names=F)




# Plot MCMCplot of parameter estimates ------------------------------------

params_plot_df <- data.frame()
for (eco in unique_ecos) {
    for (fdbk in feedback_types) {
        
        # load jags results
        jags_out_path <- file.path("data/model_results/jags-results",
                                   paste0("jags-results_Gamma_",fdbk,
                                          "-model_",eco,".RData"))
        load(jags_out_path)
        
        # Summarize parameter estimates
        mdl = regression_forms %>% 
            dplyr::filter(ecos==eco) %>% pull(mdl)
        jags_summ <- MCMCsummary(jags_output, round=2) %>%
            tibble::rownames_to_column(var="param") %>%
            mutate(param_grp = str_extract(param, "^[^\\[]+")) %>%   
            group_by(param_grp) %>%
            # extract 95% bounds and median
            reframe(med = mean(`50%`, na.rm = TRUE),
                      low_95ci = mean(`2.5%`, na.rm = TRUE),
                      hi_95ci = mean(`97.5%`, na.rm = TRUE)) %>%
            filter(param_grp %in% c("beta0","beta1","beta2","offset")) %>%
            # add identifying columns
            dplyr::mutate(eco=eco, 
                          fdbk=fdbk, 
                            mdl=mdl)
        # combine into df
        params_plot_df = rbind(params_plot_df,jags_summ)
    }
}


# --------------------------------------------------------------
# 4) Clean / factor variables ----------------------------------
# --------------------------------------------------------------

# Reorder ecoregions
eco_order <- params_plot_df %>%
    filter(param_grp == "beta1") %>%               
    group_by(eco) %>%                         
    reframe(beta0_med = mean(med, na.rm = TRUE)) %>%
    arrange(beta0_med) %>%                       
    pull(eco) 

# Nice names for the coefficients (optional)
params_plot_df <- params_plot_df %>%
    mutate(parameter = factor(param_grp,
                              levels = c("beta0",  "beta2","beta1", "offset"),
                              labels = c("Intercept",
                                         "\u03B2_ISI",
                                         "\u03B2_BUI",
                                         "Offset"))) %>%
    # Give the two models readable names and a fixed colour order
    mutate(fdbk_text = case_match(fdbk,
                          "feedback"    ~ "Fire-fuel feedback",
                          "no-feedback"  ~ "No feedback")) %>%
    mutate(fdbk_text = factor(fdbk_text,
                          levels = c("Fire-fuel feedback",
                                     "No feedback")))

# Plot
all_reg_estimates <- ggplot(params_plot_df  %>%
           mutate(eco = factor(eco, levels = eco_order)), 
       aes(y=factor(eco),colour=fdbk_text)) +
    # 95 % CI and median
    geom_errorbar(aes(xmin = low_95ci, xmax = hi_95ci),
                  orientation="y",height = 0.2,
                   position = position_dodge(width = 0.6)) +
    geom_point(aes(x=med), position = position_dodge(width = 0.6),
               size = 1.5) +
    # reference line
    geom_vline(data = filter(params_plot_df, parameter != "Offset"),
               aes(xintercept = 0), linetype = "dotted", colour = "grey60",
               inherit.aes = FALSE  ) +
    facet_grid( col=vars(parameter), scales = "free_x") +
    # facet_wrap(vars(eco),
    #            ncol = 4)) +   
    # Colours – pick any palette you like
    scale_colour_manual(values = c("No feedback" = "#d95f02",
                                   "Fire-fuel feedback"   = "#1b9e77")) +
    labs(x = "Posterior Median (95% CI)",
         y = "Ecoregion ID", colour = "Model") +
    theme_bw(base_size = 12) +
    theme(
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom",
        axis.text.y = element_text(size = 10)
    )
ggsave(filename = "figures/rq2/reg_estimates.png",
       plot = all_reg_estimates,
       width = 7, height = 6, dpi = 300, bg = "white")

# Plot - just betas + mdl predictors
betas_only <- ggplot(params_plot_df %>% 
           filter(param_grp %in% c("beta1","beta2")) %>%
           mutate(eco = factor(eco, levels = eco_order)), 
       aes(y=factor(eco),colour=fdbk_text)) +
    # 95 % CI and median
    geom_errorbar(aes(xmin = low_95ci, xmax = hi_95ci),
                  orientation="y",height = 0.2,
                  position = position_dodge(width = 0.6)) +
    geom_point(aes(x=med), position = position_dodge(width = 0.6),
               size = 1.5) +
    # reference line
    geom_vline(aes(xintercept = 0), linetype = "dotted", colour = "grey60",
               inherit.aes = FALSE  ) +
    facet_grid( col=vars(parameter))+  
    # Colours – pick any palette you like
    scale_colour_manual(values = c("No feedback" = "#d95f02",
                                   "Fire-fuel feedback"   = "#1b9e77")) +
    labs(x = "Posterior Median (95% CI)",
         y = "Ecoregion ID", colour = "Model") +
    theme_bw(base_size = 12) +
    theme(
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom",
        axis.text.y = element_text(size = 10)
    )
ggsave(filename = "figures/rq2/betas_only.png",
       plot = betas_only,
       width = 7, height = 6, dpi = 300, bg = "white")
# EDITS
# why is offset so large (shouldn't aab be smaller for feedbacks model, so offset would be negative?')
# 
# add in banner (edit facet)
# so what? what does this figure add? not much difference between feedback types in MODEL, but there is a big difference between their climate projections
# plot panel with temp/precip

ggsave(filename = "figures/ecoregion_model_summary.png",
       plot = p,
       width = 14, height = 12, dpi = 300, bg = "white")


