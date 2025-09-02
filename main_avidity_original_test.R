######################################################
## SIMULATION RECOVERY TEST -- NOROVIRUS DATA WITH ANTIBODY TITER AND BINDING AVIDITY DATA
## Author: James Hay
## Date: 10 Jun 2025
## Summary: simulates some serosurvey data with two "observation types" and fits serosolver

rm(list=ls())

## Load in libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(coda)
library(tidyverse)

## Load serosolver locally, or load as proper package after installation
#devtools::load_all("~/Documents/GitHub/serosolver")
#install.packages("~/Documents/GitHub/serosolver", repos=NULL, type="source")
#devtools::install_github("seroanalytics/serosolver")
library(serosolver)

## Set working directly locally and create a directory to store the MCMC chains in
setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/local_data/norovirus_test")

run_name <- "real_one_avidity"
main_wd <- "~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/local_data/norovirus_true"
chain_wd <- paste0(main_wd,"/chains/",run_name)
save_wd <- paste0(main_wd,"/figures/",run_name)

if(!dir.exists(save_wd)) dir.create(save_wd,recursive = TRUE)
if(!dir.exists(chain_wd)) dir.create(chain_wd,recursive = TRUE)

cart <- "Debbink" #"Sim"

## Simulation parameters
mcmc_pars <- c("save_block"=10000,"thin"=10,"thin_hist"=5,
               "iterations"=200000,
               "adaptive_period"=10000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0.05,
               "year_swap_propn"=0.8,"swap_propn"=0.5,
               "inf_propn"=0.5,"hist_sample_prob"=1,"move_size"=3, "hist_opt"=0,
               "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=TRUE)

#antibody_data <- read_csv(paste0("/Users/lsh1603970/Documents/GitHub/serosolver-norovirus-eng-serology-v2/local_data/norovirus_test/figures/sim_two_biomarker/","sim_two_biomarker_titre_data.csv"))
antibody_data <- read_csv(paste0("/Users/lsh1603970/Documents/GitHub/serosolver-norovirus-eng-serology-v2/Data/","titre_dat_norovirus_age_avidity.csv"))

## IMPORTANT to set this to a data frame, otherwise can cause some errors later
antibody_data <- as.data.frame(antibody_data) %>% select(-age,-Date)
#antibody_data <- antibody_data[order(antibody_data$individual,antibody_data$biomarker_group,antibody_data$biomarker_id),]

antibody_data_keep <- antibody_data

antibody_data <- antibody_data %>% filter(biomarker_group == 2)
antibody_data$biomarker_group <- 1
#tmp <- antibody_data
#tmp$biomarker_group <- 2
#antibody_data <- rbind(antibody_data,tmp)

# shake up the order...(does this break the sim?).. nope
#antibody_data <- antibody_data[order(antibody_data$biomarker_group,antibody_data$measurement),]

n_obs_types <- 1#2 ## Number of observation types (eg. antibody titer and binding avidity)
obs_type_dists <- 2#c(2,2) ## Vector of observation models for each observation type -- use 1 for discretized normal and 2 for continuous, truncated normal. I assume that titers are discrete and avidity is continuous
#n_indivs <- max(antibody_data$individual) ## Number of individuals for the simulation
#n_samps <- 2 ## Number of samples per person
#repeats <- 1 ## Number of repeat measurements per variant/sample combination
#samp_min <- 2009 ## First sample year
#samp_max <- 2012 ## Final sample year
#year_min <- 2000 ## First year of possible circulation (ie. time 0 of the simulation)
#year_max <- 2012 ## Final year of possible circulation
#age_min <- 2 ## Age minimum and maximum in years, simulated from a uniform distribution
#age_max <- 10

################################ 
## IMPORTANT
################################ 
## Set to FALSE to remove the measurement offset parameters
## Set to TRUE to assume that each antigen/virus has a measurement offset parameter for each observation type
use_measurement_offset_parameters <- FALSE #TRUE #

## Viruses and times for samples
#sampled_viruses <- c(2002,2006,2009,2012)
#sampling_times <- seq(samp_min, samp_max, by=1)

################################ 
## IMPORTANT
################################ 

if(cart == "Sim"){
  ## Create a fake antigenic map -- can put in what you like here
  antigenic_coords <- data.frame(Strain=sampled_viruses,X=c(0.1,0.5,3,3.5,4),Y=c(0.1,2,1,3,4))
  antigenic_map <- generate_antigenic_map_flexible(antigenic_coords,year_max=year_max+1,year_min=year_min,spar = 0.1)
  ggplot(antigenic_map) + geom_line(aes(x=x_coord,y=y_coord,col=inf_times)) +
    geom_point(data=antigenic_coords,aes(x=X,y=Y))
  
  ## This doubles the antigenic map, one entry for each biomarker_group (i.e., one for titers one for avidity)
  antigenic_map <- setup_antigenic_map(antigenic_map,n_biomarker_groups=n_obs_types,unique_biomarker_groups=1:n_obs_types)[[1]]
}
if(cart == "Debbink"){
  antigenic_map <- read_csv(paste0(save_wd,"/","antigenic_map_hand.csv"))
  p1 <- ggplot(antigenic_map,aes(x=x_coord,y=y_coord,col=inf_times,label=inf_times)) + geom_line() +
    geom_point(,cex=2) +
    geom_label() +
    #geom_point(data=antigenic_coords,aes(x=X,y=Y,col=inf_times),pch=0) +
    ylim(-3,4)
  p1
  antigenic_map <- setup_antigenic_map(antigenic_map,n_biomarker_groups=n_obs_types,unique_biomarker_groups=1:n_obs_types)[[1]]
  
}

## Get vector of times individuals can be infected
#possible_exposure_times <- unique(antigenic_map$inf_times)
#n_times <- length(possible_exposure_times)

## Set up parameter table
par_tab <- read.csv("par_tab.csv",stringsAsFactors=FALSE)



################################ 
## IMPORTANT
################################ 
## This allows you to change prior on the probability of infection per year 
## These are the parameters of a beta distribution -- see rbeta(10000, infection_model_prior_shape1,infection_model_prior_shape2)
par_tab[par_tab$names %in% c("infection_model_prior_shape1","infection_model_prior_shape2"),c("values")] <- c(1,10) ## Can also try c(1,1), or something informative. Just the parameters of a beta distribution which acts as the prior on the per-time attack rate.
################################ 

## Just some setup of the parameter table to get parameters vaguely like you showed me
par_tab$fixed <- 1
par_tab[par_tab$names %in% c("boost_long","antigenic_seniority","cr_long","obs_sd","boost_short","wane_short","cr_short"),"fixed"] <- 0
par_tab[par_tab$names =="cr_long","values"] <- 0.4
#par_tab[par_tab$names =="cr_long","lower_start"] <- 0.5. # needed for avidity
#par_tab[par_tab$names =="cr_long","upper_start"] <- 0.8 # needed for avidity
par_tab[par_tab$names =="cr_short","values"] <- 0.1

## Extends the parameter table for multiple observation types

par_tab <- extend_par_tab_biomarker_groups(par_tab,n_obs_types)
if(use_measurement_offset_parameters){
  tmp <- add_rhos_par_tab(par_tab, sampled_viruses, n_obs_types)
  par_tab <- tmp[[1]]
  measurement_indices <- tmp[[2]]
} else {
  measurement_indices <- NULL
}
par_tab[par_tab$biomarker_group == 2 & par_tab$names %in%c("boost_long","boost_short","obs_sd"),"values"] <- c(0.7,1.2,0.75)

################################ 
## IMPORTANT
################################ 
## Set the upper bound of observations for the two observation types (i.e., max observable value for titers and avidity)
par_tab[par_tab$biomarker_group == 1 & par_tab$names == "max_measurement","values"] <- 10
par_tab[par_tab$biomarker_group == 2 & par_tab$names == "max_measurement","values"] <- 10

## Simulate realistic-ish attack rates
#attack_rates <- simulate_attack_rates(infection_years=possible_exposure_times,mean_par=0.15,sd_par=0.5,n_groups =1)
## Plot attack rates
#ggplot(attack_rates) + geom_line(aes(x=time,y=prob_infection,col=as.factor(population_group))) 


################################ 
## SIMULATE DATA
################################ 
# all_simulated_data <- simulate_data(par_tab=par_tab, 
#                                     n_indiv=n_indivs, ## Number of individuals for the overall simulation
#                                     measured_biomarker_ids = sampled_viruses, 
#                                     sampling_times=sampling_times, 
#                                     age_min=age_min,age_max=age_max,
#                                     nsamps=n_samps,
#                                     antigenic_map=antigenic_map, 
#                                     attack_rates=attack_rates, 
#                                     repeats=repeats, 
#                                     missing_data = 0, ## Important -- sets a proportion of observations to be missing at random
#                                     data_type=obs_type_dists, 
#                                     demographics=NULL,
#                                     measurement_bias = measurement_indices)
## Save titre data
# write_csv(all_simulated_data$antibody_data, file=paste0(save_wd,"/",run_name,"_titre_data.csv"))
# ## Save parameter table
# write_csv(par_tab, file=paste0(save_wd,"/",run_name,"_par_tab.csv"))
# ## Save attack rates
# write_csv(all_simulated_data$attack_rates, file=paste0(save_wd,"/",run_name,"_attack_rates.csv"))
# ## Save infection histories
# write_csv(as.data.frame(all_simulated_data$infection_histories), file=paste0(save_wd,"/",run_name,"_infection_histories.csv"))


## Have a look at the simulated data
#plot_antibody_data(antibody_data,possible_exposure_times,9,study_design = "longitudinal",infection_histories = all_simulated_data$infection_histories) + facet_grid(individual~biomarker_group)
#plot_antibody_data(antibody_data %>% filter(biomarker_group == 1),possible_exposure_times,5,study_design = "cross-sectional",infection_histories = all_simulated_data$infection_histories) + facet_grid(individual~sample_time)

#antibody_data$measurement <- jitter(antibody_data$measurement,amount=5)
#antibody_data$measurement[antibody_data$measurement<0] <- 0

ggplot(antibody_data,aes(x=(measurement+1))) + 
  geom_histogram() +
  scale_x_continuous(trans='log2') +
  facet_grid(rows = vars(biomarker_group),cols = vars(biomarker_id))

ggplot(antibody_data %>% filter(biomarker_group==2),aes(x=birth,y=measurement)) + geom_point() +
  facet_wrap(~biomarker_id) + geom_smooth()

################################ 
## IMPORTANT
################################ 
## Set values for the model parameter priors.
## You would need to add your own priors here as new lines and values added to the final return value
## Note that to set a prior for a parameter, you need to get the right index in the par_tab table
prior_func <- function(par_tab){
  ## Finds any stratification coefficient parameters
  coef_pars <- which(par_tab$names %like% "coef")
  #rho_pars <- which(par_tab$names %like% "rho")
  #par_names <- par_tab$names
  pars <- par_tab$names
  f <- function(pars){
    prior_p <- 0
    ## You can add your own priors on the model parameters here, for example, this places a log-normal prior with mean 2 on the boost_short parameter:
    #p1 <- sum(dlnorm(pars[which(names(pars) == "boost_long")],log(4), 0.5,log=TRUE)) # Looser prior
    p1 <- sum(dlnorm(pars[which(names(pars) == "boost_long")],log(4), 1,log=TRUE)) # Looser prior
    #p2 <- sum(dlnorm(pars[which(names(pars) == "boost_short")],log(2), 0.5,log=TRUE)) # Looser prior
    p2 <- sum(dlnorm(pars[which(names(pars) == "boost_short")],log(2), 1,log=TRUE)) # Looser prior
    p3 <- sum(dbeta(pars[which(names(pars) == "wane_short")],2, 2,log=TRUE) )
    p4 <- sum(dbeta(pars[which(names(pars) == "antigenic_seniority")],1, 1,log=TRUE))
    #p4 <- sum(dbeta(pars[which(names(pars) == "antigenic_seniority")],2, 2,log=TRUE))
    p5 <- sum(dbeta(pars[which(names(pars) == "cr_long")],1, 5,log=TRUE) ) # Skewed toward low cross-reactivity
    p6 <- sum(dbeta(pars[which(names(pars) == "cr_short")],1, 5,log=TRUE))  # Skewed toward low cross-reactivity
    p7 <- sum(dnorm(pars[which(names(pars) == "obs_sd")],0, 100,log=TRUE)) # No change
    #p_rhos <- sum(dnorm(pars[rho_pars],0, 3,log=TRUE) )
    prior_p <- prior_p + p1 + p2 + p3 + p4 + p5 + p6 + p7 #+ p_rhos
    
    ## ==================================================
    sum(prior_p)
  }
}

################################ 
## RUN SEROSOLVER
################################ 
res <- serosolver(par_tab, 
                  antibody_data, 
                  #possible_exposure_times = possible_exposure_times,
                  antigenic_map=antigenic_map,
                  prior_func=prior_func,
                  filename=paste0(chain_wd,"/sim"),
                  n_chains=5, ## Run 3 chains
                  parallel=TRUE, ## Run in parallel
                  mcmc_pars=mcmc_pars, 
                  verbose=TRUE,
                  data_type=obs_type_dists,
                  measurement_bias= measurement_indices
) 

## Save plots from serosolver
res$all_diagnostics$p_thetas[[1]]
res$all_diagnostics$p_thetas[[2]]

## Save parameter estimates
ggsave(paste0(save_wd,"/theta_traceplot.pdf"),res$all_diagnostics$p_thetas[[1]],height=7,width=8,units="in",dpi=300)
ggsave(paste0(save_wd,"/theta_densities.pdf"),res$all_diagnostics$p_thetas[[2]],height=7,width=8,units="in",dpi=300)
ggsave(paste0(save_wd,"/rho_traceplot.pdf"),res$all_diagnostics$p_thetas[[5]],height=7,width=8,units="in",dpi=300)
ggsave(paste0(save_wd,"/rho_densities.pdf"),res$all_diagnostics$p_thetas[[6]],height=7,width=8,units="in",dpi=300)
write_csv(as.data.frame(res$all_diagnostics$theta_estimates),paste0(save_wd,"/estimated_parameters.csv"))

## Save attack rate/infection history estimates
ggsave(paste0(save_wd,"/indiv_infections.pdf"),res$all_diagnostics$p_inf_hists$indiv_infections,height=4,width=6,units="in",dpi=300)
ggsave(paste0(save_wd,"/attack_rate_traceplot.pdf"),res$all_diagnostics$p_inf_hists$indiv_infections$by_time_trace[[1]],height=4,width=6,units="in",dpi=300)
ggsave(paste0(save_wd,"/attack_rate_density.pdf"),res$all_diagnostics$p_inf_hists$indiv_infections$by_time_trace[[2]],height=4,width=6,units="in",dpi=300)

write_csv(as.data.frame(res$all_diagnostics$inf_hist_estimates$by_year),paste0(save_wd,"/estimated_attack_rates.csv"))


## Post-hoc plots
chains <- load_mcmc_chains(chain_wd,par_tab,burnin=mcmc_pars["adaptive_period"],estimated_only = TRUE)

## Plot parameter estimates vs. true values
estimated_par_names <- par_tab[par_tab$fixed == 0 & !(par_tab$names %like% "rho"),"names"]
estimated_par_names <- make.names(estimated_par_names,unique=TRUE)


################################ 
## IMPORTANT
################################
## The vertical lines are for the true values of the simulation.
## For real data, comment out the geom_vline
## Plot theta estimates by biomarker group
p_theta <- chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no)) %>%
  filter(name %in% estimated_par_names) %>%
  dplyr::rename(names_new=name,est=value) %>%
  left_join(par_tab %>% filter(fixed == 0& !(par_tab$names %like% "rho")) %>% mutate(names_new=estimated_par_names) %>% select(names,names_new,values,biomarker_group) %>%
              mutate(biomarker_group = as.factor(biomarker_group))) %>%
  ggplot() + geom_density(aes(x=est,fill=biomarker_group),alpha=0.5) + 
  geom_vline(aes(xintercept=values,linetype="True value"),col="black") +
  scale_fill_manual(name="Biomarker group",values=c("1"="blue","2"="red")) +
  scale_linetype_manual(name="",values=c("True value"="dashed")) +
  facet_wrap(~names_new,scales="free",ncol=3) +
  xlab("Value") +
  ylab("Density") +
  theme_minimal() +
  theme(strip.text=element_text(size=6))
ggsave(paste0(save_wd,"/par_estimates.pdf"),p_theta,height=8,width=7,units="in",dpi=300)

## Plot rho estimates by biomarker group
if(!is.null(measurement_indices)){
estimated_par_names <- par_tab[par_tab$fixed == 0 & par_tab$names %like% "rho","names"]
estimated_par_names <- make.names(estimated_par_names,unique=TRUE)

p_rhos <- chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no)) %>%
  filter(name %in% estimated_par_names) %>%
  dplyr::rename(names_new=name,est=value) %>%
  left_join(par_tab %>% filter(fixed == 0 & par_tab$names %like% "rho") %>% mutate(names_new=estimated_par_names) %>% select(names,names_new,values,biomarker_group) %>%
              mutate(biomarker_group = as.factor(biomarker_group))) %>%
  ggplot() + geom_density(aes(x=est,fill=biomarker_group),alpha=0.5) + 
  geom_vline(aes(xintercept=values,linetype="True value"),col="black") +
  scale_fill_manual(name="Biomarker group",values=c("1"="blue","2"="red")) +
  scale_linetype_manual(name="",values=c("True value"="dashed")) +
  facet_wrap(~names_new,scales="free",ncol=3) +
  xlab("Value") +
  ylab("Density") +
  theme_minimal() +
  theme(strip.text=element_text(size=6))
ggsave(paste0(save_wd,"/rho_estimates.pdf"),p_rhos,height=6,width=7,units="in",dpi=300)

}

## Post-hoc plots
chains <- load_mcmc_chains(chain_wd,par_tab,burnin=mcmc_pars["adaptive_period"],estimated_only = FALSE)

## Plot model predictions against observations
################################ 
## IMPORTANT
################################
## Set the known_infection_history argument to NULL for the real data
p_fits <- plot_model_fits(chains$theta_chain,chains$inf_chain,
                          settings=res$settings,
                                  individuals = 1:5,
                                  known_infection_history=as.matrix(all_simulated_data$infection_histories),
                                  orientation="cross-sectional",expand_to_all_times = FALSE)

ggsave(paste0(save_wd,"/model_fits.pdf"),p_fits,height=7,width=8,units="in",dpi=300)

################################ 
## IMPORTANT
################################
## Set the real_inf_hist argument to NULL for the real data
## Plot estimated infections against truth
p_inf_comb <- plot_cumulative_infection_histories(chains$inf_chain,indivs=1:25,
                                                  real_inf_hist=as.matrix(all_simulated_data$infection_histories),
                                                  possible_exposure_times=possible_exposure_times,number_col=5,pad_chain=TRUE,nsamp=100)
## Densities
p_inf_density <- p_inf_comb[[2]] + ylab("Posterior probability of infection") + 
  scale_x_continuous(breaks=seq(2000,2012,by=4),labels=seq(2000,2012,by=4)) +
  theme(axis.text.x=element_text(size=7,angle=45,hjust=1))+ xlab("Date (years)")
## Cumulative infections
p_inf_cumu <- p_inf_comb[[1]] + ylab("Cumulative number of infections") + 
  scale_x_continuous(breaks=seq(2000,2012,by=4),labels=seq(2000,2012,by=4)) +
  theme(axis.text.x=element_text(size=7,angle=45,hjust=1)) + xlab("Date (years)")


ggsave(paste0(save_wd,"/inf_hists_density.pdf"),p_inf_density,height=7,width=8,units="in",dpi=300)
ggsave(paste0(save_wd,"/inf_hists_cumulative.pdf"),p_inf_cumu,height=7,width=8,units="in",dpi=300)

## Look at model-predicted observations against actual observations
antibody_predictions <- plot_antibody_predictions(chains$theta_chain,chains$inf_chain,settings=res$settings)
print(antibody_predictions$proportion_correct)
p_preds <- antibody_predictions$p_pointrange
ggsave(paste0(save_wd,"/predicted_observations.pdf"),p_preds,height=4,width=7,units="in",dpi=300)

################################ 
## IMPORTANT
################################
## Set the true_ar argument to NULL for the real data
p_ar <- plot_attack_rates(chains$inf_chain,settings = res$settings,
                          true_ar=all_simulated_data$attack_rates, ## Set this to NULL for real data
                          by_group=FALSE, 
                          plot_den = TRUE,
                          pad_chain=TRUE)
ggsave(paste0(save_wd,"/attack_rates.pdf"),p_ar,height=4,width=7,units="in",dpi=300)

