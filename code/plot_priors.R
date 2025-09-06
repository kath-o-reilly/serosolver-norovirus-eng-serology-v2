library(tidyverse)
library(ggplot2)
library(patchwork)
setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/")
source("code/simulate_priors.R")

theme_use <- theme_classic() + theme(axis.text=element_text(size=6),axis.title=element_text(size=8),
                                     legend.title=element_text(size=6),legend.text=element_text(size=6),
                                     strip.text=element_text(size=8),
                                     plot.title=element_text(size=10))

## Exponential model
N <- 1000
priors <- simulate_priors(N)
priors <- priors %>% pivot_wider(names_from=name,values_from=value)
times <- seq(0,10,by=0.1)
traj <- matrix(0, ncol=length(times),nrow=N)
priors_ic50 <- priors[priors$biomarker_group == 1,]
priors_avidity <- priors[priors$biomarker_group == 2,]

simulate_and_plot_prior_traj <- function(N, priors){
for(i in 1:N){
  traj[i,] <- priors$boost_long[i] + priors$boost_short[i]*exp(-times*priors$wane_short[i])
}

traj <- as.data.frame(traj)
colnames(traj) <- times
traj$samp <- 1:N
draws <- traj %>% pivot_longer(-samp)
draws$name <- as.numeric(draws$name)
traj <- traj %>% pivot_longer(-samp) %>% group_by(name) %>% 
  dplyr::summarize(mean=median(value),
                   lower50=quantile(value,0.25),upper50=quantile(value,0.75),
                   lower=quantile(value,0.025),upper=quantile(value,0.975))
traj$name <- as.numeric(traj$name) 
list(traj,draws)
}

traj_ic50 <- simulate_and_plot_prior_traj(N, priors_ic50)
traj_avidity <- simulate_and_plot_prior_traj(N, priors_avidity)
traj <- bind_rows(traj_ic50[[1]] %>% mutate("Assay"="IC50"),traj_avidity[[1]] %>% mutate("Assay"="Avidity"))
draws <- bind_rows(traj_ic50[[2]] %>% mutate("Assay"="IC50"),traj_avidity[[2]] %>% mutate("Assay"="Avidity"))

p_lhs <- ggplot(traj) +
  geom_ribbon(aes(x=name,ymin=lower,ymax=upper,fill=Assay),alpha=0.25) +
  geom_ribbon(aes(x=name,ymin=lower50,ymax=upper50,fill=Assay),alpha=0.5) +
  geom_line(aes(x=name,y=mean,col=Assay)) +
  #geom_line(data=draws %>% filter(samp %in% 1:100), aes(x=name,y=value,group=samp,col=Assay),alpha=0.25)+
  labs(x='Time since boost (years)',y='Antibody level') +
  theme_use +
  scale_x_continuous(breaks=seq(0,10,by=2)) +
  scale_y_continuous(breaks=seq(0,16,by=2)) +
  coord_cartesian(ylim=c(0,16)) +
  facet_wrap(~Assay,ncol=1) +
  theme(legend.position="bottom")+
  ggtitle("Prior median/50% quantiles/95% quantiles")
p_rhs <- ggplot(traj) +
  #geom_ribbon(aes(x=name,ymin=lower,ymax=upper,fill=Assay),alpha=0.25) +
  #geom_ribbon(aes(x=name,ymin=lower50,ymax=upper50,fill=Assay),alpha=0.5) +
  #geom_line(aes(x=name,y=mean,col=Assay)) +
  geom_line(data=draws %>% filter(samp %in% 1:100), aes(x=name,y=value,group=samp,col=Assay),alpha=0.25)+
  labs(x='Time since boost (years)',y='Antibody level') +
  theme_use+
  scale_x_continuous(breaks=seq(0,10,by=2)) +
  scale_y_continuous(breaks=seq(0,16,by=2)) +
  coord_cartesian(ylim=c(0,16)) +
  facet_wrap(~Assay,ncol=1)+
  theme(legend.position="bottom")+
  ggtitle("Prior draws")

p_lhs | p_rhs

ggsave("figures/prior_exponential_traj.png",p_lhs | p_rhs,height=5,width=7)



## Linear model
N <- 1000
priors <- simulate_priors_original_model(N)
priors <- priors %>% pivot_wider(names_from=name,values_from=value)
times <- seq(0,10,by=0.1)
traj <- matrix(0, ncol=length(times),nrow=N)
priors_ic50 <- priors[priors$biomarker_group == 1,]
priors_avidity <- priors[priors$biomarker_group == 2,]

simulate_and_plot_prior_traj_linear <- function(N, priors){
  for(i in 1:N){
    traj[i,] <- priors$boost_long[i] + pmax(0,priors$boost_short[i]*(1-times*priors$wane_short[i]))
  }
  
  traj <- as.data.frame(traj)
  colnames(traj) <- times
  traj$samp <- 1:N
  draws <- traj %>% pivot_longer(-samp)
  draws$name <- as.numeric(draws$name)
  traj <- traj %>% pivot_longer(-samp) %>% group_by(name) %>% 
    dplyr::summarize(mean=median(value),
                     lower50=quantile(value,0.25),upper50=quantile(value,0.75),
                     lower=quantile(value,0.025),upper=quantile(value,0.975))
  traj$name <- as.numeric(traj$name) 
  list(traj,draws)
}

traj_ic50 <- simulate_and_plot_prior_traj_linear(N, priors_ic50)
traj_avidity <- simulate_and_plot_prior_traj_linear(N, priors_avidity)
traj <- bind_rows(traj_ic50[[1]] %>% mutate("Assay"="IC50"),traj_avidity[[1]] %>% mutate("Assay"="Avidity"))
draws <- bind_rows(traj_ic50[[2]] %>% mutate("Assay"="IC50"),traj_avidity[[2]] %>% mutate("Assay"="Avidity"))

p_lhs <- ggplot(traj) +
  geom_ribbon(aes(x=name,ymin=lower,ymax=upper,fill=Assay),alpha=0.25) +
  geom_ribbon(aes(x=name,ymin=lower50,ymax=upper50,fill=Assay),alpha=0.5) +
  geom_line(aes(x=name,y=mean,col=Assay)) +
  #geom_line(data=draws %>% filter(samp %in% 1:100), aes(x=name,y=value,group=samp,col=Assay),alpha=0.25)+
  labs(x='Time since boost (years)',y='Antibody level') +
  theme_use +
  scale_x_continuous(breaks=seq(0,10,by=2)) +
  scale_y_continuous(breaks=seq(0,16,by=2)) +
  coord_cartesian(ylim=c(0,16)) +
  facet_wrap(~Assay,ncol=1) +
  theme(legend.position="bottom")+
  ggtitle("Prior median/50% quantiles/95% quantiles")
p_rhs <- ggplot(traj) +
  #geom_ribbon(aes(x=name,ymin=lower,ymax=upper,fill=Assay),alpha=0.25) +
  #geom_ribbon(aes(x=name,ymin=lower50,ymax=upper50,fill=Assay),alpha=0.5) +
  #geom_line(aes(x=name,y=mean,col=Assay)) +
  geom_line(data=draws %>% filter(samp %in% 1:100), aes(x=name,y=value,group=samp,col=Assay),alpha=0.25)+
  labs(x='Time since boost (years)',y='Antibody level') +
  theme_use +
  scale_x_continuous(breaks=seq(0,10,by=2)) +
  scale_y_continuous(breaks=seq(0,16,by=2)) +
  coord_cartesian(ylim=c(0,16)) +
  facet_wrap(~Assay,ncol=1)+
  theme(legend.position="bottom") +
  ggtitle("Prior draws")

p_lhs | p_rhs

ggsave("figures/prior_linear_traj.png",p_lhs | p_rhs,height=5,width=7)
