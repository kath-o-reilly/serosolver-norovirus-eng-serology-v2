## Compare attack rates scenarios

## Plot attack rates and compare to real data
library(tidyverse)
library(ggplot2)
library(serosolver)
library(data.table)
library(patchwork)

setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/")
source("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/code/adjust_age_distribution.R", echo=TRUE)

## Read in attack rates from ic50 Debbink scenario
serosolver_ar <- read_csv("results_linear/ic50_kendra_attack_rates.csv")
inf_chain_debbink_exp <- load_infection_chains(location="local_data_exponential/norovirus_true/chains/ic50_real/",burnin = 300000)$chain
inf_chain_kendra_exp <- load_infection_chains(location="local_data_exponential/norovirus_true/chains/ic50_kendra/",burnin = 300000)$chain
inf_chain_debbink_linear <- load_infection_chains(location="local_data_linear/norovirus_true/chains/ic50_real/",burnin = 300000)$chain
inf_chain_kendra_linear <- load_infection_chains(location="local_data_linear/norovirus_true/chains/ic50_kendra/",burnin = 300000)$chain

inf_chains <- list(inf_debbink_exp=inf_chain_debbink_exp,
                   inf_kendra_exp=inf_chain_kendra_exp,
                   inf_debbink_linear=inf_chain_debbink_linear,
                   inf_kendra_linear=inf_chain_kendra_linear)

## Read in our data
antibody_data <- read_csv("Data/titre_dat_norovirus_age_avidity.csv")
age_distribution <- read_csv("Data/age_distribution_england.csv")

file_ver <- "debbink_ic50"
## Read in UKHSA reports per capita for children
real_dat <- read_csv("Data/noro_ukhsa_children.csv") %>% filter(year > 2000)

## ggplot 2 options
theme_use<- theme_classic() + theme(axis.text=element_text(size=6),
                                    axis.title=element_text(size=8),
                                    legend.title=element_text(size=6),legend.text=element_text(size=6),
                                    plot.tag=element_text(size=10,face="bold"))

colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Novovirus transition times
start_years <- c(2000,2002,2004,2006,2009,2012)
end_years <- c(2002,2004,2006,2009,2012,2013)
strains <- c("Grimsby 1996","Framington Hills\n 2002","Hunter 2004", "Den Haag 2006","New Orleans 2009","Sydney\n 2012")
noro_epi <- data.frame(strain=strains,start=start_years,end=end_years)

norovirus_key <- c("2002"="FH-2002","2006"="DH-2006","2009"="NO-2009","2012"="SY-2012")


## Get age distribution of data and ONS
## Read in age distribution data from ONS
age_distribution <- age_distribution %>% select(-Sex) %>% pivot_longer(-Age) %>%
  group_by(name) %>% dplyr::mutate(weighted=value/sum(value)) %>%
  mutate(year = str_split(name,"Mid-")[[1]][2]) %>% mutate(year=as.integer(year)) %>%
  filter(year >= 2001, year <= 2012) %>%
  rename(time=year)

data_age_dist <- antibody_data %>% 
  select(individual,birth,sample_time) %>% 
  group_by(individual,birth) %>%
  filter(sample_time == max(sample_time)) %>% ## Note latest sample time -- individuals don't contribute to population after this time
  distinct() %>% 
  expand_grid(time = seq(2001,2012,by=1)) %>% 
  mutate(age=time-birth) %>% 
  filter(age >=0) %>% 
  filter(time <= sample_time) %>%
  group_by(time,age) %>% tally() %>%
  arrange(time,age) %>%
  group_by(time) %>% mutate(prop=n/sum(n))

resampled_all <- NULL
key <- c("Debbink, exponential","Kendra, exponential","Debbink, linear","Kendra, linear")
for(i in seq_along(inf_chains)){
  ## For each read in chain, re-sample ages to get age-adjusted representation
  ## Read in the main serosolver infection history chains
  inf_chain <- pad_inf_chain(inf_chains[[i]],indivs = unique(antibody_data$individual))        
  colnames(inf_chain) <- c("individual","time","infected","samp_no","chain_no")
  times <- seq(2001,2012,by=1)
  inf_chain$time <- inf_chain$time - min(inf_chain$time) + 1
  inf_chain$time <- times[inf_chain$time]
  
  ## Merge with age data
  ages <- antibody_data %>% select(individual,birth) %>% distinct() %>% expand_grid(time=times) %>% mutate(age=time-birth)
  inf_chain <- inf_chain %>% left_join(ages) %>% filter(age >=0)
  inf_chain <- inf_chain %>% mutate(age_group_simple = if_else(age < 5,as.character(age),"5+"))
  
  ## Merge with sample data, as we don't count individuals after they've left the study
  sample_times <- antibody_data %>% select(individual,sample_time) %>% group_by(individual) %>% filter(sample_time == max(sample_time)) %>% ungroup() %>% distinct()
  inf_chain <- inf_chain %>% left_join(sample_times) %>% filter(time <= sample_time)
  
  ######################################################
  ## For each sample number/chain number, resample individuals by age for each year to get the same age distribution as the ONS census data
  data_age_dist <- data_age_dist %>% group_by(time) %>% mutate(total_n=sum(n))
  inf_chain <- inf_chain %>% left_join(data_age_dist)
  ######################################################
  inf_chains[[i]] <- inf_chain %>% mutate(key=key[i])
  resampled_all[[i]] <- resample_to_target_age_distribution(inf_chain,age_distribution)%>% mutate(key=key[i])
}
resampled_comb <- do.call("bind_rows",resampled_all)
inf_chains_comb <- do.call("bind_rows",inf_chains)

## Plot attack rates over time, weighted vs. reweighted
ars_unweighted <- inf_chains_comb %>% group_by(samp_no,chain_no,time,key) %>% summarize(ar=sum(infected)/n(),N=n()) %>% group_by(time,key) %>% 
  summarize(median=median(ar),lower=quantile(ar,0.025),upper=quantile(ar,0.975),N=median(N)) %>% mutate(ver="Original data")
ars_weighted <- resampled_comb %>% group_by(samp_no,chain_no,time,key) %>% summarize(ar=sum(infected)/n(),N=n()) %>% group_by(time,key) %>% 
  summarize(median=median(ar),lower=quantile(ar,0.025),upper=quantile(ar,0.975),N=median(N)) %>% mutate(ver="ONS census-weighted")
ars_comb <- bind_rows(ars_unweighted,ars_weighted)

#ars_comb %>% mutate("Incidence (per capita)" = paste0(signif(median,3)," (",signif(lower,3),"-",signif(upper,3),")")) %>%
#  select(time,ver,`Incidence (per capita)`,N) %>% rename(`Year`=time,`Estimate`=ver) %>%
#  pivot_wider(names_from=Estimate,values_from=`Incidence (per capita)`) %>%
#  write_csv(paste0("figures/annual_incidence_weighted_vs_unweighted_",file_ver,".csv"))

res1 <- resampled_comb %>% group_by(samp_no, chain_no,key) %>% summarize(ar = sum(infected) / n()) %>%
  group_by(key) %>% summarize(median_ar=median(ar), lower_ar=quantile(ar,0.025), upper_ar=quantile(ar,0.975)) %>% mutate("ver"="Weighted")
res2 <- inf_chains_comb %>% group_by(samp_no, chain_no,key) %>% summarize(ar = sum(infected) / n()) %>%
  group_by(key) %>% summarize(median_ar=median(ar), lower_ar=quantile(ar,0.025), upper_ar=quantile(ar,0.975)) %>%
  mutate("ver"="Unweighted")

p_ars_compare <- ggplot(ars_weighted) + 
  geom_rect(xmin=min(antibody_data$sample_time),xmax=max(antibody_data$sample_time),ymin=-Inf,ymax=0.9,fill="aliceblue",alpha=0.1) +
  
  geom_hline(data=res1,aes(yintercept=median_ar,col=key,linetype="Overall AR (posterior median)"),linewidth=0.25)+
  scale_linetype_manual(name="",values=c("Overall AR (posterior median)"="dashed"))+
  geom_errorbar(aes(x=time,ymin=lower,ymax=upper,col=key),position=position_dodge(width=0.5),width=0.2) +
  geom_point(aes(x=time,y=median,col=key),position=position_dodge(width=0.5),size=1.5) +
  geom_rect(data=noro_epi,aes(xmin=start,xmax=end,ymin=0.9,ymax=1,fill=strain),alpha=0.5,inherit.aes = FALSE) + 
  geom_text(data=noro_epi,aes(x=(start+end)/2,y=0.95,label=strain),inherit.aes = FALSE,size=1.5) +
  scale_x_continuous(breaks=seq(2001,2012,by=1),expand=c(0,0)) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  scale_fill_discrete(name="Dominant strain") +
  ylab("Attack rate (infections per capita)") +
  xlab("Norovirus year") +
  theme_use +
  theme(legend.position="bottom") +
  guides(fill="none") +
  coord_cartesian(xlim=c(2000,2013)) +
  scale_color_manual(name="Model assumptions",values=c("Debbink, exponential"=colors[1],"Kendra, exponential"=colors[2],"Debbink, linear"=colors[3],"Kendra, linear"=colors[4])) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2))

ggsave(p_ars_compare,filename=paste0("figures/attack_rates_weighted_all_scenarios.png"),width=7,height=4,dpi=300)
ggsave(p_ars_compare,filename=paste0("figures/attack_rates_weighted_all_scenarios.pdf"),width=7,height=4,dpi=300)
