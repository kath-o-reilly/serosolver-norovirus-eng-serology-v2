## Plot attack rates and compare to real data
library(tidyverse)
library(ggplot2)
library(serosolver)
library(data.table)
library(patchwork)

setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/")

## Read in attack rates from ic50 Debbink scenario
serosolver_ar <- read_csv("results_linear/ic50_kendra_attack_rates.csv")
inf_chains <- load_infection_chains(location="local_data_exponential/norovirus_true/chains/ic50_real/",burnin = 400000)

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

## Compare attack rates
p_ar_raw <- ggplot(serosolver_ar %>% group_by(time) %>% summarize(mean=mean(V1),lower=quantile(V1,0.025),upper=quantile(V1,0.975))) +
  geom_rect(data=noro_epi,aes(xmin=start,xmax=end,ymin=0.9,ymax=1,fill=strain),alpha=0.5,inherit.aes = FALSE) + 
  geom_text(data=noro_epi,aes(x=(start+end)/2,y=0.95,label=strain),inherit.aes = FALSE) +
  geom_pointrange(aes(x=time,ymin=lower,y=mean,ymax=upper,group=time,col="serosolver (per capita)")) + 
  geom_line(data=real_dat,aes(x=year,y=`1. 0-4 years`,col="0-4 years"),linewidth=0.25) +
  geom_point(data=real_dat,aes(x=year,y=`1. 0-4 years`,col="0-4 years")) +
  geom_point(data=real_dat,aes(x=year,y=`2. 5-9 years`,col="5-9 years")) +
  geom_line(data=real_dat,aes(x=year,y=`2. 5-9 years`,col="5-9 years"),linewidth=0.25) +
  geom_point(data=real_dat,aes(x=year,y=`all`,col="All")) +
  geom_line(data=real_dat,aes(x=year,y=`all`,col="All"),linewidth=0.25) +
  scale_color_manual(name="",values=c("0-4 years"=colors[2],"5-9 years"=colors[3],"All"=colors[4],"serosolver (per capita)"="grey40")) +
  theme_use +
  ylab("Attack rate (infections per capita)")  +
# Add a second y-axis to this plot which is a factor of 1000 lower
  scale_y_continuous(limits=c(0,1),sec.axis = sec_axis(~./1, name = "Reported cases (per 1000 people)")) +
  scale_fill_viridis_d(name="Norovirus strain") +
  theme(legend.position="bottom") +
  scale_x_continuous(breaks=seq(2001,2013,by=1)) +
  xlab("Norovirus year")

ggsave(p_ar_raw,filename=paste0("figures/attack_rate_comparison_",file_ver,".png"),width=6,height=4,units="in",dpi=300)
  
## Read in age distribution data from ONS
age_distribution <- age_distribution %>% select(-Sex) %>% pivot_longer(-Age) %>%
  group_by(name) %>% dplyr::mutate(weighted=value/sum(value)) %>%
  mutate(year = str_split(name,"Mid-")[[1]][2]) %>% mutate(year=as.integer(year)) %>%
  filter(year >= 2001, year <= 2012) %>%
  rename(time=year)


data_age_dist <- antibody_data %>% select(individual,birth) %>% distinct() %>% expand_grid(time = seq(2001,2012,by=1)) %>% mutate(age=time-birth) %>% filter(age >=0) %>% 
  group_by(time,age) %>% tally() %>%
  arrange(time,age) %>%
  group_by(time) %>% mutate(prop=n/sum(n))

## Age distribution of dataset vs. population
p_age_dist_compare <- ggplot(age_distribution %>% group_by(time) %>% mutate(weighted=weighted/max(weighted))) + 
  geom_bar(aes(x=Age,y=weighted,fill="ONS census"),stat="identity",alpha=0.5) +
  geom_bar(data=data_age_dist %>% group_by(time) %>% mutate(prop=prop/max(prop)),aes(x=age,y=prop,fill="Norovirus data"),stat="identity",alpha=0.5) +
  scale_fill_manual(name="",values=c("ONS census"=colors[2],"Norovirus data"=colors[3])) +
  theme_use +
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
  xlab("Age") +
  ylab("Proportion of population \n (scaled)") +
  theme(legend.position="bottom") +
  facet_wrap(~time)
ggsave(p_age_dist_compare,filename=paste0("figures/checks/age_distribution_comparison_",file_ver,".png"),width=8,height=6,units="in",dpi=300)


## Read in the main serosolver infection history chains
inf_chain <- pad_inf_chain(inf_chains$chain,indivs = unique(antibody_data$individual))        
colnames(inf_chain) <- c("individual","time","infected","samp_no","chain_no")
times <- seq(2001,2012,by=1)
inf_chain$time <- inf_chain$time - min(inf_chain$time) + 1
inf_chain$time <- times[inf_chain$time]

## Merge with age data
ages <- antibody_data %>% select(individual,birth) %>% distinct() %>% expand_grid(time=times) %>% mutate(age=time-birth)
inf_chain <- inf_chain %>% left_join(ages) %>% filter(age >=0)
inf_chain <- inf_chain %>% mutate(age_group_simple = if_else(age < 5,as.character(age),"5+"))
######################################################
## For each sample number/chain number, resample individuals by age for each year to get the same age distribution as the ONS census data
data_age_dist <- data_age_dist %>% group_by(time) %>% mutate(total_n=sum(n))
inf_chain <- inf_chain %>% left_join(data_age_dist)
######################################################

# 0) Target proportions by time & age
target_prop <- age_distribution %>%
  transmute(time, age = as.integer(Age), target_prop = weighted)

# 1) Start from the set of individuals alive at each (samp_no, chain_no, time)
#    If inf_chain can have multiple rows per individual-time (e.g., events),
#    keep one row per individual so we sample persons, not occurrences.
pop <- inf_chain %>%
  distinct(samp_no, chain_no, time, individual, age, .keep_all = TRUE)

# 2) Compute the current age distribution within each (samp_no, chain_no, time)
curr_age <- pop %>%
  count(samp_no, chain_no, time, age, name = "n_age") %>%
  group_by(samp_no, chain_no, time) %>%
  mutate(N_group = sum(n_age),
         prop_curr = n_age / N_group) %>%
  ungroup()

# 3) Join target proportions and compute re-weighting factors per age
weights_by_age <- curr_age %>%
  left_join(target_prop, by = c("time", "age")) %>%
  mutate(
    target_prop = coalesce(target_prop, 0),
    # If current prop is 0, set weight to 0 (no individuals of that age to sample)
    w_age = if_else(prop_curr > 0, target_prop / prop_curr, 0)
  ) %>%
  select(samp_no, chain_no, time, age, N_group, w_age)

# 4) Attach the age-level weight to each individual, then normalize within group
pop_weighted <- pop %>%
  left_join(weights_by_age, by = c("samp_no", "chain_no", "time", "age")) %>%
  group_by(samp_no, chain_no, time) %>%
  # If all w_age are 0 (e.g., target has zero mass on all available ages or missing),
  # fall back to equal weights; otherwise normalize w_age.
  mutate(
    w_age = replace_na(w_age, 0),
    sum_w = sum(w_age),
    prob = if_else(sum_w > 0, w_age / sum_w, 1 / n())
  ) %>%
  select(-sum_w) %>%
  ungroup()

# 5) Sample WITH replacement to keep the same total size per group
#    (rows are repeated according to how many times that individual is drawn)
resampled <- pop_weighted %>%
  group_by(samp_no, chain_no, time) %>%
  group_modify(function(df, keys) {
    size_out <- df$N_group[match(TRUE, !is.na(df$N_group))]  # same N as original group
    idx <- sample.int(n = nrow(df), size = size_out, replace = TRUE, prob = df$prob)
    # Return the selected rows; add draw index if you want
    df[idx, , drop = FALSE] %>%
      mutate(draw_id = row_number())
  }) %>%
  ungroup() %>%
  arrange(samp_no, chain_no, time, draw_id)

resampled_check <- resampled %>% group_by(samp_no,chain_no,time,age) %>% tally() %>% group_by(samp_no,chain_no,time) %>% mutate(total_n=sum(n)) %>% group_by(time,age) %>% summarize(med_n = median(n)) %>% group_by(time) %>% mutate(prop=med_n/sum(med_n)) %>%
  mutate(prop = prop/max(prop))

## Age distribution of dataset vs. population
p_age_dist_compare_updated <- ggplot(age_distribution %>% group_by(time) %>% mutate(weighted=weighted/max(weighted))) + 
  geom_bar(aes(x=Age,y=weighted,fill="ONS census"),stat="identity",alpha=0.5) +
  geom_bar(data=resampled_check,aes(x=age,y=prop,fill="Norovirus data"),stat="identity",alpha=0.5) +
  scale_fill_manual(name="",values=c("ONS census"=colors[2],"Norovirus data"=colors[3])) +
  theme_use +
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
  xlab("Age") +
  ylab("Proportion of population \n (scaled)") +
  theme(legend.position="bottom") +
  facet_wrap(~time)
ggsave(p_age_dist_compare_updated,filename=paste0("figures/checks/age_distribution_comparison_updated_",file_ver,".png"),width=8,height=6,units="in",dpi=300)


resampled <- resampled %>% left_join(resampled_check)

## Setup alternative age groups, 0-4, 5-10, and 5+
resampled <- resampled %>% mutate(age_group_wide = if_else(age <= 4, "0-4","5-10"),
                                  age_group_simple = if_else(age < 5, as.character(age), "5+"))


######################################################
## Get overall annual median attack rate from both weighted and unweighted
res1 <- resampled %>% group_by(samp_no, chain_no) %>% summarize(ar = sum(infected) / n()) %>%
  ungroup() %>% summarize(median_ar=median(ar), lower_ar=quantile(ar,0.025), upper_ar=quantile(ar,0.975)) %>% mutate("ver"="Weighted")
res2 <- inf_chain %>% group_by(samp_no, chain_no) %>% summarize(ar = sum(infected) / n()) %>%
  ungroup() %>% summarize(median_ar=median(ar), lower_ar=quantile(ar,0.025), upper_ar=quantile(ar,0.975)) %>%
  mutate("ver"="Unweighted")

overall_incidence <- bind_rows(res1,res2) %>% mutate("Cases per 1000 persons per year"=paste0(signif(median_ar*1000,3)," (",signif(lower_ar*1000,3),"-",signif(upper_ar*1000,3),")")) %>%
  select(ver,`Cases per 1000 persons per year`) %>% rename(`Estimate`=ver) %>%
  mutate(Estimate = paste0(Estimate, " cases per 1000 persons per year (posterior median and 95% CrI)")) %>%
  pivot_wider(names_from=Estimate,values_from=`Cases per 1000 persons per year`)

write_csv(overall_incidence,paste0("figures/overall_incidence_weighted_vs_unweighted_",file_ver,".csv"))

## Plot attack rates over time, weighted vs. reweighted
ars_unweighted <- inf_chain %>% group_by(samp_no,chain_no,time) %>% summarize(ar=sum(infected)/n(),N=n()) %>% group_by(time) %>% 
  summarize(median=median(ar),lower=quantile(ar,0.025),upper=quantile(ar,0.975),N=median(N)) %>% mutate(ver="Original data")
ars_weighted <- resampled %>% group_by(samp_no,chain_no,time) %>% summarize(ar=sum(infected)/n(),N=n()) %>% group_by(time) %>% 
  summarize(median=median(ar),lower=quantile(ar,0.025),upper=quantile(ar,0.975),N=median(N)) %>% mutate(ver="ONS census-weighted")
ars_comb <- bind_rows(ars_unweighted,ars_weighted)

ars_comb %>% mutate("Incidence (per capita)" = paste0(signif(median,3)," (",signif(lower,3),"-",signif(upper,3),")")) %>%
  select(time,ver,`Incidence (per capita)`,N) %>% rename(`Year`=time,`Estimate`=ver) %>%
  pivot_wider(names_from=Estimate,values_from=`Incidence (per capita)`) %>%
  write_csv(paste0("figures/annual_incidence_weighted_vs_unweighted_",file_ver,".csv"))

p_ars_compare <- ggplot(ars_comb) + 
  geom_hline(data=bind_rows(res1 %>% mutate(ver="ONS census-weighted"),res2 %>% mutate(ver="Original data")),
             aes(yintercept=median_ar,col=ver,linetype="Overall AR (posterior median)"),linewidth=0.25)+
  scale_linetype_manual(name="",values=c("Overall AR (posterior median)"="dashed")) +
  geom_errorbar(aes(x=time,ymin=lower,ymax=upper,col=ver),position=position_dodge(width=0.5),width=0.2) +
  geom_point(aes(x=time,y=median,col=ver),position=position_dodge(width=0.5),size=1.5) +
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
  scale_color_manual(name="AR estimate",values=c("Original data"=colors[2],"ONS census-weighted"=colors[3])) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2))

ggsave(p_ars_compare,filename=paste0("figures/attack_rate_comparison_weighted_vs_unweighted_",file_ver,".png"),width=6,height=4,units="in",dpi=300)

p_ars_weighted <- ggplot(ars_weighted) + 
  geom_hline(data=res1, aes(yintercept=median_ar,linetype="Overall AR (posterior median)"),linewidth=0.25,col=colors[6])+
  geom_hline(data=res1, aes(yintercept=lower_ar,linetype="Overall AR (95% CrI)"),linewidth=0.25,col=colors[6])+
  geom_hline(data=res1, aes(yintercept=upper_ar,linetype="Overall AR (95% CrI)"),linewidth=0.25,col=colors[6])+
  scale_linetype_manual(name="",values=c("Overall AR (posterior median)"="solid","Overall AR (95% CrI)"="dashed")) +
  
  geom_errorbar(aes(x=time,ymin=lower,ymax=upper),width=0.2,col="grey40") +
  geom_point(aes(x=time,y=median),size=1.5,col="grey40") +
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
  scale_color_manual(name="AR estimate",values=c("Original data"=colors[2],"ONS census-weighted"=colors[3])) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1))

ggsave(p_ars_weighted,filename=paste0("figures/attack_rate_weighted",file_ver,".png"),width=6,height=4,units="in",dpi=300)



## Plot attack rates over time and by age, weighted vs. reweighted
## Get age-stratified annual median attack rate from both weighted and unweighted
res1_age <- resampled %>% group_by(samp_no, chain_no,age_group_simple) %>% summarize(ar = sum(infected) / n()) %>%
  group_by(age_group_simple) %>% summarize(median_ar=median(ar), lower_ar=quantile(ar,0.025), upper_ar=quantile(ar,0.975)) %>% mutate("ver"="Weighted")
res2_age <- inf_chain %>% group_by(samp_no, chain_no,age_group_simple) %>% summarize(ar = sum(infected) / n()) %>%
  group_by(age_group_simple) %>% summarize(median_ar=median(ar), lower_ar=quantile(ar,0.025), upper_ar=quantile(ar,0.975)) %>%
  mutate("ver"="Unweighted")

overall_incidence_by_age <- bind_rows(res1_age,res2_age) %>% mutate("Cases per 1000 persons per year"=paste0(signif(median_ar*1000,3)," (",signif(lower_ar*1000,3),"-",signif(upper_ar*1000,3),")")) %>% 
  select(age_group_simple, ver, `Cases per 1000 persons per year`) %>% rename(`Age (years)`=age_group_simple,`Estimate`=ver) %>%
  mutate(Estimate = paste0(Estimate, " cases per 1000 persons per year (posterior median and 95% CrI)")) %>%
  pivot_wider(names_from=Estimate,values_from=`Cases per 1000 persons per year`)

write_csv(overall_incidence_by_age,paste0("figures/overall_incidence_by_age_weighted_vs_unweighted_",file_ver,".csv"))

ars_weighted_age <- resampled %>% group_by(samp_no,chain_no,time,age_group_simple) %>% summarize(ar=sum(infected)/n())
ars_weighted_summary_age <- ars_weighted_age %>% group_by(time,age_group_simple) %>% summarize(median=median(ar),lower=quantile(ar,0.025),upper=quantile(ar,0.975)) %>% mutate(ver="ONS census-weighted")

## And the same for unweighted
ars_unweighted_age <- inf_chain %>% group_by(samp_no,chain_no,time,age_group_simple) %>% summarize(ar=sum(infected)/n())
ars_unweighted_summary_age <- ars_unweighted_age %>% group_by(time,age_group_simple) %>% summarize(median=median(ar),lower=quantile(ar,0.025),upper=quantile(ar,0.975)) %>% mutate(ver="Original data")

ars_summary_age_comb <- bind_rows(ars_unweighted_summary_age,ars_weighted_summary_age)

p_ars_age_compare <- ggplot(ars_summary_age_comb) + 
  geom_hline(data=bind_rows(res1_age %>% mutate(ver="ONS census-weighted"),res2_age %>% mutate(ver="Original data")),
             aes(yintercept=median_ar,col=ver,linetype="Overall AR (posterior median)"),linewidth=0.25) +
  scale_linetype_manual(name="",values=c("Overall AR (posterior median)"="dashed")) +
  geom_errorbar(aes(x=time,ymin=lower,ymax=upper,col=ver),position=position_dodge(width=0.5),width=0.2) +
  geom_point(aes(x=time,y=median,col=ver),position=position_dodge(width=0.5),size=1.5) +
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
  scale_color_manual(name="AR estimate",values=c("Original data"=colors[2],"ONS census-weighted"=colors[3])) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2)) +
  facet_wrap(~paste0(age_group_simple, " years"),ncol=2)

ggsave(p_ars_age_compare,filename=paste0("figures/attack_rate_comparison_by_age_weighted_vs_unweighted_",file_ver,".png"),width=7,height=8,units="in",dpi=300)

p_ars_age_weighted <- ggplot(ars_weighted_summary_age) + 
  geom_hline(data=res1_age %>% filter(age_group_simple %in% unique(ars_weighted_summary_age$age_group_simple)), aes(yintercept=median_ar,linetype="Overall AR (posterior median)"),linewidth=0.25,col=colors[6])+
  geom_hline(data=res1_age %>% filter(age_group_simple %in% unique(ars_weighted_summary_age$age_group_simple)), aes(yintercept=lower_ar,linetype="Overall AR (95% CrI)"),linewidth=0.25,col=colors[6])+
  geom_hline(data=res1_age %>% filter(age_group_simple %in% unique(ars_weighted_summary_age$age_group_simple)), aes(yintercept=upper_ar,linetype="Overall AR (95% CrI)"),linewidth=0.25,col=colors[6])+
  scale_linetype_manual(name="",values=c("Overall AR (posterior median)"="solid","Overall AR (95% CrI)"="dashed")) +
  geom_errorbar(aes(x=time,ymin=lower,ymax=upper),width=0.2) +
  geom_point(aes(x=time,y=median),size=1.5) +
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
  scale_color_manual(name="AR estimate",values=c("Original data"=colors[2],"ONS census-weighted"=colors[3])) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2)) +
  facet_wrap(~paste0(age_group_simple, " years"),ncol=2)

ggsave(p_ars_age_weighted,filename=paste0("figures/attack_rate_comparison_by_age_weighted",file_ver,".png"),width=7,height=8,units="in",dpi=300)

## Plot cumulative incidence over age
## All times
# Arrange inf_chain by samp no, chain no and individual in ascending time
inf_chain_cumu <- inf_chain %>% arrange(chain_no, samp_no, individual,time) %>%
  group_by(chain_no,samp_no,individual) %>%
  mutate(cumu_infection = cumsum(infected)) %>%
  mutate(cumu_infection = if_else(cumu_infection == 0, 0, 1))

## Cumulative infections by age
inf_chain_cumulative_incidence <- inf_chain %>% arrange(chain_no, samp_no, individual,time) %>%
  group_by(chain_no,samp_no,individual) %>%
  mutate(cumu_infection = cumsum(infected)) %>% group_by(age, samp_no,chain_no) %>%
  summarize(cumu_infections = mean(cumu_infection))

ggplot(inf_chain_cumulative_incidence) + geom_line(aes(x=age,y=cumu_infections,group=interaction(samp_no,chain_no)),col="grey80",alpha=0.5) +
  scale_y_continuous(limits=c(0,3))

## Look at number/distribution of reinfections
inf_totals_summary <- inf_chain %>% arrange(chain_no, samp_no, individual,time) %>%
  group_by(chain_no,samp_no,individual) %>%
  mutate(total_inf = sum(infected)) %>%
  group_by(individual) %>%
  summarize(med_infections=median(total_inf),lower=quantile(total_inf,0.025),upper=quantile(total_inf,0.975))

p_total_infs <- ggplot(inf_totals_summary) + geom_histogram(aes(x=med_infections),binwidth=1,col="black",fill="grey70") +
  theme_use +
  xlab("Total number of infections\n(posterior median)") +
  scale_y_continuous(expand=c(0,0),limits=c(0,300),breaks=seq(0,300,by=50)) +
  ylab("Count")
write_csv(inf_totals_summary, file=paste0("figures/total_infections_per_person_",file_ver,".csv"))
ggsave(p_total_infs,filename=paste0("figures/checks/total_infections_per_person_",file_ver,".png"),width=4,height=3,units="in",dpi=300)

inf_chain_cumu_summary <- inf_chain_cumu %>% 
  #filter(time > 2007) %>%
  group_by(samp_no, chain_no, age) %>% 
  summarize(prop = sum(cumu_infection)/n(),N=n()) %>%
  group_by(age) %>% summarize(mean1=mean(prop), med=median(prop),lower=quantile(prop,0.025),upper=quantile(prop,0.975),N=median(N))

p_cumu_by_age <- ggplot(inf_chain_cumu %>% 
                          #filter(time > 2007) %>%
                          group_by(samp_no, chain_no, age) %>% 
                          summarize(prop = sum(cumu_infection)/n(),N=n())) +
  #geom_ribbon(aes(x=age,ymin=lower,ymax=upper),fill="grey80") +
  geom_line(aes(x=age,y=prop,group=interaction(chain_no,samp_no),col="Posterior draw"),alpha=0.1) +
  #geom_line(data=inf_chain_cumulative_incidence, aes(x=age,y=cumu_infections,group=interaction(chain_no,samp_no),col="Posterior draw"),alpha=0.1) +
  geom_line(data=inf_chain_cumu_summary,aes(x=age,y=mean1,col="Posterior mean"),linewidth=1) +
 # geom_line(data=inf_chain_cumulative_incidence %>% group_by(age) %>% dplyr::summarize(mean_cumu=mean(cumu_infections)),
#            aes(x=age,y=mean_cumu,col="Posterior mean (number infections)"),linewidth=1) +
  #geom_hline(yintercept=1,linetype="dashed") +
  theme_use +
  xlab("Age (years)") +
  ylab("Cumulative incidence\n (proportion seroconverted)") +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1),expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0,11,by=1))+
  scale_color_manual(name="",values=c("Posterior draw"="grey40","Posterior mean"="red")) +
  labs(tag="A") +
theme(legend.position=c(0.8,0.2))
p_cumu_by_age

p_cumu_incidence_by_age <- ggplot() +
  geom_line(data=inf_chain_cumulative_incidence, aes(x=age,y=cumu_infections,group=interaction(chain_no,samp_no),col="Posterior draw"),alpha=0.1) +
  geom_line(data=inf_chain_cumulative_incidence %>% group_by(age) %>% dplyr::summarize(mean_cumu=mean(cumu_infections)),
              aes(x=age,y=mean_cumu,col="Posterior mean"),linewidth=1) +
  geom_hline(yintercept=1,linetype="dashed") +
  theme_use +
  xlab("Age (years)") +
  ylab("Cumulative incidence\n (total infections)") +
  scale_y_continuous(limits=c(0,2.5),breaks=seq(0,2.5,by=0.2),expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0,11,by=1))+
  scale_color_manual(name="",values=c("Posterior draw"="grey40","Posterior mean"="red")) +
  theme(legend.position=c(0.8,0.2))+
  labs(tag="A")
p_cumu_incidence_by_age

write_csv(inf_chain_cumu_summary, file=paste0("figures/cumulative_incidence_by_age_all_times_",file_ver,".csv"))

## Raw incidence by age
inf_chain_summary <- inf_chain %>% 
  #filter(time > 2007) %>%
  group_by(samp_no, chain_no, age) %>% 
  summarize(prop = sum(infected)/n(),N=n())

p_incidence_by_age <- ggplot(inf_chain_summary) +
  geom_boxplot(aes(x=age,y=prop,group=age),fill="grey80") +
  theme_use +
  xlab("Age (years)") +
  ylab("Incidence per year\n (infections per capita)") +
  scale_x_continuous(breaks=seq(0,11,by=1))+
  scale_y_continuous(limits=c(0,0.4),breaks=seq(0,0.3,by=0.1)) +
  labs(tag="B") 

p_incidence_by_age
write_csv(inf_chain_summary %>% group_by(age) %>% summarize(med=median(prop),lower=quantile(prop,0.025),upper=quantile(prop,0.975),N=median(N)), file=paste0("figures/incidence_by_age_all_times_",file_ver,".csv"))

ggsave(p_cumu_by_age/p_incidence_by_age,filename=paste0("figures/cumulative_and_incidence_by_age_all_times_",file_ver,".png"),width=4,height=5,units="in",dpi=300)


## Just after 2006
inf_chain_cumu_summary <- inf_chain_cumu %>% 
  filter(time > 2007) %>%
  group_by(samp_no, chain_no, age) %>% 
  summarize(prop = sum(cumu_infection)/n(),N=n()) %>%
  group_by(age) %>% summarize(med=median(prop),lower=quantile(prop,0.025),upper=quantile(prop,0.975),N=median(N))

p_cumu_by_age <- ggplot(inf_chain_cumu_summary) +
  geom_ribbon(aes(x=age,ymin=lower,ymax=upper),fill="grey80") +
  geom_line(aes(x=age,y=med),col="grey40") +
  theme_use +
  xlab("Age (years)") +
  ylab("Cumulative incidence\n (proportion seroconverted)") +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1),expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0,11,by=1))+
  labs(tag="A")
p_cumu_by_age

write_csv(inf_chain_cumu_summary, file=paste0("figures/cumulative_incidence_by_age_2007+_,",file_ver,".csv"))

## Raw incidence by age
inf_chain_summary <- inf_chain %>% 
  filter(time > 2007) %>%
  group_by(samp_no, chain_no, age) %>% 
  summarize(prop = sum(infected)/n(),N=n())

p_incidence_by_age <- ggplot(inf_chain_summary) +
  geom_boxplot(aes(x=age,y=prop,group=age),fill="grey80") +
  theme_use +
  xlab("Age (years)") +
  ylab("Incidence per year\n (infections per capita)") +
  scale_x_continuous(breaks=seq(0,11,by=1))+
  scale_y_continuous(limits=c(0,0.4),breaks=seq(0,0.3,by=0.1)) +
  labs(tag="B")
p_incidence_by_age
write_csv(inf_chain_summary %>% group_by(age) %>% summarize(med=median(prop),lower=quantile(prop,0.025),upper=quantile(prop,0.975),N=median(N)), file=paste0("figures/incidence_by_age_2007+_",file_ver,".csv"))

ggsave(p_cumu_by_age/p_incidence_by_age,filename=paste0("figures/cumulative_and_incidence_by_age_2007+_",file_ver,".png"),width=4,height=4,units="in",dpi=300)

## Plot spline of age at isolation
antibody_data <- antibody_data %>% mutate(age_at_circulation = biomarker_id - birth)
antibody_data$strain <- norovirus_key[as.character(antibody_data$biomarker_id)]
antibody_data$strain <- factor(antibody_data$strain, levels=norovirus_key)
p_ic50_seniority <- ggplot(antibody_data %>% filter(biomarker_group == 1) %>% mutate(measurement=log2(measurement/5)) %>% filter(measurement > 0)) +
  geom_vline(xintercept=0,col="red") +
  geom_point(aes(y=measurement,x=age_at_circulation),alpha=0.75,size=0.5,col="grey40") +
  geom_smooth(aes(y=measurement,x=age_at_circulation)) +
  scale_y_continuous(limits=c(0,13),breaks=seq(0,12,by=2)) +
  scale_x_continuous(breaks=seq(-12,12,by=2)) +
  facet_wrap(~strain) +
  theme_use +
  ylab("log IC50") +
  xlab("Age (years) at first virus isolation")

p_ic50_seniority_DH2006 <- ggplot(antibody_data %>% filter(biomarker_group == 1) %>% mutate(measurement=log2(measurement/5)) %>% filter(measurement > 0,biomarker_id==2006)) +
  geom_vline(xintercept=0,col="red") +
  geom_point(aes(y=measurement,x=age_at_circulation),alpha=0.75,size=0.5,col="grey40") +
  geom_smooth(aes(y=measurement,x=age_at_circulation)) +
  scale_y_continuous(limits=c(0,13),breaks=seq(0,12,by=2)) +
  scale_x_continuous(breaks=seq(-5,5,by=1)) +
  theme_use +
  ylab("log IC50 against DH-2006") +
  xlab("Age (years) at first virus isolation")

ggsave(p_ic50_seniority,filename=paste0("figures/ic50_seniority_effect_",file_ver,".png"),width=5,height=3,units="in",dpi=300)
ggsave(p_ic50_seniority_DH2006,filename=paste0("figures/ic50_seniority_effect_dh2006",file_ver,".png"),width=5,height=3,units="in",dpi=300)


((p_ic50_seniority_DH2006 + labs(tag="A")) / plot_spacer()) | ((p_ars_weighted + labs(tag="C") + theme(legend.position=c(0.8,0.6))) / (p_incidence_by_age + labs(tag="D")))


p_incidence_all <- (p_ars_weighted + labs(tag="A") + theme(legend.position=c(0.8,0.6))) + 
  (p_cumu_by_age + labs(tag="C")) + 
  (p_incidence_by_age + labs(tag="B")) + 
  (p_cumu_incidence_by_age + labs(tag="D")) + plot_layout(ncol=2)
ggsave(p_incidence_all,filename=paste0("figures/incidence_all_",file_ver,".png"),width=8,height=6,units="in",dpi=300)                                                               
  