library(ggplot2)
library(cowplot)

project_dir <- normalizePath(getwd())

data_folder <- file.path(project_dir, "Data")

# load data
antibody_data <- read_csv(file.path(data_folder, "titre_dat_norovirus_age_avidity.csv"))
antibody_data <- as.data.frame(antibody_data)
antibody_data_keep <- antibody_data
#antibody_data <- antibody_data %>% select(-age,-Date) 
antibody_data$measurement <- log2(antibody_data$measurement/5)

summary(antibody_data$measurement)

theme_use <- theme_classic() + theme(axis.title = element_text(size=8),
                                     axis.text = element_text(size=6),
                                     legend.text = element_text(size=6),
                                     legend.title = element_text(size=6),plot.tag=element_text(face="bold",size=10),
                                     plot.margin=margin(5,5,5,5),
                                     legend.margin=margin(0,0,0,0),
                                     legend.box.margin=margin(-10,-10,-10,-10)) # Reduce spacing around plot and legend

# so the plot we want is 
# A) Sample Year: 2007, Age=1
# B) Sample Year: 2007, Age=3
# C) Sample Year: 2009, Age=1
# D) Sample Year: 2009, Age=3

# (and of responder data!)
summary(tmp2$measurement)

tmp <- antibody_data %>% filter(sample_time==2007 & round(age)==1)
un_tmp <- unique(tmp$individual)
kk <- sample(un_tmp,size=9)
tmp2 <- antibody_data %>% filter(individual %in% kk)
p1 <- ggplot(tmp2 ,aes(x=factor(biomarker_id),y=measurement)) + 
  geom_point() +
  scale_y_continuous(trans='log2') +
  facet_wrap(~individual) + 
  theme(axis.text.x = element_text(angle=45,vjust=0.5)) +
  labs(title="Sample Year: 2007, Age=1",y="log2 IC50", x="Variant response") 

tmp <- antibody_data %>% filter(sample_time==2007 & round(age)==3)
un_tmp <- unique(tmp$individual)
kk <- sample(un_tmp,size=9)
tmp2 <- antibody_data %>% filter(individual %in% kk)
p2 <- ggplot(tmp2 ,aes(x=factor(biomarker_id),y=measurement)) + 
  geom_point() +
  scale_y_continuous(trans='log2') +
  facet_wrap(~individual) + 
  theme(axis.text.x = element_text(angle=45,vjust=0.5)) +
  labs(title="Sample Year: 2007, Age=3",y="log2 IC50", x="Variant response")

tmp <- antibody_data %>% filter(sample_time==2009 & round(age)==1)
un_tmp <- unique(tmp$individual)
kk <- sample(un_tmp,size=9)
tmp2 <- antibody_data %>% filter(individual %in% kk)
p3 <- ggplot(tmp2 ,aes(x=factor(biomarker_id),y=measurement)) + 
  geom_point() +
  scale_y_continuous(trans='log2') +
  facet_wrap(~individual) + 
  theme(axis.text.x = element_text(angle=45,vjust=0.5)) +
  labs(title="Sample Year: 2009, Age=1",y="log2 IC50", x="Variant response") 

tmp <- antibody_data %>% filter(sample_time==2009 & round(age)==3)
un_tmp <- unique(tmp$individual)
kk <- sample(un_tmp,size=9)
tmp2 <- antibody_data %>% filter(individual %in% kk)
p4 <- ggplot(tmp2 ,aes(x=factor(biomarker_id),y=measurement)) + 
  geom_point() +
  scale_y_continuous(trans='log2') +
  facet_wrap(~individual) + 
  theme(axis.text.x = element_text(angle=45,vjust=0.5)) +
  labs(title="Sample Year: 2009, Age=3",y="log2 IC50", x="Variant response") 

pp <- plot_grid(p1,p2,p3,p4,labels=AUTO)

dm <- 0.8
ggsave(paste0("figures/data/data_A2_indiv_resp.pdf"),pp,height=10*dm,width=12*dm,units="in",dpi=300)

