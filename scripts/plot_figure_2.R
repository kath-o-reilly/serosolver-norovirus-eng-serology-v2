# create a barchart of the proportion of tested norovirus isolates that correspond to each 
# variant
# these data are based on PHE data, grouped by "norovirus year", ie. July-June.

# load libraries
library(janitor)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

### *** Note that the data for this figure are not publicly available (contact UKHSA for data access) *** 
dat <- read_csv("Data/filtered_metadata.csv")
head(dat)

# reshape
c1 <- colnames(dat)
dat2 <- dat[1:10,] %>% select(-"total tested") 
dat2b <- dat2
dat2b[,2:11] <- dat2[,2:11]/rowSums(dat2[,2:11],na.rm = T)

dat3 <- dat2b %>%
  melt(id.var="Row Labels") %>%
  arrange(`Row Labels`,variable) 


dat3$value[is.na(dat2$value)] <- 0
dat3$`Row Labels` <- as.integer(dat3$`Row Labels`)

# adjust order of group labelling
dat3$variable <- fct_relevel(dat3$variable,
                             "Farmington Hills (FH) 2002","Asia (AS) 2003","Hunter (HU) 2004",
                             "Den Haag (DH) 2006","Yerseke (YE) 2006","Apeldoorn (AD) 2007",
                             "Osaka (OS) 2007","New Orleans (NO) 2009","Sydney (SY) 2012",
                             "Could not assign")

# be specific about colours (use in later plots)
col <- brewer.pal(n = 9, name = 'Set1') # can use for later plots
display.brewer.pal(n = 9, name = 'Set1') # noting the final one is blank white

p1 <- ggplot(dat3,aes(x=`Row Labels`,y=value,group=variable,fill=variable)) + 
  geom_bar(stat = "identity",color="black") +
  xlab("Norovirus year (Jul-Jun)") + ylab("Proportion of Samples") + 
  scale_fill_brewer(palette = "Set1") +
  theme_bw() + # must be before specifics
  theme(legend.title = element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=6)) 

pdf("Results/Figure1_variant_proportions.pdf",height=5,width=7)
p1
dev.off()

# would be preferable to add in 'first year' of sera...
sera <- read_csv("Data/sera_noro_year.csv") %>% clean_names()
sera$noro_year <- sera$var1-0.5
sera$variable <- "Could not assign"
sera <- sera[sera$noro_year>=2003,]

p1b <- ggplot(dat3,aes(x=`Row Labels`,y=value,group=variable,fill=variable)) + 
  geom_bar(stat = "identity",color="black") +
  scale_fill_brewer(palette = "Set1") +
  geom_point(data=sera,aes(x=noro_year,y=freq/110),col="black",pch=18,size=4) +
  # Custom the Y scales:
  scale_y_continuous(
    # Features of the first axis
    name = "Proportion of Samples",
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*110, name="Number of Sera with Corresponding Birth Year")
  ) +
  xlab("Norovirus year (Jul-Jun)") + #ylab("Proportion of Samples") + 
  scale_x_continuous(breaks=seq(2003,2012,1), labels=seq(2003,2012,1), limits=c(2002.5,2012.5)) +
  theme_bw() + # must be before specifics
  theme(legend.title = element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=6)) 

pdf("Results/Figure1_variant_proportions_with_sera.pdf",height=5,width=7)
p1b
dev.off()

# table summary

data <- read_csv("Data/2008-2012_verified_IC50_23AUG2020_complete.csv")

data$NoroSampleYear <- as.numeric(format(as.Date(data$Date,"%d/%m/%Y"),"%Y"))
data$NoroSampleMonth <- as.numeric(format(as.Date(data$Date,"%d/%m/%Y"),"%m"))
data$NoroSampleYear[data$NoroSampleMonth <= 6] <- data$NoroSampleYear[data$NoroSampleMonth <= 6] -1

tmp <- table(birth = data$NoroBirthYear,sample = data$NoroSampleYear)
tmp2 <- rowSums(tmp)
tmp3 <- colSums(tmp)
tmp <- rbind(tmp,tmp3)
tmp <- cbind(tmp,c(tmp2))

write_csv(as.data.frame(tmp),"summary_sample_birth_coll.csv")

summary(data$Age)
