# simple plots of attack rates from simulations

rm(list=ls())

library(ggplot2)
library(tidyverse)


setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2")

dat <- read_csv("attack_rate_year_cart.csv")

dat$time[dat$cartography=="Kendra"] <- dat$time[dat$cartography=="Kendra"] + 0.2

p1 <- ggplot(dat,aes(x=time,y=q50,group=cartography,col=cartography)) + geom_point() +
  geom_errorbar(aes(ymin=q025, ymax=q975), width=.2) +
  theme_classic() +
  xlab("Date") + ylab("Estimated attack rate") +
  theme(legend.position='top')
p1

pdf("Annual_attack_rate_cart.pdf",height=3.5,width=5)
p1
dev.off()
