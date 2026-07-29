# compare cartography maps from Debbink and Kendra data

# libraries
library(tidyverse)
library(ggplot2)
library(lubridate)
library(dplyr)
library(cowplot)
library(RColorBrewer)

rm(list = ls())

# ---- Kendra data

dat1 <- read_csv("Data/ac-kendra-all-variants.csv") #
names(dat1) <- c("type","ref","x","y")
head(dat1)

dat1$names <- substr(dat1$ref,7,13)
dat1$variant <- substr(dat1$names,1,2)
dat1$year <- as.numeric(substr(dat1$names,4,7))
dat1 <- dat1b[order(dat1b$year),]
dat1b <- dat1 %>% filter(type == "AG")

# ---- Debbink data

dat2 <- read_csv("Data/ac-debbink-only-04Jan22.csv")
names(dat2) <- c("type","ref","x","y")
head(dat2)
dat2$names <- substr(dat2$ref,7,13)
dat2$variant <- substr(dat2$names,1,2)
dat2$year <- as.numeric(substr(dat2$names,4,8))
dat2 <- dat2[order(dat2$year),]
dat2b <- dat2 %>% filter(type == "AG")

# adjustments
# - on x coord for p2
# make the FH on each the 0,0 co-rdinate
kdxf <-  mean(dat1b$x[dat1b$names == "FH/2004"])
kdyf <-  mean(dat1b$y[dat1b$names == "FH/2004"])
dbxf <-  dat2b$x[dat2b$names == "FH/2002"]
dbyf <-  dat2b$y[dat2b$names == "FH/2002"]

# ** Kendra data - all variants
#             - centered on FH variant
#             - note y has been flipped
#             - no colors
p1 <- ggplot(dat1b,aes(x=x-kdxf,y=-(y-kdyf),group=type,shape=type,col=variant)) +
  geom_point() +
  geom_path(linetype = "dashed",col="black") +
  geom_text(data=dat1b %>% filter(type=="AG"),aes(x=x-kdxf+1.2,y=-(y-kdyf),label=names)) +
  ggtitle("Kendra data (all variants)") +
  scale_fill_brewer(palette="Set1") +
  xlab("Dimension 1") + ylab("Dimension 2") +
  theme_bw()
p1

# use scale_fill_manual() for own color palette 
# see https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
# assigning colors to specific values.... does this work?
col1 <- read_csv("Data/noro_variant_colors.csv")

# make variable for cols that corresponds to palette "col1"
dat1b$var_col <- dat1b$variant
oo <- match(dat1b$var_col,col1$Variant) # if NA then not in color scheme so make XX
dat1b$var_col[is.na(oo)] <- "XX"
# relevel required?
dat1b$var_col <- fct_relevel(dat1b$var_col,col1$Variant) 

# palette for present variants
col2 <- col1 
oo <- match(c("FH","HU","DH","YE","OS","NO","SY","XX"),col2$Variant)
col2 <- col1[oo,]

# ** Kendra data - all variants
#                - same plot with colours
p1 <- ggplot(dat1b,aes(x=x-kdxf,y=-(y-kdyf),color=var_col)) +
  geom_point(size=3) +
  scale_color_manual(name = "variant",values = col2$Color) +
  geom_path(linetype = "dashed",col="black") +
  geom_text(data=dat1b %>% filter(type=="AG"),aes(x=x-kdxf+1.2,y=-(y-kdyf),label=names)) +
  ggtitle("Kendra data (all variants)") +
  scale_x_continuous(name = "Dimension 1", breaks = seq(-8,8,4), labels = seq(-8,8,4), limits = c(-6,9.5)) +
  scale_y_continuous(name = "Dimension 2", breaks = seq(-8,8,4), labels = seq(-8,8,4), limits = c(-2,9.5)) +
  theme_bw() +
  theme(aspect.ratio=1)

pdf("Results/kendra_ac_map_all_vars.pdf",height=6,width=7)
p1
dev.off()

# for the cartography, we want co-ordinates that are an average for each year
# note the geom_path is what we want the cartography to follow

cart_ken <- dat1b %>% filter(type == "AG") %>% group_by(var_col) %>% summarize(n=n(),
                                                                               x_coord = mean(x)-kdxf,
                                                                               y_coord = -(mean(y)-kdyf),
                                                                               year=mean(year))
cart_ken <- cart_ken[order(cart_ken$year),]
write.csv(cart_ken,paste0("Data/antigenic_map_kendra.csv"))

dat1b$x_coord <- dat1b$x-kdxf
dat1b$y_coord <- -(dat1b$y-kdyf)

# check what these look like on the data
p1b <- ggplot(dat1b,aes(x=x_coord,y=y_coord,color=var_col)) +
  geom_point(size=3) +
  scale_color_manual(name = "variant",values = col2$Color) +
  #geom_path(linetype = "dashed",col="black") +
  geom_text(data=dat1b %>% filter(type=="AG"),aes(x=x-kdxf+1.2,y=-(y-kdyf),label=names)) +
  geom_path(data=cart_ken,aes(x=x_coord,y=y_coord,label=var_col),lty=2,col="black") + 
  geom_point(data=cart_ken,aes(x=x_coord,y=y_coord,label=var_col),pch=5,size=3) +
  #geom_path(linetype = "dashed",col="pink")
  ggtitle("Kendra data (all variants)") +
  scale_x_continuous(name = "Dimension 1", breaks = seq(-8,8,4), labels = seq(-8,8,4), limits = c(-6,9.5)) +
  scale_y_continuous(name = "Dimension 2", breaks = seq(-8,8,4), labels = seq(-8,8,4), limits = c(-2,9.5)) +
  theme_bw() +
  theme(aspect.ratio=1)
  
p1b

# output "antigenic_map_kendra_noro"
# fields x_coord y_coord inf_times

# --- plots to check how to use correct color scheme
# --- this took spme time to figure out!!

# tmp2 <- data.frame(x=c(1:9),y=c(1:9),cc=col1$Variant[1:9],dd=col1$Color[1:9])
# # relevel required?
# tmp2$cc <- fct_relevel(tmp2$cc,tmp2$cc) 
# 
# ggplot(tmp2, aes(x = x, y = y, colour = cc)) +       
#   geom_point(size=5) + 
#   geom_text(aes(x = x+1, y = y, label = cc)) +
#   scale_color_manual(values=tmp2$dd) # works if one each, but requires relevel
# 
# tmp3 <- rbind(tmp2,tmp2[c(3,5),])
# tmp3$x2 <- 1:11
# 
# # check it works if multiple observations
# ggplot(tmp3, aes(x = x2, y = y, colour = cc)) +       
#   geom_point(size=5) + 
#   geom_text(aes(x = x2+1, y = y, label = cc)) +
#   scale_color_manual(values=tmp2$dd) # works if one each, but requires relevel
# 
# # check it works if one group is missing
# tmp4 <- data.frame(x=c(1:10),y=c(1:10),cc=col1$Variant[1:10],dd=col1$Color[1:10])
# tmp4$cc <- fct_relevel(tmp4$cc,tmp4$cc) 
# tmp4$cc[tmp4$cc=="NO"] <- "XX"
# cp <- tmp4$dd
# cp <- cp[-8]
# 
# ggplot(tmp4, aes(x = x, y = y, colour = cc)) +       
#   geom_point(size=5) + 
#   geom_text(aes(x = x+1, y = y, label = cc)) +
#   scale_color_manual(values=cp) # pink is missing from dataset, and removed from palette.

# debbink data - all variants tested
#              - use same color mapping and scale

# palette for present variants
col3 <- col1 
oo <- match(c("XX","FH","DH","NO","SY"),col3$Variant) # CW in list but not allocated color
col3 <- col1[oo,]

dat2b$var_col <- dat2b$variant
oo <- match(dat2b$var_col,col1$Variant) # if NA then not in color scheme so make XX
dat2b$var_col[is.na(oo)] <- "XX"
# relevel required?
dat2b$var_col <- fct_relevel(dat2b$var_col,col3$Variant) 

p3 <- ggplot(dat2b,aes(x=-x+dbxf,y=y-dbyf,color=var_col)) +
  geom_point(size=3) +
  scale_color_manual(name = "variant",values = col3$Color) +
  geom_path(linetype = "dashed",col="black") +
  geom_text(data=dat2b,aes(x=-(x+dbxf-7),y=(y-dbyf),label=names)) +
  ggtitle("Debbink data (all variants)") +
  scale_x_continuous(name = "Dimension 1", breaks = seq(-8,8,4), labels = seq(-8,8,4), limits = c(-6,9.5)) +
  scale_y_continuous(name = "Dimension 2", breaks = seq(-8,8,4), labels = seq(-8,8,4), limits = c(-2,9.5)) +
  theme_bw() +
  theme(aspect.ratio=1)

pdf("Results/debbink_ac_map_all_vars.pdf",height=6,width=7)
p3
dev.off()

p4 <- plot_grid(p3,p1b)

pdf("Results/Figure2_both_data_ac_update.pdf",height=6,width=12)
p4
dev.off()



# and what happens if I overlay?
#update x and y
dat1b$x_m <- dat1b$x-kdxf# -dat1b$x+kdxf # kendra 
dat1b$y_m <- dat1b$y-kdyf
dat2b$x_m <- -dat2b$x+dbxf # debbink
dat2b$y_m <- dat2b$y-dbyf
# merge
dat1b$source <- "Kendra"
dat2b$source <- "Debbink"
data <- bind_rows(dat1b,dat2b)

p3 <- ggplot(data,aes(x=x_m,y=y_m,label=names,col=source,group=source)) +
  geom_point() +
  geom_path(linetype = "dashed") +
  geom_text(aes(x=x_m+1,y=y_m)) +
  ggtitle("Kendra & Debbink Data") +
  scale_color_brewer(palette="Dark2") +
  xlab("x") + ylab("y")

pdf("compare_kendra_debbink_cartography.pdf",height=4,width=7)
p3
dev.off()

pdf("compare_kendra_debbink_cartography_noHu.pdf",height=4,width=7)
p3
dev.off()

# end