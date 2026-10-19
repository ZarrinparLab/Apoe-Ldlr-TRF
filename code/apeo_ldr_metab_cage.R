setwd("/mnt/zarrinpar/scratch/sfloresr/maria_apoe/")

library(tidyverse)
library(data.table)
library(ggvenn)
library(gplots)
library(ggpubfigs)
library(viridis)
library(RColorBrewer)
library(BioVenn) 
library(VennDiagram)
library("qiime2R")
library(ggpubr)
library(car)
library(agricolae)

###########################################################
#metabolic cages

df<- read.csv("2026_reprocess_16s/data/2021/MetabolicCage/merged.csv", na.strings = '.',
              colClasses = c(Animal="factor")) %>%
  filter(!is.na(Animal))


df_alive<-filter(df,Animal!=2)%>%
  mutate(Group = gsub("_", "", Group))%>%
  mutate(zt_time_binned = case_when(
    between(DateTime, 1, 4) ~ "4",
    between(DateTime, 5, 8) ~ "8",
    between(DateTime, 9, 12) ~ "12",
    between(DateTime, 13, 16) ~ "16",
    between(DateTime, 17, 20) ~ "20",
    TRUE ~ "0" ))%>%
  mutate(zt_time_binned=factor(zt_time_binned,levels=c("0","4","8","12","16","20")))


for (i in 1:nrow(df_alive)){
  if (i-1 == 0){
    df_alive$WaterCons[i] <- 0
    df_alive$FoodCons[i] <- 0
    df$food_cumulative[i] <- 0
  } else if (df_alive$Animal[i] != df_alive$Animal[i-1]){
    df_alive$WaterCons[i] <- 0
    df_alive$FoodCons[i] <- 0
    df_alive$food_cumulative[i] <- 0
  }else {
    water_diff <- df_alive$WaterInA_M[i] - df_alive$WaterInA_M[i-1]
    food_diff <- df_alive$FoodInA_M[i] - df_alive$FoodInA_M[i-1]
    food_cumulative <- df_alive$FoodInA_M[i] + df_alive$FoodInA_M[i-1]
    df_alive$WaterCons[i] <- water_diff
    df_alive$FoodCons[i] <- food_diff
    df_alive$food_cumulative[i] <- food_cumulative
  }
}

df_apoe<-df_alive%>%filter(Genotype=="ApoeKO")
df_ldlr<-df_alive%>%filter(Genotype=="LdlrKO")

###########################################################
#RER

RER <- ggplot(data=df_alive, aes(x=DateTime, y=RER_M, colour=Group, fill=Group)) +
  facet_wrap(~Genotype, scales="free") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  labs(title="Mean respiratory exchange ratio (RER)") +
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))

rermod1 <- aov(RER_M~ Group * zt_time_binned, data = df_apoe)
summary(rermod1)

LSD_Test_rer1<- (LSD.test(rermod1, c("Group", "zt_time_binned")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Group","DateTime"))
write.table(LSD_Test_rer1,"SFR26_0205_apoerer_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)

rermod2 <- aov(RER_M~ Group * zt_time_binned, data = df_ldlr)
summary(rermod2)

LSD_Test_rer2<- (LSD.test(rermod2, c("Group", "zt_time_binned")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Group","DateTime"))
write.table(LSD_Test_rer2,"SFR26_0205_ldlrrer_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)

###########################################################
#mean energy expenditure

kcalh <- ggplot(data=df_alive, aes(x=DateTime, y=kcal_hr_M, colour=Group, fill=Group)) +
  facet_wrap(~Genotype, scales="free") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  labs(title="Mean energy expenditure (kcal/hr)") +
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))

kcalmod1 <- aov(kcal_hr_M~ Group * zt_time_binned, data = df_apoe)
summary(kcalmod1)

LSD_Test_kcal1<- (LSD.test(kcalmod1, c("Group", "zt_time_binned")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Group","DateTime"))
write.table(LSD_Test_kcal1,"SFR26_0205_apoekcalh_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)

kcalmod2 <- aov(kcal_hr_M~ Group * zt_time_binned, data = df_ldlr)
summary(kcalmod2)

LSD_Test_kcal2<- (LSD.test(kcalmod2, c("Group", "zt_time_binned")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Group","DateTime"))
write.table(LSD_Test_kcal2,"SFR26_0205_ldlrkcalh_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)

###########################################################
#wheel activity

wheel <- ggplot(data=df_alive, aes(x=DateTime, y=wheel_activity, colour=Group, fill=Group)) +
  facet_wrap(~Genotype, scales="free") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  labs(title="Wheel Activity (m)") +
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))

wmod1 <- aov(wheel_activity~ Group * zt_time_binned, data = df_apoe)
summary(wmod1)

LSD_Test_w1<- (LSD.test(wmod1, c("Group", "zt_time_binned")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Group","DateTime"))
write.table(LSD_Test_w1,"SFR26_0205_apoewheel_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)

wmod2 <- aov(wheel_activity~ Group * zt_time_binned, data = df_ldlr)
summary(wmod2)

LSD_Test_w2<- (LSD.test(wmod2, c("Group", "zt_time_binned")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Group","DateTime"))
write.table(LSD_Test_w2,"SFR26_0205_ldlrwheel_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)

###########################################################
#food consumption --done separately

foodcons<-fread("2026_reprocess_16s/data/2021/MetabolicCage/foodconsumption_formatted.csv")

foodcons_cln<-foodcons%>%
  mutate(Group = paste(genotype,feeding,sep=""))%>%
  mutate(zt_time_binned = case_when(
    between(ZT_bin, 1, 4) ~ "4",
    between(ZT_bin, 5, 8) ~ "8",
    between(ZT_bin, 9, 12) ~ "12",
    between(ZT_bin, 13, 16) ~ "16",
    between(ZT_bin, 17, 20) ~ "20",
    TRUE ~ "0" ))%>%
  mutate(zt_time_binned=factor(zt_time_binned,levels=c("0","4","8","12","16","20")))

fc_apoe<-foodcons_cln%>%filter(genotype=="ApoeKO")
fc_ldlr<-foodcons_cln%>%filter(genotype=="LdlrKO")

foodc <- ggplot(data=foodcons_cln, aes(x=ZT_bin, y=food_change_g, colour=Group, fill=Group)) +
  facet_wrap(~genotype, scales="free") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  labs(title="Food Consumption (g)") +
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))

fcmod1 <- aov(food_change_g~ Group * zt_time_binned, data = fc_apoe)
summary(fcmod1)

LSD_Test_fc1<- (LSD.test(fcmod1, c("Group", "zt_time_binned")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Group","ZT"))
write.table(LSD_Test_fc1,"SFR26_0205_apoefoodconsump_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)

fcmod2 <- aov(food_change_g~ Group * zt_time_binned, data = fc_ldlr)
summary(fcmod2)

LSD_Test_fc2<- (LSD.test(fcmod2, c("Group", "zt_time_binned")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Group","ZT"))
write.table(LSD_Test_fc2,"SFR26_0205_ldlrfoodconsump_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)

