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
#cumul food intake

foodint<-fread("2026_reprocess_16s/data/2021/ldlr_apoe_foodintake.csv")%>%
  gather(mouse,kcal,-Week)%>%
  mutate(group = case_when(
    str_starts(mouse, "AA") ~ "Apoeadlib",
    str_starts(mouse, "AT") ~ "ApoeTRF",
    str_starts(mouse, "LA") ~ "Ldlradlib",
    str_starts(mouse, "LT") ~ "LdlrTRF",
    TRUE ~ NA_character_
  ))%>%
  arrange(mouse, Week) %>%
  group_by(mouse) %>%
  mutate(cum_kcal = cumsum(kcal)) %>%
  ungroup() %>% filter(mouse!="AA_3")%>%
  mutate(genotype = ifelse(str_starts(group, "Apoe"), "ApoeKO", "LdlrKO"),
         condition = ifelse(str_ends(group,"TRF"),"TRF","adlib"))

food_summarized <- foodint %>%group_by(group, Week) %>%
  summarise(mean_cum = mean(cum_kcal, na.rm = TRUE),
            SEM = sd(cum_kcal, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")%>%
  mutate(genotype = ifelse(str_starts(group, "Apoe"), "ApoeKO", "LdlrKO"),
         condition = ifelse(str_ends(group,"TRF"),"TRF","adlib"))

foodint_a<-foodint%>%filter(genotype=="ApoeKO")
foodint_l<-foodint%>%filter(genotype=="LdlrKO")


#cumulative food intake
food <- ggplot(food_summarized, aes(x= Week, y= mean_cum))+
  geom_line(aes(color=condition), size=0.4)+
  geom_errorbar(aes(ymin=mean_cum-SEM, ymax=mean_cum+SEM, color=condition),
                width=.2, padding=0.2)+
  geom_point(aes(fill=condition,shape=genotype), colour="black",size=2) +
  facet_wrap(~genotype, scales="free") +
  labs(title="Average Cumulative Food\nIntake vs Time") +
  ylab("Food Intake (kcal/mouse)")+
  xlab("Weeks")+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0),limits=c(0,2300),breaks = seq(0, 2300, by = 500))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,25),breaks = seq(0, 25, by = 5))+
  scale_fill_manual(values = c("#8c471f","#278c43")) +
  scale_color_manual(values = c("#8c471f","#278c43")) +
  scale_shape_manual(values = c(22,21))+
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top")

ggsave("cumul_food_intake_noaa3.pdf", plot=food,height=4, width=6.5)

foodint_m1 <- aov(cum_kcal ~condition * Week, data = foodint_a)
summary(foodint_m1)

LSD_Test<- (LSD.test(foodint_m1, c("condition", "Week")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition","Week"))
write.table(LSD_Test,"SFR26_0323_foodintvweek_apoe_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)

foodint_m2 <- aov(cum_kcal ~condition * Week, data = foodint_l)
summary(foodint_m2)

LSD_Test<- (LSD.test(foodint_m2, c("condition", "Week")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition","Week"))
write.table(LSD_Test,"SFR26_0323_foodintvweek_ldlr_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)

###########################################################
#cholesterol
serum<-fread("2026_reprocess_16s/data/2021/Lesions/apoe_ldlr_serum_chol_tg.csv")

chol<-serum%>%filter(metric=="CHOL")
summary_chol <- chol %>%
  group_by(genotype,condition) %>%
  summarise(mean = mean(dilution),SEM = sd(dilution)/sqrt(n()),.groups = "drop")

chol_a<-chol%>%filter(genotype=="ApoeKO")
chol_l<-chol%>%filter(genotype=="LdlrKO")

p <- ggplot(chol, aes(x = genotype, y = dilution,color=condition,group=condition)) +
  geom_jitter(aes(color = condition, shape = genotype),width = 0.12,size = 3) +
  geom_errorbar(data = summary_chol,aes(y = mean, ymin = mean - SEM, ymax = mean + SEM,color = condition),
                width = 0.15,size = 0.7,position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("#8c471f","#278c43")) +
  scale_shape_manual(values = c(15, 16)) +
  scale_y_continuous(expand = c(0,0),limits=c(0,2500))+
  labs(y = "mg/dl", x = NULL) +
  theme_classic()

pairwise.t.test(chol_a$dilution,chol_a$condition,p.adjust.method = "fdr",pool.sd = TRUE)
# adlib
# TRF 0.15
pairwise.t.test(chol_l$dilution,chol_l$condition,p.adjust.method = "fdr",pool.sd = TRUE)
# adlib
# TRF 0.017
###########################################################
#trilglycerides

tg<-serum%>%filter(metric=="TG")
summary_tg <- tg %>%
  group_by(genotype,condition) %>%
  summarise(mean = mean(dilution),SEM = sd(dilution)/sqrt(n()),.groups = "drop")

tg_a<-tg%>%filter(genotype=="ApoeKO")
tg_l<-tg%>%filter(genotype=="LdlrKO")

p <- ggplot(tg, aes(x = genotype, y = dilution,color=condition,group=condition)) +
  geom_jitter(aes(color = condition, shape = genotype),width = 0.12,size = 3) +
  geom_errorbar(data = summary_tg,aes(y = mean, ymin = mean - SEM, ymax = mean + SEM,color = condition),
                width = 0.15,size = 0.7,position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("#8c471f","#278c43")) +
  scale_shape_manual(values = c(15, 16)) +
  scale_y_continuous(expand = c(0,0),limits=c(0,600))+
  labs(y = "mg/dl", x = NULL) +
  theme_classic()

pairwise.t.test(tg_a$dilution,tg_a$condition,p.adjust.method = "fdr",pool.sd = TRUE)
# adlib
# TRF 0.92
pairwise.t.test(tg_l$dilution,tg_l$condition,p.adjust.method = "fdr",pool.sd = TRUE)
# adlib
# TRF 0.29
###########################################################
#lesion
les<-fread("2026_reprocess_16s/data/2021/Lesions/apoe_ldlr_lesion_data.csv")

summary_les <- les %>%
  group_by(genotype,region, condition) %>%
  summarise(mean = mean(percent),SEM = sd(percent)/sqrt(n()),.groups = "drop")

les_aa<-les%>%filter(genotype=="ApoeKO"&region=="Aorta")
les_la<-les%>%filter(genotype=="LdlrKO"&region=="Aorta")
les_ar<-les%>%filter(genotype=="ApoeKO"&region=="Aortic_arches")
les_lr<-les%>%filter(genotype=="LdlrKO"&region=="Aortic_arches")

p <- ggplot(les, aes(x = region, y = percent)) +
  geom_jitter(aes(color = condition, shape = genotype),width = 0.12,size = 3) +
  geom_errorbar(data = summary_les,aes(y = mean, ymin = mean - SEM, ymax = mean + SEM,color = condition),
    width = 0.15,size = 0.7,position = position_dodge(width = 0.4)) +
  geom_point(data = summary_les,aes(y = mean),size = 3,position = position_dodge(width = 0.4)) +
  facet_wrap(~genotype, scales = "free_x") +
  scale_color_manual(values = c("#8c471f","#278c43")) +
  scale_shape_manual(values = c(15, 16)) +
  scale_y_continuous(expand = c(0,0),limits=c(0,100))+
  labs(y = "Lesion %", x = NULL) +
  theme_classic()

pairwise.t.test(les_aa$percent,les_aa$condition,p.adjust.method = "fdr")
# adlib
# TRF 0.71  
pairwise.t.test(les_ar$percent,les_ar$condition,p.adjust.method = "fdr")
# adlib
# TRF 0.37
pairwise.t.test(les_la$percent,les_la$condition,p.adjust.method = "fdr")
# adlib
# TRF 0.04 
pairwise.t.test(les_lr$percent,les_lr$condition,p.adjust.method = "fdr")
# adlib
# TRF 0.012


