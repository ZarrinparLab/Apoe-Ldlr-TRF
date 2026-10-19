setwd("/mnt/zarrinpar/scratch/sfloresr/maria_apoe/")

library(MetaCycle)
library(tidyverse)
library(data.table)
library(ggvenn)
library(VennDiagram) 
library(grDevices)
library(gplots)
library(ggpubr)
library(BioVenn) 

######################################################################
metadata<-fread("2026_reprocess_16s/data/metadata_2026.tsv")%>%
  mutate(zt = parse_number(ZT_time))%>%
  filter(host_age==11.5)%>%
  arrange(zt)

feattab<-fread("2026_reprocess_16s/data/gg2_feature-table_filt11400_rclr/rclr_table/feature-table.tsv")%>%
  dplyr::rename(FeatureID=`#OTU ID`)
######################################################################

#apoe-adlib
mdAA<-metadata%>%filter(title2=="Apoe_adlib_11.5")

mtbAA<-feattab%>%dplyr::select(FeatureID, any_of(mdAA$sample_name))

mtbAA_filtered <- mtbAA %>%
  pivot_longer(-FeatureID,
               names_to = "sample_name",
               values_to = "value") %>%
  left_join(mdAA %>% select(sample_name, zt),
            by = "sample_name") %>%
  group_by(FeatureID, zt) %>%
  summarise(has_value = any(!is.na(value)), .groups = "drop") %>%
  group_by(FeatureID) %>%
  summarise(keep = all(has_value), .groups = "drop") %>%
  filter(keep) %>%
  select(FeatureID) %>%
  inner_join(mtbAA, by = "FeatureID")
write.table(mtbAA,"metacycle/apoe_adlib_otufeaturetable_rclr.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#apoe-TRF
mdAT<-metadata%>%filter(title2=="Apoe_TRF_11.5")

mtbAT<-feattab%>%dplyr::select(FeatureID, any_of(mdAT$sample_name))

mtbAT_filtered <- mtbAT %>%
  pivot_longer(-FeatureID,
               names_to = "sample_name",
               values_to = "value") %>%
  left_join(mdAT %>% select(sample_name, zt),
            by = "sample_name") %>%
  group_by(FeatureID, zt) %>%
  summarise(has_value = any(!is.na(value)), .groups = "drop") %>%
  group_by(FeatureID) %>%
  summarise(keep = all(has_value), .groups = "drop") %>%
  filter(keep) %>%
  select(FeatureID) %>%
  inner_join(mtbAT, by = "FeatureID")
write.table(mtbAT,"metacycle/apoe_TRF_otufeaturetable_rclr.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#ldlr-adlib
mdLA<-metadata%>%filter(title2=="Ldlr_adlib_11.5")

mtbLA<-feattab%>%dplyr::select(FeatureID, any_of(mdLA$sample_name))

mtbLA_filtered <- mtbLA %>%
  pivot_longer(-FeatureID,
               names_to = "sample_name",
               values_to = "value") %>%
  left_join(mdLA %>% select(sample_name, zt),
            by = "sample_name") %>%
  group_by(FeatureID, zt) %>%
  summarise(has_value = any(!is.na(value)), .groups = "drop") %>%
  group_by(FeatureID) %>%
  summarise(keep = all(has_value), .groups = "drop") %>%
  filter(keep) %>%
  select(FeatureID) %>%
  inner_join(mtbLA, by = "FeatureID")
write.table(mtbLA,"metacycle/ldlr_adlib_otufeaturetable_rclr.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#ldlr-TRF
mdLT<-metadata%>%filter(title2=="Ldlr_TRF_11.5")

mtbLT<-feattab%>%dplyr::select(FeatureID, any_of(mdLT$sample_name))

mtbLT_filtered <- mtbLT %>%
  pivot_longer(-FeatureID,
               names_to = "sample_name",
               values_to = "value") %>%
  left_join(mdAA %>% select(sample_name, zt),
            by = "sample_name") %>%
  group_by(FeatureID, zt) %>%
  summarise(has_value = any(!is.na(value)), .groups = "drop") %>%
  group_by(FeatureID) %>%
  summarise(keep = all(has_value), .groups = "drop") %>%
  filter(keep) %>%
  select(FeatureID) %>%
  inner_join(mtbLT, by = "FeatureID")
write.table(mtbLT,"metacycle/ldlr_TRF_otufeaturetable_rclr.txt",sep = "\t",row.names = FALSE,quote=FALSE)
#################################################
#run metacycle

#on apoe_adlib
meta2d(infile="metacycle/apoe_adlib_otufeaturetable_rclr.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="metacycle/apoe_adlib_metacycle",timepoints=rep(seq(4, 24, by=4),times=c(4,5,5,5,4,4)),minper=20,maxper=24)

#on apoe_TRF
meta2d(infile="metacycle/apoe_TRF_otufeaturetable_rclr.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="metacycle/apoe_TRF_metacycle",timepoints=rep(seq(4, 24, by=4),times=c(4,5,5,5,4,3)),minper=20,maxper=24)

#on ldlr_adlib
meta2d(infile="metacycle/ldlr_adlib_otufeaturetable_rclr.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="metacycle/ldlr_adlib_metacycle",timepoints=rep(seq(4, 24, by=4),times=c(4,4,4,4,4,4)),minper=20,maxper=24)

#on ldlr_TRF
meta2d(infile="metacycle/ldlr_TRF_otufeaturetable_rclr.txt", cycMethod=c("LS"),filestyle="txt", 
       outdir="metacycle/ldlr_TRF_metacycle",timepoints=rep(seq(4, 24, by=4),times=c(4,4,4,4,4,4)),minper=20,maxper=24)

#too many zeroes, will need to try cosinor

# 
# library(MetaCycle)
# 
# dat <- read.table(
#   "metacycle/apoe_adlib_otufeaturetable_rclr.txt",
#   header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
# )
# 
# rownames(dat) <- dat$FeatureID
# dat$FeatureID <- NULL
# 
# # force numeric (keeps NA as NA)
# dat[] <- lapply(dat, function(x) as.numeric(as.character(x)))
# 
# # KEEP only features with enough observed values across the 30 samples
# # (try 12–18; pick based on how strict you want to be)
# keep <- rowSums(!is.na(dat)) >= 3
# dat2 <- dat[keep, , drop = FALSE]
# 
# outf <- "metacycle/apoe_adlib_otufeaturetable_rclr.filtered.txt"
# write.table(
#   cbind(FeatureID = rownames(dat2), dat2),
#   file = outf, sep = "\t", quote = FALSE, row.names = FALSE
# )
# 
# meta2d(
#   infile = outf,
#   filestyle = "txt",
#   outdir = "metacycle/apoe_adlib_metacycle",
#   cycMethod = c("JTK","LS"),
#   timepoints = rep(seq(4, 24, by = 4), times = rep(5, 6)),
#   minper = 20, maxper = 24
# )

######################################################################
#read metacycle hits (early)

apoeadlib_metacyc<-fread("meta2d_11.5 Apoe_adlib_input_metacycle.csv")
apoetrf_metacyc<-fread("meta2d_11.5 Apoe_TRF_input_metacycle.csv")
ldlradlib_metacyc<-fread("meta2d_11.5 Ldlr_adlib_input_metacycle.csv")
ldlrtrf_metacyc<-fread("meta2d_11.5 Ldlr_TRF_input_metacycle.csv")
adlibw20_metacyc<-fread("meta2d_29.5 Ldlr_adlib_input_metacycle.csv")
trfw20_metacyc<-fread("meta2d_29.5 Ldlr_TRF_input_metacycle.csv")
######################################################################
#make a bar plot summarizing the number of cycling and non-cycling hits

sigAA<-apoeadlib_metacyc%>%filter(JTK_pvalue<0.05) #1
sigAT<-apoetrf_metacyc%>%filter(JTK_pvalue<0.05) #4
sigLA<-ldlradlib_metacyc%>%filter(JTK_pvalue<0.05) #3
sigLT<-ldlrtrf_metacyc%>%filter(JTK_pvalue<0.05) #8

cyc_summ<-data.frame(condition=c("Apoe_adlib","Apoe_TRF","Ldlr_adlib","Ldlr_TRF"),
                     cycling_n=c(1,4,3,8),
                     total_n=c(483,464,388,412))%>%
  mutate(not_cycling_n=total_n-cycling_n,
         cycling_perc=(cycling_n/total_n)*100,
         not_cycling_perc=(not_cycling_n/total_n)*100)
write.table(cyc_summ,"SFR26_0206_apoeldr_w1_cyclingsumm.txt",sep = "\t",row.names = FALSE, quote=FALSE)

cyc_summ_plot<-cyc_summ%>%gather(cycling_grp, value, -condition)%>%
  filter(cycling_grp=="cycling_perc"|cycling_grp=="not_cycling_perc")%>%
  mutate(condition=factor(condition,levels=c("Apoe_adlib","Apoe_TRF","Ldlr_adlib","Ldlr_TRF")),
         genotype = sub("_.*$", "", condition),
         feed = sub("^.*_", "", condition))%>%
  filter(cycling_grp=="cycling_perc")

p<-ggplot(data = cyc_summ_plot,aes(x = genotype, y = value, fill = feed))+ 
   geom_bar(stat='identity', position='dodge')+ theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0),limits=c(0,2),breaks = seq(0, 2, by = 0.5))+
  labs(y="% cycling ASV") +
  scale_fill_manual(values=c("#8c471f","#278c43","#8c471f","#278c43"))

ggsave("SFR26_0206_apoeldr_w1_cycling_summ_barplot.pdf", plot=p,height=3, width=2.2)

apoe_chi<-cyc_summ%>%dplyr::select(1:2,4)%>%
  #filter(condition=="Apoe_adlib"|condition=="Apoe_TRF")%>%
  #filter(condition=="Ldlr_adlib"|condition=="Ldlr_TRF")%>%
  #filter(condition=="Ldlr_adlib"|condition=="Apoe_adlib")%>%
  filter(condition=="Ldlr_TRF"|condition=="Apoe_TRF")%>%
  column_to_rownames("condition")%>%as.matrix()

chisq.test(apoe_chi,simulate.p.value=TRUE, B=2000) 
#X-squared = 1.9334, df = NA, p-value = 0.2199
#X-squared = 2.0121, df = NA, p-value = 0.2369
#X-squared = 1.5087, df = NA, p-value = 0.3198
#X-squared = 1.8828, df = NA, p-value = 0.2514

fisher.test(apoe_chi)
#p-value = 0.2087
#p-value = 0.2257
#p-value = 0.329
#p-value = 0.2445


######################################################################
#make a bar plot summarizing the number of cycling and non-cycling hits (reads)

summ_tab<-fread("OTU Cycling.csv")%>%
  filter(grepl("reads", cycling))%>%
  filter(grepl("early", allinfo))%>%
  mutate(not_cycling_n=totals-counts,
         not_cycling_perc=(not_cycling_n/totals)*100)%>%
  dplyr::rename(condition=allinfo,cycling_n=counts,total_n=totals,cycling_perc=percent)%>%
  dplyr::select(condition,cycling_n,total_n,not_cycling_n,cycling_perc,not_cycling_perc)%>%
  mutate(condition = sub("_early$", "", condition),
    condition = paste0(toupper(substr(condition, 1, 1)),substr(condition, 2, nchar(condition))))

cyc_summ_plot<-summ_tab%>%gather(cycling_grp, value, -condition)%>%
  filter(cycling_grp=="cycling_perc"|cycling_grp=="not_cycling_perc")%>%
  mutate(condition=factor(condition,levels=c("Apoe_adlib","Apoe_TRF","Ldlr_adlib","Ldlr_TRF")),
         genotype = sub("_.*$", "", condition),
         feed = sub("^.*_", "", condition))%>%
  filter(cycling_grp=="cycling_perc")

p<-ggplot(data = cyc_summ_plot,aes(x = genotype, y = value, fill = feed))+ 
  geom_bar(stat='identity', position='dodge')+ theme_classic() +
  theme(legend.position = "top")+
  #scale_y_continuous(expand=c(0,0),limits=c(0,2),breaks = seq(0, 2, by = 0.5))+
  labs(y="% cycling reads") +
  scale_fill_manual(values=c("#8c471f","#278c43","#8c471f","#278c43"))

ggsave("SFR26_0206_apoeldr_w1_cycling_summ_barplot_reads.pdf", plot=p,height=3, width=2.2)

######################################################################

taxa<-fread("taxonomy.tsv")%>%dplyr::rename(CycID=`Feature ID`)

sigLA_annot<-sigLA%>%
  left_join(.,taxa,by="CycID")

sigLT_annot<-sigLT%>%
  left_join(.,taxa,by="CycID")

######################################################################
#make a bar plot summarizing the number of cycling and non-cycling hits (reads)

sigLAe<-ldlradlib_metacyc%>%filter(JTK_pvalue<0.05) #3
sigLTe<-ldlrtrf_metacyc%>%filter(JTK_pvalue<0.05) #8
sigLAl<-adlibw20_metacyc%>%filter(JTK_pvalue<0.05) #2
sigLTl<-trfw20_metacyc%>%filter(JTK_pvalue<0.05) #15

cyc_summ<-data.frame(condition=c("adlib_early","adlib_late","TRF_early","TRF_late"),
                     cycling_n=c(3,8,2,15),
                     total_n=c(388,412,360,324))%>%
  mutate(not_cycling_n=total_n-cycling_n,
         cycling_perc=(cycling_n/total_n)*100,
         not_cycling_perc=(not_cycling_n/total_n)*100)
write.table(cyc_summ,"SFR26_0206_ldr_w1w20_cyclingsumm.txt",sep = "\t",row.names = FALSE, quote=FALSE)

cyc_summ_plot<-cyc_summ%>%gather(cycling_grp, value, -condition)%>%
  filter(cycling_grp=="cycling_perc"|cycling_grp=="not_cycling_perc")%>%
  mutate(condition=factor(condition,levels=c("adlib_early","adlib_late","TRF_early","TRF_late")),
         feed = sub("_.*$", "", condition),
         time = sub("^.*_", "", condition))%>%
  filter(cycling_grp=="cycling_perc")

p<-ggplot(data = cyc_summ_plot,aes(x = feed, y = value, fill = time))+ 
  geom_bar(stat='identity', position='dodge')+ theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0),limits=c(0,5),breaks = seq(0, 5, by = 1))+
  labs(y="% cycling ASV") +
  scale_fill_manual(values=c("#8c471f","#278c43","#8c471f","#278c43"))

ggsave("SFR26_0206_apoeldr_w1_cycling_summ_barplot.pdf", plot=p,height=3, width=2.2)

apoe_chi<-cyc_summ%>%dplyr::select(1:2,4)%>%
  #filter(condition=="adlib_early"|condition=="adlib_late")%>%
  #filter(condition=="TRF_early"|condition=="TRF_late")%>%
  #filter(condition=="adlib_early"|condition=="TRF_early")%>%
  filter(condition=="adlib_late"|condition=="TRF_late")%>%
  column_to_rownames("condition")%>%as.matrix()

chisq.test(apoe_chi,simulate.p.value=TRUE, B=2000) 
#X-squared = 2.0121, df = NA, p-value = 0.2414
#X-squared = 11.678, df = NA, p-value = 0.0009995
#X-squared = 0.13322, df = NA, p-value = 1
#X-squared = 4.3284, df = NA, p-value = 0.05347
