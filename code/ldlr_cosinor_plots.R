setwd("/mnt/zarrinpar/scratch/sfloresr/maria_apoe/")

library(tidyverse)
library(data.table)
library(ggvenn)
library(VennDiagram) 
library(grDevices)
library(gplots)
library(ggpubr)
library(BioVenn)
library(ggdendro)

######################################################################
#late apoe vs ldlr
#note that only the patten of the node rclr makes sense here, otherwise the cycling for apoe ad lib is highest which makes no sense
metadata<-fread("2026_reprocess_16s/data/metadata_2026.tsv")%>%
  mutate(zt = parse_number(ZT_time))

# cosAAl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_AAl.csv")
# cosATl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_ATl.csv")
# cosLAl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LAl.csv")
# cosLTl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LTl.csv")

cosAAl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_AAl_nn.csv")
cosATl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_ATl_nn.csv")
cosLAl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LAl_nn.csv")
cosLTl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LTl_nn.csv")

# cosAAl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_AAl_notnorm.csv")
# cosATl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_ATl_notnorm.csv")
# cosLAl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LAl_notnorm.csv")
# cosLTl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LTl_notnorm.csv")

sigcosAAl<-cosAAl%>%filter(q<0.05) #38
sigcosATl<-cosATl%>%filter(q<0.05) #28
sigcosLAl<-cosLAl%>%filter(q<0.05) #43
sigcosLTl<-cosLTl%>%filter(q<0.05) #97
write.table(sigcosAAl,"2026_reprocess_16s/data/cosinor/node_analysis/cosinor_AAl_qless0.05.csv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(sigcosATl,"2026_reprocess_16s/data/cosinor/node_analysis/cosinor_ATl_qless0.05.csv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(sigcosLAl,"2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LAl_qless0.05.csv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(sigcosLTl,"2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LTl_qless0.05.csv",sep = "\t",row.names = FALSE, quote=FALSE)

sigcosAAl<-cosAAl%>%filter(p<0.05) #157 #0 #0
sigcosATl<-cosATl%>%filter(p<0.05) #41 #0 #0
sigcosLAl<-cosLAl%>%filter(p<0.05) #56 #3 #0
sigcosLTl<-cosLTl%>%filter(p<0.05) #142 #72 #15

sigcosAAl<-cosAAl%>%filter(p<0.05) #94 #0 #0
sigcosATl<-cosATl%>%filter(p<0.05) #10 #0 #0
sigcosLAl<-cosLAl%>%filter(p<0.05) #27 #1 #0
sigcosLTl<-cosLTl%>%filter(p<0.05) #54 #6 #0

cyc_summ<-data.frame(condition=c("Apoe_adlib","Apoe_TRF","Ldlr_adlib","Ldlr_TRF"),
                     # cycling_n=c(94,10,27,54),
                     # total_n=c(450,459,404,379))%>%
                     cycling_n=c(157,41,56,142),
                     total_n=c(822,917,793,773))%>%
                     # cycling_n=c(38,28,43,97),
                     # total_n=c(1100,1254,1020,994))%>%
  mutate(not_cycling_n=total_n-cycling_n,
         cycling_perc=(cycling_n/total_n)*100,
         not_cycling_perc=(not_cycling_n/total_n)*100)
write.table(cyc_summ,"2026_reprocess_16s/data/cosinor/SFR26_0225_apoeldr_w20_cyclingsumm.txt",sep = "\t",row.names = FALSE, quote=FALSE)

cyc_summ_plot<-cyc_summ%>%gather(cycling_grp, value, -condition)%>%
  filter(cycling_grp=="cycling_perc"|cycling_grp=="not_cycling_perc")%>%
  mutate(condition=factor(condition,levels=c("Apoe_adlib","Apoe_TRF","Ldlr_adlib","Ldlr_TRF")),
         genotype = sub("_.*$", "", condition),
         feed = sub("^.*_", "", condition))%>%
  filter(cycling_grp=="cycling_perc")

p<-ggplot(data = cyc_summ_plot,aes(x = genotype, y = value, fill = feed))+ 
  geom_bar(stat='identity', position='dodge')+ theme_classic() +
  theme(legend.position = "top")+
  #scale_y_continuous(expand=c(0,0),limits=c(0,12),breaks = seq(0, 15, by = 2))+
  labs(y="% cycling ASV") +
  scale_fill_manual(values=c("#8c471f","#278c43","#8c471f","#278c43"))

ggsave("2026_reprocess_16s/data/cosinor/SFR26_0225_apoeldr_w20_cycling_summ_barplot.pdf", plot=p,height=3, width=2.2)


apoe_chi<-cyc_summ%>%dplyr::select(1:2,4)%>%
  filter(condition=="Apoe_adlib"|condition=="Apoe_TRF")%>%
  #filter(condition=="Ldlr_adlib"|condition=="Ldlr_TRF")%>%
  #filter(condition=="Ldlr_adlib"|condition=="Apoe_adlib")%>%
  #filter(condition=="Ldlr_TRF"|condition=="Apoe_TRF")%>%
  column_to_rownames("condition")%>%as.matrix()

chisq.test(apoe_chi,simulate.p.value=TRUE, B=2000) 
#X-squared = 3.2094, df = NA, p-value = 0.07696
#X-squared = 23.912, df = NA, p-value = 0.0004998
#X-squared = 0.83437, df = NA, p-value = 0.3738
#X-squared = 59.802, df = NA, p-value = 0.0004998

######################################################################
#late venn diagrams
list_venn <- list(ldlr_adlib = sigcosLAl$test,
                  ldlr_trf = sigcosLTl$test,
                  apoe_adlib = sigcosAAl$test,
                  apoe_trf = sigcosATl$test)

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

shared_cyc<-plyr::ldply (all, data.frame)%>%
  dplyr::rename(pattern=".id", test="X..i..")
write.table(shared_cyc,"2026_reprocess_16s/data/cosinor/late_apoeldlr_overlap.txt",sep = "\t",row.names = FALSE,quote=FALSE)

venn_plot <- venn.diagram(
  x = list_venn,
  filename = NULL,   # <- critical
  fill = c("#cd863a","#8c471f","#9bcb3e","#278c43"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.2,
  cat.col = c("#cd863a","#8c471f","#9bcb3e","#278c43"),
  main = ""
)

pdf("2026_reprocess_16s/data/cosinor/SFR26_0226_latevenn.pdf",
    width = 8, height = 8)
grid.draw(venn_plot)
dev.off()

venn_plot <- venn.diagram(
  x = list(
    # adlib = sigcosLAl$test,
    # TRF   = sigcosLTl$test
    adlib = sigcosAAl$test,
    TRF   = sigcosATl$test
  ),
  filename = NULL,
  fill = c("#8c471f", "#278c43"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("#8c471f", "#278c43"),
  main = ""
)

pdf("2026_reprocess_16s/data/cosinor/SFR26_0226_apoelate.pdf",
    width = 6, height = 6)
# pdf("2026_reprocess_16s/data/cosinor/SFR26_0226_ldlrlate.pdf",
#     width = 6, height = 6)

grid.draw(venn_plot)
dev.off()

######################################################################
metadata<-fread("2026_reprocess_16s/data/metadata_2026.tsv")%>%
  mutate(zt = parse_number(ZT_time))%>%
  filter(host_age==29.5)
ft<-fread("2026_reprocess_16s/data/gg2_feature-table_filt11400_rclr_phylo//rclr_table/feature-table.tsv")%>%
  dplyr::rename(FeatureID=`#OTU ID`)%>%dplyr::select(FeatureID,any_of(metadata$sample_name))

sigcosLTl<-cosLTl%>%filter(q<0.05)%>%
  dplyr::rename(FeatureID=test)
sigcosATl<-cosLTl%>%filter(q<0.05)%>%
  dplyr::rename(FeatureID=test)

#just the TRF overlaps
node_name<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LTl_qless0.05_nodeannot_species.csv")%>%
  dplyr::rename(FeatureID=node)
trfs_ldlr<-shared_cyc%>%filter(pattern=="ldlr_trf:apoe_trf")%>%
  dplyr::rename(FeatureID=test)%>%
  left_join(.,node_name,by="FeatureID")%>%
  left_join(.,sigcosLTl,by="FeatureID")%>%
  mutate(genotype="Ldlr")

ldlrtrf_heatmap<-ft %>%
  #filter(FeatureID %in% sigcosLTl$FeatureID)%>%
  filter(FeatureID %in% trfs_ldlr$FeatureID)%>%
  gather(sample_name,rclr,-FeatureID)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  left_join(.,metadata,by="sample_name")%>%
  group_by(FeatureID,title2,zt)%>%
  dplyr::summarise(mn_rclr=mean(rclr,na.rm=TRUE))%>%
  group_by(FeatureID,title2)%>%mutate(Zscore=(mn_rclr - mean(mn_rclr))/sd(mn_rclr))%>%
  mutate(zt=factor(zt, levels = c("4","8","12","16","20","24")))%>%
  #left_join(sigcosLTl %>% dplyr::select(FeatureID, `acrophase[h]`),by = "FeatureID")
  left_join(trfs_ldlr %>% dplyr::select(FeatureID, `acrophase[h]`),by = "FeatureID")

# Create ordered factor by acrophase
feature_order <- trfs_ldlr %>%
#feature_order <- sigcosLTl %>%
  arrange(`acrophase[h]`) %>%
  pull(FeatureID)

ldlrtrf_heatmap$FeatureID <- factor(
  ldlrtrf_heatmap$FeatureID,
  levels = feature_order,
  ordered = TRUE
)

plt<-ggplot(ldlrtrf_heatmap,aes(x=zt, y=FeatureID)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(~title2,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),strip.text.y = element_text(angle = 0))+
  # scale_fill_gradient2(low = "#075AFF",
  #                      mid = "#FFFFCC",
  #                      high = "#FF0000")
  scale_fill_gradient(low = "#00008B",
                       high = "#FFD700")
ggsave("2026_reprocess_16s/data/cosinor/node_analysis/ldlr_trf_97_heatmap.pdf", plot=plt,height=10, width=20)
ggsave("2026_reprocess_16s/data/cosinor/node_analysis/ldlrapoe_trf_11_heatmap.pdf", plot=plt,height=2, width=7)

##########################################################
#phase

trfs_ldlr<-shared_cyc%>%filter(pattern=="ldlr_trf:apoe_trf")%>%
  dplyr::rename(FeatureID=test)%>%
  left_join(.,node_name,by="FeatureID")%>%
  left_join(.,sigcosLTl,by="FeatureID")%>%
  mutate(genotype="Ldlr")

trfs_aploe<-shared_cyc%>%filter(pattern=="ldlr_trf:apoe_trf")%>%
  dplyr::rename(FeatureID=test)%>%
  left_join(.,node_name,by="FeatureID")%>%
  left_join(.,sigcosATl,by="FeatureID")%>%
  mutate(genotype="Apoe")

trfs<-rbind(trfs_ldlr,trfs_aploe)

p<-ggplot(trfs, aes(x=`acrophase[h]`,fill=genotype)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0), limits = c(0,2.5)) +
  scale_x_continuous(expand=c(0,0),limits = c(0,24)) +
  scale_fill_manual(values=c("#ffffff","#59595b"))
ggsave("2026_reprocess_16s/data/cosinor/node_analysis/ldlrapoe_trf_11_acro_density.pdf",p, width = 4.5, height = 4)

trfs_ldlr_97<-sigcosLTl%>%
  left_join(.,node_name,by="FeatureID")%>%
  mutate(title2="Ldlr_TRF")

p<-ggplot(trfs_ldlr_97, aes(x=`acrophase[h]`,fill=title2)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0), limits = c(0,0.5)) +
  scale_x_continuous(expand=c(0,0),limits = c(0,24)) +
  scale_fill_manual(values=c("#278c43"))
ggsave("2026_reprocess_16s/data/cosinor/node_analysis/ldlrTRF_trf_97_acro_density.pdf",p, width = 4.5, height = 4)

##########################################################
#plot unique over time
sem <- function(x) sd(x)/sqrt(length(x))
metadata<-fread("2026_reprocess_16s/data/metadata_2026.tsv")%>%
  mutate(zt = parse_number(ZT_time))%>%
  filter(host_age==29.5)

sigcosLTl<-cosLTl%>%filter(q<0.05)%>%
  dplyr::rename(FeatureID=test)

mtbZT<-fread("2026_reprocess_16s/data/gg2_feature-table_filt11400_rclr_phylo//rclr_table/feature-table.tsv")%>%
  dplyr::rename(FeatureID=`#OTU ID`)%>%dplyr::select(FeatureID,any_of(metadata$sample_name))%>%
  gather(sample_name,rclr,-FeatureID)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  filter(FeatureID %in% trfs_aploe$FeatureID)%>% 
  left_join(.,metadata,by="sample_name")# %>%
  #group_by(title2,genotype,feeding_type,zt)%>%summarise(mn_count=mean(rclr,na.rm=TRUE),sem=sem(rclr))

p<-mtbZT%>%
  ggplot(aes(x=zt, y=rclr, color=feeding_type,fill=feeding_type)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point", size = 2) + 
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  #ggplot(aes(x=zt, y=mn_count, color=feeding_type)) +
  #geom_point(alpha=1.0) + geom_line() +
  theme_pubr() +
  facet_wrap(~genotype)+
  scale_color_manual(values=c("#8c471f", "#278c43"))+
  scale_fill_manual(values=c("#8c471f", "#278c43"))+
  scale_x_continuous(breaks=c(4,8,12,16,20,24))+
  #geom_ribbon(aes(ymin = mn_count-sem, ymax = mn_count+sem, fill=feeding_type),alpha=0.3,colour = NA)+
  ggtitle("shared nodes apoe and ldlr TRF")
ggsave("2026_reprocess_16s/data/cosinor/node_analysis/ldlrapoe_trf_11_overZT.pdf",p, width = 6.5, height = 4)

mtbZT<-fread("2026_reprocess_16s/data/gg2_feature-table_filt11400_rclr_phylo//rclr_table/feature-table.tsv")%>%
  dplyr::rename(FeatureID=`#OTU ID`)%>%dplyr::select(FeatureID,any_of(metadata$sample_name))%>%
  gather(sample_name,rclr,-FeatureID)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  #filter(FeatureID %in% sigcosLTl$FeatureID)%>% 
  filter(FeatureID %in% trfs_aploe$FeatureID)%>% 
  left_join(.,metadata,by="sample_name") #%>%
  #group_by(FeatureID,title2,genotype,feeding_type,zt)%>%summarise(mn_count=mean(rclr,na.rm=TRUE),sem=sem(rclr))


ZT_dist <- function(data,node) {
  
  data%>%
    ggplot(aes(x=zt, y=rclr, color=feeding_type,fill=feeding_type)) +
    stat_summary(fun = mean, geom = "line") +
    stat_summary(fun = mean, geom = "point", size = 2) + 
    stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
    #ggplot(aes(x=zt, y=mn_count, color=feeding_type)) +
    #geom_point(alpha=1.0) + geom_line() +
    theme_pubr() +
    #facet_wrap(~genotype)+
    facet_wrap(~genotype,scales="free")+
    scale_color_manual(values=c("#8c471f", "#278c43"))+
    scale_fill_manual(values=c("#8c471f", "#278c43"))+
    scale_x_continuous(breaks=c(4,8,12,16,20,24))+
    #geom_ribbon(aes(ymin = mn_count-sem, ymax = mn_count+sem, fill=feeding_type),alpha=0.3,colour = NA)+
    labs(color="Feeding Type",
         fill="Feeding Typee",
         y ="RCLR",
         x ="ZT",
         title = node)+
    theme(plot.title = element_text(face = "bold",hjust = 0.5))
  ggsave(paste("2026_reprocess_16s/data/cosinor/node_analysis/apoe_ldlr_TRF_zt/",node,".pdf",sep=""), height=4, width=6.5)
  #ggsave(paste("2026_reprocess_16s/data/cosinor/node_analysis/ldlr_TRF_zt/",node,".pdf",sep=""), height=4, width=6.5)
}


mtbZT_nested <- mtbZT %>% 
  group_by(FeatureID) %>% 
  nest()

mtbZT_plots <- 
  mtbZT_nested %>% 
  mutate(plot = map2(data, FeatureID ,  ~ ZT_dist(.x,.y)))
##########################################################
#early vs late ldlr

cosLA<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LA.csv")
cosLT<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LT.csv")
cosLAl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LAl.csv")
cosLTl<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LTl.csv")

sigcosLA<-cosLA%>%filter(q<0.05) #13
sigcosLT<-cosLT%>%filter(q<0.05) #11
sigcosLAl<-cosLAl%>%filter(q<0.05) #43
sigcosLTl<-cosLTl%>%filter(q<0.05) #97
write.table(sigcosLA,"2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LA_qless0.05.csv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(sigcosLT,"2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LT_qless0.05.csv",sep = "\t",row.names = FALSE, quote=FALSE)


cyc_summ<-data.frame(condition=c("adlib_early","adlib_late","TRF_early","TRF_late"),
                     cycling_n=c(13,43,11,97),
                     total_n=c(1038,1020,1068,994))%>%
  mutate(not_cycling_n=total_n-cycling_n,
         cycling_perc=(cycling_n/total_n)*100,
         not_cycling_perc=(not_cycling_n/total_n)*100)
write.table(cyc_summ,"2026_reprocess_16s/data/cosinor/SFR26_0225_ldr_w1w20_cyclingsumm.txt",sep = "\t",row.names = FALSE, quote=FALSE)

cyc_summ_plot<-cyc_summ%>%gather(cycling_grp, value, -condition)%>%
  filter(cycling_grp=="cycling_perc"|cycling_grp=="not_cycling_perc")%>%
  mutate(condition=factor(condition,levels=c("adlib_early","adlib_late","TRF_early","TRF_late")),
         feed = sub("_.*$", "", condition),
         time = sub("^.*_", "", condition))%>%
  filter(cycling_grp=="cycling_perc")

p<-ggplot(data = cyc_summ_plot,aes(x = feed, y = value, fill = time))+ 
  geom_bar(stat='identity', position='dodge')+ theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0),limits=c(0,12),breaks = seq(0, 15, by = 2))+
  labs(y="% cycling ASV") +
  scale_fill_manual(values=c("#8c471f","#278c43","#8c471f","#278c43"))

ggsave("2026_reprocess_16s/data/cosinor/SFR26_0225_apoeldr_w1w20_cycling_summ_barplot.pdf", plot=p,height=3, width=2.2)

apoe_chi<-cyc_summ%>%dplyr::select(1:2,4)%>%
  #filter(condition=="adlib_early"|condition=="adlib_late")%>%
  #filter(condition=="TRF_early"|condition=="TRF_late")%>%
  #filter(condition=="adlib_early"|condition=="TRF_early")%>%
  filter(condition=="adlib_late"|condition=="TRF_late")%>%
  column_to_rownames("condition")%>%as.matrix()

chisq.test(apoe_chi,simulate.p.value=TRUE, B=2000)
#X-squared = 17.066, df = NA, p-value = 0.0004998
#X-squared = 79.029, df = NA, p-value = 0.0004998
#X-squared = 0.2312, df = NA, p-value = 0.6722
#X-squared = 23.912, df = NA, p-value = 0.0004998
######################################################################
#ldlr venn diagrams

list_venn <- list(ldlr_adlib_w1 = sigcosLA$test,
                  ldlr_trf_w1 = sigcosLT$test,
                  ldlr_adlib_w20 = sigcosLAl$test,
                  ldlr_trf_w20 = sigcosLTl$test)

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

shared_cyc<-plyr::ldply (all, data.frame)%>%
  dplyr::rename(pattern=".id", test="X..i..")
write.table(shared_cyc,"2026_reprocess_16s/data/cosinor/node_analysis/ldlr_w1w20_overlap.txt",sep = "\t",row.names = FALSE,quote=FALSE)

venn_plot <- venn.diagram(
  x = list_venn,
  filename = NULL,   # <- critical
  fill = c("#cd863a","#9bcb3e","#8c471f","#278c43"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.2,
  cat.col = c("#cd863a","#9bcb3e","#8c471f","#278c43"),
  main = ""
)

pdf("2026_reprocess_16s/data/cosinor/node_analysis/SFR26_0226_ldlrvenn.pdf",
    width = 8, height = 8)
grid.draw(venn_plot)
dev.off()
######################################################################
metadata<-fread("2026_reprocess_16s/data/metadata_2026.tsv")%>%
  mutate(zt = parse_number(ZT_time))%>%
  filter(genotype=="Ldlrtm1Her/J")
ft<-fread("2026_reprocess_16s/data/gg2_feature-table_filt11400_rclr_phylo/rclr_table/feature-table.tsv")%>%
  dplyr::rename(FeatureID=`#OTU ID`)%>%dplyr::select(FeatureID,any_of(metadata$sample_name))

sigcosLTl<-cosLTl%>%filter(q<0.05)%>%
  dplyr::rename(FeatureID=test)
sigcosLT<-cosLT%>%filter(q<0.05)%>%
  dplyr::rename(FeatureID=test)

#just the TRF overlaps
# node_name<-fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LTl_qless0.05_nodeannot_species.csv")%>%
#   dplyr::rename(FeatureID=node)
trfs_ldlr<-shared_cyc%>%filter(pattern=="ldlr_trf_w1:ldlr_trf_w20")%>%
  dplyr::rename(FeatureID=test)%>%
  #left_join(.,node_name,by="FeatureID")%>%
  left_join(.,sigcosLTl,by="FeatureID")%>%
  mutate(genotype="Ldlr") #since n=2 heatmap wont work

ldlrtrf_heatmap<-ft %>%
  filter(FeatureID %in% sigcosLTl$FeatureID)%>%
  #filter(FeatureID %in% trfs_ldlr$FeatureID)%>%
  gather(sample_name,rclr,-FeatureID)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  left_join(.,metadata,by="sample_name")%>%
  group_by(FeatureID,title2,zt)%>%
  dplyr::summarise(mn_rclr=mean(rclr,na.rm=TRUE))%>%
  group_by(FeatureID,title2)%>%mutate(Zscore=(mn_rclr - mean(mn_rclr))/sd(mn_rclr))%>%
  mutate(zt=factor(zt, levels = c("4","8","12","16","20","24")))%>%
  left_join(
    sigcosLTl %>% dplyr::select(FeatureID, `acrophase[h]`),by = "FeatureID")

# Create ordered factor by acrophase
feature_order <- sigcosLTl %>%
  arrange(`acrophase[h]`) %>%
  pull(FeatureID)

ldlrtrf_heatmap$FeatureID <- factor(
  ldlrtrf_heatmap$FeatureID,
  levels = feature_order,
  ordered = TRUE
)

plt<-ggplot(ldlrtrf_heatmap,aes(x=zt, y=FeatureID)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(~title2,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),strip.text.y = element_text(angle = 0))+
  # scale_fill_gradient2(low = "#075AFF",
  #                      mid = "#FFFFCC",
  #                      high = "#FF0000") 
  scale_fill_gradient(low = "#00008B",
                      high = "#FFD700")
ggsave("2026_reprocess_16s/data/cosinor/node_analysis/ldlr_w1w20_97_heatmap.pdf", plot=plt,height=10, width=20)
######################################################################
sem <- function(x) sd(x)/sqrt(length(x))
metadata<-fread("2026_reprocess_16s/data/metadata_2026.tsv")%>%
  mutate(zt = parse_number(ZT_time))%>%
  filter(genotype=="Ldlrtm1Her/J")

sigcosLTl<-cosLTl%>%filter(q<0.05)%>%
  dplyr::rename(FeatureID=test)

sigcosLT<-cosLT%>%filter(q<0.05)%>%
  dplyr::rename(FeatureID=test)

trfs_ldlr<-shared_cyc%>%filter(pattern=="ldlr_trf_w1:ldlr_trf_w20")%>%
  dplyr::rename(FeatureID=test)

mtbZT<-fread("2026_reprocess_16s/data/gg2_feature-table_filt11400_rclr_phylo//rclr_table/feature-table.tsv")%>%
  dplyr::rename(FeatureID=`#OTU ID`)%>%dplyr::select(FeatureID,any_of(metadata$sample_name))%>%
  gather(sample_name,rclr,-FeatureID)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  filter(FeatureID %in% trfs_ldlr$FeatureID)%>% 
  left_join(.,metadata,by="sample_name") #%>%
  #group_by(title2,host_age,feeding_type,zt)%>%summarise(mn_count=mean(rclr,na.rm=TRUE),sem=sem(rclr))

p<-mtbZT%>%
  ggplot(aes(x=zt, y=rclr, color=feeding_type,fill=feeding_type)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point", size = 2) + 
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  #ggplot(aes(x=zt, y=mn_count, color=feeding_type)) +
  #geom_point(alpha=1.0) + geom_line() +
  theme_pubr() +
  facet_wrap(~host_age)+
  scale_color_manual(values=c("#8c471f", "#278c43"))+
  scale_fill_manual(values=c("#8c471f", "#278c43"))+
  scale_x_continuous(breaks=c(4,8,12,16,20,24))+
  #geom_ribbon(aes(ymin = mn_count-sem, ymax = mn_count+sem, fill=feeding_type),alpha=0.3,colour = NA)+
  ggtitle("shared nodes ldlr TRF early late")
ggsave("2026_reprocess_16s/data/cosinor/node_analysis/ldlrtrf_w1w20_2_overZT.pdf",p, width = 6.5, height = 4)

mtbZT<-fread("2026_reprocess_16s/data/gg2_feature-table_filt11400_rclr_phylo//rclr_table/feature-table.tsv")%>%
  dplyr::rename(FeatureID=`#OTU ID`)%>%dplyr::select(FeatureID,any_of(metadata$sample_name))%>%
  gather(sample_name,rclr,-FeatureID)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  filter(FeatureID %in% sigcosLTl$FeatureID)%>% 
  #filter(FeatureID %in% sigcosLT$FeatureID)%>% 
  left_join(.,metadata,by="sample_name") #%>%
  #group_by(FeatureID,title2,host_age,feeding_type,zt)%>%summarise(mn_count=mean(rclr,na.rm=TRUE),sem=sem(rclr))


ZT_dist <- function(data,node) {
  
  data%>%
    #ggplot(aes(x=zt, y=mn_count, color=feeding_type)) +
    ggplot(aes(x=zt, y=rclr, color=feeding_type,fill=feeding_type)) +
    stat_summary(fun = mean, geom = "line") +
    stat_summary(fun = mean, geom = "point", size = 2) + 
    stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
    #geom_point(alpha=1.0) + geom_line() +
    theme_pubr() +
    #facet_wrap(~genotype)+
    facet_wrap(~host_age,scales="free")+
    scale_color_manual(values=c("#8c471f", "#278c43"))+
    scale_fill_manual(values=c("#8c471f", "#278c43"))+
    scale_x_continuous(breaks=c(4,8,12,16,20,24))+
    #geom_ribbon(aes(ymin = mn_count-sem, ymax = mn_count+sem, fill=feeding_type),alpha=0.3,colour = NA)+
    labs(color="Feeding Type",
         fill="Feeding Typee",
         y ="RCLR",
         x ="ZT",
         title = node)+
    theme(plot.title = element_text(face = "bold",hjust = 0.5))
ggsave(paste("2026_reprocess_16s/data/cosinor/node_analysis/ldlrTRF_w1w20_zt/",node,".pdf",sep=""), height=4, width=6.5)
#ggsave(paste("2026_reprocess_16s/data/cosinor/node_analysis/ldlrTRF_w1_zt/",node,".pdf",sep=""), height=4, width=6.5)
}


mtbZT_nested <- mtbZT %>% 
  group_by(FeatureID) %>% 
  nest()

mtbZT_plots <- 
  mtbZT_nested %>% 
  mutate(plot = map2(data, FeatureID ,  ~ ZT_dist(.x,.y)))
######################################################################
mtbZT<-fread("2026_reprocess_16s/data/gg2_feature-table_filt11400_rclr_phylo//rclr_table/feature-table.tsv")%>%
  dplyr::rename(FeatureID=`#OTU ID`)%>%dplyr::select(FeatureID,any_of(metadata$sample_name))%>%
  gather(sample_name,rclr,-FeatureID)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  filter(FeatureID=="n1674")%>% 
  left_join(.,metadata,by="sample_name")%>%
  mutate(ZT_time=factor(ZT_time,levels=c("ZT4","ZT8","ZT12","ZT16","ZT20","ZT24")))

p<-mtbZT%>%
  ggplot(aes(x=zt, y=rclr,fill=feeding_type)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  #geom_boxplot(alpha=0.3) +
  #geom_point(position = position_jitterdodge(dodge.width = 0.8,jitter.width = 0.15),size = 3,alpha = 0.8,shape=21) +
  theme_pubr() +
  facet_wrap(~host_age)+
  scale_color_manual(values=c("#8c471f", "#278c43"))+
  scale_fill_manual(values=c("#8c471f", "#278c43"))+
  ggtitle("n1674")
ggsave("2026_reprocess_16s/data/cosinor/node_analysis/ldlrtrf_w1w20_2_overZT.pdf",p, width = 6.5, height = 4)

##########################################################
#phase

trfs_ldlr_11<-sigcosLT%>%
  #left_join(.,node_name,by="FeatureID")%>%
  mutate(title2="Ldlr_TRF")

p<-ggplot(trfs_ldlr_11, aes(x=`acrophase[h]`,fill=title2)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0), limits = c(0,8)) +
  scale_x_continuous(expand=c(0,0),limits = c(0,24)) +
  scale_fill_manual(values=c("#9bcb3e"))
ggsave("2026_reprocess_16s/data/cosinor/node_analysis/ldlrTRF_trf_w1_11_acro_density.pdf",p, width = 4.5, height = 4)

