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
metadata<-fread("2026_reprocess_16s/data/metadata_2026.tsv")%>%
  mutate(zt = parse_number(ZT_time))

cosAAl<-fread("2026_reprocess_16s/data/cosinor/asv_analysis/cosinor_AAl.csv")
cosATl<-fread("2026_reprocess_16s/data/cosinor/asv_analysis/cosinor_ATl.csv")
cosLAl<-fread("2026_reprocess_16s/data/cosinor/asv_analysis/cosinor_LAl.csv")
cosLTl<-fread("2026_reprocess_16s/data/cosinor/asv_analysis/cosinor_LTl.csv")

# cosAAl<-fread("2026_reprocess_16s/data/cosinor/asv_analysis/cosinor_AAl_clr.csv")
# cosATl<-fread("2026_reprocess_16s/data/cosinor/asv_analysis/cosinor_ATl_clr.csv")
# cosLAl<-fread("2026_reprocess_16s/data/cosinor/asv_analysis/cosinor_LAl_clr.csv")
# cosLTl<-fread("2026_reprocess_16s/data/cosinor/asv_analysis/cosinor_LTl_clr.csv")

sigcosAAl<-cosAAl%>%filter(p<0.05) #4 #11 #11
sigcosATl<-cosATl%>%filter(p<0.05) #3 #6 #6
sigcosLAl<-cosLAl%>%filter(p<0.05) #11 #17 #16
sigcosLTl<-cosLTl%>%filter(p<0.05) #16 #20 #20

sigcosAAl<-cosAAl%>%filter(p<0.05) #1752 #1765 #1745
sigcosATl<-cosATl%>%filter(p<0.05) #10 #0 #0
sigcosLAl<-cosLAl%>%filter(p<0.05) #8 #0 #0
sigcosLTl<-cosLTl%>%filter(p<0.05) #28 #1855 #1

cyc_summ<-data.frame(condition=c("Apoe_adlib","Apoe_TRF","Ldlr_adlib","Ldlr_TRF"),
                     cycling_n=c(4,3,11,16),
                     #cycling_n=c(11,6,16,20),
                     total_n=c(458,541,399,364))%>%
  mutate(not_cycling_n=total_n-cycling_n,
         cycling_perc=(cycling_n/total_n)*100,
         not_cycling_perc=(not_cycling_n/total_n)*100)
write.table(cyc_summ,"2026_reprocess_16s/data/cosinor/asv_analysis/SFR26_0225_apoeldr_w20_cyclingsumm.txt",sep = "\t",row.names = FALSE, quote=FALSE)

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

ggsave("2026_reprocess_16s/data/cosinor/asv_analysis/SFR26_0225_apoeldr_w20_cycling_summ_barplot.pdf", plot=p,height=3, width=2.2)


apoe_chi<-cyc_summ%>%dplyr::select(1:2,4)%>%
  #filter(condition=="Apoe_adlib"|condition=="Apoe_TRF")%>%
  filter(condition=="Ldlr_adlib"|condition=="Ldlr_TRF")%>%
  #filter(condition=="Ldlr_adlib"|condition=="Apoe_adlib")%>%
  #filter(condition=="Ldlr_TRF"|condition=="Apoe_TRF")%>%
  column_to_rownames("condition")%>%as.matrix()

chisq.test(apoe_chi,simulate.p.value=TRUE, B=2000) 
#X-squared = 2.4777, df = NA, p-value = 0.1369
#X-squared = 0.93306, df = NA, p-value = 0.3883
#X-squared = 1.8076, df = NA, p-value = 0.2374
#X-squared = 14.997, df = NA, p-value = 0.0004998
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
# taxa<-fread("2026_reprocess_16s/data/modern_beta/phyloRPCA-taxonomy/taxonomy.tsv")%>%
#   dplyr::rename(FeatureID=`Feature ID`)
metadata<-fread("2026_reprocess_16s/data/metadata_2026.tsv")%>%
  mutate(zt = parse_number(ZT_time))%>%
  filter(host_age==29.5)
ft<-fread("2026_reprocess_16s/data/gg2_feature-table_filt11400_rclr/rclr_table/feature-table.tsv")%>%
  dplyr::rename(FeatureID=`#OTU ID`)%>%dplyr::select(FeatureID,any_of(metadata$sample_name))

sigcosLTl<-cosLTl%>%filter(q<0.05) #97

ldlrtrf_heatmap<-ft %>%
  filter(FeatureID %in% sigcosLTl$test)%>%
  gather(sample_name,rclr,-FeatureID)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  left_join(.,metadata,by="sample_name")%>%
  group_by(FeatureID,title2,zt)%>%
  dplyr::summarise(mn_rclr=mean(rclr,na.rm=TRUE))%>%
  group_by(FeatureID,title2)%>%mutate(Zscore=(mn_rclr - mean(mn_rclr))/sd(mn_rclr))%>%
  mutate(zt=factor(zt, levels = c("4","8","12","16","20","24")))

# Run clustering
ft_df<-ft%>%
  gather(sample_name,rclr,-FeatureID)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  filter(FeatureID %in% sigcosLTl$test)%>%
  spread(sample_name,rclr)

ft_matrix<-ft_df%>%column_to_rownames("FeatureID") %>%as.matrix()
ft_dendro <- as.dendrogram(hclust(d = dist(x = ft_matrix)))
dendro_plot <- ggdendrogram(data = ft_dendro, rotate = TRUE)

ft_order <- order.dendrogram(ft_dendro)
ldlrtrf_heatmap$FeatureID <- factor(x = ldlrtrf_heatmap$FeatureID,
                                       levels = ft_df$FeatureID[ft_order], 
                                       ordered = TRUE)

plt<-ggplot(ldlrtrf_heatmap,aes(x=zt, y=FeatureID)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(~title2,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),strip.text.y = element_text(angle = 0))+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") 
ggsave("2026_reprocess_16s/data/cosinor/ldlr_trf_97_heatmap.pdf", plot=plt,height=10, width=20)

#just the TRF overlaps

trfs<-shared_cyc%>%filter(pattern=="ldlr_trf:apoe_trf")%>%
  dplyr::rename(FeatureID=test)%>%
  left_join(.,taxa,by="FeatureID")


##########################################################
#early vs late ldlr

cosLA<-fread("2026_reprocess_16s/data/cosinor/cosinor_LA.csv")
cosLT<-fread("2026_reprocess_16s/data/cosinor/cosinor_LT.csv")
cosLAl<-fread("2026_reprocess_16s/data/cosinor/cosinor_LAl.csv")
cosLTl<-fread("2026_reprocess_16s/data/cosinor/cosinor_LTl.csv")

sigcosLA<-cosLA%>%filter(q<0.05) #13
sigcosLT<-cosLT%>%filter(q<0.05) #11
sigcosLAl<-cosLAl%>%filter(q<0.05) #43
sigcosLTl<-cosLTl%>%filter(q<0.05) #97

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


