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
library(vegan)

#functions

make_rpcaw1<-function(qfile,titlenam,shapenum){
  
  
  wk1a_df<-qfile$data$Vectors%>%
    dplyr::rename(sample_name=SampleID)%>%
    left_join(.,metadata,by="sample_name")%>%
    mutate(feed_zt=paste(feeding_type,zt,sep="-"))%>%
    mutate(feed_zt=factor(feed_zt,levels=c("adlib-4","adlib-8","adlib-12","adlib-16","adlib-20","adlib-24",
                                           "TRF-4","TRF-8","TRF-12","TRF-16","TRF-20","TRF-24")))
  
  p<-wk1a_df%>%
    ggplot(aes(x=PC1, y=PC2, fill=feed_zt, group=feed_zt)) +
    geom_point(alpha=1.0, shape=shapenum, size=4)+
    theme_pubr() +
    scale_fill_manual(values=c("#FFB347","#E69934","#CC8833","#B36B29","#995C20","#8B4513",
                               "#8CCB90","#6ABC79","#49A961","#359151","#367347","#35533E"))+
    scale_color_manual(values=c("#FFB347","#E69934","#CC8833","#B36B29","#995C20","#8B4513",
                                "#8CCB90","#6ABC79","#49A961","#359151","#367347","#35533E"))+
    theme(legend.position = "right")+
    labs(x =paste("PC1 (",round(qfile$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
         y =paste("PC2 (",round(qfile$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
         title=titlenam)
  return(p)
}
#############################################################################
metadata<-fread("2026_reprocess_16s/data/metadata_2026.tsv")%>%
  mutate(zt = parse_number(ZT_time))

w20md<-metadata%>%filter(host_age==29.5)
#############################################################################
#alpha diversity--faith PD (just week1)
faithpd<-fread("2026_reprocess_16s/data/core_diversity_r11400/faith_pd_vector/alpha-diversity.tsv")%>%
  dplyr::rename(sample_name=V1)%>%
  left_join(.,metadata,by="sample_name")%>%
  filter(host_age==29.5)%>%
  mutate(cond_dayphase=paste(feeding_type,day_phase,sep=""),
         day_phase=factor(day_phase,levels=c("daylight","dark")),
         zt = parse_number(ZT_time))

apoe_fpd<-faithpd%>%filter(genotype=="Apoetm1Unc/J")
ldlr_fpd<-faithpd%>%filter(genotype=="Ldlrtm1Her/J")

p<-ggplot(faithpd, aes(x=feeding_type, y=faith_pd, fill=day_phase)) +
  geom_boxplot(alpha=0.3) +
  geom_point(aes(shape = genotype),position = position_jitterdodge(dodge.width = 0.8,jitter.width = 0.15),size = 3,alpha = 0.8) +
  scale_fill_manual(values=c("#fbb040","#2e3192"))+
  facet_wrap(~genotype)+
  scale_shape_manual(values = c("Apoetm1Unc/J" = 22,"Ldlrtm1Her/J" = 21))+
  theme_pubr()+
  labs(x="condition",y="Faith's PD", title="alpha diversity (week 20)")+
  theme(legend.position = "right")

ggsave("alphafaithpd_wk20.pdf", plot=p,height=3.5, width=6.5)

pairwise.wilcox.test(apoe_fpd$faith_pd, apoe_fpd$cond_dayphase,p.adjust.method="fdr")

# adlibdark adlibdaylight TRFdark
# adlibdaylight 0.170     -             -      
#   TRFdark       0.353     0.019         -      
#   TRFdaylight   0.487     0.139         0.755  

pairwise.wilcox.test(ldlr_fpd$faith_pd, ldlr_fpd$cond_dayphase,p.adjust.method="fdr")

# adlibdark adlibdaylight TRFdark
# adlibdaylight 0.97      -             -      
#   TRFdark       0.68      0.72          -      
#   TRFdaylight   0.72      0.68          0.68   

faithpd <- faithpd %>%
  bind_rows(faithpd %>%filter(zt == 24) %>%mutate(zt = 0))

apoe_fpd<-faithpd%>%filter(genotype=="Apoetm1Unc/J")
ldlr_fpd<-faithpd%>%filter(genotype=="Ldlrtm1Her/J")


p<-ggplot(data=faithpd, aes(x=zt, y=faith_pd, colour=feeding_type, fill=feeding_type)) +
  facet_wrap(~genotype, scales="free") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values=c("#8c471f","#278c43"))+
  scale_color_manual(values=c("#8c471f","#278c43"))+
  scale_x_continuous(breaks = seq(0, 24, by = 4)) + 
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("alphafaithpd_byzt_week20.pdf", plot=p,height=4, width=6.5)


mod1 <- aov(faith_pd~ feeding_type * zt, data = apoe_fpd)
summary(mod1)

LSD_Test_fd1<- (LSD.test(mod1, c("feeding_type", "zt")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Group","DateTime"))
write.table(LSD_Test_fd1,"apoezt_faithpd_pvals_anovaLSD_wk20.txt",sep = "\t",row.names = FALSE,quote=FALSE)

mod2 <- aov(faith_pd~ feeding_type * zt, data = ldlr_fpd)
summary(mod2)

LSD_Test_fd2<- (LSD.test(mod2, c("feeding_type", "zt")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Group","DateTime"))
write.table(LSD_Test_fd2,"ldlrzt_faithpd_pvals_anovaLSD_wk20.txt",sep = "\t",row.names = FALSE,quote=FALSE)
#############################################################################
#alpha diversity--faith PD (ldlr w1 vs w20)
faithpd<-fread("2026_reprocess_16s/data/core_diversity_r11400/faith_pd_vector/alpha-diversity.tsv")%>%
  dplyr::rename(sample_name=V1)%>%
  left_join(.,metadata,by="sample_name")%>%
  filter(genotype=="Ldlrtm1Her/J")%>%
  mutate(cond_dayphase=paste(feeding_type,day_phase,sep=""),
         day_phase=factor(day_phase,levels=c("daylight","dark")),
         zt = parse_number(ZT_time))

early_fpd<-faithpd%>%filter(host_age==11.5)
late_fpd<-faithpd%>%filter(host_age==29.5)

p<-ggplot(faithpd, aes(x=feeding_type, y=faith_pd, fill=day_phase)) +
  geom_boxplot(alpha=0.3) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,jitter.width = 0.15),size = 3,alpha = 0.8,shape=21) +
  scale_fill_manual(values=c("#fbb040","#2e3192"))+
  facet_wrap(~host_age)+
  theme_pubr()+
  labs(x="condition",y="Faith's PD", title="alpha diversity Ldlr-/-")+
  theme(legend.position = "right")

ggsave("alphafaithpd_ldlr.pdf", plot=p,height=3.5, width=6.5)

pairwise.wilcox.test(early_fpd$faith_pd, early_fpd$cond_dayphase,p.adjust.method="fdr")

# adlibdark adlibdaylight TRFdark
# adlibdaylight 0.98      -             -      
#   TRFdark       0.98      0.98          -      
#   TRFdaylight   0.89      0.89          0.89   

pairwise.wilcox.test(late_fpd$faith_pd, late_fpd$cond_dayphase,p.adjust.method="fdr")

# adlibdark adlibdaylight TRFdark
# adlibdaylight 0.97      -             -      
#   TRFdark       0.68      0.72          -      
#   TRFdaylight   0.72      0.68          0.68 

faithpd <- faithpd %>%
  bind_rows(faithpd %>%filter(zt == 24) %>%mutate(zt = 0))

early_fpd<-faithpd%>%filter(host_age==11.5)
late_fpd<-faithpd%>%filter(host_age==29.5)


p<-ggplot(data=faithpd, aes(x=zt, y=faith_pd, colour=feeding_type, fill=feeding_type)) +
  facet_wrap(~host_age, scales="free") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values=c("#8c471f","#278c43"))+
  scale_color_manual(values=c("#8c471f","#278c43"))+
  scale_x_continuous(breaks = seq(0, 24, by = 4)) + 
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("alphafaithpd_byzt_wldlr.pdf", plot=p,height=4, width=6.5)


mod1 <- aov(faith_pd~ feeding_type * zt, data = early_fpd)
summary(mod1)

LSD_Test_fd1<- (LSD.test(mod1, c("feeding_type", "zt")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Group","DateTime"))
write.table(LSD_Test_fd1,"ldlr_early_zt_faithpd_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)

mod2 <- aov(faith_pd~ feeding_type * zt, data = late_fpd)
summary(mod2)

LSD_Test_fd2<- (LSD.test(mod2, c("feeding_type", "zt")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Group","DateTime"))
write.table(LSD_Test_fd2,"ldlr_late_zt_faithpd_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)
#############################################################################
set.seed(1234)
wk20_pc<-read_qza("2026_reprocess_16s/data/modern_beta/late_phyloRPCA-ordination.qza")

wk20_df<-wk20_pc$data$Vectors%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(.,metadata,by="sample_name")

p<-wk20_df%>%
  ggplot(aes(x=PC1, y=PC2, fill=genotype, shape=genotype, group=genotype)) +
  geom_point(alpha=1.0, size=4)+
  theme_pubr() +
  scale_fill_manual(values=c("#ffffff","#59595b"))+
  scale_color_manual(values=c("#ffffff","#59595b"))+
  scale_shape_manual(values=c(22,21))+
  theme(legend.position = "right")+
  labs(x =paste("PC1 (",round(wk1_pc$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(wk1_pc$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       title="Apoe-/- vs. Ldlr-/-")
ggsave("phyloRPCA_apoeldlr_wk20.pdf", plot=p,height=4, width=6)

dismat<-fread("2026_reprocess_16s/data/modern_beta/late_phyloRPCA-dm/distance-matrix.tsv")
md_w20<-metadata%>%filter(sample_name %in% dismat$V1)
dismat<-dismat%>%column_to_rownames("V1")%>%as.matrix()
permanova_result <- adonis2(dismat ~ genotype, data = md_w20, permutations = 999, by="terms")
# adonis2(formula = dismat ~ genotype, data = md_w20, permutations = 999, by = "terms")
# Df SumOfSqs      R2      F Pr(>F)    
# genotype  1   29.798 0.15644 17.247  0.001 ***
#   Residual 93  160.674 0.84356                  
# Total    94  190.472 1.00000
#############################################################################

wk20a_pc<-read_qza("2026_reprocess_16s/data/modern_beta/ApoE-late_phyloRPCA-ordination.qza")

p<-make_rpcaw1(wk20a_pc,"Apoe-/-",22)
ggsave("phyloRPCA_apoezt_wk20.pdf", plot=p,height=4, width=6)

dismat<-fread("2026_reprocess_16s/data/modern_beta/ApoE-late_phyloRPCA-dm/distance-matrix.tsv")
md_apoew20<-metadata%>%filter(sample_name %in% dismat$V1)
dismat<-dismat%>%column_to_rownames("V1")%>%as.matrix()
permanova_result <- adonis2(dismat ~ feeding_type*zt, data = md_apoew20, permutations = 999, by="terms")
# adonis2(formula = dismat ~ feeding_type * zt, data = md_apoew20, permutations = 999, by = "terms")
# Df SumOfSqs      R2      F Pr(>F)  
# feeding_type     1    4.369 0.04042 2.1872  0.095 .
# zt               1    0.353 0.00326 0.1766  0.850  
# feeding_type:zt  1    3.493 0.03231 1.7486  0.172  
# Residual        50   99.882 0.92400                
# Total           53  108.098 1.00000

wk20l_pc<-read_qza("2026_reprocess_16s/data/modern_beta/LdlR-late_phyloRPCA-ordination.qza")

p<-make_rpcaw1(wk20l_pc,"Ldlr-/-",21)
ggsave("phyloRPCA_ldlrzt_wk20.pdf", plot=p,height=4, width=6)
dismat<-fread("2026_reprocess_16s/data/modern_beta/LdlR-late_phyloRPCA-dm/distance-matrix.tsv")
md_ldlrw20<-metadata%>%filter(sample_name %in% dismat$V1)
dismat<-dismat%>%column_to_rownames("V1")%>%as.matrix()
permanova_result <- adonis2(dismat ~ feeding_type*zt, data = md_ldlrw20, permutations = 999, by="terms")
# adonis2(formula = dismat ~ feeding_type * zt, data = md_ldlrw20, permutations = 999, by = "terms")
# Df SumOfSqs      R2      F Pr(>F)   
# feeding_type     1   11.935 0.12598 6.3656  0.004 **
#   zt               1    1.972 0.02082 1.0520  0.367   
# feeding_type:zt  1    0.211 0.00222 0.1123  0.905   
# Residual        43   80.626 0.85098                 
# Total           46   94.744 1.00000

wk1l_pc<-read_qza("2026_reprocess_16s/data/modern_beta/LdlR-early_phyloRPCA-ordination.qza")

p<-make_rpcaw1(wk1l_pc,"Ldlr-/-",21)
ggsave("phyloRPCA_ldlrzt_wk1.pdf", plot=p,height=4, width=6)
dismat<-fread("2026_reprocess_16s/data/modern_beta/LdlR-early_phyloRPCA-dm/distance-matrix.tsv")
md_ldlrw1<-metadata%>%filter(sample_name %in% dismat$V1)
dismat<-dismat%>%column_to_rownames("V1")%>%as.matrix()
permanova_result <- adonis2(dismat ~ feeding_type*zt, data = md_ldlrw1, permutations = 999, by="terms")
# adonis2(formula = dismat ~ feeding_type * zt, data = md_ldlrw1, permutations = 999, by = "terms")
# Df SumOfSqs      R2       F Pr(>F)    
# feeding_type     1   23.858 0.24808 15.6326  0.001 ***
#   zt               1    2.570 0.02673  1.6841  0.171    
# feeding_type:zt  1    2.590 0.02694  1.6973  0.190    
# Residual        44   67.150 0.69826                   
# Total           47   96.169 1.00000  

############################################################################
#qurro ldlr late
metadata<-fread("2026_reprocess_16s/data/metadata_2026.tsv")%>%
  mutate(zt = parse_number(ZT_time))
ldlr_lr<-fread("2026_reprocess_16s/data/modern_beta/qurro/LdlR-late_qurro_res-10perc.tsv")%>%
  dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,metadata,by="sample_name")

p<-ggplot(data=ldlr_lr, aes(x=zt, y=Current_Natural_Log_Ratio, colour=title2, fill=title2)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point", size = 2) + 
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values=c("#8c471f","#278c43"))+
  scale_color_manual(values=c("#8c471f","#278c43"))+
  scale_x_continuous(breaks = seq(0, 24, by = 4)) + 
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("tpbtm10perc_nat_log_ldlr-late_overzt.pdf", plot=p,height=4, width=4)

ldlr_lr_forstat<-ldlr_lr%>%mutate(title2=gsub("_","",title2))

mod <- aov(Current_Natural_Log_Ratio~ title2 * zt, data = ldlr_lr_forstat)
summary(mod)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# title2       1  21.79  21.794  12.372 0.00104 **
#   zt           1   3.60   3.603   2.045 0.15992   
# title2:zt    1   0.52   0.517   0.294 0.59065   
# Residuals   43  75.75   1.762    
LSD_Test_fd<- (LSD.test(mod, c("title2", "zt")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  tidyr::separate_wider_delim(group,delim = ":",names = c("Group", "zt"))
write.table(LSD_Test_fd,"ldlr-late_zt_logratio10perc_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#qurro apoe late
metadata<-fread("2026_reprocess_16s/data/metadata_2026.tsv")%>%
  mutate(zt = parse_number(ZT_time))
apoe_lr<-fread("2026_reprocess_16s/data/modern_beta/qurro/ApoE-late_qurro_res-10perc.tsv")%>%
  dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,metadata,by="sample_name")

p<-ggplot(data=apoe_lr, aes(x=zt, y=Current_Natural_Log_Ratio, colour=title2, fill=title2)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point", size = 2) + 
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values=c("#8c471f","#278c43"))+
  scale_color_manual(values=c("#8c471f","#278c43"))+
  scale_x_continuous(breaks = seq(0, 24, by = 4)) + 
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("tpbtm10perc_nat_log_apoe-late_overzt.pdf", plot=p,height=4, width=4)

apoe_lr_forstat<-apoe_lr%>%mutate(title2=gsub("_","",title2))

mod <- aov(Current_Natural_Log_Ratio~ title2 * zt, data = apoe_lr_forstat)
summary(mod)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# title2       1  13.70  13.701   8.183 0.00616 **
#   zt           1   6.03   6.031   3.602 0.06349 . 
# title2:zt    1   0.17   0.175   0.104 0.74807   
# Residuals   50  83.72   1.674   
LSD_Test_fda<- (LSD.test(mod, c("title2", "zt")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  tidyr::separate_wider_delim(group,delim = ":",names = c("Group", "zt"))
write.table(LSD_Test_fda,"apoe-late_zt_logratio10perc_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)


#############################################################################
ldlr_pc<-read_qza("2026_reprocess_16s/data/modern_beta/LdlR_phyloRPCA-ordination.qza")

ldlr_df<-ldlr_pc$data$Vectors%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(.,metadata,by="sample_name")

p<-ldlr_df%>%
  ggplot(aes(x=PC1, y=PC2, fill=title2, group=title2)) +
  geom_point(alpha=1.0, shape=21, size=4)+
  theme_pubr() +
  scale_fill_manual(values=c("#cd863c","#8c471f","#9bcb3d","#278c43"))+
  scale_color_manual(values=c("#cd863c","#8c471f","#9bcb3d","#278c43"))+
  
  theme(legend.position = "right")+
  labs(x =paste("PC1 (52.0%)",sep=""),
       y =paste("PC2 (39.2%)",sep=""),
       title="Ldlr-/-")
ggsave("phyloRPCA_ldlr_wks.pdf", plot=p,height=4, width=6)

dismat<-fread("2026_reprocess_16s/data/modern_beta/LdlR_phyloRPCA-dm/distance-matrix.tsv")
md_ldlr<-metadata%>%filter(sample_name %in% dismat$V1)
dismat<-dismat%>%column_to_rownames("V1")%>%as.matrix()
permanova_result <- adonis2(dismat ~ title2, data = md_ldlr, permutations = 999, by="terms")

#run pairwise permanova
groups <- unique(md_ldlr$title2)

pairwise_comparisons <- combn(groups, 2, simplify = FALSE)

pairwise_results <- map_dfr(pairwise_comparisons, function(pair) {
  
  md_sub <- md_ldlr %>%
    filter(title2 %in% pair)
  
  dist_sub <- dismat[rownames(dismat) %in% md_sub$sample_name,
                     rownames(dismat) %in% md_sub$sample_name]
  
  md_sub <- md_sub[match(rownames(dist_sub), md_sub$sample_name), ]
  
  res <- adonis2(dist_sub ~ title2,
                 data = md_sub,
                 permutations = 999,
                 by = "terms")
  
  tibble(
    group1 = pair[1],
    group2 = pair[2],
    F = res$F[1],
    R2 = res$R2[1],
    p_value = res$`Pr(>F)`[1]
  )
})

pairwise_results <- pairwise_results %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

#############################################################################
#qurro ldlr early vs. late
metadata<-fread("2026_reprocess_16s/data/metadata_2026.tsv")%>%
  mutate(zt = parse_number(ZT_time))
ldlr_lr<-fread("2026_reprocess_16s/data/modern_beta/qurro/LdlR_qurro_res-10perc.tsv")%>%
  dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,metadata,by="sample_name")

p<-ggplot(data=ldlr_lr, aes(x=zt, y=Current_Natural_Log_Ratio, colour=title2, fill=title2)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point", size = 2) + 
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values=c("#cd863a","#8c471f","#9bcb3e","#278c43"))+
  scale_color_manual(values=c("#cd863a","#8c471f","#9bcb3e","#278c43"))+
  scale_x_continuous(breaks = seq(0, 24, by = 4)) + 
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("tpbtm10perc_nat_log_ldlr_overzt.pdf", plot=p,height=4, width=4)

ldlr_lr_forstat<-ldlr_lr%>%mutate(title2=gsub("_","",title2))

mod <- aov(Current_Natural_Log_Ratio~ title2 * zt, data = ldlr_lr_forstat)
summary(mod)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# title2       3  74.07  24.690  33.533 1.69e-14 ***
#   zt           1  10.72  10.722  14.562 0.000253 ***
#   title2:zt    3   2.43   0.811   1.102 0.352994    
# Residuals   87  64.06   0.736 
LSD_Test_fd<- (LSD.test(mod, c("title2", "zt")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  tidyr::separate_wider_delim(group,delim = ":",names = c("Group", "zt"))
write.table(LSD_Test_fd,"ldl_zt_logratio10perc_pvals_anovaLSD.txt",sep = "\t",row.names = FALSE,quote=FALSE)
