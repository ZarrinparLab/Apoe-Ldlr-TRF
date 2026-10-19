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

bdm<-fread("2026_reprocess_16s/data/BIRDMAn/metacycle_birdman_results.csv")
clk<-fread("2026_reprocess_16s/data/BIRDMAn/top10_phase_shifters.csv")
