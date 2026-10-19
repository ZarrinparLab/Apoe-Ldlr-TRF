#in qiime2-amplicon-2026.1, run in terminal

##get rclr table by node
qiime gemelli phylogenetic-rclr-transformation \
  --i-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/gg2_feature-table_filt11400.qza \
  --i-phylogeny /mnt/zarrinpar/Pynchon/Databases/greengenes-2024.9/2024.09.phylogeny.asv.nwk.qza \
  --output-dir /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/gg2_feature-table_filt11400_rclr_phylo

##get rclr table by ASV
qiime gemelli rclr-transformation \
  --i-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/gg2_feature-table_filt11400.qza \
  --output-dir /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/gg2_feature-table_filt11400_rclr

##get clr table by ASV
qiime gemelli clr-transformation \
  --i-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/gg2_feature-table_filt11400.qza \
  --output-dir /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/gg2_feature-table_filt11400_clr


#get wk1 subset
qiime feature-table filter-samples \
    --i-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/gg2_feature-table_filt11400.qza \
    --m-metadata-file /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/metadata_2026.tsv \
    --p-where '[host_age]= 11.5' \
    --o-filtered-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/gg2_ftable11400_early.qza
  
#get wk20 subset
qiime feature-table filter-samples \
    --i-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/gg2_feature-table_filt11400.qza \
    --m-metadata-file /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/metadata_2026.tsv \
    --p-where '[host_age]= 29.5' \
    --o-filtered-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/gg2_ftable11400_late.qza
    
    
#get wk1 phyloRPCA
qiime gemelli phylogenetic-rpca-with-taxonomy \
    --i-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/gg2_ftable11400_early.qza \
    --i-phylogeny /mnt/zarrinpar/Pynchon/Databases/greengenes-2024.9/2024.09.phylogeny.asv.nwk.qza \
    --m-taxonomy-file /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/ApoeLdlr_taxonomy.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot  /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/early_phyloRPCA-ordination.qza \
    --o-distance-matrix /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/early_phyloRPCA-dm.qza \
    --o-counts-by-node-tree /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/early_phyloRPCA-tree.qza \
    --o-counts-by-node /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/early_phyloRPCA-table.qza \
    --o-t2t-taxonomy /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/early_phyloRPCA-taxonomy.qza
    
#get wk20 phyloRPCA
qiime gemelli phylogenetic-rpca-with-taxonomy \
    --i-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/gg2_ftable11400_late.qza \
    --i-phylogeny /mnt/zarrinpar/Pynchon/Databases/greengenes-2024.9/2024.09.phylogeny.asv.nwk.qza \
    --m-taxonomy-file /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/ApoeLdlr_taxonomy.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot  /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/late_phyloRPCA-ordination.qza \
    --o-distance-matrix /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/late_phyloRPCA-dm.qza \
    --o-counts-by-node-tree /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/late_phyloRPCA-tree.qza \
    --o-counts-by-node /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/late_phyloRPCA-table.qza \
    --o-t2t-taxonomy /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/late_phyloRPCA-taxonomy.qza
    
#qurro need to do this in older qiime version (qurro-qiime2-2023.2)

#ldlr
qiime gemelli phylogenetic-rpca-with-taxonomy \
    --i-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/gg2_ftable11400_LdlR.qza \
    --i-phylogeny /mnt/zarrinpar/Pynchon/Databases/greengenes-2024.9/2024.09.phylogeny.asv.nwk.qza \
    --m-taxonomy-file /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/phyloRPCA-taxonomy.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot  /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR_phyloRPCA-ordination.qza \
    --o-distance-matrix /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR_phyloRPCA-dm.qza \
    --o-counts-by-node-tree /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR_phyloRPCA-tree.qza \
    --o-counts-by-node /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR_phyloRPCA-table.qza \
    --o-t2t-taxonomy /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR_phyloRPCA-taxonomy.qza
    
qiime qurro loading-plot \
    --i-ranks /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR_phyloRPCA-ordination.qza \
    --i-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR_phyloRPCA-table.qza\
    --m-sample-metadata-file /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/metadata_2026.tsv \
    --m-feature-metadata-file /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR_phyloRPCA-taxonomy.qza \
    --o-visualization /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR_phyloRPCA-qurro.qzv
    
#ldlr w20
qiime gemelli phylogenetic-rpca-with-taxonomy \
    --i-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/gg2_ftable11400_LdlR_late.qza \
    --i-phylogeny /mnt/zarrinpar/Pynchon/Databases/greengenes-2024.9/2024.09.phylogeny.asv.nwk.qza \
    --m-taxonomy-file /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/phyloRPCA-taxonomy.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot  /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR-late_phyloRPCA-ordination.qza \
    --o-distance-matrix /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR-late_phyloRPCA-dm.qza \
    --o-counts-by-node-tree /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR-late_phyloRPCA-tree.qza \
    --o-counts-by-node /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR-late_phyloRPCA-table.qza \
    --o-t2t-taxonomy /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR-late_phyloRPCA-taxonomy.qza
  
qiime qurro loading-plot \
    --i-ranks /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR-late_phyloRPCA-ordination.qza \
    --i-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR-late_phyloRPCA-table.qza\
    --m-sample-metadata-file /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/metadata_2026.tsv \
    --m-feature-metadata-file /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/phyloRPCA-taxonomy.qza \
    --o-visualization /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/LdlR-late_phyloRPCA-qurro.qzv
    
#apoe w20
qiime gemelli phylogenetic-rpca-with-taxonomy \
    --i-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/gg2_ftable11400_ApoE_late.qza \
    --i-phylogeny /mnt/zarrinpar/Pynchon/Databases/greengenes-2024.9/2024.09.phylogeny.asv.nwk.qza \
    --m-taxonomy-file /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/phyloRPCA-taxonomy.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot  /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/ApoE-late_phyloRPCA-ordination.qza \
    --o-distance-matrix /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/ApoE-late_phyloRPCA-dm.qza \
    --o-counts-by-node-tree /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/ApoE-late_phyloRPCA-tree.qza \
    --o-counts-by-node /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/ApoE-late_phyloRPCA-table.qza \
    --o-t2t-taxonomy /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/ApoE-late_phyloRPCA-taxonomy.qza
  
qiime qurro loading-plot \
    --i-ranks /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/ApoE-late_phyloRPCA-ordination.qza \
    --i-table /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/ApoE-late_phyloRPCA-table.qza\
    --m-sample-metadata-file /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/metadata_2026.tsv \
    --m-feature-metadata-file /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/phyloRPCA-taxonomy.qza \
    --o-visualization /mnt/zarrinpar/scratch/sfloresr/maria_apoe/2026_reprocess_16s/data/modern_beta/qurro/ApoE-late_phyloRPCA-qurro.qzv