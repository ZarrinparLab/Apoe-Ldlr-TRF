setwd("~/scratch/maria_apoe")

library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(ape)

tree <- read.tree("2026_reprocess_16s/data/gg2_feature-table_filt11400_rclr_phylo/counts_by_node_tree/tree.nwk")
tax <- fread("2026_reprocess_16s/data/ApoeLdlr_taxonomy/taxonomy.tsv")%>%dplyr::rename(Feature.ID=`Feature ID`)
#sig_nodes <- fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LTl.csv")%>%filter(q<0.05)%>%dplyr::rename(node=test)
sig_nodes <- fread("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LT.csv")%>%filter(q<0.05)%>%dplyr::rename(node=test)

tax <- tax %>%
  mutate(
    Genus = str_extract(Taxon, "g__[^;]+"),
    Genus = str_remove(Genus, "g__"),
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)
  ) %>%
  mutate(
    Family = str_extract(Taxon, "f__[^;]+"),
    Family = str_remove(Family, "f__"),
    Family = ifelse(is.na(Family) | Family == "", "Unclassified", Family)
  )%>%
  mutate(
    Species = str_extract(Taxon, "s__[^;]+"),
    Species = str_remove(Species, "s__"),
    Species = ifelse(is.na(Species) | Species == "", "Unclassified", Species)
  )

get_node_index <- function(tree, node_label) {
  idx <- which(tree$node.label == node_label)
  if (length(idx) == 0) return(NA)
  idx + Ntip(tree)
}

get_side_tips <- function(child, tree) {
  if (child <= Ntip(tree)) {
    tree$tip.label[child]
  } else {
    extract.clade(tree, child)$tip.label
  }
}

get_children <- function(tree, node_index) {
  tree$edge[tree$edge[,1] == node_index, 2]
}

sig_nodes <- read_csv("2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LT.csv") %>%
  filter(q<0.05)%>%
  rename(node = test) %>%
  mutate(node = trimws(node)) %>%
  distinct(node, .keep_all = TRUE)

sig_nodes_asv <- sig_nodes%>%filter(!str_detect(node, "n"))%>%
  dplyr::rename(Feature.ID=node)%>%
  left_join(.,tax,by="Feature.ID")


node_species_summary <- map_df(sig_nodes$node, function(nlabel) {
  
  node_index <- get_node_index(tree, nlabel)
  if (is.na(node_index)) return(NULL)
  
  children <- get_children(tree, node_index)
  
  side1 <- get_side_tips(children[1], tree)
  side2 <- get_side_tips(children[2], tree)
  
  side1_species <- tax %>%
    filter(Feature.ID %in% side1) %>%
    count(Species, sort = TRUE)
  
  side2_species <- tax %>%
    filter(Feature.ID %in% side2) %>%
    count(Species, sort = TRUE)
  
  tibble(
    node = nlabel,
    side1_top_species = side1_species$Species[1],
    side2_top_species = side2_species$Species[1],
    side1_nASVs = length(side1),
    side2_nASVs = length(side2)
  )
})

write.table(node_species_summary,"2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LTl_qless0.05_nodeannot_species.csv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(node_species_summary,"2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LT_qless0.05_nodeannot_species.csv",sep = "\t",row.names = FALSE, quote=FALSE)


node_genus_summary <- map_df(sig_nodes$node, function(nlabel) {
  
  node_index <- get_node_index(tree, nlabel)
  if (is.na(node_index)) return(NULL)
  
  children <- get_children(tree, node_index)
  
  side1 <- get_side_tips(children[1], tree)
  side2 <- get_side_tips(children[2], tree)
  
  side1_genera <- tax %>%
    filter(Feature.ID %in% side1) %>%
    count(Genus, sort = TRUE)
  
  side2_genera <- tax %>%
    filter(Feature.ID %in% side2) %>%
    count(Genus, sort = TRUE)
  
  tibble(
    node = nlabel,
    side1_top_genus = side1_genera$Genus[1],
    side2_top_genus = side2_genera$Genus[1],
    side1_nASVs = length(side1),
    side2_nASVs = length(side2)
  )
})

write.table(node_genus_summary,"2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LTl_qless0.05_nodeannot_genus.csv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(node_genus_summary,"2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LT_qless0.05_nodeannot_genus.csv",sep = "\t",row.names = FALSE, quote=FALSE)


sig_nodes_asv <- sig_nodes%>%filter(!str_detect(node, "n"))%>%
  dplyr::rename(Feature.ID=node)%>%
  left_join(.,tax,by="Feature.ID")


node_family_summary <- map_df(sig_nodes$node, function(nlabel) {
  
  node_index <- get_node_index(tree, nlabel)
  if (is.na(node_index)) return(NULL)
  
  children <- get_children(tree, node_index)
  
  side1 <- get_side_tips(children[1], tree)
  side2 <- get_side_tips(children[2], tree)
  
  side1_family <- tax %>%
    filter(Feature.ID %in% side1) %>%
    count(Family, sort = TRUE)
  
  side2_family <- tax %>%
    filter(Feature.ID %in% side2) %>%
    count(Family, sort = TRUE)
  
  tibble(
    node = nlabel,
    side1_top_family = side1_family$Family[1],
    side2_top_family = side2_family$Family[1],
    side1_nASVs = length(side1),
    side2_nASVs = length(side2)
  )
})

write.table(node_family_summary,"2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LTl_qless0.05_nodeannot_family.csv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(node_family_summary,"2026_reprocess_16s/data/cosinor/node_analysis/cosinor_LT_qless0.05_nodeannot_family.csv",sep = "\t",row.names = FALSE, quote=FALSE)


