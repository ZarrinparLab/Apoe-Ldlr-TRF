setwd("~/scratch/maria_apoe")

library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(ape)
library(data.table)

# ── Load inputs ────────────────────────────────────────────────────────────────

tree <- read.tree("2026_reprocess_16s/data/gg2_feature-table_filt11400_rclr_phylo/counts_by_node_tree/tree.nwk")

# All internal nodes in the tree
all_nodes <- tibble(node = tree$node.label) %>%
  filter(!is.na(node), node != "") %>%
  mutate(node = trimws(node)) %>%
  distinct(node)

# Raw counts feature table: rows = ASVs (Feature ID), cols = samples
counts <- fread("2026_reprocess_16s/data/gg2_feature-table_filt11400.tsv") %>%   # <-- update path
  as.data.frame()

# Make Feature.ID the rowname index for easy lookup
rownames(counts) <- counts[[1]]
counts <- counts[, -1]   # drop the Feature.ID column; keep only sample columns
sample_cols <- colnames(counts)

# ── Tree helper functions ──────────────────────────────────────────────────────

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
  tree$edge[tree$edge[, 1] == node_index, 2]
}

# ── Compute log-ratios ─────────────────────────────────────────────────────────
# For each internal node in the tree and each sample:
#   ln( sum(counts_side1) ) - ln( sum(counts_side2) )
# Mirrors: np.log(df[list_a].sum(axis=1)) - np.log(df[list_b].sum(axis=1))

logratio_rows <- map(all_nodes$node, function(nlabel) {
  
  node_index <- get_node_index(tree, nlabel)
  if (is.na(node_index)) return(NULL)
  
  children <- get_children(tree, node_index)
  if (length(children) < 2) return(NULL)
  
  side1_tips <- get_side_tips(children[1], tree)
  side2_tips <- get_side_tips(children[2], tree)
  
  # Subset to ASVs actually present in the feature table
  side1_tips <- side1_tips[side1_tips %in% rownames(counts)]
  side2_tips <- side2_tips[side2_tips %in% rownames(counts)]
  
  # Sum counts across ASVs (rows) per sample (col) — mirrors pandas .sum(axis=1)
  # after transposing: side1/2_tips are the column list equivalents
  sum_side1 <- colSums(counts[side1_tips, , drop = FALSE])
  sum_side2 <- colSums(counts[side2_tips, , drop = FALSE])
  
  # 0.5 pseudocount avoids log(0) NaNs — standard in compositional 16S analysis
  log_ratio <- log(sum_side1 + 0.5) - log(sum_side2 + 0.5)
  
  # Return as a named vector with node label as the row identifier
  c(node = nlabel, setNames(log_ratio, sample_cols))
})

# ── Build output feature table ─────────────────────────────────────────────────
# Rows = nodes, columns = samples  (mirrors QIIME2 feature table orientation)

logratio_table <- bind_rows(logratio_rows) %>%
  column_to_rownames("node")   # optional: keeps node as rowname rather than column

# Convert sample columns back to numeric (bind_rows coerces to character)
logratio_table <- logratio_table %>%
  mutate(across(everything(), as.numeric))

# ── Write output ───────────────────────────────────────────────────────────────

out_path <- "2026_reprocess_16s/data/gg2_feature-table_filt11400_nodenorm.tsv"

write.table(logratio_table,out_path,sep= "\t",row.names = TRUE,quote= FALSE)

message("Done. Output written to: ", out_path)
message("Dimensions: ", nrow(logratio_table), " nodes x ", ncol(logratio_table), " samples")
