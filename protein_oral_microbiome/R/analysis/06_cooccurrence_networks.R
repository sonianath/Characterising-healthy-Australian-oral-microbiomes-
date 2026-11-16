# Genus co-occurrence networks for the high animal and high plant groups.
# Edges: Spearman correlations of CLR abundances, BH-adjusted p < 0.05 and
# |rho| >= 0.3. Both networks use the same genera (present in >=50% of
# samples in each group). Exploratory, given the group sizes.

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(phyloseq); library(igraph); library(psych)
})

ps <- readRDS(file.path(PATHS$results, "phyloseq.rds"))$ps

ps_genus <- tax_glom(ps, taxrank = "Genus")
taxa_names(ps_genus) <- make.unique(as.character(tax_table(ps_genus)[, "Genus"]))
ps_net <- prune_samples(
  sample_data(ps_genus)$protein_group %in% c("High animal", "High plant"),
  ps_genus)

otu <- as(otu_table(ps_net), "matrix")
if (taxa_are_rows(ps_net)) otu <- t(otu)
grp <- as.character(sample_data(ps_net)$protein_group)

keep <- colMeans(otu[grp == "High animal", ] > 0) >= 0.5 &
        colMeans(otu[grp == "High plant",  ] > 0) >= 0.5
otu <- otu[, keep]
cat("Genera:", ncol(otu), "\n"); print(table(grp))

build_net <- function(counts) {
  ct <- corr.test(clr_transform(counts), method = "spearman",
                  adjust = "BH", ci = FALSE)
  A <- ct$r
  A[ct$p > 0.05] <- 0          # ct$p holds BH-adjusted p above the diagonal
  A[abs(A) < 0.3] <- 0
  diag(A) <- 0
  graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE,
                              diag = FALSE)
}

g_animal <- build_net(otu[grp == "High animal", ])
g_plant  <- build_net(otu[grp == "High plant",  ])

net_props <- function(g, label) {
  gw <- g
  E(gw)$weight <- abs(E(gw)$weight)
  comp <- components(gw)
  data.frame(
    network      = label,
    n_nodes      = vcount(gw),
    n_edges      = ecount(gw),
    density      = edge_density(gw),
    mean_degree  = mean(degree(gw)),
    modularity   = modularity(cluster_fast_greedy(gw)),
    n_components = comp$no,
    largest_component = max(comp$csize),
    pct_negative = 100 * mean(E(g)$weight < 0))
}
global <- rbind(net_props(g_animal, "High animal"),
                net_props(g_plant,  "High plant"))
print(global)

taxa <- colnames(otu)
centrality <- data.frame(
  taxon  = taxa,
  animal = eigen_centrality(g_animal)$vector[taxa],
  plant  = eigen_centrality(g_plant)$vector[taxa])
centrality$delta <- centrality$animal - centrality$plant

saveRDS(list(global = global, centrality = centrality,
             n = table(grp), n_genera = ncol(otu)),
        file.path(PATHS$results, "networks.rds"))
