# Supplementary Table 2. Eigenvector centrality of genera in the high
# animal and high plant co-occurrence networks, ordered by the difference.

source("R/00_setup.R")

net <- readRDS(file.path(PATHS$results, "networks.rds"))
cen <- net$centrality[order(-net$centrality$delta), ]

supp2 <- data.frame(
  Genus         = cen$taxon,
  `More central in` = ifelse(cen$delta > 0, "High animal", "High plant"),
  `High animal` = fmt_num(cen$animal),
  `High plant`  = fmt_num(cen$plant),
  Difference    = fmt_num(cen$delta),
  check.names = FALSE)

write_table(supp2, "supp_table2_network_centrality")
write_table(net$global, "network_global_properties")
