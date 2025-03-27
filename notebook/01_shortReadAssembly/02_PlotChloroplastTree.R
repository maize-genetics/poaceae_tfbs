# Plots a newick tree from RaxML (matK locus)
# Charlie Hale, 2024.04.22
library(ape)
library(Cairo)

# Read tree file
par(mar = c(0,0,0,0))

# Plot the tree
plot(tree, type = "cladogram", edge.color = "blue", edge.width = 2, cex = 0.5)
tree <- read.tree(file = "/workdir/coh22/poaceae_tfbs/output/chloroplastTree/RAxML_bestTree.matk_bold_tree")

library(ggtree)
ggtree(tree, layout = "rectangular") + geom_tree(linewidth = 0.05) + 
  geom_tiplab(size = 0.1) + 
  theme_tree2() + 
  scale_color_manual(values=c('black', '#976391'))
ggsave("/workdir/coh22/poaceae_tfbs/output/large_phylogenetic_tree.png", 
       width = 32, height = 24, units = "in", dpi = 1000, type = "cairo-png")

# Display the interactive plot
p_plotly

plot(tree, type = "cladogram", edge.color = "blue", edge.width = 2, cex = 0.2)
png("/workdir/coh22/poaceae_tfbs/output/chloroplastTree/large_phylo_tree.png", width = 1600, height = 1200, res = 300)
plot(tree, type = "cladogram", edge.color = "blue", edge.width = 2, cex = 0.1)
dev.off()